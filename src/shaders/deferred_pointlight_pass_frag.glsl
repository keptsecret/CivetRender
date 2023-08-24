#version 420 core
out vec4 FragColor;

struct DirLight {
    bool valid;
    vec3 direction;
    vec3 color;

    mat4 light_space_mat;
    sampler2D shadow_map;
};

struct PointLight {
    bool valid;
    vec3 position;
    vec3 color;

    float constant;
    float linear;
    float quadratic;

    float radius;
    samplerCube shadow_map;
};

uniform sampler2D PositionMap;
uniform sampler2D AlbedoMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D NormalMap;
uniform vec2 screenSize;

uniform vec3 viewPos;
uniform PointLight light;

const float PI = 3.14159265359;
const uint NUM_PCF_SAMPLES = 32;
const vec2 Poisson[32] = vec2[](
    vec2(-0.975402, -0.0711386),
    vec2(-0.920347, -0.41142),
    vec2(-0.883908, 0.217872),
    vec2(-0.884518, 0.568041),
    vec2(-0.811945, 0.90521),
    vec2(-0.792474, -0.779962),
    vec2(-0.614856, 0.386578),
    vec2(-0.580859, -0.208777),
    vec2(-0.53795, 0.716666),
    vec2(-0.515427, 0.0899991),
    vec2(-0.454634, -0.707938),
    vec2(-0.420942, 0.991272),
    vec2(-0.261147, 0.588488),
    vec2(-0.211219, 0.114841),
    vec2(-0.146336, -0.259194),
    vec2(-0.139439, -0.888668),
    vec2(0.0116886, 0.326395),
    vec2(0.0380566, 0.625477),
    vec2(0.0625935, -0.50853),
    vec2(0.125584, 0.0469069),
    vec2(0.169469, -0.997253),
    vec2(0.320597, 0.291055),
    vec2(0.359172, -0.633717),
    vec2(0.435713, -0.250832),
    vec2(0.507797, -0.916562),
    vec2(0.545763, 0.730216),
    vec2(0.56859, 0.11655),
    vec2(0.743156, -0.505173),
    vec2(0.736442, -0.189734),
    vec2(0.843562, 0.357036),
    vec2(0.865413, 0.763726),
    vec2(0.872005, -0.927));

float distributionGGX(vec3 N, vec3 H, float roughness);
float geometrySchlickGGX(float NdotV, float roughness);
float geometrySmith(vec3 N, vec3 V, vec3 L, float roughness);
vec3 fresnelSchlick(float cosTheta, vec3 F0);

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

float calcShadowCube(vec3 fragPos, float cosTheta, vec3 normal) {
    vec3 fragToLight = fragPos - light.position;

    float currentDepth = length(fragToLight);
    float shadow = 0.0;
    float bias = 0.001;//max(0.05 * (1.0 - dot(normal, -fragToLight)), 0.005);
    float viewDistance = length(viewPos - fragPos);
    float radius = (1.0 + (viewDistance / light.radius)) / 25.0;

    vec3 T = vec3(1.0);
    vec3 t1 = cross(normal, vec3(0, 1, 0));
    vec3 t2 = cross(normal, vec3(0, 0, 1));
    if (length(t1) > length(t2)) {
        T = t1;
    } else {
        T = t2;
    }
    T = normalize(T);
    vec3 B = cross(T, normal);

    for (int i = 0; i < NUM_PCF_SAMPLES; i++) {
        vec2 P = Poisson[i];
        float closestDepth = texture(light.shadow_map, fragToLight + (T * P.x + B * P.y) * radius).r;
        closestDepth *= light.radius;
        if (currentDepth - bias > closestDepth) {
            shadow += 1.0;
        }
    }
    shadow /= float(NUM_PCF_SAMPLES);

    return shadow;
}

vec3 calcPointLight(vec3 worldPos, vec3 normal, vec2 texCoords) {
    vec3 N = normal;
    vec3 V = normalize(viewPos - worldPos);

    vec3 albedo = texture(AlbedoMap, texCoords).xyz;
    float ao = texture(AORoughMetallicMap, texCoords).r;
    float roughness = texture(AORoughMetallicMap, texCoords).g;
    float metallic = texture(AORoughMetallicMap, texCoords).b;

    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);

    // reflectance equation
    vec3 L = normalize(light.position - worldPos);
    vec3 H = normalize(V + L);
    float distance = length(light.position - worldPos);
    float dmr = distance / light.radius;
    float s = clamp(1 - (dmr * dmr * dmr * dmr), 0.0, 1.0);
    float s2 = s * s;
    float attenuation = s2 / (distance * distance + 1);
    vec3 radiance = light.color * attenuation;

    // cook-torrance brdf
    float NDF = distributionGGX(N, H, roughness);
    float G = geometrySmith(N, V, L, roughness);
    vec3 F = fresnelSchlick(max(dot(H, V), 0.0), F0);

    vec3 kS = F;
    vec3 kD = vec3(1.0) - kS;
    kD *= 1.0 - metallic;

    vec3 numerator = NDF * G * F;
    float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.0001;
    vec3 specular = numerator / denominator;

    // add to outgoing radiance Lo
    float NdotL = max(dot(N, L), 0.0);
    vec3 Lo = (kD * albedo / PI + specular) * radiance * NdotL;

    // shadow
    float shadow = calcShadowCube(worldPos, NdotL, normal);

    vec3 color = (1.0 - shadow) * Lo;

    return color;
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;
    vec3 normal = normalize(texture(NormalMap, texCoords).xyz);

    FragColor.rgb = calcPointLight(worldPos, normal, texCoords);
}

float distributionGGX(vec3 N, vec3 H, float roughness) {
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;

    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return num / denom;
}

float geometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return num / denom;
}

float geometrySmith(vec3 N, vec3 V, vec3 L, float roughness) {
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = geometrySchlickGGX(NdotV, roughness);
    float ggx1  = geometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}

vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}