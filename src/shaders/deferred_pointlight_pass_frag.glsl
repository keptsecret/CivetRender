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

float distributionGGX(vec3 N, vec3 H, float roughness);
float geometrySchlickGGX(float NdotV, float roughness);
float geometrySmith(vec3 N, vec3 V, vec3 L, float roughness);
vec3 fresnelSchlick(float cosTheta, vec3 F0);

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

vec3 sampleOffsetDirections[20] = vec3[]
(
vec3( 1,  1,  1), vec3( 1, -1,  1), vec3(-1, -1,  1), vec3(-1,  1,  1),
vec3( 1,  1, -1), vec3( 1, -1, -1), vec3(-1, -1, -1), vec3(-1,  1, -1),
vec3( 1,  1,  0), vec3( 1, -1,  0), vec3(-1, -1,  0), vec3(-1,  1,  0),
vec3( 1,  0,  1), vec3(-1,  0,  1), vec3( 1,  0, -1), vec3(-1,  0, -1),
vec3( 0,  1,  1), vec3( 0, -1,  1), vec3( 0, -1, -1), vec3( 0,  1, -1)
);

float calcShadowCube(vec3 fragPos, float cosTheta) {
    vec3 fragToLight = fragPos - light.position;

    float currentDepth = length(fragToLight);
    float shadow = 0.0;
    float bias = 0.005 * tan(acos(cosTheta));
    bias = clamp(bias, 0.0, 0.2);
    int samples = 20;
    float viewDistance = length(viewPos - fragPos);
    float radius = (1.0 + (viewDistance / light.radius)) / 25.0;;

    for (int i = 0; i < samples; i++) {
        float closestDepth = texture(light.shadow_map, fragToLight + sampleOffsetDirections[i] * radius).r;
        closestDepth *= light.radius;
        if (currentDepth - bias > closestDepth) {
            shadow += 1.0;
        }
    }
    shadow /= float(samples);

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
    float s = min(distance / light.radius, 1);
    float s2 = s * s;
    float attenuation = (1.0 - s2) * (1.0 - s2) / (distance * distance);
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
    float shadow = calcShadowCube(worldPos, NdotL);

    vec3 ambient = vec3(0.03) * albedo * ao * attenuation; // should be 0.03, will revert when indirect lighting
    vec3 color = ambient + (1.0 - shadow) * Lo;

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