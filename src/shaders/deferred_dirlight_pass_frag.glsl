#version 420 core
out vec4 FragColor;

layout (std140, binding = 1) uniform LightSpaceMatrices {
    mat4 lightSpaceMatrices[4];
};

struct DirLight {
    bool valid;
    vec3 direction;
    vec3 color;

    bool use_cascaded_shadows;
    float far_plane;
    int cascade_count;
    float cascade_distances[4]; // limit to 4 cascade levels
    sampler2DArray shadow_cascades;

    mat4 light_space_mat;   // only if shadow not cascade
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

uniform mat4 cam_view;
uniform vec3 viewPos;
uniform DirLight light;

const float PI = 3.14159265359;

float distributionGGX(vec3 N, vec3 H, float roughness);
float geometrySchlickGGX(float NdotV, float roughness);
float geometrySmith(vec3 N, vec3 V, vec3 L, float roughness);
vec3 fresnelSchlick(float cosTheta, vec3 F0);

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

float calcShadow(vec3 worldPos, vec3 normal) {
    vec4 fragPosLightSpace = vec4(1.0);
    int layer = -1;

    if (light.use_cascaded_shadows) {
        vec4 fragPosViewSpace = cam_view * vec4(worldPos, 1.0);
        float depthValue = abs(fragPosViewSpace.z);

        vec4 res = step(depthValue, vec4(light.cascade_distances[0], light.cascade_distances[1], light.cascade_distances[2], light.cascade_distances[3]));
        layer = light.cascade_count - int(res.x + res.y + res.z + res.w);

        fragPosLightSpace = lightSpaceMatrices[layer] * vec4(worldPos, 1.0);
    } else {
        fragPosLightSpace = light.light_space_mat * vec4(worldPos, 1.0);
    }

    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;
    float currentDepth = projCoords.z;

    float bias = max(0.05 * (1.0 - dot(normal, -light.direction)), 0.005);
    if (light.use_cascaded_shadows) {
        if (layer == light.cascade_count) {
            bias *= 1.0 / (light.far_plane * 0.5);
        } else {
            bias *= 1.0 / (light.cascade_distances[layer] * 0.5);
        }
    }

    float shadow = 0.0;
    if (!(currentDepth > 1.0)) {
        vec2 texelSize = vec2(1.0);
        if (light.use_cascaded_shadows) {
            texelSize /= vec2(textureSize(light.shadow_cascades, 0));
        } else {
            texelSize /= textureSize(light.shadow_map, 0);
        }

        // PCF, 4x4 filter
        for (int x = -1; x <= 2; x++) {
            for (int y = -1; y <= 2; y++) {
                float pcfDepth = 1.0;
                if (light.use_cascaded_shadows) {
                    pcfDepth = texture(light.shadow_cascades, vec3(projCoords.xy + vec2(x, y) * texelSize, layer)).r;
                } else {
                    pcfDepth = texture(light.shadow_map, projCoords.xy + vec2(x, y) * texelSize).r;
                }
                shadow += currentDepth - bias > pcfDepth ? 1.0 : 0.0;
            }
        }
        shadow /= 16.0;
    }

    return shadow;
}

vec3 calcDirLight(vec3 worldPos, vec3 normal, vec2 texCoords) {
    vec3 N = normal;
    vec3 V = normalize(viewPos - worldPos);

    vec3 albedo = texture(AlbedoMap, texCoords).xyz;
    float ao = texture(AORoughMetallicMap, texCoords).r;
    float roughness = texture(AORoughMetallicMap, texCoords).g;
    float metallic = texture(AORoughMetallicMap, texCoords).b;

    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);

    // reflectance equation
    vec3 Lo = vec3(0.0);
    vec3 L = normalize(-light.direction);
    vec3 H = normalize(V + L);
    // skip attenuation calculation because directional light is global
    vec3 radiance = light.color;

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
    Lo += (kD * albedo / PI + specular) * radiance * NdotL;

    // shadow
    float shadow = calcShadow(worldPos, normal);

    vec3 ambient = vec3(0.03) * albedo * ao;
    vec3 color = ambient + (1.0 - shadow) * Lo;

    return color;
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;
    vec3 normal = normalize(texture(NormalMap, texCoords).xyz);

    FragColor.rgb = calcDirLight(worldPos, normal, texCoords);
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