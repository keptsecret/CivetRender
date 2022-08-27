#version 420 core
out vec4 FragColor;

struct Material {
    sampler2D texture_diffuse1;
    sampler2D texture_specular1;
    float shininess;
};

struct DirLight {
    bool valid;
    vec3 direction;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    mat4 light_space_mat;
    sampler2D shadow_map;
};

struct PointLight {
    bool valid;
    vec3 position;

    float constant;
    float linear;
    float quadratic;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float far_plane;
    samplerCube shadow_map;
};

#define NUM_DIR_LIGHTS 4
#define NUM_POINT_LIGHTS 32

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;

uniform vec3 viewPos;
uniform DirLight dirLights[NUM_DIR_LIGHTS];
uniform PointLight pointLights[NUM_POINT_LIGHTS];
uniform Material material;

vec3 calcDirLight(DirLight light, vec3 normal, vec3 viewDir);
vec3 calcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir);
float calcShadow(DirLight light, vec4 fragPosLightSpace, vec3 normal, vec3 lightDir);
float calcShadowCube(PointLight light, vec3 fragPos, vec3 lightPos);

void main() {
    vec3 normal = normalize(Normal);
    vec3 viewDir = normalize(viewPos - FragPos);

    vec3 result = vec3(0,0,0);
    for (int i = 0; i < NUM_DIR_LIGHTS; i++) {
         if (dirLights[i].valid) {
            result += calcDirLight(dirLights[i], normal, viewDir);
        }
    }

    for (int i = 0; i < NUM_POINT_LIGHTS; i++) {
        if (pointLights[i].valid) {
            result += calcPointLight(pointLights[i], normal, FragPos, viewDir);
        }
    }

    vec4 fragColor = vec4(result, 1.0);

    // manual gamma correction
    float gamma = 2.2;
    FragColor.rgb = pow(fragColor.rgb, vec3(1.0 / gamma));

    // for debugging shadow maps
    //vec3 fragToLight = FragPos - pointLights[0].position;
    //float closestDepth = texture(pointLights[0].shadow_map, fragToLight).r;
    //FragColor = vec4(vec3(closestDepth), 1.0);
}

vec3 calcDirLight(DirLight light, vec3 normal, vec3 viewDir) {
    vec3 ambient = light.ambient * vec3(texture(material.texture_diffuse1, TexCoords));

    // diffuse shading
    vec3 lightDir = normalize(-light.direction);
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * vec3(texture(material.texture_diffuse1, TexCoords));

    // specular shading - using blinn-phong halfway vector
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);
    vec3 specular = light.specular * spec * vec3(texture(material.texture_specular1, TexCoords));

    // shadow
    vec4 fragPosLightSpace = light.light_space_mat * vec4(FragPos, 1.0);
    float shadow = calcShadow(light, fragPosLightSpace, normal, lightDir);

    return (ambient + (1.0 - shadow) * (diffuse + specular));
}

vec3 calcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir) {
    vec3 ambient = light.ambient * vec3(texture(material.texture_diffuse1, TexCoords));

    // diffuse shading
    vec3 lightDir = normalize(light.position - fragPos);
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * vec3(texture(material.texture_diffuse1, TexCoords));

    // specular shading - using blinn-phong halfway vector
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);
    vec3 specular = light.specular * spec * vec3(texture(material.texture_specular1, TexCoords));

    // attenuation
    float distance = length(light.position - fragPos);
    float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * (distance * distance));

    ambient *= attenuation;
    diffuse *= attenuation;
    specular *= attenuation;

    // shadow
    float shadow = calcShadowCube(light, fragPos, light.position);

    return (ambient + (1.0 - shadow) * (diffuse + specular));
}

float calcShadow(DirLight light, vec4 fragPosLightSpace, vec3 normal, vec3 lightDir) {
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;

    float shadow = 0.0;
    if (!(projCoords.z > 1.0)) {
        float bias = max(0.05 * (1.0 - dot(normal, lightDir)), 0.005);
        vec2 texelSize = 1.0 / textureSize(light.shadow_map, 0);
        float currentDepth = projCoords.z;

        // PCF, offset by 1 on each side
        for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <= 1; y++) {
                float pcfDepth = texture(light.shadow_map, projCoords.xy + vec2(x, y) * texelSize).r;
                shadow += currentDepth - bias > pcfDepth ? 1.0 : 0.0f;
            }
        }
        shadow /= 9.0;
    }

    return shadow;
}

vec3 sampleOffsetDirections[20] = vec3[]
(
vec3( 1,  1,  1), vec3( 1, -1,  1), vec3(-1, -1,  1), vec3(-1,  1,  1),
vec3( 1,  1, -1), vec3( 1, -1, -1), vec3(-1, -1, -1), vec3(-1,  1, -1),
vec3( 1,  1,  0), vec3( 1, -1,  0), vec3(-1, -1,  0), vec3(-1,  1,  0),
vec3( 1,  0,  1), vec3(-1,  0,  1), vec3( 1,  0, -1), vec3(-1,  0, -1),
vec3( 0,  1,  1), vec3( 0, -1,  1), vec3( 0, -1, -1), vec3( 0,  1, -1)
);

float calcShadowCube(PointLight light, vec3 fragPos, vec3 lightPos) {
    vec3 fragToLight = fragPos - lightPos;

    float currentDepth = length(fragToLight);
    float shadow = 0.0;
    float bias = 0.15;
    int samples = 20;
    float viewDistance = length(viewPos - fragPos);
    float radius = 0.05;

    for (int i = 0; i < samples; i++) {
        float closestDepth = texture(light.shadow_map, fragToLight + sampleOffsetDirections[i] * radius).r;
        closestDepth *= light.far_plane;
        if (currentDepth - bias > closestDepth) {
            shadow += 1.0;
        }
    }
    shadow /= float(samples);

    return shadow;
}