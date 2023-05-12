#version 420 core
out vec4 FragColor;

struct Material {
    sampler2D texture_diffuse1;
    sampler2D texture_specular1;
    float shininess;

    sampler2D texture_normal;
    bool use_normal_map;
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

uniform sampler2D PositionMap;
uniform sampler2D DiffuseMap;
uniform sampler2D SpecularMap;
uniform sampler2D NormalMap;
uniform vec2 screenSize;

uniform vec3 viewPos;
uniform DirLight light;
uniform Material material;

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

float calcShadow(vec4 fragPosLightSpace, vec3 normal) {
    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;

    float shadow = 0.0;
    if (!(projCoords.z > 1.0)) {
        float bias = max(0.05 * (1.0 - dot(normal, -light.direction)), 0.005);
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

vec3 calcDirLight(vec3 worldPos, vec3 normal, vec2 texCoords) {
    vec3 ambient = light.ambient * texture(DiffuseMap, texCoords).xyz;

    // diffuse shading
    vec3 lightDir = normalize(-light.direction);
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * texture(DiffuseMap, texCoords).xyz;

    // specular shading - using blinn-phong halfway vector
    vec3 viewDir = normalize(viewPos - worldPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);
    vec3 specular = light.specular * spec * texture(SpecularMap, texCoords).xyz;

    // shadow
    vec4 fragPosLightSpace = light.light_space_mat * vec4(worldPos, 1.0);
    float shadow = calcShadow(fragPosLightSpace, normal);

    return (ambient + (1.0 - shadow) * (diffuse + specular));
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;
    vec3 normal = normalize(texture(NormalMap, texCoords).xyz);

    FragColor.rgb = calcDirLight(worldPos, normal, texCoords);
}