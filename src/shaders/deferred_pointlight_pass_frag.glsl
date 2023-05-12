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

    float radius;
    samplerCube shadow_map;
};

uniform sampler2D PositionMap;
uniform sampler2D DiffuseMap;
uniform sampler2D SpecularMap;
uniform sampler2D NormalMap;
uniform vec2 screenSize;

uniform vec3 viewPos;
uniform PointLight light;
uniform Material material;

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
    vec3 ambient = light.ambient * texture(DiffuseMap, texCoords).xyz;

    // diffuse shading
    vec3 lightDir = light.position - worldPos;
    float cosTheta = dot(normal, lightDir);
    float diff = max(cosTheta, 0.0);
    vec3 diffuse = light.diffuse * diff * texture(DiffuseMap, texCoords).xyz;

    // specular shading - using blinn-phong halfway vector
    vec3 viewDir = normalize(viewPos - worldPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);
    vec3 specular = light.specular * spec * texture(SpecularMap, texCoords).xyz;

    // attenuation
    float distance = length(lightDir);
    float s = min(distance / light.radius, 1);
    float s2 = s * s;
    float attenuation = (1.0 - s2) * (1.0 - s2) / (light.constant + light.linear * distance + light.quadratic * (distance * distance));

    ambient *= attenuation;
    diffuse *= attenuation;
    specular *= attenuation;

    // shadow
    float shadow = calcShadowCube(worldPos, cosTheta);

    return (ambient + (1.0 - shadow) * (diffuse + specular));
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;
    vec3 normal = normalize(texture(NormalMap, texCoords).xyz);

    FragColor.rgb = calcPointLight(worldPos, normal, texCoords);
}