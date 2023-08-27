#version 420 core
out vec4 FragColor;

uniform sampler2D PositionMap;
uniform sampler2D AlbedoMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D NormalMap;
uniform sampler2D ReflectedMap;

uniform vec2 screenSize;
uniform vec3 viewPos;

vec3 fresnelSchlick(float cosTheta, vec3 F0);
vec3 blur(sampler2D image, vec2 uv, vec2 direction);

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;

    vec3 albedo = texture(AlbedoMap, texCoords).rgb;
    float roughness = texture(AORoughMetallicMap, texCoords).g;
    float metallic = texture(AORoughMetallicMap, texCoords).b;

    vec3 reflection = texture(ReflectedMap, texCoords).rgb;
    vec3 blurred = blur(ReflectedMap, texCoords, vec2(5, 0));

    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);

    vec3 N = normalize(texture(NormalMap, texCoords).xyz);
    vec3 V = normalize(viewPos - worldPos);
    vec3 F = fresnelSchlick(max(dot(N, V), 0.0), F0);

    FragColor = vec4(mix(reflection, blurred, roughness) * F, metallic);
}

vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

vec3 blur(sampler2D image, vec2 texCoord, vec2 direction) {
    const float TwoPi = 6.28318530718;
    const float Directions = 12.0;
    const float Quality = 4.0;
    const float Size = 8.0;

    vec2 radius = Size / textureSize(image, 0);
    vec3 color = texture(image, texCoord).rgb;

    for(float d = 0.0; d < TwoPi; d += TwoPi / Directions) {
        for(float i = 1.0 / Quality; i < 1.001; i += 1.0 / Quality) {
            color += texture(image, texCoord + vec2(cos(d), sin(d)) * radius * i).rgb;
        }
    }

    color /= Quality * Directions + 1.0;
    return color;
}