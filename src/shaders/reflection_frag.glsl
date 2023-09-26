#version 420 core
out vec4 FragColor;

uniform sampler2D PositionMap;
uniform sampler2D AlbedoMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D NormalMap;
uniform sampler2D ReflectedMap;

uniform vec2 screenSize;
uniform vec3 viewPos;

float fresnelSchlickRoughness(float cosTheta, float F0, float roughness);
vec3 blur(sampler2D image, vec2 uv, vec2 direction);

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 worldPos = texture(PositionMap, texCoords).xyz;

    vec3 albedo = texture(AlbedoMap, texCoords).rgb;
    float roughness = texture(AORoughMetallicMap, texCoords).g;
    roughness = max(roughness, 0.025);  // going to 0 causes flickering and errors
    float metallic = texture(AORoughMetallicMap, texCoords).b;

    vec4 reflectedColor = texture(ReflectedMap, texCoords);
    vec3 reflection = reflectedColor.rgb;
    float fade = reflectedColor.a;
    vec3 blurred = blur(ReflectedMap, texCoords, vec2(5, 0));

    float F0 = mix(0.04, 1.0, metallic);
    float roughnessFade = smoothstep(0.4, 0.7, 1.0 - roughness);

    vec3 N = normalize(texture(NormalMap, texCoords).xyz);
    vec3 V = normalize(viewPos - worldPos);
    float F = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, roughness);

    FragColor.rgb = mix(reflection, blurred, roughness);
    FragColor.a = F * roughnessFade * fade;
}

float fresnelSchlickRoughness(float cosTheta, float F0, float roughness) {
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

vec3 blur(sampler2D image, vec2 texCoord, vec2 direction) {
    const float TwoPi = 6.28318530718;
    const float Directions = 16.0;
    const float Quality = 4.0;
    const float Size = 16.0;

    vec2 radius = Size / textureSize(image, 0);
    vec3 color = texture(image, texCoord).rgb;

    for(float d = 0.0; d < TwoPi; d += TwoPi / Directions) {
        for(float i = 1.001 / Quality; i < 1.001; i += 1.0 / Quality) {
            color += texture(image, texCoord + vec2(cos(d), sin(d)) * radius * i).rgb;
        }
    }

    color /= Quality * Directions + 1.0;
    return color;
}