#version 420 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;

uniform vec3 viewPos;

void main() {
    vec3 lightColor = vec3(1.0, 1.0, 1.0);
    vec3 objectColor = vec3(0.4, 0.4, 0.4);

    vec3 ambient = 0.1 * lightColor;

    vec3 normal = normalize(Normal);
    vec3 viewDir = normalize(viewPos - FragPos);
    float diff = max(dot(normal, viewDir), 0.0);
    vec3 diffuse = diff * lightColor;

    vec3 result = (ambient + diffuse) * objectColor;
    vec4 fragColor = vec4(result, 1.0);

    // manual gamma correction
    float gamma = 2.2;
    FragColor.rgb = pow(fragColor.rgb, vec3(1.0 / gamma));
}