#version 420 core
out vec4 FragColor;

uniform sampler2D RawFinalImage;

uniform vec2 screenSize;

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

vec3 aces_input_matrix[3] = vec3[]
(
vec3(0.59719, 0.35458, 0.04823),
vec3(0.07600, 0.90834, 0.01566),
vec3(0.02840, 0.13383, 0.83777)
);

vec3 aces_output_matrix[3] = vec3[]
(
vec3( 1.60475, -0.53108, -0.07367),
vec3(-0.10208,  1.10813, -0.00605),
vec3(-0.00327, -0.07276,  1.07602)
);

vec3 mul(vec3[3] m, vec3 v) {
    float x = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
    float y = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
    float z = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
    return vec3(x, y, z);
}

vec3 rtt_and_odt_fit(vec3 color)
{
    vec3 a = color * (color + 0.0245786f) - 0.000090537f;
    vec3 b = color * (0.983729f * color + 0.4329510f) + 0.238081f;
    return a / b;
}

vec3 aces_fitted(vec3 color) {
    // ACES filmic tonemapping curve adapted to GLSL from: https://64.github.io/tonemapping/
    color = mul(aces_input_matrix, color);
    color = rtt_and_odt_fit(color);
    return mul(aces_output_matrix, color);
}

vec3 aces_approx(vec3 color) {
    // Approximated ACES filmic tonemapping curve adapted to GLSL from: https://64.github.io/tonemapping/
    color *= 0.6;
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    return clamp((color * (a * color + b)) / (color * (c * color + d) + e), 0.0, 1.0);
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 color = texture(RawFinalImage, texCoords).xyz;

    // apply tonemapping
    vec3 tonemapped = aces_approx(color);

    // apply gamma correction
    float gamma = 2.2;
    FragColor.rgb = pow(tonemapped, vec3(1.0 / gamma));
}