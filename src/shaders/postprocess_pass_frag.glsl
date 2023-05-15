#version 420 core
out vec4 FragColor;

uniform sampler2D RawFinalImage;

uniform vec2 screenSize;

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

// ACES filmic tonemapping curve adapted to GLSL from: https://github.com/godotengine/godot/blob/master/drivers/gles3/shaders/tonemap_inc.glsl
vec3 aces_fitted(vec3 color) {
    const float exposure_bias = 1.8f;
    const float A = 0.0245786f;
    const float B = 0.000090537f;
    const float C = 0.983729f;
    const float D = 0.432951f;
    const float E = 0.238081f;

    // Exposure bias baked into transform to save shader instructions. Equivalent to `color *= exposure_bias`
    const mat3 rgb_to_rrt = mat3(
        vec3(0.59719f * exposure_bias, 0.35458f * exposure_bias, 0.04823f * exposure_bias),
        vec3(0.07600f * exposure_bias, 0.90834f * exposure_bias, 0.01566f * exposure_bias),
        vec3(0.02840f * exposure_bias, 0.13383f * exposure_bias, 0.83777f * exposure_bias));

    const mat3 odt_to_rgb = mat3(
        vec3(1.60475f, -0.53108f, -0.07367f),
        vec3(-0.10208f, 1.10813f, -0.00605f),
        vec3(-0.00327f, -0.07276f, 1.07602f));

    color *= rgb_to_rrt;
    vec3 color_tonemapped = (color * (color + A) - B) / (color * (C * color + D) + E);
    color_tonemapped *= odt_to_rgb;

    return clamp(color_tonemapped, 0.0, 1.0);
}

// Approximated ACES filmic tonemapping curve adapted to GLSL from: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 color) {
    color *= 0.6;
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    return clamp((color * (a * color + b)) / (color * (c * color + d) + e), 0.0, 1.0);
}

// Approximate conversion to srgb from http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
// Use instead of gamma correction
vec3 linear_to_srgb(vec3 color) {
    return max(vec3(1.055) * pow(color, vec3(0.416666667)) - vec3(0.055), vec3(0.0));
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 color = texture(RawFinalImage, texCoords).xyz;

    vec3 tonemapped = aces_fitted(color);

    FragColor.rgb = linear_to_srgb(tonemapped);
}