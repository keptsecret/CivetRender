#version 420 core
out vec4 FragColor;

in vec3 TexCoords;

uniform samplerCube skytexture;

vec3 calcExposedColor(vec3 color) {
    const float exposure = -15.0;
    return exp2(exposure) * color;
}

void main() {
    FragColor.rgb = calcExposedColor(texture(skytexture, TexCoords.rgb).rgb);
}