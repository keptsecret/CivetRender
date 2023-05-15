#version 420 core
out vec4 FragColor;

in vec3 TexCoords;

uniform samplerCube skytexture;

void main() {
    FragColor = texture(skytexture, TexCoords.rgb);
}