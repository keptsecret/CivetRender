#version 420 core
layout (location = 0) in vec3 aPos;

out vec2 TexCoord;

void main() {
    TexCoord = aPos.xy * vec2(0.5) + vec2(0.5);
    gl_Position = vec4(aPos.xy, 1.0, 1.0);
}