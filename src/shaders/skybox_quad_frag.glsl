#version 420 core
layout (location = 0) out vec4 colorOut;

in vec2 TexCoord;

uniform float distanceValue;

void main() {
    colorOut = vec4(distanceValue, distanceValue * distanceValue, 0.0, 0.0);
}