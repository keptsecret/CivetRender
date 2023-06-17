#version 420 core
in vec4 FragPos;

layout (location = 0) out vec4 distanceOut;

void main() {
    float distanceToFrag = length(FragPos.xyz);

    distanceOut = vec4(distanceToFrag, distanceToFrag * distanceToFrag, 0.0, 0.0);
}