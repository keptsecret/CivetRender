#version 420 core
layout (triangles) in;
layout (triangle_strip, max_vertices=18) out;

uniform mat4 viewMatrices[6];

out vec4 FragPos;

void main() {
    for (int face = 0; face < 6; face++) {
        gl_Layer = face;
        for (int i = 0; i < 3; i++) {
            FragPos = viewMatrices[face] * gl_in[i].gl_Position;
            gl_Position = FragPos;
            EmitVertex();
        }
        EndPrimitive();
    }
}