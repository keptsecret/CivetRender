#version 420 core

struct Material {
    sampler2D texture_diffuse1;
    sampler2D texture_specular1;
    float shininess;

    sampler2D texture_normal;
    bool use_normal_map;
};

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;
in mat3 TBN;

layout (location = 0) out vec3 FragPosOut;
layout (location = 1) out vec3 DiffuseOut;
layout (location = 2) out vec3 SpecularOut;
layout (location = 3) out vec3 NormalOut;

uniform Material material;

void main() {
    FragPosOut = FragPos;
    DiffuseOut = texture(material.texture_diffuse1, TexCoords).xyz;
    SpecularOut = texture(material.texture_specular1, TexCoords).xyz;

    NormalOut = Normal;
    if (material.use_normal_map) {
        NormalOut = texture(material.texture_normal, TexCoords).rgb;
        NormalOut = NormalOut * 2.0 - 1.0;
    }
    NormalOut = normalize(TBN * NormalOut);     // normals in world space
}