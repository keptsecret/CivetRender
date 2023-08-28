#version 420 core

struct Material {
    vec3 albedo;
    float metallic;
    float roughness;
    float ambient;

    sampler2D texture_albedo;
    bool use_albedo_map;

    sampler2D texture_metallic;
    bool use_metallic_map;

    sampler2D texture_roughness;
    bool use_roughness_map;

    sampler2D texture_ao;
    bool use_ao_map;

    bool is_glossy_rough;

    sampler2D texture_normal;
    bool use_normal_map;

    sampler2D texture_bump;
    float bump_scale;
    bool use_bump_map;
};

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoords;
in mat3 TBN;

layout (location = 0) out vec3 FragPosOut;
layout (location = 1) out vec3 AlbedoOut;
layout (location = 2) out vec3 MetallicRoughAOOut;
layout (location = 3) out vec3 NormalOut;

uniform Material material;

vec3 srgb_to_linear(vec3 color) {
    // Approximation from http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
    return color * (color * (color * 0.305306011 + 0.682171111) + 0.012522878);
}

void main() {
    vec4 diffuse = material.use_albedo_map ? texture(material.texture_albedo, TexCoords) : vec4(material.albedo, 1.0);
    if (diffuse.a < 0.1) {
        discard;
    }

    FragPosOut = FragPos;
    AlbedoOut = srgb_to_linear(diffuse.rgb);
    MetallicRoughAOOut.r = material.use_ao_map ? texture(material.texture_ao, TexCoords).r : material.ambient;
    MetallicRoughAOOut.g = material.use_roughness_map ? texture(material.texture_roughness, TexCoords).g : material.roughness;
    MetallicRoughAOOut.g = material.is_glossy_rough ? 1 - MetallicRoughAOOut.g : MetallicRoughAOOut.g;
    MetallicRoughAOOut.b = material.use_metallic_map ? texture(material.texture_metallic, TexCoords).b : material.metallic;

    NormalOut = normalize(Normal);
    if (material.use_normal_map) {
        NormalOut = texture(material.texture_normal, TexCoords).rgb;
        NormalOut = NormalOut * 2.0 - 1.0;
        NormalOut = normalize(TBN * NormalOut); // normals in world space
    }
    if (material.use_bump_map) {
        // adapted from https://www.shadertoy.com/view/MsScRt
        vec2 texelSize = 1. / textureSize(material.texture_bump, 0);
        float H = texture(material.texture_bump, TexCoords).r;
        float Hx = texture(material.texture_bump, TexCoords + dFdx(TexCoords.xy)).r;
        float Hy = texture(material.texture_bump, TexCoords + dFdy(TexCoords.xy)).r;
        vec2 dxy = H - vec2(Hx, Hy);

        vec3 bump = normalize(vec3(dxy * material.bump_scale / texelSize, 1.0));
        bump = normalize(bump * 0.5 + 0.5);
        NormalOut = normalize(TBN * bump);
    }
}