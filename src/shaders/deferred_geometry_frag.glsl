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
    MetallicRoughAOOut.r = material.use_roughness_map ? texture(material.texture_roughness, TexCoords).r : material.roughness;
    MetallicRoughAOOut.r = material.is_glossy_rough ? 1 - MetallicRoughAOOut.r : MetallicRoughAOOut.r;
    MetallicRoughAOOut.g = material.use_metallic_map ? texture(material.texture_metallic, TexCoords).g : material.metallic;
    MetallicRoughAOOut.b = material.use_ao_map ? texture(material.texture_ao, TexCoords).b : material.ambient;

    NormalOut = Normal;
    if (material.use_normal_map) {
        NormalOut = texture(material.texture_normal, TexCoords).rgb;
        NormalOut = NormalOut * 2.0 - 1.0;
        NormalOut = normalize(TBN * NormalOut); // normals in world space
    }
    if (material.use_bump_map) {
        // adapted from "Bump Mapping Unparametrized Surfaces on the GPU"
        vec3 vn = Normal;
        vec3 posDX = dFdx(FragPos.xyz);// choose dFdx (#version 420) or dFdxFine (#version 450) here
        vec3 posDY = dFdy(FragPos.xyz);
        vec3 r1 = cross(posDY, vn);
        vec3 r2 = cross(vn, posDX);
        float det = dot(posDX, r1);
        float Hll = texture(material.texture_bump, TexCoords).x;//-- height from bump map texture, tc=texture coordinates
        float Hlr = texture(material.texture_bump, TexCoords + dFdx(TexCoords.xy)).x;
        float Hul = texture(material.texture_bump, TexCoords + dFdy(TexCoords.xy)).x;
        // float dBs = ddx_fine ( height );     //-- optional explicit height
        // float dBt = ddy_fine ( height );

        // gradient of surface texture. dBs=Hlr-Hll, dBt=Hul-Hll
        vec3 surf_grad = sign(det) * ((Hlr - Hll) * r1 + (Hul - Hll)* r2);
        float bump_amt = 0.1;// bump_amt = adjustable bump amount
        vec3 vbumpnorm = vn*(1.0-bump_amt) + bump_amt * normalize (abs(det)*vn - surf_grad);// bump normal
        NormalOut = vbumpnorm;
    }
}