#include <core/material.h>

#include <core/primitive.h>
#include <core/spectrum.h>
#include <core/texture.h>
#include <utils/reflection.h>

namespace civet {

Material::~Material() {}

void Material::bump(const std::shared_ptr<Texture<float>>& map, SurfaceInteraction* si) {
	// Compute offset positions and evaluate displacement texture
	SurfaceInteraction si_eval = *si;

	// Shift _si_eval_ _du_ in the $u$ direction
	float du = .5f * (std::abs(si->dudx) + std::abs(si->dudy));
	// The most common reason for du to be zero is for ray that start from
	// light sources, where no differentials are available. In this case,
	// we try to choose a small enough du so that we still get a decently
	// accurate bump value.
	if (du == 0) {
		du = .0005f;
	}
	si_eval.p = si->p + du * si->shading.dpdu;
	si_eval.uv = si->uv + Vector2f(du, 0.f);
	si_eval.n = normalize((Normal3f)cross(si->shading.dpdu, si->shading.dpdv) + du * si->dndu);
	float u_displace = map->evaluate(si_eval);

	// Shift _si_eval_ _dv_ in the $v$ direction
	float dv = .5f * (std::abs(si->dvdx) + std::abs(si->dvdy));
	if (dv == 0) {
		dv = .0005f;
	}
	si_eval.p = si->p + dv * si->shading.dpdv;
	si_eval.uv = si->uv + Vector2f(0.f, dv);
	si_eval.n = normalize((Normal3f)cross(si->shading.dpdu, si->shading.dpdv) + dv * si->dndv);
	float v_displace = map->evaluate(si_eval);
	float displace = map->evaluate(*si);

	// Compute bump-mapped differential geometry
	Vector3f dpdu = si->shading.dpdu +
			(u_displace - displace) / du * Vector3f(si->shading.n) + displace * Vector3f(si->shading.dndu);
	Vector3f dpdv = si->shading.dpdv +
			(v_displace - displace) / dv * Vector3f(si->shading.n) + displace * Vector3f(si->shading.dndv);
	si->setShadingGeometry(dpdu, dpdv, si->shading.dndu, si->shading.dndv, false);
}

void GLMaterial::bind(Shader& shader, unsigned int tex_offset) {
	// unsigned int diffuse_no = 1;
	// unsigned int specular_no = 1;
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + (i + tex_offset));
		std::string name = textures[i]->type;

		// std::string number;
		// if (name == "texture_diffuse") {
		//	number = std::to_string(diffuse_no++);
		//	} else if (name == "texture_specular") {
		//		number = std::to_string(specular_no++);
		//	} else if (name == "texture_normal") {
		//		number = "";
		//	} else if (name == "texture_bump") {
		//		number = "";
		//	}

		shader.setInt(("material." + name).c_str(), i + tex_offset);
		glBindTexture(GL_TEXTURE_2D, textures[i]->id);
	}
	glActiveTexture(GL_TEXTURE0);

	if (!use_albedo_map) {
		shader.setVec3("material.albedo", albedo);
	}
	if (!use_metallic_map) {
		shader.setFloat("material.metallic", metallic);
	}
	if (!use_roughness_map) {
		shader.setFloat("material.roughness", roughness);
	}
	if (!use_ao_map) {
		shader.setFloat("material.ambient", ambient);
	}

	shader.setBool("material.use_albedo_map", use_albedo_map);
	shader.setBool("material.use_metallic_map", use_metallic_map);
	shader.setBool("material.use_roughness_map", use_roughness_map);
	shader.setBool("material.use_ao_map", use_ao_map);
	shader.setBool("material.is_glossy_rough", is_glossy_rough);
	shader.setBool("material.use_normal_map", use_normal_map);
	shader.setBool("material.use_bump_map", use_bump_map);
	shader.setFloat("material.bump_scale", bump_scale);
}

void GLMaterial::unbind(unsigned int tex_offset) {
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + (i + tex_offset));
		glBindTexture(GL_TEXTURE_2D, 0);
	}
}

} // namespace civet