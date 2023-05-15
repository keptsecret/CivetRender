#include <core/material.h>

namespace civet {

Material::~Material() {}

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
}

void GLMaterial::unbind(unsigned int tex_offset) {
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + (i + tex_offset));
		glBindTexture(GL_TEXTURE_2D, 0);
	}
}

} // namespace civet