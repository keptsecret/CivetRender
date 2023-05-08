#include <core/mesh.h>

#include <shapes/triangle.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>

namespace civet {

unsigned int loadTextureFromFile(const char* path, const std::string& directory, bool gamma) {
	std::string filename = std::string(path);
	filename = directory + '/' + filename;

	unsigned int texture_id;
	glGenTextures(1, &texture_id);

	int width, height, n_channels;
	unsigned char* data = stbi_load(filename.c_str(), &width, &height, &n_channels, 0);
	if (data) {
		GLenum format_in;
		GLenum format_out;
		if (n_channels == 1) {
			format_in = format_out = GL_RED;
		} else if (n_channels == 3) {
			format_in = gamma ? GL_SRGB : GL_RGB;
			format_out = GL_RGB;
		} else if (n_channels == 4) {
			format_in = gamma ? GL_SRGB_ALPHA : GL_RGBA;
			format_out = GL_RGBA;
		}

		glBindTexture(GL_TEXTURE_2D, texture_id);
		glTexImage2D(GL_TEXTURE_2D, 0, format_in, width, height, 0, format_out, GL_UNSIGNED_BYTE, data);
		glGenerateMipmap(GL_TEXTURE_2D);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	} else {
		std::cout << "ERROR::loadTextureFromFile(): Texture failed to load at path: " << path << '\n';
	}

	stbi_image_free(data);

	return texture_id;
}

GLMesh::GLMesh(std::vector<GLVertex> _vertices, std::vector<unsigned int> _indices, std::vector<std::shared_ptr<GLTexture>> _textures, const std::string& name, bool _use_indices) :
		vertices(_vertices), indices(_indices), textures(_textures), use_indices(_use_indices), VAO(0), VBO(0), EBO(0), Node(name, Mesh) {
	for (const auto& v : vertices) {
		bounds = bUnion(bounds, v.position);
	}

	setupMesh();
}

void GLMesh::draw(Shader& shader, unsigned int tex_offset) {
	///< will probably have to add shader specific code, for different shader structures
	// currently follows naming convention of material.texture_diffuse0 and so on or material.texture_specular0
	///< shadow maps always bind to GL_TEXTURE0, i = 0
	// TODO: make checks that tex_offset and drawing active textures stay within texture limit, until alternative solution found

	unsigned int diffuse_no = 1;
	unsigned int specular_no = 1;
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + (i + tex_offset));
		std::string number;
		std::string name = textures[i]->type;
		if (name == "texture_diffuse") {
			number = std::to_string(diffuse_no++);
		} else if (name == "texture_specular") {
			number = std::to_string(specular_no++);
		} else if (name == "texture_normal") {
			// TODO: might change depending on how many normal maps provided
			number = "";
		}

		shader.setInt(("material." + name + number).c_str(), i + tex_offset);
		glBindTexture(GL_TEXTURE_2D, textures[i]->id);
	}
	glActiveTexture(GL_TEXTURE0);

	shader.setMat4("model", getWorldTransform().m);

	glBindVertexArray(VAO);
	if (use_indices) {
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, nullptr);
	} else {
		glDrawArrays(GL_TRIANGLES, 0, vertices.size());
	}

	glBindVertexArray(0);
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + (i + tex_offset));
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	glCheckError("ERROR::GLMesh::draw: OpenGL error code");
}

void GLMesh::updateBounds() {
	bounds = Bounds3f();
	for (const auto& v : vertices) {
		bounds = bUnion(bounds, transform_data.transform(v.position));
	}
}

void GLMesh::setupMesh() {
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);

	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLVertex), &vertices[0], GL_STATIC_DRAW);

	if (use_indices) {
		glGenBuffers(1, &EBO);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
	}

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (void*)offsetof(GLVertex, position));
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (void*)offsetof(GLVertex, normal));
	glEnableVertexAttribArray(1);

	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (void*)offsetof(GLVertex, uv));
	glEnableVertexAttribArray(2);

	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (void*)offsetof(GLVertex, tangent));
	glEnableVertexAttribArray(3);

	glBindVertexArray(0);
	glCheckError("ERROR::GLMesh::setupMesh: OpenGL error code");
}

GLModel::GLModel(const std::string& name, bool gamma) :
		gamma_correction(gamma), Node(name, Model) {
}

void GLModel::loadModel(std::string path) {
	Assimp::Importer importer;
	unsigned int importer_flags = aiProcess_Triangulate | aiProcess_GenNormals | aiProcess_CalcTangentSpace | aiProcess_FlipUVs;
	const aiScene* scene = importer.ReadFile(path, importer_flags);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
		std::cout << "ERROR::GLModel: assimp error: " << importer.GetErrorString() << '\n';
		return;
	}

	loadModel(scene, path);
}

void GLModel::loadModel(const aiScene* scene, std::string path) {
	directory = path.substr(0, path.find_last_of('/'));

	processNode(scene->mRootNode, scene);

	for (const auto& m : meshes) {
		bounds = bUnion(bounds, m->bounds);
	}
}

void GLModel::draw(Shader& shader, unsigned int tex_offset) {
	shader.setBool("material.use_normal_map", use_normal_map);
	for (auto& mesh : meshes) {
		mesh->draw(shader, tex_offset);
	}
}

void GLModel::updateBounds() {
	bounds = Bounds3f();
	for (const auto& m : meshes) {
		bounds = bUnion(bounds, m->bounds);
	}
}

void GLModel::setTransform(Transform t) {
	transform_data.transform = t;
	updateBounds();
}

std::vector<std::shared_ptr<GLMesh>> GLModel::getMeshes() {
	return meshes;
}

void GLModel::processNode(aiNode* node, const aiScene* scene) {
	for (unsigned int i = 0; i < node->mNumMeshes; i++) {
		aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
		std::shared_ptr<GLMesh> newmesh = processMesh(mesh, scene);
		newmesh->parent = this;
		newmesh->name = mesh->mName.C_Str();
		meshes.push_back(newmesh);
	}
	for (unsigned int i = 0; i < node->mNumChildren; i++) {
		processNode(node->mChildren[i], scene);
	}
}

std::shared_ptr<GLMesh> GLModel::processMesh(aiMesh* mesh, const aiScene* scene) {
	std::vector<GLVertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<std::shared_ptr<GLTexture>> textures;

	for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
		GLVertex vertex;
		Point3f position(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);
		vertex.position = position;

		Normal3f normal(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
		vertex.normal = normal;

		Normal3f tangent(mesh->mTangents[i].x, mesh->mTangents[i].y, mesh->mTangents[i].z);
		vertex.tangent = tangent;

		if (mesh->mTextureCoords[0]) {
			Point2f uv(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);
			vertex.uv = uv;
		} else {
			vertex.uv = Point2f(0, 0);
		}

		vertices.push_back(vertex);
	}

	for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
		aiFace face = mesh->mFaces[i];
		for (unsigned int j = 0; j < face.mNumIndices; j++) {
			indices.push_back(face.mIndices[j]);
		}
	}

	///< handle loading textures
	if (mesh->mMaterialIndex >= 0) {
		aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
		std::vector<std::shared_ptr<GLTexture>> diffuse_maps = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
		textures.insert(textures.end(), diffuse_maps.begin(), diffuse_maps.end());
		std::vector<std::shared_ptr<GLTexture>> specular_maps = loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
		textures.insert(textures.end(), specular_maps.begin(), specular_maps.end());

		// TODO: make checks for file types; obj stores normals in HEIGHT
		std::vector<std::shared_ptr<GLTexture>> normal_maps = loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
		use_normal_map = normal_maps.size() > 0;
		textures.insert(textures.end(), normal_maps.begin(), normal_maps.end());
	}

	return std::make_shared<GLMesh>(vertices, indices, textures, "Mesh");
}

std::vector<std::shared_ptr<GLTexture>> GLModel::loadMaterialTextures(aiMaterial* mat, aiTextureType type, std::string type_name) {
	std::vector<std::shared_ptr<GLTexture>> textures;
	for (unsigned int i = 0; i < mat->GetTextureCount(type); i++) {
		aiString path;
		mat->GetTexture(type, i, &path);
		bool skip = false;
		for (unsigned int j = 0; j < loaded_textures.size(); j++) {
			if (std::strcmp(loaded_textures[j]->path.data(), path.C_Str()) == 0) {
				textures.push_back(loaded_textures[j]);
				skip = true;
				break;
			}
		}

		if (!skip) {
			// texture hasn't been loaded before
			bool gamma_correct = type == aiTextureType_DIFFUSE;
			std::shared_ptr<GLTexture> texture = std::make_shared<GLTexture>();
			texture->id = loadTextureFromFile(path.C_Str(), this->directory, gamma_correct);
			texture->type = type_name;
			texture->path = path.C_Str();
			textures.push_back(texture);
			loaded_textures.push_back(texture);
		}
	}

	return textures;
}

TriangleMesh::TriangleMesh(const Transform& otw, int n_tris, const int* v_idx, int n_verts,
		const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
		std::shared_ptr<Texture<float>> am) :
		n_triangles(n_tris), n_vertices(n_verts), vertex_indices(v_idx, v_idx + 3 * n_tris), alpha_mask(am) {
	p.reset(new Point3f[n_verts]);
	for (int i = 0; i < n_verts; i++) {
		p[i] = otw(_p[i]);
	}

	if (_uv) {
		uv.reset(new Point2f[n_verts]);
		memcpy(uv.get(), _uv, n_verts * sizeof(Point2f));
	}
	if (_n) {
		n.reset(new Normal3f[n_verts]);
		for (int i = 0; i < n_verts; i++) {
			n[i] = otw(_n[i]);
		}
	}
	if (_s) {
		s.reset(new Vector3f[n_verts]);
		for (int i = 0; i < n_verts; i++) {
			s[i] = otw(_s[i]);
		}
	}
}

std::vector<std::shared_ptr<Shape>> createTriangleMesh(const Transform* otw, const Transform* wto, bool _reverse_orientation,
		int n_tris, const int* v_idx, int n_verts,
		const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
		std::shared_ptr<Texture<float>> am) {
	auto _mesh = std::make_shared<TriangleMesh>(*otw, n_tris, v_idx, n_verts, _p, _s, _n, _uv, am);
	std::vector<std::shared_ptr<Shape>> tris;
	tris.reserve(n_tris);
	for (int i = 0; i < n_tris; i++) {
		tris.push_back(std::make_shared<Triangle>(otw, wto, _reverse_orientation, _mesh, i));
	}
	return tris;
}

} // namespace civet