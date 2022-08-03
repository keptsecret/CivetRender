#include <core/mesh.h>

#include <assimp/postprocess.h>
#include <shapes/triangle.h>
#include <stb/stb_image.h>
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
		GLenum format;
		if (n_channels == 1) {
			format = GL_RED;
		} else if (n_channels == 3) {
			format = GL_RGB;
		} else if (n_channels == 4) {
			format = GL_RGBA;
		}

		glBindTexture(GL_TEXTURE_2D, texture_id);
		glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
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

GLMesh::GLMesh(std::vector<GLVertex> _vertices, std::vector<unsigned int> _indices, std::vector<GLTexture> _textures, bool _use_indices) :
		vertices(_vertices), indices(_indices), textures(_textures), use_indices(_use_indices), VAO(0), VBO(0), EBO(0) {
	setupMesh();
}

void GLMesh::draw(Shader& shader) {
	///< will probably have to add shader specific code, for different shader structures
	// currently follows naming convention of material.texture_diffuse0 and so on or material.texture_specular0

	unsigned int diffuse_no = 1;
	unsigned int specular_no = 1;
	for (unsigned int i = 0; i < textures.size(); i++) {
		glActiveTexture(GL_TEXTURE0 + i);
		std::string number;
		std::string name = textures[i].type;
		if (name == "texture_diffuse") {
			number = std::to_string(diffuse_no++);
		} else if (name == "texture_specular") {
			number = std::to_string(specular_no++);
		}

		shader.setInt(("material." + name + number).c_str(), i);
		glBindTexture(GL_TEXTURE_2D, textures[i].id);
	}
	glActiveTexture(GL_TEXTURE0);

	glBindVertexArray(VAO);
	if (use_indices) {
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, nullptr);
	} else {
		glDrawArrays(GL_TRIANGLES, 0, vertices.size());
	}
	glBindVertexArray(0);
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

	glBindVertexArray(0);
}

void GLModel::draw(Shader& shader) {
	for (auto& mesh : meshes) {
		mesh.draw(shader);
	}
}

void GLModel::loadModel(std::string path) {
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
		std::cout << "ERROR::GLModel: Assimp error: " << importer.GetErrorString() << '\n';
		return;
	}
	directory = path.substr(0, path.find_last_of('/'));

	processNode(scene->mRootNode, scene);
}

void GLModel::processNode(aiNode* node, const aiScene* scene) {
	for (unsigned int i = 0; i < node->mNumMeshes; i++) {
		aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
		meshes.push_back(processMesh(mesh, scene));
	}
	for (unsigned int i = 0; i < node->mNumChildren; i++) {
		processNode(node->mChildren[i], scene);
	}
}

GLMesh GLModel::processMesh(aiMesh* mesh, const aiScene* scene) {
	std::vector<GLVertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<GLTexture> textures;

	for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
		GLVertex vertex;
		Point3f position(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);
		vertex.position = position;

		Normal3f normal(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
		vertex.normal = normal;

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
		std::vector<GLTexture> diffuse_maps = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
		textures.insert(textures.end(), diffuse_maps.begin(), diffuse_maps.end());
		std::vector<GLTexture> specular_maps = loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
		textures.insert(textures.end(), specular_maps.begin(), specular_maps.end());
	}

	return GLMesh(vertices, indices, textures);
}

std::vector<GLTexture> GLModel::loadMaterialTextures(aiMaterial* mat, aiTextureType type, std::string type_name) {
	std::vector<GLTexture> textures;
	for (unsigned int i = 0; i < mat->GetTextureCount(type); i++) {
		aiString path;
		mat->GetTexture(type, i, &path);
		bool skip = false;
		for (unsigned int j = 0; j < loaded_textures.size(); j++) {
			if (std::strcmp(loaded_textures[j].path.data(), path.C_Str()) == 0) {
				textures.push_back(loaded_textures[j]);
				skip = true;
				break;
			}
		}

		if (!skip) {
			// texture hasn't been loaded before
			GLTexture texture;
			texture.id = loadTextureFromFile(path.C_Str(), this->directory);
			texture.type = type_name;
			texture.path = path.C_Str();
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