#include <core/mesh.h>

#include <assimp/postprocess.h>
#include <shapes/triangle.h>
#include <assimp/Importer.hpp>

namespace civet {

GLMesh::GLMesh(std::vector<GLVertex> _vertices, std::vector<unsigned int> _indices, bool _use_indices) :
		vertices(_vertices), indices(_indices), use_indices(_use_indices) {
	setupMesh();
}

void GLMesh::draw(Shader& shader) {
	///< will probably have to add shader specific code, for different shader structures

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
	///< std::vector<Texture> textures;

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

	return GLMesh(vertices, indices);
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