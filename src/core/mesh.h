#ifndef CIVET_MESH_H
#define CIVET_MESH_H

#include <assimp/scene.h>
#include <core/civet.h>
#include <core/shader.h>
#include <core/node.h>

namespace civet {

unsigned int loadTextureFromFile(const char* path, const std::string& directory, bool gamma = false);

struct GLVertex {
	Point3f position;
	Normal3f normal;
	Point2f uv;
	Normal3f tangent;
};

struct GLTexture {
	unsigned int id;
	std::string type;
	std::string path;
};

class GLMesh : public Node {
public:
	// TODO: implement this maybe when material class
	GLMesh(std::vector<GLVertex> _vertices, std::vector<unsigned int> _indices, std::vector<GLTexture> _textures, const std::string& name, bool _use_indices = true);

	void draw(Shader& shader, unsigned int tex_offset);

	std::vector<GLVertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<GLTexture> textures; ///< TODO: review when image textures implemented

private:
	void setupMesh();

	unsigned int VAO, VBO, EBO;
	bool use_indices;
};

class GLModel : public Node {
public:
	GLModel() :
			Node("GLModel", Model) {}

	GLModel(const char* path, const std::string& name, bool gamma = false);

	GLModel(const aiScene* scene, const std::string& name, const char* path, bool gamma = false) :
			gamma_correction(gamma), Node(name, Model) {
		loadModel(scene, path);
	}

	void draw(Shader& shader, unsigned int tex_offset);

	std::vector<GLMesh> getMeshes();

private:
	void loadModel(const aiScene* scene, std::string path);
	void processNode(aiNode* node, const aiScene* scene);
	GLMesh processMesh(aiMesh* mesh, const aiScene* scene);
	///< TODO: needs method to load textures here
	std::vector<GLTexture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, std::string type_name);

	std::vector<GLTexture> loaded_textures;
	std::vector<GLMesh> meshes;
	std::string directory;
	bool gamma_correction;
	bool use_normal_map;
};

class TriangleMesh {
public:
	/**
	 * Initializes a triangle mesh
	 * @param otw object-to-world transformation for the mesh
	 * @param n_tris total number of triangles in the mesh
	 * @param v_idx pointer to array of vertex indices, where i-th triangle has vertices at elements (3i, 3i+1, 3i+2)
	 * @param n_verts total number of vertices in the mesh
	 * @param _p array of n_verts vertex positions
	 * @param _s (optional) array of tangent vectors, one per vertex for computing shading tangents
	 * @param _n (optional) array of normals, one per vertex, interpolated across triangle faces as shading normals
	 * @param _uv (optional) array of parametric u,v-values, one per vertex
	 * @param am (optional) alpha mask texture
	 */
	TriangleMesh(const Transform& otw, int n_tris, const int* v_idx, int n_verts,
			const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
			std::shared_ptr<Texture<float>> am);

	const int n_triangles;
	const int n_vertices;
	std::vector<int> vertex_indices;
	std::unique_ptr<Point3f[]> p;
	std::unique_ptr<Normal3f[]> n;
	std::unique_ptr<Vector3f[]> s;
	std::unique_ptr<Point2f[]> uv;
	std::shared_ptr<Texture<float>> alpha_mask;
};

std::vector<std::shared_ptr<Shape>> createTriangleMesh(const Transform* otw, const Transform* wto, bool _reverse_orientation,
		int n_tris, const int* v_idx, int n_verts,
		const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
		std::shared_ptr<Texture<float>> am);

} // namespace civet

#endif // CIVET_MESH_H
