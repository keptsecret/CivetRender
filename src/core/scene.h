#ifndef CIVET_SCENE_H
#define CIVET_SCENE_H

#include <core/civet.h>
#include <core/skybox.h>
#include <map>

namespace civet {

class Editor;

class Scene {
public:
	Scene() {}

	Scene(const char* path);

	void loadScene(const char* path);
	void addNode(std::shared_ptr<Node> node, NodeType type);

	void clearSelectedNode() {
		selected_node = nullptr;
		selected_self = false;
	}

	// ray tracing functions
	const Bounds3f& worldBound() const { return world_bound; }
	void buildScene();
	bool isBuildForRT() const { return is_built; }

	bool intersect(const Ray& ray, SurfaceInteraction* isect) const;
	bool intersectP(const Ray& ray) const;
	bool intersectTr(Ray ray, Sampler& sampler, SurfaceInteraction* isect, Spectrum* transmittance) const;

public:
	std::vector<std::shared_ptr<GLModel>> models;
	std::vector<std::shared_ptr<GLDirectionalLight>> dir_lights;
	std::vector<std::shared_ptr<GLPointLight>> point_lights;

	std::vector<std::shared_ptr<Node>> nodes;
	std::shared_ptr<Skybox> skybox;

	std::vector<std::shared_ptr<Light>> lights;

private:
	friend Editor;

	bool selected_self = false;
	std::shared_ptr<Node> selected_node = nullptr;

	// ray tracing variables
	std::shared_ptr<Primitive> aggregate;
	// std::map<std::string, std::shared_ptr<Texture<float>>> float_textures;
	// std::map<std::string, std::shared_ptr<Texture<Spectrum>>> spectrum_textures;
	Bounds3f world_bound;
	bool is_built = false;
};

} // namespace civet

#endif // CIVET_SCENE_H
