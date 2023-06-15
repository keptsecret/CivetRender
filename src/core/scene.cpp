#include <core/scene.h>

#include <lights/distant.h>
#include <lights/point.h>
#include <materials/disney.h>
#include <shapes/triangle.h>
#include <textures/constant.h>
#include <textures/imagemap.h>
#include <utils/bvh.h>

#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>

namespace civet {

Scene::Scene(const char* path) {
	SkyboxParameters params{ Vector3f(0.f, 0.5f, 0.5f), 128 };
	skybox = std::make_shared<Skybox>();
	skybox->init(params);
	nodes.push_back(skybox);

	probe_grid = std::make_shared<IlluminanceField>();

	loadScene(path);
}

void Scene::loadScene(const char* path) {
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_GenNormals | aiProcess_CalcTangentSpace | aiProcess_FlipUVs);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
		std::cout << "ERROR::GLModel: assimp error: " << importer.GetErrorString() << '\n';
		return;
	}

	auto model = std::make_shared<GLModel>("Model");
	model->loadModel(scene, path);
	models.push_back(model);
	nodes.push_back(model);

	for (int i = 0; i < scene->mNumLights; i++) {
		aiLight* light = scene->mLights[i];
		std::string light_name = ("Light_" + std::to_string(i));

		switch (light->mType) {
			case aiLightSource_DIRECTIONAL: {
				auto new_light = std::make_shared<GLDirectionalLight>(light_name, Vector3f(light->mDirection.x, light->mDirection.y, light->mDirection.z));
				new_light->color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light->init();
				dir_lights.push_back(new_light);
				nodes.push_back(new_light);
				break;
			}
			case aiLightSource_POINT: {
				auto new_light = std::make_shared<GLPointLight>(light_name, Point3f(light->mPosition.x, light->mPosition.y, light->mPosition.z));
				new_light->color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light->attenuation.constant = light->mAttenuationConstant;
				new_light->attenuation.linear = light->mAttenuationLinear;
				new_light->attenuation.quadratic = light->mAttenuationQuadratic;
				new_light->init();
				point_lights.push_back(new_light);
				nodes.push_back(new_light);
				break;
			}
			default:
				break;
		}
	}
}

void Scene::addNode(std::shared_ptr<Node> node, NodeType type) {
	switch (type) {
		case NodeType::Model:
			models.push_back(std::static_pointer_cast<GLModel>(node));
			break;
		case NodeType::DirectionalLight:
			dir_lights.push_back(std::static_pointer_cast<GLDirectionalLight>(node));
			break;
		case NodeType::PointLight:
			point_lights.push_back(std::static_pointer_cast<GLPointLight>(node));
			break;
	}

	nodes.push_back(node);
}

void Scene::buildScene() {
	printf("Building scene for ray tracing...\n");
	std::vector<std::shared_ptr<Primitive>> prims;
	for (const auto model : models) {
		auto meshes = model->getMeshes();
		for (const auto mesh : meshes) {
			// Convert mesh data into shapes
			std::vector<int> indices(mesh->indices.begin(), mesh->indices.end());

			std::vector<Point3f> vertices;
			std::vector<Normal3f> normals;
			std::vector<Vector3f> tangents;
			std::vector<Point2f> uvs;
			vertices.reserve(mesh->vertices.size());
			normals.reserve(mesh->vertices.size());
			tangents.reserve(mesh->vertices.size());
			uvs.reserve(mesh->vertices.size());

			for (const auto tri : mesh->vertices) {
				vertices.push_back(tri.position);
				normals.push_back(tri.normal);
				tangents.push_back(tri.tangent);
				uvs.push_back(tri.uv);
			}

			mesh->setStaticWorldTransforms();
			int n_tris = indices.size() / 3;
			std::vector<std::shared_ptr<Shape>> tris = createTriangleMesh(&mesh->transform_data.static_transform, &mesh->transform_data.static_inv_transform,
					false, n_tris, indices.data(), vertices.size(), vertices.data(), tangents.data(), normals.data(), uvs.data(), nullptr);

			// Convert material
			///< we only have Disney BRDF so far in glsl
			///< doesn't map over well and is missing lots of textures, including normal maps
			auto material = mesh->material;
			DisneyBSDFParams params;
			for (const auto tex : material->textures) {
				std::unique_ptr<UVMapping2D> map;
				map.reset(new UVMapping2D(1., 1., 0., 0.));
				std::string filepath = model->getDirectory() + '/' + tex->path;

				if (tex->type == "texture_albedo") {
					auto texture = new ImageTexture<RGBSpectrum, Spectrum>(std::move(map), filepath, false, 8.f, ImageWrap::Repeat, 1.f, true);
					params.color = std::shared_ptr<Texture<Spectrum>>(texture);
				} else if (tex->type == "texture_metallic") {
					auto texture = new ImageTexture<RGBSpectrum, float>(std::move(map), filepath, false, 8.f, ImageWrap::Repeat, 1.f, true);
					params.metallic = std::shared_ptr<Texture<float>>(texture);
				} else if (tex->type == "texture_roughness") {
					auto texture = new ImageTexture<RGBSpectrum, float>(std::move(map), filepath, false, 8.f, ImageWrap::Repeat, 1.f, true);
					params.roughness = std::shared_ptr<Texture<float>>(texture);
				} else if (tex->type == "texture_bump") {
					auto texture = new ImageTexture<RGBSpectrum, float>(std::move(map), filepath, false, 8.f, ImageWrap::Repeat, 1.f, true);
					params.bump = std::shared_ptr<Texture<float>>(texture);
				}
			}

			float rgb[3] = {material->albedo.x, material->albedo.y, material->albedo.z};
			std::shared_ptr<Texture<Spectrum>> color = params.color ? params.color : std::make_shared<ConstantTexture<Spectrum>>(Spectrum::fromRGB(rgb));
			std::shared_ptr<Texture<float>> metallic = params.metallic ? params.metallic : std::make_shared<ConstantTexture<float>>(material->metallic);
			std::shared_ptr<Texture<float>> eta = params.eta ? params.eta : std::make_shared<ConstantTexture<float>>(1.5f);
			std::shared_ptr<Texture<float>> roughness = params.roughness ? params.roughness : std::make_shared<ConstantTexture<float>>(material->roughness);
			std::shared_ptr<Texture<float>> spec_tint = params.spec_tint ? params.spec_tint : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<float>> anisotropic = params.anis ? params.anis : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<float>> sheen = params.sheen ? params.sheen : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<float>> sheen_tint = params.sheen_tint ? params.sheen_tint : std::make_shared<ConstantTexture<float>>(0.5f);
			std::shared_ptr<Texture<float>> clearcoat = params.cc ? params.cc : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<float>> clearcoat_gloss = params.cc_gloss ? params.cc_gloss : std::make_shared<ConstantTexture<float>>(1.f);
			std::shared_ptr<Texture<float>> spec_trans = params.spec_trans ? params.spec_trans : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<Spectrum>> scatter_dist = params.scatter_dist ? params.scatter_dist : std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.f));
			bool thin = params.t;
			std::shared_ptr<Texture<float>> flatness = params.flatness ? params.flatness : std::make_shared<ConstantTexture<float>>(0.f);
			std::shared_ptr<Texture<float>> diff_trans = params.diff_trans ? params.diff_trans : std::make_shared<ConstantTexture<float>>(1.f);
			std::shared_ptr<Texture<float>> bump = params.bump;
			auto mapped_mat = std::make_shared<DisneyMaterial>(color, metallic, eta, roughness, spec_tint,
					anisotropic, sheen, sheen_tint, clearcoat,
					clearcoat_gloss, spec_trans, scatter_dist, thin,
					flatness, diff_trans, bump, material->is_glossy_rough);

			for (auto& t : tris) {
				prims.push_back(std::make_shared<GeometricPrimitive>(t, mapped_mat, nullptr, nullptr));
			}
		}
	}

	printf("Building BVH...\n");
	aggregate = std::make_shared<BVH>(prims, 4, BVH::SplitMethod::HLBVH);
	world_bound = aggregate->worldBound();

	// Convert lights
	printf("Building lights...\n");
	for (const auto& dir : dir_lights) {
		if (!dir->active) {
			continue;
		}

		RGBSpectrum rgb;
		rgb[0] = dir->color.x * dir->power;
		rgb[1] = dir->color.y * dir->power;
		rgb[2] = dir->color.z * dir->power;
		lights.push_back(std::make_shared<DistantLight>(dir->transform_data.transform, rgb, -dir->direction));
	}
	for (const auto& point : point_lights) {
		if (!point->active) {
			continue;
		}

		RGBSpectrum rgb;
		rgb[0] = point->color.x * point->power;
		rgb[1] = point->color.y * point->power;
		rgb[2] = point->color.z * point->power;
		lights.push_back(std::make_shared<PointLight>(point->transform_data.transform, MediumInterface(), point->position, rgb));
	}

	for (const auto& light : lights) {
		light->preprocess(*this);
	}

	is_built = true;
}

bool Scene::intersect(const Ray& ray, SurfaceInteraction* isect) const {
	return aggregate->intersect(ray, isect);
}

bool Scene::intersectP(const Ray& ray) const {
	return aggregate->intersectP(ray);
}

bool Scene::intersectTr(Ray ray, Sampler& sampler, SurfaceInteraction* isect, Spectrum* transmittance) const {
	*transmittance = Spectrum(1.f);
	while (true) {
		bool hit = intersect(ray, isect);
		if (ray.medium) {
			*transmittance = ray.medium->Tr(ray, sampler);
		}
		if (!hit) {
			return false;
		}
		if (isect->primitive->getMaterial() != nullptr) {
			return true;
		}
		ray = isect->spawnRay(ray.d);
	}
}

} // namespace civet