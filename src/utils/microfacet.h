#ifndef CIVET_MICROFACET_H
#define CIVET_MICROFACET_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>

namespace civet {

class MicrofacetDistribution {
public:
	/**
	 * Distribution function D returns the differential area of microfacets oriented with given normal vector
	 * @param wh surface normal vector
	 * @return differential area of microfacets with given surface normal
	 */
	virtual float D(const Vector3f& wh) const = 0;

	/**
	 * Auxiliary function Lambda for shadowing-masking that measures invisible masked microfacet area per visible microfacet area
	 * @param w viewing direction
	 * @return invisible masked microfacet area per visible microfacet area
	 */
	virtual float lambda(const Vector3f& w) const = 0;

	/**
	 * Shadowing masking function
	 * @param w viewing direction
	 * @return fraction of the total microfacet area over differential area of surface dA that is visible in the given direction
	 */
	float G1(const Vector3f& w) const {
		return 1.0f / (1 + lambda(w));
	}

	/**
	 * Another shadowing masking function
	 * @param wo viewing direction
	 * @param wi another viewing direction
	 * @return fraction of the total microfacet area in a differential area that is visible in both given directions
	 */
	float G(const Vector3f& wo, const Vector3f& wi) const {
		return 1.0f / (1 + lambda(wo) + lambda(wi));
	}

protected:
	MicrofacetDistribution(bool sample) :
			sample_visible_area(sample) {}

	const bool sample_visible_area;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
	BeckmannDistribution(const float x, const float y, bool sample_visible_area = true) :
			MicrofacetDistribution(sample_visible_area), alphax(x), alphay(y) {}

	static float roughnessToAlpha(float roughness) {
		roughness = std::max(roughness, 1e-3f);
		float x = std::log(roughness);
		return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
	}

	float D(const Vector3f& wh) const override;

	float lambda(const civet::Vector3f &w) const override;

private:
	const float alphax, alphay;
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
	TrowbridgeReitzDistribution(const float x, const float y, bool sample_visible_area = true) :
			MicrofacetDistribution(sample_visible_area), alphax(x), alphay(y) {}

	static float roughnessToAlpha(float roughness) {
		roughness = std::max(roughness, 1e-3f);
		float x = std::log(roughness);
		return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
	}

	float D(const civet::Vector3f &wh) const override;

	float lambda(const civet::Vector3f &w) const override;

private:
	const float alphax, alphay;
};

} // namespace civet

#endif // CIVET_MICROFACET_H
