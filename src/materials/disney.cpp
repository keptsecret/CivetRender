#include <materials/disney.h>

#include <utils/reflection.h>

/**
 * Based on descriptions of Disney BSDF from:
 * https://media.disneyanimation.com/uploads/production/publication_asset/48/asset/s2012_pbs_disney_brdf_notes_v3.pdf
 * https://blog.selfshadow.com/publications/s2015-shading-course/burley/s2015_pbs_disney_bsdf_notes.pdf
 * and lots of other places
 *
 * TODO: implement clearcoat and disney bssrdf, computescatteringfunction
 */

namespace civet {

class DisneyDiffuse : public BxDF {
public:
	DisneyDiffuse(const Spectrum& R) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f& wo, int n_samples, const Point2f* samples) const override { return R; }
	Spectrum rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const override { return R; }

private:
	Spectrum R;
};

Spectrum DisneyDiffuse::f(const Vector3f& wo, const Vector3f& wi) const {
	float Fo = schlickWeight(absCosTheta(wo));
	float Fi = schlickWeight(absCosTheta(wi));

	return R * InvPi * (1 - Fo / 2) * (1 - Fi / 2);
}

class DisneyFakeSS : public BxDF {
public:
	DisneyFakeSS(const Spectrum& R, float roughness) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
			R(R),
			roughness(roughness) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f& wo, int n_samples, const Point2f* samples) const override { return R; }
	Spectrum rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const override { return R; }

private:
	Spectrum R;
	float roughness;
};

Spectrum DisneyFakeSS::f(const Vector3f& wo, const Vector3f& wi) const {
	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.f);
	}
	wh = normalize(wh);
	float cos_theta_d = dot(wi, wh);

	float Fss90 = cos_theta_d * cos_theta_d * roughness;
	float Fo = schlickWeight(absCosTheta(wo));
	float Fi = schlickWeight(absCosTheta(wi));
	float Fss = lerp(Fo, 1.0f, Fss90) * lerp(Fi, 1.0f, Fss90);
	float ss = 1.25f * (Fss * (1 / (absCosTheta(wo) + absCosTheta(wi)) - 0.5f) + 0.5f);

	return R * InvPi * ss;
}

class DisneyRetro : public BxDF {
public:
	DisneyRetro(const Spectrum& R, float roughness) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
			R(R),
			roughness(roughness) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f& wo, int n_samples, const Point2f* samples) const override { return R; }
	Spectrum rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const override { return R; }

private:
	Spectrum R;
	float roughness;
};

Spectrum DisneyRetro::f(const Vector3f& wo, const Vector3f& wi) const {
	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.f);
	}
	wh = normalize(wh);
	float cos_theta_d = dot(wi, wh);

	float Fo = schlickWeight(absCosTheta(wo));
	float Fi = schlickWeight(absCosTheta(wi));
	float Rr = 2 * roughness * cos_theta_d * cos_theta_d;

	return R * InvPi * Rr * (Fo + Fi + Fo * Fi * (Rr - 1));
}

class DisneySheen : public BxDF {
public:
	DisneySheen(const Spectrum& R) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f& wo, int n_samples, const Point2f* samples) const override { return R; }
	Spectrum rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const override { return R; }

private:
	Spectrum R;
};

Spectrum DisneySheen::f(const Vector3f& wo, const Vector3f& wi) const {
	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.f);
	}
	wh = normalize(wh);
	float cos_theta_d = dot(wi, wh);

	return R * schlickWeight(cos_theta_d);
}

class DisneyFresnel : public Fresnel {
public:
	DisneyFresnel(const Spectrum& R0, float metallic, float eta) :
			R0(R0), metallic(metallic), eta(eta) {}

	Spectrum evaluate(float cos_theta_I) const override {
		return lerp(metallic, Spectrum(frDielectric(cos_theta_I, 1, eta)), frSchlick(R0, cos_theta_I));
	}

private:
	const Spectrum R0;
	const float metallic, eta;
};

class DisneyMicrofacetDistribution : public TrowbridgeReitzDistribution {
public:
	DisneyMicrofacetDistribution(float a_x, float a_y) :
			TrowbridgeReitzDistribution(a_x, a_y) {}

	float G(const Vector3f& wo, const Vector3f& wi) const {
		return G1(wo) * G1(wi);
	}
};

void DisneyMaterial::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
}

} // namespace civet