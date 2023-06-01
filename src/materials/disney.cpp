#include <materials/disney.h>

#include <core/interaction.h>
#include <core/texture.h>
#include <utils/reflection.h>

/**
 * Based on descriptions of Disney BSDF from:
 * https://media.disneyanimation.com/uploads/production/publication_asset/48/asset/s2012_pbs_disney_brdf_notes_v3.pdf
 * https://blog.selfshadow.com/publications/s2015-shading-course/burley/s2015_pbs_disney_bsdf_notes.pdf
 * and lots of other places
 *
 * TODO: implement disney bssrdf
 */

namespace civet {

inline float sqr(float x) {
	return x * x;
}

// For a dielectric, R(0) = (eta - 1)^2 / (eta + 1)^2, assuming we're
// coming from air..
inline float schlickR0FromEta(float eta) {
	return sqr(eta - 1) / sqr(eta + 1);
}

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

// Implements Generalized Trowbridge-Reitz, using fixed gamma in its normalized form
inline float GTR1(float cos_theta, float alpha) {
	float alpha2 = alpha * alpha;
	return (alpha2 - 1) / (Pi * std::log(alpha2) * (1 + (alpha2 - 1) * cos_theta * cos_theta));
}

// Smith masking-shadowing term for GGX
inline float smithGGXG1(float cos_theta, float alpha) {
	float cos_theta2 = cos_theta * cos_theta;
	float alpha2 = alpha * alpha;
	return 1 / (cos_theta + std::sqrt(alpha2 + cos_theta2 - alpha2 * cos_theta2));
}

class DisneyClearcoat : public BxDF {
public:
	DisneyClearcoat(float weight, float gloss) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
			weight(weight),
			gloss(gloss) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type) const override;
	float pdf(const Vector3f& wo, const Vector3f& wi) const;

private:
	float weight, gloss;
};

Spectrum DisneyClearcoat::f(const Vector3f& wo, const Vector3f& wi) const {
	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.f);
	}
	wh = normalize(wh);

	float Dr = GTR1(absCosTheta(wh), gloss);
	// Fresnel term uses Schlick approximation with hardcoded IOR of 1.5 -> F0 = 0.04
	float Fr = frSchlick(0.04f, dot(wo, wh));
	// Fixed roughness of 0.25 for masking-shadowing
	float Gr = smithGGXG1(absCosTheta(wo), 0.25f) * smithGGXG1(absCosTheta(wi), 0.25f);
	return weight * Gr * Fr * Dr * 0.25f;
}

Spectrum DisneyClearcoat::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	if (wo.z == 0) {
		return Spectrum(0.f);
	}

	float alpha2 = gloss * gloss;
	float cos_theta = std::sqrt(std::max(float(0), (1 - std::pow(alpha2, 1 - sample[0])) / (1 - alpha2)));
	float sin_theta = std::sqrt(std::max((float)0, 1 - cos_theta * cos_theta));
	float phi = 2 * Pi * sample[1];
	Vector3f wh = sphericalDirection(sin_theta, cos_theta, phi);
	if (!sameHemisphere(wo, wh)) {
		wh = -wh;
	}

	*wi = reflect(wo, wh);
	if (!sameHemisphere(wo, *wi)) {
		return Spectrum(0.f);
	}

	*pdf_val = pdf(wo, *wi);
	return f(wo, *wi);
}

float DisneyClearcoat::pdf(const Vector3f& wo, const Vector3f& wi) const {
	if (!sameHemisphere(wo, wi)) {
		return 0.f;
	}

	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return 0.f;
	}
	wh = normalize(wh);

	// The sampling routine samples wh exactly from the GTR1 distribution.
	// Thus, the final value of the PDF is just the value of the
	// distribution for wh converted to a measure with respect to the
	// surface normal.
	float Dr = GTR1(absCosTheta(wh), gloss);
	return Dr * absCosTheta(wh) / (4 * dot(wo, wh));
}

void DisneyMaterial::computeScatteringFunctions(SurfaceInteraction* si, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	// Perform bump mapping with _bumpMap_, if present
	if (bump_map) {
		bump(bump_map, si);
	}

	// Evaluate textures for _DisneyMaterial_ material and allocate BRDF
	si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

	// Diffuse
	Spectrum c = color->evaluate(*si).clamp();
	float metallic_weight = metallic->evaluate(*si, 2);
	float e = eta->evaluate(*si);
	float strans = specular_transmission->evaluate(*si);
	float diffuse_weight = (1 - metallic_weight) * (1 - strans);
	float dt = diffuse_transmission->evaluate(*si) / 2; // 0: all diffuse is reflected -> 1, transmitted
	float rough = roughness->evaluate(*si, 1);
	float lum = c.y();
	// normalize lum. to isolate hue+sat
	Spectrum Ctint = lum > 0 ? (c / lum) : Spectrum(1.);

	float sheen_weight = sheen->evaluate(*si);
	Spectrum Csheen;
	if (sheen_weight > 0) {
		float stint = sheen_tint->evaluate(*si);
		Csheen = lerp(stint, Spectrum(1.), Ctint);
	}

	if (diffuse_weight > 0) {
		if (thin) {
			float flat = flatness->evaluate(*si);
			// Blend between DisneyDiffuse and fake subsurface based on
			// flatness.  Additionally, weight using diffTrans.
			si->bsdf->add(ARENA_ALLOC(arena, DisneyDiffuse)(diffuse_weight * (1 - flat) * (1 - dt) * c));
			si->bsdf->add(ARENA_ALLOC(arena, DisneyFakeSS)(diffuse_weight * flat * (1 - dt) * c, rough));
		} else {
			Spectrum sd = scatter_distance->evaluate(*si);
			if (sd.isBlack()) { // No subsurface scattering; use regular (Fresnel modified)
				// diffuse.
				si->bsdf->add(ARENA_ALLOC(arena, DisneyDiffuse)(diffuse_weight * c));
			} else {
				// Use a BSSRDF instead.
				// TODO: implement BSSRDFs
				// si->bsdf->add(ARENA_ALLOC(arena, SpecularTransmission)(1.f, 1.f, e, mode));
				// si->bssrdf = ARENA_ALLOC(arena, DisneyBSSRDF)(c * diffuse_weight, sd, *si, e, this, mode);
			}
		}

		// Retro-reflection.
		si->bsdf->add(
				ARENA_ALLOC(arena, DisneyRetro)(diffuse_weight * c, rough));

		// Sheen (if enabled)
		if (sheen_weight > 0) {
			si->bsdf->add(ARENA_ALLOC(arena, DisneySheen)(diffuse_weight * sheen_weight * Csheen));
		}
	}

	// Create the microfacet distribution for metallic and/or specular
	// transmission.
	float aspect = std::sqrt(1 - anisotropic->evaluate(*si) * .9);
	float ax = std::max(float(.001), sqr(rough) / aspect);
	float ay = std::max(float(.001), sqr(rough) * aspect);
	MicrofacetDistribution* distrib = ARENA_ALLOC(arena, DisneyMicrofacetDistribution)(ax, ay);

	// Specular is Trowbridge-Reitz with a modified Fresnel function.
	float specTint = specular_tint->evaluate(*si);
	Spectrum Cspec0 = lerp(metallic_weight, schlickR0FromEta(e) * lerp(specTint, Spectrum(1.), Ctint), c);
	Fresnel* fresnel = ARENA_ALLOC(arena, DisneyFresnel)(Cspec0, metallic_weight, e);
	si->bsdf->add(ARENA_ALLOC(arena, MicrofacetReflection)(Spectrum(1.), distrib, fresnel));

	// Clearcoat
	float cc = clearcoat->evaluate(*si);
	if (cc > 0) {
		si->bsdf->add(ARENA_ALLOC(arena, DisneyClearcoat)(cc, lerp(clearcoat_gloss->evaluate(*si), .1, .001)));
	}

	// BTDF
	if (strans > 0) {
		// Walter et al's model, with the provided transmissive term scaled
		// by sqrt(color), so that after two refractions, we're back to the
		// provided color.
		Spectrum T = strans * sqrt(c);
		if (thin) {
			// Scale roughness based on IOR (Burley 2015, Figure 15).
			float rscaled = (0.65f * e - 0.35f) * rough;
			float ax = std::max(float(.001), sqr(rscaled) / aspect);
			float ay = std::max(float(.001), sqr(rscaled) * aspect);
			MicrofacetDistribution* scaledDistrib =
					ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(ax, ay);
			si->bsdf->add(ARENA_ALLOC(arena, MicrofacetTransmission)(T, scaledDistrib, 1., e, mode));
		} else {
			si->bsdf->add(ARENA_ALLOC(arena, MicrofacetTransmission)(T, distrib, 1., e, mode));
		}
	}
	if (thin) {
		// Lambertian, weighted by (1 - diffTrans)
		si->bsdf->add(ARENA_ALLOC(arena, LambertianTransmission)(dt * c));
	}
}

} // namespace civet