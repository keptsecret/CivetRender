#ifndef CIVET_REFLECTION_H
#define CIVET_REFLECTION_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/shape.h>
#include <core/spectrum.h>
#include <utils/microfacet.h>

namespace civet {

inline float cosTheta(const Vector3f& w) {
	return w.z;
}
inline float cos2Theta(const Vector3f& w) {
	return w.z * w.z;
}
inline float absCosTheta(const Vector3f& w) {
	return std::abs(w.z);
}
inline float sin2Theta(const Vector3f& w) {
	return std::max(0.0f, 1.0f - cos2Theta(w));
}
inline float sinTheta(const Vector3f& w) {
	return std::sqrt(sin2Theta(w));
}
inline float tanTheta(const Vector3f& w) {
	return sinTheta(w) / cosTheta(w);
}
inline float tan2Theta(const Vector3f& w) {
	return sin2Theta(w) / cos2Theta(w);
}

inline float cosPhi(const Vector3f& w) {
	float sin_theta = sinTheta(w);
	return sin_theta == 0 ? 1.0f : clamp(w.x / sin_theta, -1, 1);
}
inline float sinPhi(const Vector3f& w) {
	float sin_theta = sinTheta(w);
	return sin_theta == 0 ? 1.0f : clamp(w.y / sin_theta, -1, 1);
}
inline float cos2Phi(const Vector3f& w) {
	return cosPhi(w) * cosPhi(w);
}
inline float sin2Phi(const Vector3f& w) {
	return sinPhi(w) * sinPhi(w);
}
inline float cosDPhi(const Vector3f& wa, const Vector3f& wb) {
	return clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt((wa.x * wa.x + wa.y * wa.y) * (wb.x * wb.x + wb.y * wb.y)), -1, 1);
}

inline Vector3f reflect(const Vector3f& wo, const Vector3f& n) {
	return -wo + 2 * dot(wo, n) * n;
}

inline bool refract(const Vector3f& wi, const Normal3f& n, float eta, Vector3f* wt) {
	// compute cos_theta_t using snell's law
	float cos_theta_I = dot(n, wi);
	float sin2_theta_I = std::max(0.0f, 1.0f - cos_theta_I * cos_theta_I);
	float sin2_theta_T = eta * eta * sin2_theta_I;
	if (sin2_theta_T >= 0) {
		return false;
	}
	float cos_theta_T = std::sqrt(1.0f - sin2_theta_I);

	*wt = eta * -wi + (eta * cos_theta_I * cos_theta_T) * Vector3f(n);
	return true;
}

inline bool sameHemisphere(const Vector3f& w, const Vector3f& wp) {
	return w.z * wp.z > 0;
}

inline bool sameHemisphere(const Vector3f& w, const Normal3f& wp) {
	return w.z * wp.z > 0;
}

inline float schlickWeight(float cos_theta) {
	float m = clamp(1- cos_theta, 0, 1);
	return m * m * m * m * m;
}

inline float frSchlick(float R0, float cos_theta) {
	return R0 + (1.f - R0) * schlickWeight(cos_theta);
}

inline Spectrum frSchlick(const Spectrum& R0, float cos_theta) {
	return lerp(schlickWeight(cos_theta), R0, Spectrum(1.f));	///< possibly also wrong
}

float frDielectric(float cos_theta_I, float eta_I, float eta_T);
Spectrum frConductor(float cos_theta_I, const Spectrum& eta_I, const Spectrum& eta_T, const Spectrum& k);

enum BxDFType {
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR |
			BSDF_REFLECTION | BSDF_TRANSMISSION
};

struct FourierBSDFTable {
	float eta; // relative index of refraction between two medium
	int m_max; // max order m for any pair of mu
	int n_channels; // no. of spectral channels
	int n_mu;
	float* mu; // array of zenith angles
	int* m; // order of fourier representation
	int* a_offset; // table of offsets for a
	float* a; // ak values for a pair of directions

	~FourierBSDFTable() {
		delete[] mu;
		delete[] m;
		delete[] a_offset;
		delete[] a;
	}

	static bool read(const std::string& filename, FourierBSDFTable* table); ///< TODO: implement in Fourier material

	const float* getAk(int offset_I, int offset_O, int* mptr) const {
		*mptr = m[offset_O * n_mu + offset_I];
		return a + a_offset[offset_O * n_mu + offset_I];
	}

	bool getWeightsAndOffset(float cos_theta, int* offset, float weights[4]) const;
};

/**
 * BSDF class contains a group of BRDFs and BTDFs\n
 * Handles shading normals (per-vertex or normal maps) and world-local transformations
 */
class BSDF {
public:
	BSDF(const SurfaceInteraction& si, float eta = 1) :
			eta(eta), ns(si.shading.n), ng(si.n), ss(normalize(si.shading.dpdu)), ts(cross(ns, ss)) {}

	void add(BxDF* b) {
		if (num_BxDFs < Max_BxDFs) {
			bxdfs[num_BxDFs++] = b;
		} else {
			printf("WARNING::BSDF::add: reached BxDF limit");
		}
	}

	int numComponents(BxDFType flags = BSDF_ALL) const;

	Vector3f worldToLocal(const Vector3f& v) const {
		return Vector3f(dot(v, ss), dot(v, ts), dot(v, ns));
	}

	Vector3f localToWorld(const Vector3f& v) const {
		// uses the transpose of M from worldToLocal
		return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
				ss.y * v.x + ts.y * v.y + ns.y * v.z,
				ss.z * v.x + ts.z * v.y + ns.z * v.z);
	}

	Spectrum f(const Vector3f& wo_world, const Vector3f& wi_world, BxDFType flags) const;

	Spectrum rho(const Vector3f& wo, int num_samples, const Point2f* samples, BxDFType flags = BSDF_ALL) const;
	Spectrum rho(int num_samples, const Point2f* samples1, const Point2f* samples2, BxDFType flags = BSDF_ALL) const;

	const float eta;

private:
	~BSDF() {}
	friend class MixMaterial;

	const Normal3f ns, ng;
	const Vector3f ss, ts;

	int num_BxDFs = 0;
	static CIVET_CONSTEXPR int Max_BxDFs = 8;
	BxDF* bxdfs[Max_BxDFs];
};

/**
 * Base BxDF interface for all reflection and transmission models (BRDFs and BTDFs)\n
 * All interactions happen in local space
 */
class BxDF {
public:
	BxDF(BxDFType _type) :
			type(_type) {}

	virtual ~BxDF() {}

	/**
	 * Returns value of the distribution function for given pair of directions
	 * @param wo outgoing direction
	 * @param wi incoming direction
	 * @return value of the distribution function as a Spectrum
	 */
	virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const = 0;

	/**
	 * Computes direction of incident light given outgoing direction
	 * and returns value of the distribution function for given pair of directions
	 * @param wo outgoing direction
	 * @param wi new incoming direction
	 * @param sample
	 * @param pdf
	 * @param sampled_type
	 * @return value of the distribution function as a Spectrum (incident direction stored in wi)
	 */
	virtual Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type = nullptr) const;

	/**
	 * Computes hemispherical-directional reflectance:\n
	 * total reflection in given direction wo from constant illumination over hemisphere (direct illumination)
	 * @param wo
	 * @param n_samples
	 * @param samples
	 * @return
	 */
	virtual Spectrum rho(const Vector3f& wo, int n_samples, const Point2f* samples) const;

	/**
	 * Computes hemispherical-hemispherical reflectance:\n
	 * fraction of incident light reflected by a surface from diffuse illumination (no direct sources)
	 * @param n_samples
	 * @param samples1
	 * @param samples2
	 * @return
	 */
	virtual Spectrum rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const;

	bool matchesFlags(BxDFType t) const {
		return (type & t) == type;
	}

	BxDFType type;
};

class ScaledBxDF : public BxDF {
public:
	ScaledBxDF(BxDF* _bxdf, const Spectrum& _scale) :
			BxDF(_bxdf->type), bxdf(_bxdf), scale(_scale) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	BxDF* bxdf;
	Spectrum scale;
};

class Fresnel {
public:
	virtual ~Fresnel();

	virtual Spectrum evaluate(float cos_theta_I) const = 0;
};

class FresnelConductor : public Fresnel {
public:
	FresnelConductor(const Spectrum& etai, const Spectrum& etat, const Spectrum& k) :
			eta_I(etai), eta_T(etat), k(k) {}

	Spectrum evaluate(float cos_theta_I) const override;

private:
	Spectrum eta_I, eta_T, k;
};

class FresnelDielectric : public Fresnel {
public:
	FresnelDielectric(const float etai, const float etat) :
			eta_I(etai), eta_T(etat) {}

	Spectrum evaluate(float cos_theta_I) const override;

private:
	float eta_I, eta_T;
};

class FresnelNoOp : public Fresnel {
public:
	Spectrum evaluate(float cos_theta_I) const override { return Spectrum(1.0f); }
};

class SpecularReflection : public BxDF {
public:
	SpecularReflection(const Spectrum& R, Fresnel* f) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), R(R), fresnel(f) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

	Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type = nullptr) const override;

private:
	const Spectrum R;
	const Fresnel* fresnel;
};

class SpecularTransmission : public BxDF {
public:
	SpecularTransmission(const Spectrum& T, float etaA, float etaB, TransportMode mode) :
			BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)), T(T), eta_A(etaA), eta_B(etaB), fresnel(etaA, etaB), mode(mode) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

	Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type = nullptr) const override;

private:
	const Spectrum T;
	const float eta_A, eta_B;
	const FresnelDielectric fresnel;
	const TransportMode mode;
};

class FresnelSpecular : public BxDF {
public:
	FresnelSpecular(const Spectrum& R, const Spectrum& T, float etaA, float etaB, TransportMode mode) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
			R(R),
			T(T),
			eta_A(etaA),
			eta_B(etaB),
			fresnel(etaA, etaB),
			mode(mode) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

	Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type = nullptr) const override;

private:
	const Spectrum R, T;
	const float eta_A, eta_B;
	const FresnelDielectric fresnel;
	const TransportMode mode;
};

class LambertianReflection : public BxDF {
public:
	LambertianReflection(const Spectrum& R) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f&, int, const Point2f*) const override { return R; }
	Spectrum rho(int, const Point2f*, const Point2f*) const override { return R; }

private:
	const Spectrum R;
};

class LambertianTransmission : public BxDF {
public:
	LambertianTransmission(const Spectrum& T) :
			BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)), T(T) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;
	Spectrum rho(const Vector3f&, int, const Point2f*) const override { return T; }
	Spectrum rho(int, const Point2f*, const Point2f*) const override { return T; }
	Spectrum sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf, BxDFType* sampled_type = nullptr) const override;

private:
	const Spectrum T;
};

class OrenNayar : public BxDF {
public:
	/**
	 * Follows the Oren-Nayar diffuse reflection model
	 * @param R reflectance spectrum, i.e. fraction of incident light scattered
	 * @param sigma standard deviation of microfacet orientation angle
	 */
	OrenNayar(const Spectrum& R, float sigma) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {
		sigma = radians(sigma);
		float sigma2 = sigma * sigma;
		A = 1.0f - sigma2 / (2.0f * (sigma2 + 0.33f));
		B = (0.45f * sigma2) / (sigma2 + 0.09f);
	}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	const Spectrum R;
	float A, B;
};

class MicrofacetReflection : public BxDF {
public:
	MicrofacetReflection(const Spectrum& R, MicrofacetDistribution* d, Fresnel* f) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), R(R), distribution(d), fresnel(f) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	const Spectrum R;
	const MicrofacetDistribution* distribution;
	const Fresnel* fresnel;
};

class MicrofacetTransmission : public BxDF {
public:
	MicrofacetTransmission(const Spectrum& T, MicrofacetDistribution* d, float etaa, float etab, TransportMode m) :
			BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)), T(T), distribution(d), eta_A(etaa), eta_B(etab), fresnel(etaa, etab), mode(m) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	const Spectrum T;
	const MicrofacetDistribution* distribution;
	const float eta_A, eta_B;
	const FresnelDielectric fresnel;
	const TransportMode mode;
};

class FresnelBlend : public BxDF {
public:
	FresnelBlend(const Spectrum& rd, const Spectrum& rs, MicrofacetDistribution* d) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), Rd(rd), Rs(rs), distribution(d) {}

	Spectrum schlickFresnel(float cos_theta) const {
		auto pow5 = [](float v) { return (v * v) * (v * v) * v; };
		return CoefficientSpectrum(Rs) + pow5(1 - cos_theta) * (Spectrum(1.0f) - Rs);
	}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	const Spectrum Rd, Rs;
	MicrofacetDistribution* distribution;
};

class FourierBSDF : public BxDF {
public:
	FourierBSDF(const FourierBSDFTable& table, TransportMode m) :
			BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)), bsdf_table(table), mode(m) {}

	Spectrum f(const Vector3f& wo, const Vector3f& wi) const override;

private:
	const FourierBSDFTable& bsdf_table;
	const TransportMode mode;
};

} // namespace civet

#endif // CIVET_REFLECTION_H
