#include <utils/reflection.h>

#include <core/interaction.h>
#include <core/sampler.h>
#include <core/spectrum.h>
#include <utils/interpolation.h>
#include <utils/sampling.h>

namespace civet {

float frDielectric(float cos_theta_I, float eta_I, float eta_T) {
	cos_theta_I = clamp(cos_theta_I, -1, 1);

	bool entering = cos_theta_I > 0.0f;
	if (!entering) {
		swapElem(eta_I, eta_T);
		cos_theta_I = std::abs(cos_theta_I);
	}

	float sin_theta_I = std::sqrt(std::max(0.0f, 1 - cos_theta_I * cos_theta_I));
	float sin_theta_T = eta_I / eta_T * sin_theta_I;
	if (sin_theta_T >= 1) {
		return 1;
	}
	float cos_theta_T = std::sqrt(std::max(0.0f, 1 - sin_theta_T * sin_theta_T));

	float r_para = ((eta_T * cos_theta_I) - (eta_I * cos_theta_T)) / ((eta_T * cos_theta_I) + (eta_I * cos_theta_T));
	float r_perp = ((eta_I * cos_theta_I) - (eta_T * cos_theta_T)) / ((eta_I * cos_theta_I) + (eta_T * cos_theta_T));
	return (r_para * r_para + r_perp * r_perp) / 2;
}

Spectrum frConductor(float cos_theta_I, const Spectrum& eta_I, const Spectrum& eta_T, const Spectrum& k) {
	cos_theta_I = clamp(cos_theta_I, -1, 1);
	Spectrum eta = eta_T / eta_I;
	Spectrum eta_k = k / eta_I;

	float cos_theta_I2 = cos_theta_I * cos_theta_I;
	float sin_theta_I2 = 1.0f - cos_theta_I2;
	Spectrum eta2 = eta * eta;
	Spectrum etak2 = eta_k * eta_k;

	Spectrum t0 = eta2 - etak2 - sin_theta_I2;
	Spectrum a2pb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
	Spectrum t1 = a2pb2 + cos_theta_I2;
	Spectrum a = sqrt(0.5f * (a2pb2 + t0));
	Spectrum t2 = 2.0f * cos_theta_I * a;
	Spectrum r_perp = (t1 - t2) / (t1 + t2);

	Spectrum t3 = cos_theta_I2 * a2pb2 + sin_theta_I2 * sin_theta_I2;
	Spectrum t4 = t2 * sin_theta_I2;
	Spectrum r_para = r_perp * (t3 - t4) / (t3 + t4);

	return 0.5f * (r_perp + r_para);
}

bool FourierBSDFTable::getWeightsAndOffset(float cos_theta, int* offset, float weights[4]) const {
	return catmullRomWeights(n_mu, mu, cos_theta, offset, weights);
}

Spectrum BSDF::f(const Vector3f& wo_world, const Vector3f& wi_world, BxDFType flags) const {
	Vector3f wi = worldToLocal(wi_world), wo = worldToLocal(wo_world);
	bool reflect = dot(wi_world, ng) * dot(wo_world, ng) > 0;

	Spectrum f(0.0f);
	for (int i = 0; i < num_BxDFs; i++) {
		if (bxdfs[i]->matchesFlags(flags) &&
				((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
						(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION)))) {
			f += bxdfs[i]->f(wo, wi);
		}
	}
	return f;
}

Spectrum BSDF::rho(const Vector3f& wo, int num_samples, const Point2f* samples, BxDFType flags) const {
	Spectrum val(0.0f);
	for (int i = 0; i < num_BxDFs; i++) {
		if (bxdfs[i]->matchesFlags(flags)) {
			val += bxdfs[i]->rho(wo, num_samples, samples);
		}
	}
	return val;
}

Spectrum BSDF::rho(int num_samples, const Point2f* samples1, const Point2f* samples2, BxDFType flags) const {
	Spectrum val(0.0f);
	for (int i = 0; i < num_BxDFs; i++) {
		if (bxdfs[i]->matchesFlags(flags)) {
			val += bxdfs[i]->rho(num_samples, samples1, samples2);
		}
	}
	return val;
}

Spectrum BSDF::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& u, float* pdf_val, BxDFType flags, BxDFType* sampled_type) const {
	int matching_comps = numComponents(flags);
	if (matching_comps == 0) {
		*pdf_val = 0;
		return Spectrum(0.f);
	}
	int comp = std::min((int)std::floor(u[0] * matching_comps), matching_comps - 1);
	BxDF* bxdf = nullptr;
	int count = comp;
	for (int i = 0; i < num_BxDFs; i++) {
		if (bxdfs[i]->matchesFlags(flags) && count-- == 0) {
			bxdf = bxdfs[i];
			break;
		}
	}

	Point2f u_remapped(u[0] * matching_comps - comp, u[1]);
	Vector3f wil, wol = worldToLocal(wo);
	*pdf_val = 0;
	if (sampled_type) {
		*sampled_type = bxdf->type;
	}
	Spectrum f = bxdf->sample_f(wol, &wil, u_remapped, pdf_val, sampled_type);
	if (*pdf_val == 0) {
		return 0;
	}
	*wi = localToWorld(wil);

	if (!(bxdf->type & BSDF_SPECULAR) && matching_comps > 1) {
		for (int i = 0; i < num_BxDFs; i++) {
			if (bxdfs[i] != bxdf && bxdfs[i]->matchesFlags(flags)) {
				*pdf_val += bxdfs[i]->pdf(wol, wil);
			}
		}
	}
	if (matching_comps > 1) {
		*pdf_val /= matching_comps;
	}

	if (!(bxdf->type & BSDF_SPECULAR) && matching_comps > 1) {
		bool reflect = dot(*wi, ng) * dot(wo, ng) > 0;
		f = 0;
		for (int i = 0; i < num_BxDFs; i++) {
			if (bxdfs[i]->matchesFlags(flags) &&
					((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
							(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION)))) {
				f += bxdfs[i]->f(wol, wil);
			}
		}
	}
	return f;
}

float BSDF::pdf(const Vector3f& wo, const Vector3f& wi, BxDFType flags) const {
	if (num_BxDFs == 0.f) {
		return 0.f;
	}
	Vector3f wol = worldToLocal(wo), wil = worldToLocal(wi);
	if (wol.z == 0) {
		return 0.;
	}
	float pdf_val = 0.f;
	int matching_comps = 0;
	for (int i = 0; i < num_BxDFs; i++) {
		if (bxdfs[i]->matchesFlags(flags)) {
			++matching_comps;
			pdf_val += bxdfs[i]->pdf(wol, wil);
		}
	}
	float v = matching_comps > 0 ? pdf_val / matching_comps : 0.f;
	return v;
}

Spectrum BxDF::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	*wi = cosineSampleHemisphere(sample);
	if (wo.z < 0) {
		wi->z *= -1;
	}
	*pdf_val = pdf(wo, *wi);
	return f(wo, *wi);
}

Spectrum BxDF::rho(const Vector3f& wo, int n_samples, const Point2f* samples) const {
	Spectrum r(0.f);
	for (int i = 0; i < n_samples; ++i) {
		// Estimate one term of rho_hd
		Vector3f wi;
		float pdf_val = 0;
		Spectrum f = sample_f(wo, &wi, samples[i], &pdf_val);
		if (pdf_val > 0) {
			r += f * absCosTheta(wi) / pdf_val;
		}
	}
	return r / n_samples;
}

Spectrum BxDF::rho(int n_samples, const Point2f* samples1, const Point2f* samples2) const {
	Spectrum r(0.f);
	for (int i = 0; i < n_samples; ++i) {
		// Estimate one term of rho_hh
		Vector3f wo, wi;
		wo = uniformSampleHemisphere(samples1[i]);
		float pdfo = uniformHemispherePdf(), pdfi = 0;
		Spectrum f = sample_f(wo, &wi, samples2[i], &pdfi);
		if (pdfi > 0) {
			r += f * absCosTheta(wi) * absCosTheta(wo) / (pdfo * pdfi);
		}
	}
	return r / (Pi * n_samples);
}

float BxDF::pdf(const Vector3f& wo, const Vector3f& wi) const {
	return sameHemisphere(wo, wi) ? absCosTheta(wi) * InvPi : 0;
}

Spectrum ScaledBxDF::f(const Vector3f& wo, const Vector3f& wi) const {
	return scale * bxdf->f(wo, wi);
}

Fresnel::~Fresnel() {}

Spectrum FresnelConductor::evaluate(float cos_theta_I) const {
	return frConductor(std::abs(cos_theta_I), eta_I, eta_T, k);
}

Spectrum FresnelDielectric::evaluate(float cos_theta_I) const {
	return frDielectric(cos_theta_I, eta_I, eta_T);
}

Spectrum SpecularReflection::f(const Vector3f& wo, const Vector3f& wi) const {
	return Spectrum(0.0f);
}

Spectrum SpecularReflection::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	*wi = Vector3f(-wo.x, -wo.y, wo.z);
	*pdf_val = 1;
	return fresnel->evaluate(cosTheta(*wi)) * R / absCosTheta(*wi);
}

Spectrum SpecularTransmission::f(const Vector3f& wo, const Vector3f& wi) const {
	return Spectrum(0.0f);
}

Spectrum SpecularTransmission::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	bool entering = cosTheta(wo) > 0;
	float etai = entering ? eta_A : eta_B;
	float etat = entering ? eta_B : eta_A;

	if (!refract(wo, faceforward(Normal3f(0, 0, 1), wo), etai / etat, wi)) {
		return 0;
	}

	*pdf_val = 1;
	Spectrum ft = T * Spectrum(1.0f) - fresnel.evaluate(cosTheta(*wi));

	if (mode == TransportMode::Radiance) {
		ft *= (etai * etai) / (etat * etat);
	}

	return ft / absCosTheta(*wi);
}

Spectrum FresnelSpecular::f(const Vector3f& wo, const Vector3f& wi) const {
	return Spectrum(0.0f);
}

Spectrum FresnelSpecular::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	float F = frDielectric(cosTheta(wo), eta_A, eta_B);
	if (sample[0] < F) {
		// Compute specular reflection for _FresnelSpecular_

		// Compute perfect specular reflection direction
		*wi = Vector3f(-wo.x, -wo.y, wo.z);
		if (sampled_type) {
			*sampled_type = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
		}
		*pdf_val = F;
		return F * R / absCosTheta(*wi);
	} else {
		// Compute specular transmission for _FresnelSpecular_

		// Figure out which $\eta$ is incident and which is transmitted
		bool entering = cosTheta(wo) > 0;
		float etaI = entering ? eta_A : eta_B;
		float etaT = entering ? eta_B : eta_A;

		// Compute ray direction for specular transmission
		if (!refract(wo, faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi)) {
			return 0;
		}
		Spectrum ft = T * (1 - F);

		// Account for non-symmetry with transmission to different medium
		if (mode == TransportMode::Radiance) {
			ft *= (etaI * etaI) / (etaT * etaT);
		}
		if (sampled_type) {
			*sampled_type = BxDFType(BSDF_SPECULAR | BSDF_TRANSMISSION);
		}
		*pdf_val = 1 - F;
		return ft / absCosTheta(*wi);
	}
}

Spectrum LambertianReflection::f(const Vector3f& wo, const Vector3f& wi) const {
	return R * InvPi;
}

Spectrum LambertianTransmission::f(const Vector3f& wo, const Vector3f& wi) const {
	return T * InvPi;
}

Spectrum LambertianTransmission::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	*wi = cosineSampleHemisphere(sample);
	if (wo.z > 0) {
		wi->z *= -1;
	}
	*pdf_val = pdf(wo, *wi);
	return f(wo, *wi);
}

Spectrum OrenNayar::f(const Vector3f& wo, const Vector3f& wi) const {
	float sin_theta_I = sinTheta(wi);
	float sin_theta_O = sinTheta(wo);

	float max_cos = 0;
	if (sin_theta_I > 1e-4f && sin_theta_O > 1e-4f) {
		float sin_phi_I = sinPhi(wi), cos_phi_I = cosPhi(wi);
		float sin_phi_O = sinPhi(wo), cos_phi_O = cosPhi(wo);
		float d_cos = cos_phi_I * cos_phi_O + sin_phi_I * sin_phi_O;
		max_cos = std::max(0.0f, d_cos);
	}

	float sin_alpha, tan_beta;
	if (absCosTheta(wi) > absCosTheta(wo)) {
		sin_alpha = sin_theta_O;
		tan_beta = sin_theta_I / absCosTheta(wi);
	} else {
		sin_alpha = sin_theta_I;
		tan_beta = sin_theta_O / absCosTheta(wo);
	}

	return R * InvPi * (A + B * max_cos * sin_alpha * tan_beta);
}

Spectrum MicrofacetReflection::f(const Vector3f& wo, const Vector3f& wi) const {
	float cos_theta_O = absCosTheta(wo), cos_theta_I = absCosTheta(wi);
	Vector3f wh = wi + wo;
	if (cos_theta_I == 0 || cos_theta_O == 0) {
		return Spectrum(0.0f);
	}
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.0f);
	}
	wh = normalize(wh);
	Spectrum F = fresnel->evaluate(dot(wi, wh));
	return R * distribution->D(wh) * distribution->G(wo, wi) * F / (4 * cos_theta_I * cos_theta_O);
}

Spectrum MicrofacetReflection::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	Vector3f wh = distribution->sample_wh(wo, sample);
	*wi = reflect(wo, wh);
	if (!sameHemisphere(wo, *wi)) {
		return Spectrum(0.f);
	}
	*pdf_val = distribution->pdf(wo, wh) / (4 * dot(wo, wh));
	return f(wo, *wi);
}

float MicrofacetReflection::pdf(const Vector3f& wo, const Vector3f& wi) const {
	if (!sameHemisphere(wo, wi)) {
		return 0.f;
	}
	Vector3f wh = normalize(wo + wi);
	return distribution->pdf(wo, wh) / (4 * dot(wo, wh));
}

Spectrum MicrofacetTransmission::f(const Vector3f& wo, const Vector3f& wi) const {
	if (sameHemisphere(wo, wi)) {
		return 0;
	}

	float cos_theta_O = cosTheta(wo), cos_theta_I = cosTheta(wi);
	if (cos_theta_I == 0 || cos_theta_O == 0) {
		return Spectrum(0.0f);
	}

	float eta = cos_theta_O > 0 ? (eta_B / eta_A) : (eta_A / eta_B);
	Vector3f wh = normalize(wo + wi * eta);
	if (wh.z < 0) {
		wh = -wh;
	}

	if (dot(wo, wh) * dot(wi, wh) > 0) {
		// on same side
		return Spectrum(0.0f);
	}

	Spectrum F = fresnel.evaluate(dot(wo, wh));
	float sqrt_denom = dot(wo, wh) + eta * dot(wi, wh);
	float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;

	return (Spectrum(1.0f) - F) * T * std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta * absDot(wi, wh) * absDot(wo, wh) * factor * factor / (cos_theta_I * cos_theta_O * sqrt_denom * sqrt_denom));
}

Spectrum MicrofacetTransmission::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	Vector3f wh = distribution->sample_wh(wo, sample);
	float eta = cosTheta(wo) > 0 ? (eta_A / eta_B) : (eta_B / eta_A);
	if (!refract(wo, Normal3f(wh), eta, wi)) {
		return 0;
	}
	*pdf_val = pdf(wo, *wi);
	return f(wo, *wi);
}

float MicrofacetTransmission::pdf(const Vector3f& wo, const Vector3f& wi) const {
	if (sameHemisphere(wo, wi)) {
		return 0;
	}
	float eta = cosTheta(wo) > 0 ? (eta_B / eta_A) : (eta_A / eta_B);
	Vector3f wh = normalize(wo + wi * eta);
	float sqrt_denom = dot(wo, wh) + eta * dot(wi, wh);
	float dwh_dwi = std::abs((eta * eta * dot(wi, wh)) / (sqrt_denom * sqrt_denom));
	return distribution->pdf(wo, wh) * dwh_dwi;
}

Spectrum FresnelBlend::f(const Vector3f& wo, const Vector3f& wi) const {
	auto pow5 = [](float v) { return (v * v) * (v * v) * v; };
	Spectrum diffuse = (28.f / (23.f * Pi)) * Rd * (Spectrum(1.f) - Rs) *
			(1 - pow5(1 - .5f * absCosTheta(wi))) *
			(1 - pow5(1 - .5f * absCosTheta(wo)));
	Vector3f wh = wi + wo;
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
		return Spectrum(0.0f);
	}
	wh = normalize(wh);
	Spectrum specular = distribution->D(wh) / (4 * absDot(wi, wh) * std::max(absCosTheta(wi), absCosTheta(wo))) * schlickFresnel(dot(wi, wh));
	return diffuse + specular;
}

Spectrum FresnelBlend::sample_f(const Vector3f& wo, Vector3f* wi, const Point2f& sample, float* pdf_val, BxDFType* sampled_type) const {
	Point2f u = sample;
	if (u[0] < 0.5f) {
		u[0] = 2 * u[0];
		*wi = cosineSampleHemisphere(u);
		if (wo.z < 0) {
			wi->z *= -1;
		}
	} else {
		u[0] = 2 * (u[0] - 0.5f);
		Vector3f wh = distribution->sample_wh(wo, u);
		*wi = reflect(wo, wh);
		if (!sameHemisphere(wo, *wi)) {
			return Spectrum(0.f);
		}
	}
	*pdf_val = pdf(wo, *wi);
	return f(wo, *wi);
}

float FresnelBlend::pdf(const Vector3f& wo, const Vector3f& wi) const {
	if (!sameHemisphere(wo, wi)) {
		return 0;
	}
	Vector3f wh = normalize(wo + wi);
	float pdf_wh = distribution->pdf(wo, wh);
	return 0.5f * (absCosTheta(wi) * InvPi + pdf_wh / (4 * dot(wo, wh)));
}

Spectrum FourierBSDF::f(const Vector3f& wo, const Vector3f& wi) const {
	// find zenith angle cosines and azimuth difference angle
	float mu_I = cosTheta(-wi), mu_O = cosTheta(wo);
	float cos_phi = cosDPhi(-wi, wo);

	// compute fourier coefficients
	int offseti, offseto;
	float weightsi[4], weightso[4];
	if (!bsdf_table.getWeightsAndOffset(mu_I, &offseti, weightsi) || !bsdf_table.getWeightsAndOffset(mu_O, &offseto, weightso)) {
		return Spectrum(0.0f);
	}

	float* ak = ALLOCA(float, bsdf_table.m_max* bsdf_table.n_channels);
	memset(ak, 0, bsdf_table.m_max * bsdf_table.n_channels * sizeof(float));

	int m_max = 0;
	for (int b = 0; b < 4; b++) {
		for (int a = 0; a < 4; a++) {
			float weight = weightsi[a] * weightso[b];
			if (weight != 0) {
				int m;
				const float* ap = bsdf_table.getAk(offseti + a, offseto + b, &m);
				m_max = std::max(m_max, m);
				for (int c = 0; c < bsdf_table.n_channels; c++) {
					for (int k = 0; k < m; k++) {
						ak[c * bsdf_table.m_max + k] += weight * ap[c * m + k];
					}
				}
			}
		}
	}

	// evaluate fourier expansion for angle phi
	float Y = std::max(0.0f, fourier(ak, m_max, cos_phi));
	float scale = mu_I != 0 ? (1 / std::abs(mu_I)) : 0.0f;
	if (mode == TransportMode::Radiance && mu_I * mu_O > 0) {
		float eta = mu_I > 0 ? 1 / bsdf_table.eta : bsdf_table.eta;
		scale *= eta * eta;
	}

	if (bsdf_table.n_channels == 1) {
		return Spectrum(Y * scale);
	} else {
		float R = fourier(ak + 1 * bsdf_table.m_max, m_max, cos_phi);
		float B = fourier(ak + 2 * bsdf_table.m_max, m_max, cos_phi);
		float G = 1.39829f * Y - 0.100913f * B - 0.297375f * R;
		float rgb[3] = { R * scale, G * scale, B * scale };
		return Spectrum::fromRGB(rgb).clamp();
	}
}

} // namespace civet