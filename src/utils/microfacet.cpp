#include <utils/microfacet.h>

#include <utils/reflection.h>

namespace civet {

float BeckmannDistribution::D(const Vector3f& wh) const {
	float tan2_theta = tan2Theta(wh);
	if (std::isinf(tan2_theta)) {
		return 0;
	}
	float cos4_theta = cos2Theta(wh) * cos2Theta(wh);
	return std::exp(-tan2_theta * (cosPhi(wh) / (alphax * alphay) + sin2Phi(wh) / (alphax * alphay))) / (Pi * alphax * alphay * cos4_theta);
}

float BeckmannDistribution::lambda(const Vector3f& w) const {
	float abs_tan_theta = std::abs(tanTheta(w));
	if (std::isinf(abs_tan_theta)) {
		return 0;
	}

	float alpha = std::sqrt(cos2Phi(w) * alphax * alphax + sin2Phi(w) * alphay * alphay);
	float a = 1.0f / (alpha * abs_tan_theta);
	if (a >= 1.6f) {
		return 0;
	}
	return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

static void BeckmannSample11(float cos_theta_I, float U1, float U2, float* slope_x, float* slope_y) {
	/* Special case (normal incidence) */
	if (cos_theta_I > .9999) {
		float r = std::sqrt(-std::log(1.0f - U1));
		float sin_phi = std::sin(2 * Pi * U2);
		float cos_phi = std::cos(2 * Pi * U2);
		*slope_x = r * cos_phi;
		*slope_y = r * sin_phi;
		return;
	}

	/* The original inversion routine from the paper contained
	   discontinuities, which causes issues for QMC integration
	   and techniques like Kelemen-style MLT. The following code
	   performs a numerical inversion with better behavior */
	float sin_theta_I = std::sqrt(std::max(0.f, 1.f - cos_theta_I * cos_theta_I));
	float tan_theta_I = sin_theta_I / cos_theta_I;
	float cot_theta_I = 1 / tan_theta_I;

	/* Search interval -- everything is parameterized
	   in the Erf() domain */
	float a = -1, c = erf(cot_theta_I);
	float sample_x = std::max(U1, (float)1e-6f);

	/* Start with a good initial guess */
	// float b = (1-sample_x) * a + sample_x * c;

	/* We can do better (inverse of an approximation computed in
	 * Mathematica) */
	float thetaI = std::acos(cos_theta_I);
	float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
	float b = c - (1 + c) * std::pow(1 - sample_x, fit);

	/* Normalization factor for the CDF */
	static const float SQRT_PI_INV = 1.f / std::sqrt(Pi);
	float normalization = 1 / (1 + c + SQRT_PI_INV * tan_theta_I * std::exp(-cot_theta_I * cot_theta_I));

	int it = 0;
	while (++it < 10) {
		/* Bisection criterion -- the oddly-looking
		   Boolean expression are intentional to check
		   for NaNs at little additional cost */
		if (!(b >= a && b <= c)) {
			b = 0.5f * (a + c);
		}

		/* Evaluate the CDF and its derivative
		   (i.e. the density function) */
		float inv_erf = erfInv(b);
		float value = normalization * (1 + b + SQRT_PI_INV * tan_theta_I * std::exp(-inv_erf * inv_erf)) - sample_x;
		float derivative = normalization * (1 - inv_erf * tan_theta_I);

		if (std::abs(value) < 1e-5f) {
			break;
		}

		/* Update bisection intervals */
		if (value > 0) {
			c = b;
		} else {
			a = b;
		}

		b -= value / derivative;
	}

	/* Now convert back into a slope value */
	*slope_x = erfInv(b);

	/* Simulate Y component */
	*slope_y = erfInv(2.0f * std::max(U2, 1e-6f) - 1.0f);
}

static Vector3f BeckmannSample(const Vector3f& wi, float alpha_x, float alpha_y,
		float U1, float U2) {
	// 1. stretch wi
	Vector3f wiStretched = normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
	float slope_x, slope_y;
	BeckmannSample11(cosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

	// 3. rotate
	float tmp = cosPhi(wiStretched) * slope_x - sinPhi(wiStretched) * slope_y;
	slope_y = sinPhi(wiStretched) * slope_x + cosPhi(wiStretched) * slope_y;
	slope_x = tmp;

	// 4. unstretch
	slope_x = alpha_x * slope_x;
	slope_y = alpha_y * slope_y;

	// 5. compute normal
	return normalize(Vector3f(-slope_x, -slope_y, 1.f));
}

Vector3f BeckmannDistribution::sample_wh(const Vector3f& wo, const Point2f& u) const {
	if (!sample_visible_area) {
		// sample full distribution of normals
		float tan2_theta, phi;
		if (alphax == alphay) {
			float log_sample = std::log(1 - u[0]);
			tan2_theta = -alphax * alphax * log_sample;
			phi = u[1] * 2 * Pi;
		} else {
			float log_sample = std::log(1 - u[0]);
			phi = std::atan(alphay / alphax * std::tan(2 * Pi * u[1] + 0.5f * Pi));
			if (u[1] > 0.5f) {
				phi += Pi;
			}
			float sin_phi = std::sin(phi), cos_phi = std::cos(phi);
			float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			tan2_theta = -log_sample / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
		}
		float cos_theta = 1 / std::sqrt(1 + tan2_theta);
		float sin_theta = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
		Vector3f wh = sphericalDirection(sin_theta, cos_theta, phi);
		if (!sameHemisphere(wo, wh)) {
			wh = -wh;
		}
		return wh;
	} else {
		// sample visible area of normals
		Vector3f wh;
		bool flip = wo.z < 0;
		wh = BeckmannSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) {
			wh = -wh;
		}
		return wh;
	}
}

float TrowbridgeReitzDistribution::D(const Vector3f& wh) const {
	float tan2_theta = tan2Theta(wh);
	if (std::isinf(tan2_theta)) {
		return 0;
	}
	float cos4_theta = cos2Theta(wh) * cos2Theta(wh);
	float e = (cos2Phi(wh) / (alphax * alphax) + sin2Phi(wh) / (alphay * alphay)) * tan2_theta;
	return 1 / (Pi * alphax * alphay * cos4_theta * (1 + e) * (1 + e));
}

static void TrowbridgeReitzSample11(float cos_theta, float U1, float U2, float* slope_x, float* slope_y) {
	// special case (normal incidence)
	if (cos_theta > .9999) {
		float r = std::sqrt(U1 / (1 - U1));
		float phi = 6.28318530718 * U2;
		*slope_x = r * cos(phi);
		*slope_y = r * sin(phi);
		return;
	}

	float sin_theta = std::sqrt(std::max(0.f, 1.f - cos_theta * cos_theta));
	float tan_theta = sin_theta / cos_theta;
	float a = 1 / tan_theta;
	float G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));

	// sample slope_x
	float A = 2 * U1 / G1 - 1;
	float tmp = 1.f / (A * A - 1.f);
	if (tmp > 1e10) {
		tmp = 1e10;
	}
	float B = tan_theta;
	float D = std::sqrt(std::max(float(B * B * tmp * tmp - (A * A - B * B) * tmp), 0.f));
	float slope_x_1 = B * tmp - D;
	float slope_x_2 = B * tmp + D;
	*slope_x = (A < 0 || slope_x_2 > 1.f / tan_theta) ? slope_x_1 : slope_x_2;

	// sample slope_y
	float S;
	if (U2 > 0.5f) {
		S = 1.f;
		U2 = 2.f * (U2 - .5f);
	} else {
		S = -1.f;
		U2 = 2.f * (.5f - U2);
	}
	float z =
			(U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
			(U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
	*slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);
}

static Vector3f TrowbridgeReitzSample(const Vector3f& wi, float alpha_x,
		float alpha_y, float U1, float U2) {
	// 1. stretch wi
	Vector3f wi_stretched = normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
	float slope_x, slope_y;
	TrowbridgeReitzSample11(cosTheta(wi_stretched), U1, U2, &slope_x, &slope_y);

	// 3. rotate
	float tmp = cosPhi(wi_stretched) * slope_x - sinPhi(wi_stretched) * slope_y;
	slope_y = sinPhi(wi_stretched) * slope_x + cosPhi(wi_stretched) * slope_y;
	slope_x = tmp;

	// 4. unstretch
	slope_x = alpha_x * slope_x;
	slope_y = alpha_y * slope_y;

	// 5. compute normal
	return normalize(Vector3f(-slope_x, -slope_y, 1.));
}

Vector3f TrowbridgeReitzDistribution::sample_wh(const Vector3f& wo, const Point2f& u) const {
	if (!sample_visible_area) {
		// sample full distribution of normals
		Vector3f wh;
		float cos_theta = 0, phi = (2 * Pi) * u[1];
		if (alphax == alphay) {
			float tan_theta2 = alphax * alphax * u[0] / (1.0f - u[0]);
			cos_theta = 1 / std::sqrt(1 + tan_theta2);
		} else {
			phi = std::atan(alphay / alphax * std::tan(2 * Pi * u[1] + .5f * Pi));
			if (u[1] > .5f) {
				phi += Pi;
			}
			float sin_phi = std::sin(phi), cos_phi = std::cos(phi);
			const float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			const float alpha2 = 1 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
			float tan_theta2 = alpha2 * u[0] / (1 - u[0]);
			cos_theta = 1 / std::sqrt(1 + tan_theta2);
		}
		float sin_theta = std::sqrt(std::max(0.f, 1.f - cos_theta * cos_theta));
		wh = sphericalDirection(sin_theta, cos_theta, phi);
		if (!sameHemisphere(wo, wh)) {
			wh = -wh;
		}
		return wh;
	} else {
		// sample visible area of normals
		Vector3f wh;
		bool flip = wo.z < 0;
		wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) {
			wh = -wh;
		}
		return wh;
	}
}

float TrowbridgeReitzDistribution::lambda(const Vector3f& w) const {
	float abs_tan_theta = std::abs(tanTheta(w));
	if (std::isinf(abs_tan_theta)) {
		return 0;
	}

	float alpha = std::sqrt(cos2Phi(w) * alphax * alphax + sin2Phi(w) * alphay * alphay);
	float alpha2_tan2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
	return (-1 + std::sqrt(1.0f + alpha2_tan2_theta)) / 2;
}

float MicrofacetDistribution::pdf(const Vector3f& wo, const Vector3f& wh) const {
	if (sample_visible_area) {
		return D(wh) * G1(wo) * absDot(wo, wh) / absCosTheta(wo);
	} else {
		return D(wh) * absCosTheta(wh);
	}
}

} // namespace civet