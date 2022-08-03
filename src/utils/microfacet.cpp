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

float TrowbridgeReitzDistribution::D(const Vector3f& wh) const {
	float tan2_theta = tan2Theta(wh);
	if (std::isinf(tan2_theta)) {
		return 0;
	}
	float cos4_theta = cos2Theta(wh) * cos2Theta(wh);
	float e = (cos2Phi(wh) / (alphax * alphax) + sin2Phi(wh) / (alphay * alphay)) * tan2_theta;
	return 1 / (Pi * alphax * alphay * cos4_theta * (1 + e) * (1 + e));
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

} // namespace civet