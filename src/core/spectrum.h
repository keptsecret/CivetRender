#ifndef CIVET_SPECTRUM_H
#define CIVET_SPECTRUM_H

#include <core/civet.h>

namespace civet {

static const int sampled_wavelength_start = 400;
static const int sampled_wavelength_end = 700;
static const int N_spectral_samples = 60;

enum class SpectrumType { Reflectance,
	Illuminant };

// Spectral XYZ data from CIE
static const int N_CIE_samples = 471;
extern const float CIE_X[N_CIE_samples];
extern const float CIE_Y[N_CIE_samples];
extern const float CIE_Z[N_CIE_samples];
extern const float CIE_lambda[N_CIE_samples];
static const float CIE_Y_integral = 106.856895f;

extern bool spectrumSamplesSorted(const float* lambda, const float* vals, int n);
extern void sortSpectrumSamples(float* lambda, float* vals, int n);
extern float averageSpectrumSamples(const float* lambda, const float* vals, int n, float l_start, float l_end);
extern float interpolateSpectrumSamples(const float* lambda, const float* vals, int n, float l);

extern void blackbody(const float* lambda, int n, float T, float* Le);
extern void blackbodyNormalized(const float* lambda, int n, float T, float* Le);

// Spectral data for RGB to XYZ
static const int N_RGB2Spect_samples = 32;
extern const float RGB2SpectLambda[N_RGB2Spect_samples];
extern const float RGBRefl2SpectWhite[N_RGB2Spect_samples];
extern const float RGBRefl2SpectCyan[N_RGB2Spect_samples];
extern const float RGBRefl2SpectMagenta[N_RGB2Spect_samples];
extern const float RGBRefl2SpectYellow[N_RGB2Spect_samples];
extern const float RGBRefl2SpectRed[N_RGB2Spect_samples];
extern const float RGBRefl2SpectGreen[N_RGB2Spect_samples];
extern const float RGBRefl2SpectBlue[N_RGB2Spect_samples];

extern const float RGBIllum2SpectWhite[N_RGB2Spect_samples];
extern const float RGBIllum2SpectCyan[N_RGB2Spect_samples];
extern const float RGBIllum2SpectMagenta[N_RGB2Spect_samples];
extern const float RGBIllum2SpectYellow[N_RGB2Spect_samples];
extern const float RGBIllum2SpectRed[N_RGB2Spect_samples];
extern const float RGBIllum2SpectGreen[N_RGB2Spect_samples];
extern const float RGBIllum2SpectBlue[N_RGB2Spect_samples];

// Utility functions
inline void XYZtoRGB(const float xyz[3], float rgb[3]) {
	rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
	rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
	rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}

inline void RGBtoXYZ(const float rgb[3], float xyz[3]) {
	xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
	xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
	xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
}

template <int n_spectrum_samples>
class CoefficientSpectrum {
public:
	CoefficientSpectrum(float v = 0.0f) {
		for (int i = 0; i < n_spectrum_samples; i++) {
			c[i] = v;
		}
	}

	CoefficientSpectrum& operator+=(const CoefficientSpectrum& s2) {
		for (int i = 0; i < n_spectrum_samples; i++) {
			c[i] += s2.c[i];
		}
		return *this;
	}

	CoefficientSpectrum operator+(const CoefficientSpectrum& s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; i++) {
			ret.c[i] += s2.c[i];
		}
		return ret;
	}

	CoefficientSpectrum operator-() const {
		CoefficientSpectrum ret;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] = -c[i];
		}
		return ret;
	}

	CoefficientSpectrum operator-(const CoefficientSpectrum& s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] -= s2.c[i];
		}
		return ret;
	}

	CoefficientSpectrum operator/(const CoefficientSpectrum& s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] /= s2.c[i];
		}
		return ret;
	}

	CoefficientSpectrum operator*(const CoefficientSpectrum& sp) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] *= sp.c[i];
		}
		return ret;
	}

	CoefficientSpectrum& operator*=(const CoefficientSpectrum& sp) {
		for (int i = 0; i < n_spectrum_samples; ++i) {
			c[i] *= sp.c[i];
		}
		return *this;
	}

	CoefficientSpectrum operator*(float a) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] *= a;
		}
		return ret;
	}

	CoefficientSpectrum& operator*=(float a) {
		for (int i = 0; i < n_spectrum_samples; ++i) {
			c[i] *= a;
		}
		return *this;
	}

	friend inline CoefficientSpectrum operator*(float a, const CoefficientSpectrum& s) {
		return s * a;
	}

	CoefficientSpectrum operator/(float a) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] /= a;
		}
		return ret;
	}

	CoefficientSpectrum& operator/=(float a) {
		for (int i = 0; i < n_spectrum_samples; ++i) {
			c[i] /= a;
		}
		return *this;
	}

	bool operator==(const CoefficientSpectrum& sp) const {
		for (int i = 0; i < n_spectrum_samples; ++i) {
			if (c[i] != sp.c[i]) {
				return false;
			}
		}
		return true;
	}

	bool operator!=(const CoefficientSpectrum& sp) const {
		return *this != sp;
	}

	bool isBlack() const {
		for (int i = 0; i < n_spectrum_samples; i++) {
			if (c[i] != 0) {
				return false;
			}
		}
		return true;
	}

	friend CoefficientSpectrum sqrt(const CoefficientSpectrum& s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < n_spectrum_samples; i++) {
			ret.c[i] = std::sqrt(s.c[i]);
		}
		return ret;
	}

	template <int n>
	friend inline CoefficientSpectrum<n> pow(const CoefficientSpectrum<n>& s, float e);

	friend CoefficientSpectrum exp(const CoefficientSpectrum& s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < n_spectrum_samples; ++i) {
			ret.c[i] = std::exp(s.c[i]);
		}
		return ret;
	}

	CoefficientSpectrum clamp(float low = 0, float high = Infinity) const {
		CoefficientSpectrum ret;
		for (int i = 0; i < n_spectrum_samples; i++) {
			ret.c[i] = civet::clamp(c[i], low, high);
		}
		return ret;
	}

	float maxComponentValue() const {
		float m = c[0];
		for (int i = 1; i < N_spectral_samples; ++i) {
			m = std::max(m, c[i]);
		}
		return m;
	}

	bool hasNaNs() const {
		for (int i = 0; i < n_spectrum_samples; i++) {
			if (isNaN(c[i])) {
				return true;
			}
		}
		return false;
	}

	float& operator[](int i) {
		return c[i];
	}

	float operator[](int i) const {
		return c[i];
	}

	static const int n_samples = n_spectrum_samples;

protected:
	float c[n_spectrum_samples];
};

class SampledSpectrum : public CoefficientSpectrum<N_spectral_samples> {
public:
	SampledSpectrum(float v = 0.0f) :
			CoefficientSpectrum(v) {}

	SampledSpectrum(const CoefficientSpectrum<N_spectral_samples>& v) :
			CoefficientSpectrum<N_spectral_samples>(v) {}

	SampledSpectrum(const RGBSpectrum& r, SpectrumType type);

	static SampledSpectrum fromSampled(const float* lambda, const float* v, int n);

	// TODO: call on Engine start up if used
	static void init();

	void toXYZ(float xyz[3]) const {
		xyz[0] = xyz[1] = xyz[2] = 0.0f;
		for (int i = 0; i < N_spectral_samples; i++) {
			xyz[0] += X.c[i] * c[i];
			xyz[1] += Y.c[i] * c[i];
			xyz[2] += Z.c[i] * c[i];
		}
		float scale = float(sampled_wavelength_end - sampled_wavelength_start) / float(CIE_Y_integral * N_spectral_samples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
	}

	void toRGB(float rgb[3]) const {
		float xyz[3];
		toXYZ(xyz);
		XYZtoRGB(xyz, rgb);
	}

	SampledSpectrum fromRGB(const float rgb[3], SpectrumType type);

	float y() const {
		///< y coeff correlates with luminance and accessed often, so provide one for it
		float yy = 0.0f;
		for (int i = 0; i < N_spectral_samples; i++) {
			yy += Y.c[i] * c[i];
		}
		return yy * float(sampled_wavelength_end - sampled_wavelength_start) / float(N_spectral_samples);
	}

	RGBSpectrum toRGBSpectrum() const;

private:
	static SampledSpectrum X, Y, Z;

	static SampledSpectrum rgb_refl2spect_white, rgb_refl2spect_cyan;
	static SampledSpectrum rgb_refl2spect_magenta, rgb_refl2spect_yellow;
	static SampledSpectrum rgb_refl2spect_red, rgb_refl2spect_green;
	static SampledSpectrum rgb_refl2spect_blue;

	static SampledSpectrum rgb_illum2spect_white, rgb_illum2spect_cyan;
	static SampledSpectrum rgb_illum2spect_magenta, rgb_illum2spect_yellow;
	static SampledSpectrum rgb_illum2spect_red, rgb_illum2spect_green;
	static SampledSpectrum rgb_illum2spect_blue;
};

class RGBSpectrum : public CoefficientSpectrum<3> {
	using CoefficientSpectrum<3>::c;

public:
	RGBSpectrum(float v = 0.0f) :
			CoefficientSpectrum<3>(v) {}

	RGBSpectrum(const CoefficientSpectrum<3>& v) :
			CoefficientSpectrum<3>(v) {}

	RGBSpectrum(const RGBSpectrum& s, SpectrumType type = SpectrumType::Reflectance) {
		*this = s;
	}

	static RGBSpectrum fromRGB(const float rgb[3], SpectrumType type = SpectrumType::Reflectance) {
		RGBSpectrum s;
		s.c[0] = rgb[0];
		s.c[1] = rgb[1];
		s.c[2] = rgb[2];

		return s;
	}

	void toRGB(float* rgb) const {
		rgb[0] = c[0];
		rgb[1] = c[1];
		rgb[2] = c[2];
	}

	void toXYZ(float xyz[3]) const { RGBtoXYZ(c, xyz); }

	static RGBSpectrum fromXYZ(const float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
		RGBSpectrum r;
		XYZtoRGB(xyz, r.c);
		return r;
	}

	float y() const {
		const float y_weight[3] = { 0.212671f, 0.715160f, 0.072169f };
		return y_weight[0] * c[0] + y_weight[1] * c[1] + y_weight[2] * c[2];
	}

	static RGBSpectrum fromSampled(const float* lambda, const float* v, int n);
};

// Some inline functions
template <int n_spectrum_samples>
inline CoefficientSpectrum<n_spectrum_samples> pow(const CoefficientSpectrum<n_spectrum_samples>& s, float e) {
	const CoefficientSpectrum<n_spectrum_samples> ret;
	for (int i = 0; i < n_spectrum_samples; ++i) {
		ret.c[i] = std::pow(s.c[i], e);
	}
	return ret;
}

inline RGBSpectrum lerp(float t, const RGBSpectrum& s1, const RGBSpectrum& s2) {
	return (1 - t) * s1 + t * s2;
}

inline SampledSpectrum lerp(float t, const SampledSpectrum& s1, const SampledSpectrum& s2) {
	return (1 - t) * s1 + t * s2;
}

} // namespace civet

#endif // CIVET_SPECTRUM_H
