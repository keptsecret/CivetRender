#ifndef CIVET_MIPMAP_H
#define CIVET_MIPMAP_H

#include <core/civet.h>
#include <core/spectrum.h>
#include <core/texture.h>
#include <utils/parallel.h>

namespace civet {

enum class ImageWrap {
	Repeat,
	Black,
	Clamp
};

struct ResampleWeight {
	int first_texel;
	float weight[4];
};

template <typename T>
class MIPMap {
public:
	MIPMap(const Point2i& res, const T* img, bool do_trilinear = false, float max_aniso = 8.f, ImageWrap wrap_mode = ImageWrap::Repeat);

	int width() const { return resolution[0]; }
	int height() const { return resolution[1]; }
	int levels() const { return pyramid.size(); }

	const T& texel(int level, int s, int t) const;

	T lookup(const Point2f& st, float width = 0) const;
	T lookup(const Point2f& st, Vector2f dstdx, Vector2f dstdy) const;

private:
	std::unique_ptr<ResampleWeight[]> resampleWeights(int old_res, int new_res) {
		std::unique_ptr<ResampleWeight[]> wt(new ResampleWeight[new_res]);
		float filter_width = 0.2f;
		for (int i = 0; i < new_res; i++) {
			float center = (i + 0.5f) * old_res / new_res;
			wt[i].first_texel = std::floor((center - filter_width) + 0.5f);
			for (int j = 0; j < 4; j++) {
				float pos = wt[i].first_texel + j + 0.5f;
				wt[i].weight[j] = lanczos((pos - center) / filter_width);
			}
			float inv_wt_sum = 1.f / (wt[i].weight[0] + wt[i].weight[1] + wt[i].weight[2] + wt[i].weight[3]);
			for (int j = 0; j < 4; j++) {
				wt[i].weight[j] *= inv_wt_sum;
			}
		}
		return wt;
	}

	float clamp(float v) { return civet::clamp(v, 0.f, Infinity); }
	RGBSpectrum clamp(const RGBSpectrum& v) { return v.clamp(0.f, Infinity); }
	SampledSpectrum clamp(const SampledSpectrum& v) { return v.clamp(0.f, Infinity); }

	T triangle(int level, const Point2f& st) const;
	T EWA(int level, Point2f st, Vector2f dst0, Vector2f dst1) const;

	const bool do_trilinear;
	const float max_anisotropy;
	const ImageWrap wrap_mode;
	Point2i resolution;

	std::vector<std::unique_ptr<BlockedArray<T>>> pyramid;
	static CIVET_CONSTEXPR int WeightLUTSize = 128;
	static float weight_Lut[WeightLUTSize];
};

template <typename T>
MIPMap<T>::MIPMap(const Point2i& res, const T* img, bool do_trilinear, float max_aniso, ImageWrap wrap_mode) :
		do_trilinear(do_trilinear), max_anisotropy(max_aniso), wrap_mode(wrap_mode), resolution(res) {
	std::unique_ptr<T[]> resampled_image = nullptr;
	if (!isPowerOf2(resolution[0]) || !isPowerOf2(resolution[1])) {
		Point2i res_pow2(roundUpPow2(resolution[0]), roundUpPow2(resolution[1]));

		// resample in s direction
		auto s_weights = resampleWeights(resolution[0], res_pow2[0]);
		resampled_image.reset(new T[res_pow2[0] * res_pow2[1]]);
		parallelFor(
				[&](int t) {
					for (int s = 0; s < res_pow2[0]; s++) {
						resampled_image[t * res_pow2[0] + s] = 0.f;
						for (int j = 0; j < 4; j++) {
							int orig_s = s_weights[s].first_texel + j;
							if (wrap_mode == ImageWrap::Repeat) {
								orig_s = mod(orig_s, resolution[0]);
							} else if (wrap_mode == ImageWrap::Clamp) {
								orig_s = civet::clamp(orig_s, 0, resolution[0] - 1);
							}
							if (orig_s >= 0 && orig_s < (int)resolution[0]) {
								resampled_image[t * res_pow2[0] + s] += s_weights[s].weight[j] * img[t * resolution[0] + orig_s];
							}
						}
					}
				},
				resolution[1], 16);

		// resample in t direction
		auto t_weights = resampleWeights(resolution[1], res_pow2[1]);
		std::vector<T*> resample_bufs;
		int n_threads = maxThreadIndex(numSystemCores());
		for (int i = 0; i < n_threads; ++i) {
			resample_bufs.push_back(new T[res_pow2[1]]);
		}
		parallelFor(
				[&](int s) {
					T* work_data = resample_bufs[thread_index];
					for (int t = 0; t < res_pow2[1]; ++t) {
						work_data[t] = 0.f;
						for (int j = 0; j < 4; ++j) {
							int offset = t_weights[t].first_texel + j;
							if (wrap_mode == ImageWrap::Repeat) {
								offset = mod(offset, resolution[1]);
							} else if (wrap_mode == ImageWrap::Clamp) {
								offset = civet::clamp(offset, 0, (int)resolution[1] - 1);
							}
							if (offset >= 0 && offset < (int)resolution[1]) {
								work_data[t] += t_weights[t].weight[j] * resampled_image[offset * res_pow2[0] + s];
							}
						}
					}
					for (int t = 0; t < res_pow2[1]; ++t) {
						resampled_image[t * res_pow2[0] + s] = clamp(work_data[t]);
					}
				},
				res_pow2[0], 32);

		for (auto ptr : resample_bufs) {
			delete[] ptr;
		}

		resolution = res_pow2;

		// initialize levels of mipmap
		int n_levels = 1 + log2Int(std::max(resolution[0], resolution[1]));
		pyramid.resize(n_levels);

		// most detailed level
		pyramid[0].reset(new BlockedArray<T>(resolution[0], resolution[1], resampled_image ? resampled_image.get() : img));

		for (int i = 1; i < n_levels; i++) {
			int sres = std::max(1, pyramid[i - 1]->uSize() / 2);
			int tres = std::max(1, pyramid[i - 1]->vSize() / 2);
			pyramid[i].reset(new BlockedArray<T>(sres, tres));

			parallelFor(
					[&](int t) {
						for (int s = 0; s < sres; ++s) {
							(*pyramid[i])(s, t) = .25f *
									(texel(i - 1, 2 * s, 2 * t) + texel(i - 1, 2 * s + 1, 2 * t) +
											texel(i - 1, 2 * s, 2 * t + 1) + texel(i - 1, 2 * s + 1, 2 * t + 1));
						}
					},
					tres, 16);
		}
	}
}

template <typename T>
const T& MIPMap<T>::texel(int level, int s, int t) const {
	const BlockedArray<T>& l = *pyramid[level];
	switch (wrap_mode) {
		case ImageWrap::Repeat:
			s = mod(s, l.uSize());
			t = mod(t, l.vSize());
			break;
		case ImageWrap::Clamp:
			s = civet::clamp(s, 0, l.uSize() - 1);
			t = civet::clamp(t, 0, l.vSize() - 1);
			break;
		case ImageWrap::Black: {
			static const T black = 0.f;
			if (s < 0 || s >= (int)l.uSize() || t < 0 || t >= (int)l.vSize()) {
				return black;
			}
			break;
		}
	}
	return l(s, t);
}

template <typename T>
T MIPMap<T>::lookup(const Point2f& st, float width) const {
	float level = levels() - 1 + log2(std::fmax(width, 1e-8));
	if (level < 0) {
		return triangle(0, st);
	} else if (level >= levels() - 1) {
		return texel(levels() - 1, 0, 0);
	} else {
		int iLevel = std::floor(level);
		float delta = level - iLevel;
		return lerp(delta, triangle(iLevel, st), triangle(iLevel + 1, st));
	}
}

template <typename T>
T MIPMap<T>::triangle(int level, const Point2f& st) const {
	level = civet::clamp(level, 0, levels() - 1);
	float s = st[0] * pyramid[level]->uSize() - 0.5f;
	float t = st[1] * pyramid[level]->vSize() - 0.5f;
	int s0 = std::floor(s), t0 = std::floor(t);
	float ds = s - s0, dt = t - t0;
	return (1 - ds) * (1 - dt) * texel(level, s0, t0) +
			(1 - ds) * dt * texel(level, s0, t0 + 1) +
			ds * (1 - dt) * texel(level, s0 + 1, t0) +
			ds * dt * texel(level, s0 + 1, t0 + 1);
}

template <typename T>
T MIPMap<T>::lookup(const Point2f& st, Vector2f dst0, Vector2f dst1) const {
	if (do_trilinear) {
		float width = std::max(std::max(std::abs(dst0[0]), std::abs(dst0[1])),
				std::max(std::abs(dst1[0]), std::abs(dst1[1])));
		return lookup(st, width);
	}

	// Compute ellipse minor and major axes
	if (dst0.lengthSquared() < dst1.lengthSquared()) {
		std::swap(dst0, dst1);
	}
	float major_length = dst0.length();
	float minor_length = dst1.length();

	// Clamp ellipse eccentricity if too large
	if (minor_length * max_anisotropy < major_length && minor_length > 0) {
		float scale = major_length / (minor_length * max_anisotropy);
		dst1 *= scale;
		minor_length *= scale;
	}
	if (minor_length == 0) {
		return triangle(0, st);
	}

	// Choose level of detail for EWA lookup and perform EWA filtering
	float lod = std::max((float)0, levels() - 1.f + log2(minor_length));
	int ilod = std::floor(lod);
	return lerp(lod - ilod, EWA(ilod, st, dst0, dst1), EWA(ilod + 1, st, dst0, dst1));
}

template <typename T>
T MIPMap<T>::EWA(int level, Point2f st, Vector2f dst0, Vector2f dst1) const {
	if (level >= levels()) {
		return texel(levels() - 1, 0, 0);
	}
	// Convert EWA coordinates to appropriate scale for level
	st[0] = st[0] * pyramid[level]->uSize() - 0.5f;
	st[1] = st[1] * pyramid[level]->vSize() - 0.5f;
	dst0[0] *= pyramid[level]->uSize();
	dst0[1] *= pyramid[level]->vSize();
	dst1[0] *= pyramid[level]->uSize();
	dst1[1] *= pyramid[level]->vSize();

	// Compute ellipse coefficients to bound EWA filter region
	float A = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1;
	float B = -2 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
	float C = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1;
	float invF = 1 / (A * C - B * B * 0.25f);
	A *= invF;
	B *= invF;
	C *= invF;

	// Compute the ellipse's $(s,t)$ bounding box in texture space
	float det = -B * B + 4 * A * C;
	float inv_det = 1 / det;
	float usqrt = std::sqrt(det * C), vsqrt = std::sqrt(A * det);
	int s0 = std::ceil(st[0] - 2 * inv_det * usqrt);
	int s1 = std::floor(st[0] + 2 * inv_det * usqrt);
	int t0 = std::ceil(st[1] - 2 * inv_det * vsqrt);
	int t1 = std::floor(st[1] + 2 * inv_det * vsqrt);

	// Scan over ellipse bound and compute quadratic equation
	T sum(0.f);
	float sum_wts = 0;
	for (int it = t0; it <= t1; ++it) {
		float tt = it - st[1];
		for (int is = s0; is <= s1; ++is) {
			float ss = is - st[0];
			// Compute squared radius and filter texel if inside ellipse
			float r2 = A * ss * ss + B * ss * tt + C * tt * tt;
			if (r2 < 1) {
				int index =
						std::min((int)(r2 * WeightLUTSize), WeightLUTSize - 1);
				float weight = weight_Lut[index];
				sum += texel(level, is, it) * weight;
				sum_wts += weight;
			}
		}
	}
	return sum / sum_wts;
}

template <typename T>
float MIPMap<T>::weight_Lut[WeightLUTSize];

} // namespace civet

#endif // CIVET_MIPMAP_H
