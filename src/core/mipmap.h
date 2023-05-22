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
	MIPMap(const Point2i& res, const T* img, bool do_trilinear, float max_aniso, ImageWrap wrap_mode);

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

	T triangle(int level, const Point2f &st) const;
	T EWA(int level, Point2f st, Vector2f dst0, Vector2f dst1) const;

	const bool do_trilinear;
	const float max_anisotropy;
	const ImageWrap wrap_mode;
	Point2i resolution;

	std::vector<std::unique_ptr<BlockedArray<T>>> pyramid;
	static CIVET_CONSTEXPR int WeightLUTSize = 128;
	static float weight_Lut[WeightLUTSize];
};

} // namespace civet

#endif // CIVET_MIPMAP_H
