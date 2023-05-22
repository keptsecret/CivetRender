#ifndef CIVET_IMAGEIO_H
#define CIVET_IMAGEIO_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <cctype>

namespace civet {

std::unique_ptr<RGBSpectrum[]> readImage(const std::string& name, Point2i* resolution);
RGBSpectrum* readImageEXR(const std::string& name, int* width,
		int* height, Bounds2i* dataWindow = nullptr,
		Bounds2i* displayWindow = nullptr);

void writeImage(const std::string& name, const float* rgb,
		const Bounds2i& outputBounds, const Point2i& totalResolution);

} // namespace civet

#endif // CIVET_IMAGEIO_H
