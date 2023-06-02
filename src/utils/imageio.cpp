#include <utils/imageio.h>

#include <stb/stb_image.h>
#include <core/spectrum.h>

namespace civet {

static void writeImageEXR(const std::string& name, const float* pixels,
		int xRes, int yRes, int totalXRes, int totalYRes,
		int xOffset, int yOffset);
static void writeImageTGA(const std::string& name, const uint8_t* pixels,
		int xRes, int yRes, int totalXRes, int totalYRes,
		int xOffset, int yOffset);
static RGBSpectrum* readImageTGA(const std::string& name, int* w, int* h);
static RGBSpectrum* readImagePNG(const std::string& name, int* w, int* h);
static bool writeImagePFM(const std::string& filename, const float* rgb,
		int xres, int yres);
static RGBSpectrum* readImagePFM(const std::string& filename, int* xres,
		int* yres);

std::unique_ptr<RGBSpectrum[]> readImage(const std::string& name, Point2i* resolution) {
	// TODO: handle other file extensions not supported by stb_image
	return std::unique_ptr<RGBSpectrum[]>(readImagePNG(name, &resolution->x, &resolution->y));
}

static RGBSpectrum* readImagePNG(const std::string& name, int* w, int* h) {
	int width, height, n_channels;
	unsigned char* data = stbi_load(name.c_str(), &width, &height, &n_channels, 3);	// follow pbrt loading only 24-bits from images
	if (!data) {
		fprintf(stderr, "Error::readImagePNG: error reading PNG %s\n", name.c_str());
		return nullptr;
	}
	*w = width;
	*h = height;

	RGBSpectrum *ret = new RGBSpectrum[*w * *h];
	unsigned char *src = data;
	for (unsigned int y = 0; y < height; ++y) {
		for (unsigned int x = 0; x < width; ++x, src += 3) {
			float c[3];
			c[0] = src[0] / 255.f;
			c[1] = src[1] / 255.f;
			c[2] = src[2] / 255.f;
			ret[y * *w + x] = RGBSpectrum::fromRGB(c);
		}
	}

	free(data);

	return ret;
}

} // namespace civet