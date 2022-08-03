#ifndef CIVET_INTERPOLATION_H
#define CIVET_INTERPOLATION_H

#include <core/civet.h>

namespace civet {

bool catmullRomWeights(int size, const float* nodes, float x, int* offset, float* weights);

/**
 * Evaluates a weighted sum of cosines
 * @param a array of coefficients ak
 * @param m maximum order
 * @param cos_phi cosine of the angle
 * @return weighted sum of cosines
 */
float fourier(const float* a, int m, double cos_phi);

} // namespace civet

#endif // CIVET_INTERPOLATION_H
