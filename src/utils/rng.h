#ifndef CIVET_RNG_H
#define CIVET_RNG_H

#include <core/civet.h>

namespace civet {

static const float OneMinusEpsilon = 0x1.fffffep-1;

#define PCG32_DEFAULT_STATE 0x853c49e6748fea9bULL
#define PCG32_DEFAULT_STREAM 0xda3e39cb94b95bdbULL
#define PCG32_MULT 0x5851f42d4c957f2dULL
class RNG {
public:
	RNG();
	RNG(uint64_t sequenceIndex) { setSequence(sequenceIndex); }
	void setSequence(uint64_t sequenceIndex);
	uint32_t uniformUInt32();
	uint32_t uniformUInt32(uint32_t b) {
		uint32_t threshold = (~b + 1u) % b;
		while (true) {
			uint32_t r = uniformUInt32();
			if (r >= threshold) {
				return r % b;
			}
		}
	}

	float uniformFloat() {
		return std::min(OneMinusEpsilon, float(uniformUInt32() * 2.3283064365386963e-10f));
	}

	template <typename Iterator>
	void shuffle(Iterator begin, Iterator end) {
		for (Iterator it = end - 1; it > begin; --it) {
			std::iter_swap(it, begin + uniformUInt32((uint32_t)(it - begin + 1)));
		}
	}

	void advance(int64_t idelta) {
		uint64_t cur_mult = PCG32_MULT, cur_plus = inc, acc_mult = 1u,
				 acc_plus = 0u, delta = (uint64_t)idelta;
		while (delta > 0) {
			if (delta & 1) {
				acc_mult *= cur_mult;
				acc_plus = acc_plus * cur_mult + cur_plus;
			}
			cur_plus = (cur_mult + 1) * cur_plus;
			cur_mult *= cur_mult;
			delta /= 2;
		}
		state = acc_mult * state + acc_plus;
	}

	int64_t operator-(const RNG& other) const {
		uint64_t cur_mult = PCG32_MULT, cur_plus = inc, cur_state = other.state,
				 the_bit = 1u, distance = 0u;
		while (state != cur_state) {
			if ((state & the_bit) != (cur_state & the_bit)) {
				cur_state = cur_state * cur_mult + cur_plus;
				distance |= the_bit;
			}
			the_bit <<= 1;
			cur_plus = (cur_mult + 1ULL) * cur_plus;
			cur_mult *= cur_mult;
		}
		return (int64_t)distance;
	}

private:
	uint64_t state, inc;
};

inline RNG::RNG() :
		state(PCG32_DEFAULT_STATE), inc(PCG32_DEFAULT_STREAM) {}

inline void RNG::setSequence(uint64_t sequenceIndex) {
	state = 0u;
	inc = (sequenceIndex << 1u) | 1u;
	uniformUInt32();
	state += PCG32_DEFAULT_STATE;
	uniformUInt32();
}

inline uint32_t RNG::uniformUInt32() {
	uint64_t oldstate = state;
	state = oldstate * PCG32_MULT + inc;
	uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
	uint32_t rot = (uint32_t)(oldstate >> 59u);
	return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

} // namespace civet

#endif // CIVET_RNG_H
