/*
 * Implementation adapted (largely taken) from: https://github.com/mmp/pbrt-v3/blob/master/src/core/parallel.h
 */

#ifndef CIVET_PARALLEL_H
#define CIVET_PARALLEL_H

#include <core/geometry/vecmath.h>
#include <core/civet.h>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>

namespace civet {

class AtomicFloat {
public:
	// AtomicFloat Public Methods
	explicit AtomicFloat(float v = 0) { bits = floatToBits(v); }

	operator float() const { return bitsToFloat(bits); }

	float operator=(float v) {
		bits = floatToBits(v);
		return v;
	}

	void add(float v) {
		uint32_t old_bits = bits, new_bits;
		do {
			new_bits = floatToBits(bitsToFloat(old_bits) + v);
		} while (!bits.compare_exchange_weak(old_bits, new_bits));
	}

private:
	std::atomic<uint32_t> bits;
};

// Simple one-use barrier; ensures that multiple threads all reach a
// particular point of execution before allowing any of them to proceed
// past it.
//
// Note: this should be heap allocated and managed with a shared_ptr, where
// all threads that use it are passed the shared_ptr. This ensures that
// memory for the Barrier won't be freed until all threads have
// successfully cleared it.
class Barrier {
public:
	Barrier(int count) :
			count(count) {}

	~Barrier() {}

	void wait();

private:
	std::mutex mutex;
	std::condition_variable cv;
	int count;
};

void parallelFor(std::function<void(int64_t)> func, int64_t count,
		int chunk_size = 1);
extern CIVET_THREAD_LOCAL int thread_index;
void parallelFor2D(std::function<void(Point2i)> func, const Point2i& count);
int maxThreadIndex(int num_threads);
int numSystemCores();

void parallelInit(int num_threads = 0);
void parallelCleanup();
void mergeWorkerThreadStats();

} // namespace civet

#endif // CIVET_PARALLEL_H
