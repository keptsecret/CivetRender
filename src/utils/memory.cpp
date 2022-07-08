#include "memory.h"

namespace civet {

void* allocAligned(size_t size) {
#if defined(CIVET_HAVE__ALIGNED_MALLOC)
	return _aligned_malloc(size, CIVET_L1_CACHE_LINE_SIZE);
#elif defined(CiVET_HAVE_POSIX_MEMALIGN)
	void* ptr;
	if (posix_memalign(&ptr, CIVET_L1_CACHE_LINE_SIZE, size) != 0)
		ptr = nullptr;
	return ptr;
#else
	return memalign(CIVET_L1_CACHE_LINE_SIZE, size);
#endif
}

void freeAligned(void* ptr) {
	if (!ptr)
		return;
#if defined(PBRT_HAVE__ALIGNED_MALLOC)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

} // namespace civet