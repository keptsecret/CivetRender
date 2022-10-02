/*
 * Implementation adapted (largely taken) from: https://github.com/mmp/pbrt-v3/blob/master/src/core/memory.h
 */

#ifndef CIVET_MEMORY_H
#define CIVET_MEMORY_H

#include "core/civet.h"
#include <stddef.h>
#include <list>

namespace civet {

// Memory Declarations
#define ARENA_ALLOC(arena, Type) new ((arena).alloc(sizeof(Type))) Type
void* allocAligned(size_t size);
template <typename T>
T* allocAligned(size_t size) {
	return (T*)allocAligned(size * sizeof(T));
}

void freeAligned(void*);

class
#ifdef CIVET_HAVE_ALIGNAS
		alignas(CIVET_L1_CACHE_LINE_SIZE)
#endif // CIVET_HAVE_ALIGNAS
		MemoryArena {
public:
	MemoryArena(size_t blockSize = 262144) :
			block_size(blockSize) {}

	~MemoryArena() {
		freeAligned(current_block);
		for (auto& block : used_blocks) {
			freeAligned(block.second);
		}
		for (auto& block : available_blocks) {
			freeAligned(block.second);
		}
	}

	void* alloc(size_t n_bytes) {
		// Round up _nBytes_ to minimum machine alignment
#if __GNUC__ == 4 && __GNUC_MINOR__ < 9
		// gcc bug: max_align_t wasn't in std:: until 4.9.0
		const int align = alignof(::max_align_t);
#elif !defined(CIVET_HAVE_ALIGNOF)
		const int align = 16;
#else
		const int align = alignof(std::max_align_t);
#endif
#ifdef CIVET_HAVE_CONSTEXPR
		static_assert(isPowerOf2(align), "Minimum alignment not a power of two");
#endif
		n_bytes = (n_bytes + align - 1) & ~(align - 1);
		if (current_block_pos + n_bytes > current_alloc_size) {
			// Add current block to _usedBlocks_ list
			if (current_block) {
				used_blocks.push_back(std::make_pair(current_alloc_size, current_block));
				current_block = nullptr;
				current_alloc_size = 0;
			}

			// Get new block of memory for _MemoryArena_

			// Try to get memory block from _availableBlocks_
			for (auto iter = available_blocks.begin();
					iter != available_blocks.end(); ++iter) {
				if (iter->first >= n_bytes) {
					current_alloc_size = iter->first;
					current_block = iter->second;
					available_blocks.erase(iter);
					break;
				}
			}
			if (!current_block) {
				current_alloc_size = std::max(n_bytes, block_size);
				current_block = allocAligned<uint8_t>(current_alloc_size);
			}
			current_block_pos = 0;
		}
		void* ret = current_block + current_block_pos;
		current_block_pos += n_bytes;
		return ret;
	}

	template <typename T>
	T* alloc(size_t n = 1, bool run_constructor = true) {
		T* ret = (T*)alloc(n * sizeof(T));
		if (run_constructor) {
			for (size_t i = 0; i < n; ++i) {
				new (&ret[i]) T();
			}
		}
		return ret;
	}

	void reset() {
		current_block_pos = 0;
		available_blocks.splice(available_blocks.begin(), used_blocks);
	}

	size_t totalAllocated() const {
		size_t total = current_alloc_size;
		for (const auto& alloc : used_blocks) {
			total += alloc.first;
		}
		for (const auto& alloc : available_blocks) {
			total += alloc.first;
		}
		return total;
	}

private:
	MemoryArena(const MemoryArena&) = delete;
	MemoryArena& operator=(const MemoryArena&) = delete;

	const size_t block_size;
	size_t current_block_pos = 0, current_alloc_size = 0;
	uint8_t* current_block = nullptr;
	std::list<std::pair<size_t, uint8_t*>> used_blocks, available_blocks;
};

template <typename T, int logBlockSize>
class BlockedArray {
public:
	BlockedArray(int _u_res, int _v_res, const T* d = nullptr) :
			u_res(_u_res), v_res(_v_res), u_blocks(roundUp(_u_res) >> logBlockSize) {
		int n_alloc = roundUp(_u_res) * roundUp(_v_res);
		data = allocAligned<T>(n_alloc);
		for (int i = 0; i < n_alloc; ++i) {
			new (&data[i]) T();
		}
		if (d) {
			for (int v = 0; v < _v_res; ++v) {
				for (int u = 0; u < _u_res; ++u) {
					(*this)(u, v) = d[v * _u_res + u];
				}
			}
		}
	}

	constexpr int blockSize() const { return 1 << logBlockSize; }

	int roundUp(int x) const {
		return (x + blockSize() - 1) & ~(blockSize() - 1);
	}

	int uSize() const { return u_res; }
	int vSize() const { return v_res; }

	~BlockedArray() {
		for (int i = 0; i < u_res * v_res; ++i) {
			data[i].~T();
		}
		freeAligned(data);
	}

	int block(int a) const { return a >> logBlockSize; }
	int offset(int a) const { return (a & (blockSize() - 1)); }

	T& operator()(int u, int v) {
		int bu = block(u), bv = block(v);
		int ou = offset(u), ov = offset(v);
		int offset = blockSize() * blockSize() * (u_blocks * bv + bu);
		offset += blockSize() * ov + ou;
		return data[offset];
	}

	const T& operator()(int u, int v) const {
		int bu = block(u), bv = block(v);
		int ou = offset(u), ov = offset(v);
		int offset = blockSize() * blockSize() * (u_blocks * bv + bu);
		offset += blockSize() * ov + ou;
		return data[offset];
	}

	void getLinearArray(T* a) const {
		for (int v = 0; v < v_res; ++v) {
			for (int u = 0; u < u_res; ++u) {
				*a++ = (*this)(u, v);
			}
		}
	}

private:
	T* data;
	const int u_res, v_res, u_blocks;
};

} // namespace civet

#endif // CIVET_MEMORY_H
