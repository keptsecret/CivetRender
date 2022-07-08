#include <utils/parallel.h>

namespace civet {

static std::vector<std::thread> threads;
static bool shutdown_threads = false;
class ParallelForLoop;
static ParallelForLoop* workList = nullptr;
static std::mutex work_list_mutex;

// Bookkeeping variables to help with the implementation of
// mergeWorkerThreadStats().
static std::atomic<bool> report_worker_stats{ false };
// Number of workers that still need to report their stats.
static std::atomic<int> reporter_count;
// After kicking the workers to report their stats, the main thread waits
// on this condition variable until they've all done so.
static std::condition_variable report_done_condition;
static std::mutex report_done_mutex;

class ParallelForLoop {
public:
	ParallelForLoop(std::function<void(int64_t)> f1d, int64_t max_i, int chk_sz) :
			func1D(std::move(f1d)),
			max_index(max_i),
			chunk_size(chk_sz) {}

	ParallelForLoop(const std::function<void(Point2i)>& f, const Point2i& count) :
			func2D(f),
			max_index(count.x * count.y),
			chunk_size(1) {
		nx = count.x;
	}

public:
	std::function<void(int64_t)> func1D;
	std::function<void(Point2i)> func2D;
	const int64_t max_index;
	const int chunk_size;
	int64_t next_index = 0;
	int active_workers = 0;
	ParallelForLoop* next = nullptr;
	int nx = -1;

	bool finished() const {
		return next_index >= max_index && active_workers == 0;
	}
};

void Barrier::wait() {
	std::unique_lock<std::mutex> lock(mutex);
	if (--count == 0) {
		// This is the last thread to reach the barrier; wake up all of the
		// other ones before exiting.
		cv.notify_all();
	} else {
		// Otherwise there are still threads that haven't reached it. Give
		// up the lock and wait to be notified.
		cv.wait(lock, [this] { return count == 0; });
	}
}

static std::condition_variable work_list_condition;

static void workerThreadFunc(int t_idx, std::shared_ptr<Barrier> barrier) {
	printf("Started execution in worker thread %d", t_idx);
	thread_index = t_idx;

	// The main thread sets up a barrier so that it can be sure that all
	// workers have called ProfilerWorkerThreadInit() before it continues
	// (and actually starts the profiling system).
	barrier->wait();

	// Release our reference to the Barrier so that it's freed once all of
	// the threads have cleared it.
	barrier.reset();

	std::unique_lock<std::mutex> lock(work_list_mutex);
	while (!shutdown_threads) {
		if (report_worker_stats) {
			if (--reporter_count == 0) {
				// Once all worker threads have merged their stats, wake up
				// the main thread.
				report_done_condition.notify_one();
			}
			// Now sleep again.
			work_list_condition.wait(lock);
		} else if (!workList) {
			// Sleep until there are more tasks to run
			work_list_condition.wait(lock);
		} else {
			// Get work from _workList_ and run loop iterations
			ParallelForLoop& loop = *workList;

			// Run a chunk of loop iterations for _loop_

			// Find the set of loop iterations to run next
			int64_t indexStart = loop.next_index;
			int64_t indexEnd =
					std::min(indexStart + loop.chunk_size, loop.max_index);

			// Update _loop_ to reflect iterations this thread will run
			loop.next_index = indexEnd;
			if (loop.next_index == loop.max_index) {
				workList = loop.next;
			}
			loop.active_workers++;

			// Run loop indices in _[indexStart, indexEnd)_
			lock.unlock();
			for (int64_t index = indexStart; index < indexEnd; ++index) {
				if (loop.func1D) {
					loop.func1D(index);
				}
				// Handle other types of loops
				else {
					loop.func2D(Point2i(index % loop.nx, index / loop.nx));
				}
			}
			lock.lock();

			// Update _loop_ to reflect completion of iterations
			loop.active_workers--;
			if (loop.finished()) {
				work_list_condition.notify_all();
			}
		}
	}
	printf("Exiting worker thread %d", t_idx);
}

// Parallel Definitions
void parallelFor(std::function<void(int64_t)> func, int64_t count,
		int chunk_size) {
	// Run iterations immediately if not using threads or if _count_ is small
	if (threads.empty() || count < chunk_size) {
		for (int64_t i = 0; i < count; ++i) {
			func(i);
		}
		return;
	}

	// Create and enqueue _ParallelForLoop_ for this loop
	ParallelForLoop loop(std::move(func), count, chunk_size);
	work_list_mutex.lock();
	loop.next = workList;
	workList = &loop;
	work_list_mutex.unlock();

	// Notify worker threads of work to be done
	std::unique_lock<std::mutex> lock(work_list_mutex);
	work_list_condition.notify_all();

	// Help out with parallel loop iterations in the current thread
	while (!loop.finished()) {
		// Run a chunk of loop iterations for _loop_

		// Find the set of loop iterations to run next
		int64_t indexStart = loop.next_index;
		int64_t indexEnd = std::min(indexStart + loop.chunk_size, loop.max_index);

		// Update _loop_ to reflect iterations this thread will run
		loop.next_index = indexEnd;
		if (loop.next_index == loop.max_index)
			workList = loop.next;
		loop.active_workers++;

		// Run loop indices in _[indexStart, indexEnd)_
		lock.unlock();
		for (int64_t index = indexStart; index < indexEnd; ++index) {
			if (loop.func1D) {
				loop.func1D(index);
			}
			// Handle other types of loops
			else {
				loop.func2D(Point2i(index % loop.nx, index / loop.nx));
			}
		}
		lock.lock();

		// Update _loop_ to reflect completion of iterations
		loop.active_workers--;
	}
}

CIVET_THREAD_LOCAL int thread_index;

int maxThreadIndex() {
	// TODO: get option for nthreads from some config
	//return nThreads == 0 ? numSystemCores() : nThreads;
	return numSystemCores();
}

void parallelFor2D(std::function<void(Point2i)> func, const Point2i& count) {
	if (threads.empty() || count.x * count.y <= 1) {
		for (int y = 0; y < count.y; ++y)
			for (int x = 0; x < count.x; ++x)
				func(Point2i(x, y));
		return;
	}

	ParallelForLoop loop(std::move(func), count);
	{
		std::lock_guard<std::mutex> lock(work_list_mutex);
		loop.next = workList;
		workList = &loop;
	}

	std::unique_lock<std::mutex> lock(work_list_mutex);
	work_list_condition.notify_all();

	// Help out with parallel loop iterations in the current thread
	while (!loop.finished()) {
		// Run a chunk of loop iterations for _loop_

		// Find the set of loop iterations to run next
		int64_t idx_start = loop.next_index;
		int64_t idx_end = std::min(idx_start + loop.chunk_size, loop.max_index);

		// Update _loop_ to reflect iterations this thread will run
		loop.next_index = idx_end;
		if (loop.next_index == loop.max_index)
			workList = loop.next;
		loop.active_workers++;

		// Run loop indices in _[idx_start, idx_end)_
		lock.unlock();
		for (int64_t index = idx_start; index < idx_end; ++index) {
			if (loop.func1D) {
				loop.func1D(index);
			}
			// Handle other types of loops
			else {
				loop.func2D(Point2i(index % loop.nx, index / loop.nx));
			}
		}
		lock.lock();

		// Update _loop_ to reflect completion of iterations
		loop.active_workers--;
	}
}

int numSystemCores() {
	return std::max(1u, std::thread::hardware_concurrency());
}

void parallelInit() {
	int nThreads = maxThreadIndex();
	thread_index = 0;

	// Create a barrier so that we can be sure all worker threads get past
	// their call to ProfilerWorkerThreadInit() before we return from this
	// function.  In turn, we can be sure that the profiling system isn't
	// started until after all worker threads have done that.
	std::shared_ptr<Barrier> barrier = std::make_shared<Barrier>(nThreads);

	// Launch one fewer worker thread than the total number we want doing
	// work, since the main thread helps out, too.
	for (int i = 0; i < nThreads - 1; ++i) {
		threads.push_back(std::thread(workerThreadFunc, i + 1, barrier));
	}

	barrier->wait();
}

void parallelCleanup() {
	if (threads.empty()) {
		return;
	}

	{
		std::lock_guard<std::mutex> lock(work_list_mutex);
		shutdown_threads = true;
		work_list_condition.notify_all();
	}

	for (std::thread& thread : threads) {
		thread.join();
	}
	threads.erase(threads.begin(), threads.end());
	shutdown_threads = false;
}

void mergeWorkerThreadStats() {
	std::unique_lock<std::mutex> lock(work_list_mutex);
	std::unique_lock<std::mutex> doneLock(report_done_mutex);
	// Set up state so that the worker threads will know that we would like
	// them to report their thread-specific stats when they wake up.
	report_worker_stats = true;
	reporter_count = threads.size();

	// Wake up the worker threads.
	work_list_condition.notify_all();

	// wait for all of them to merge their stats.
	report_done_condition.wait(lock, []() { return reporter_count == 0; });

	report_worker_stats = false;
}

} // namespace civet