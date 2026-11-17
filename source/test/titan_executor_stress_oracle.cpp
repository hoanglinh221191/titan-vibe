// CORE-08 oracle: TitanExecutor with several threads submitting batches.
//
// The old executor parked two kinds of waiters on one condition variable
// (done_cv_): a submitter waiting for the executor to become idle and the
// owner of the running batch waiting for its workers. The last worker woke
// only one of them (notify_one), and the end of a batch woke nobody, so a
// second submitter could sleep forever (or take the owner's wake-up).
//
//   A  two threads, each submitting parallel_for and parallel_do batches in a
//      loop; every result is checked; a watchdog fails the run when the
//      batches stop progressing (a hang is the old failure)
//   B  as A, with every 37th batch throwing: the exception reaches its own
//      submitter and the executor keeps working
//   C  a batch submitted from inside a worker of the same executor runs
//      inline on that worker (the old code hit an assert there)
//
// usage: titan_executor_stress_oracle [rounds] [timeout_s]
#include "titan_executor.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

std::atomic<long> g_progress{0};
std::atomic<int> g_failures{0};

void fail(const std::string &what) {
  std::printf("FAIL %s\n", what.c_str());
  std::fflush(stdout);
  g_failures.fetch_add(1);
}

void submitter(TitanExecutor &exec, int id, int rounds, bool with_throw) {
  std::vector<long long> partial(exec.capacity(), 0);
  for (int r = 0; r < rounds; ++r) {
    const std::size_t n = 64 + static_cast<std::size_t>((r * 7 + id * 13) % 97);
    std::fill(partial.begin(), partial.end(), 0);
    const bool throw_now = with_throw && (r % 37 == id);
    bool caught = false;
    try {
      exec.parallel_for(n, exec.capacity(), 1,
                        [&](std::size_t worker, std::size_t begin, std::size_t end) {
                          long long s = 0;
                          for (std::size_t i = begin; i < end; ++i)
                            s += static_cast<long long>(i);
                          partial[worker] += s;
                          if (throw_now && begin == 0)
                            throw std::runtime_error("planned");
                        });
    } catch (const std::runtime_error &) {
      caught = true;
    }
    if (caught != throw_now)
      fail("thread " + std::to_string(id) + " round " + std::to_string(r) +
           ": exception delivered to the wrong submitter");
    if (!throw_now) {
      long long total = 0;
      for (long long v : partial)
        total += v;
      const long long expected =
          static_cast<long long>(n) * static_cast<long long>(n - 1) / 2;
      if (total != expected)
        fail("thread " + std::to_string(id) + " round " + std::to_string(r) +
             ": parallel_for sum " + std::to_string(total) + " != " +
             std::to_string(expected));
    }
    std::atomic<int> hits{0};
    exec.parallel_do(3, [&](std::size_t worker) {
      hits.fetch_add(static_cast<int>(worker) + 1);
    });
    if (hits.load() != 6)
      fail("thread " + std::to_string(id) + " round " + std::to_string(r) +
           ": parallel_do ran " + std::to_string(hits.load()) + " != 6");
    g_progress.fetch_add(1);
  }
}

bool run_two_submitters(TitanExecutor &exec, int rounds, int timeout_s,
                        bool with_throw, const char *label) {
  g_progress.store(0);
  std::thread a(submitter, std::ref(exec), 0, rounds, with_throw);
  std::thread b(submitter, std::ref(exec), 1, rounds, with_throw);
  const long target = 2L * rounds;
  long last = -1;
  auto last_change = std::chrono::steady_clock::now();
  while (g_progress.load() < target) {
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const long now = g_progress.load();
    if (now != last) {
      last = now;
      last_change = std::chrono::steady_clock::now();
    } else if (std::chrono::steady_clock::now() - last_change >
               std::chrono::seconds(timeout_s)) {
      std::printf("FAIL %s: no progress for %d s after %ld of %ld rounds "
                  "(lost wake-up: a submitter never woke)\n",
                  label, timeout_s, now, target);
      std::fflush(stdout);
      std::_Exit(2); // the hung threads cannot be joined
    }
  }
  a.join();
  b.join();
  std::printf("PASS %s: %ld rounds of parallel_for + parallel_do from 2 threads\n",
              label, target);
  return true;
}

} // namespace

int main(int argc, char **argv) {
  const int rounds = argc > 1 ? std::atoi(argv[1]) : 4000;
  const int timeout_s = argc > 2 ? std::atoi(argv[2]) : 20;
  TitanExecutor exec(4);

  run_two_submitters(exec, rounds, timeout_s, false, "A two submitters");
  run_two_submitters(exec, rounds, timeout_s, true,
                     "B two submitters, every 37th batch throws");

  // C: nested submission from a worker of the same executor runs inline
  std::atomic<long long> nested{0};
  exec.parallel_for(8, 4, 1, [&](std::size_t, std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i)
      exec.parallel_for(10, 4, 1, [&](std::size_t, std::size_t b, std::size_t e) {
        for (std::size_t j = b; j < e; ++j)
          nested.fetch_add(static_cast<long long>(j));
      });
  });
  std::atomic<int> nested_do{0};
  exec.parallel_do(2, [&](std::size_t) {
    exec.parallel_do(3, [&](std::size_t w) { nested_do.fetch_add(1 + static_cast<int>(w)); });
  });
  if (nested.load() != 8 * 45 || nested_do.load() != 2)
    fail("C nested batches: sum " + std::to_string(nested.load()) + " (expected 360), do " +
         std::to_string(nested_do.load()) + " (expected 2: one inline call per outer worker)");
  else
    std::printf("PASS C nested batch runs inline on the calling worker\n");

  if (g_failures.load() != 0) {
    std::printf("FAIL %d check(s)\n", g_failures.load());
    return 1;
  }
  std::printf("PASS executor stress\n");
  return 0;
}
