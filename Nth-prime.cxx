#include "sys.h"
#include "utils/ctz.h"
#include "fastprimes/Primes.h"
#include "threadpool/AIThreadPool.h"
#include "cwds/benchmark.h"
#include <fstream>
#include <iostream>

using integer_t = fastprimes::integer_t;
using prime_t = fastprimes::prime_t;

std::vector<uint32_t> debug_primes;
void debug_init_primes()
{
  std::ifstream ifs("primes_till_4000000000", std::ios::binary);
  if (!ifs.is_open()) {
    DoutFatal(dc::fatal, "Failed to open file primes_till_4000000000 for reading.");
  }

  // Read the size of the vector first.
  size_t size;
  ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

  // Resize the vector and read the data.
  debug_primes.resize(size);
  ifs.read(reinterpret_cast<char*>(debug_primes.data()), size * sizeof(uint32_t));

  ifs.close();
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  debug_init_primes();

  constexpr int capacity = 90;

  AIThreadPool thread_pool(30);
  AIQueueHandle handler = thread_pool.new_queue(capacity);

  integer_t n = 1000000000UL; // 1000000000000UL;
  fastprimes::Primes generator(n, handler);

  benchmark::Stopwatch stopwatch(0);
  double const cpu_frequency = 3612059050.0;        // In cycles per second.

  std::cout << "Starting reading out sieve..." << std::endl;
  stopwatch.start();

  uint64_t cnt = 0;
  prime_t last_prime;
  try
  {
    for (;;)
    {
      auto p = generator.next_prime();
      if (p > n)
        break;
      last_prime = p;
      if (!(cnt > uint64_t{189961811} || p == debug_primes[cnt]))
        throw std::runtime_error("next_prime() returned the wrong prime!");
      ++cnt;
    }
  }
  catch (std::out_of_range const&)
  {
  }

  stopwatch.stop();
  uint64_t cycles = stopwatch.diff_cycles() - benchmark::Stopwatch::s_stopwatch_overhead;
  float delta = cycles / cpu_frequency;
  std::cout << "Time spent calling next_prime: " << delta << " seconds." << std::endl;

  std::cout << "The " << cnt << "-th prime is " << last_prime << std::endl;
}
