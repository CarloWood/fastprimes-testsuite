#include "sys.h"
#include "fastprimes/Primes.h"
#include "threadpool/AIThreadPool.h"
#include <iostream>

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  constexpr int num_threads = 8;
  constexpr int capacity = num_threads + 2;

  AIThreadPool thread_pool(num_threads);
  AIQueueHandle handler = thread_pool.new_queue(capacity);

  fastprimes::integer_t n = 100000000UL;
  fastprimes::Primes generator(n, num_threads, handler);

  uint64_t cnt = 0;
  uint64_t sum = 0;
  uint64_t N = std::atoi(argv[1]);
  uint64_t last_prime;

  try
  {
    while (cnt < N)
    {
      auto p = generator.next_prime();
      last_prime = p;
      ++cnt;
      sum += p;
    }
  }
  catch (std::out_of_range const&)
  {
  }

  std::cout << "sum of first " << N << " primes numbers is " << sum << std::endl;
  std::cout << "The last added prime was " << last_prime << std::endl;
}
