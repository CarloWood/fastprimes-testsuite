#include "sys.h"
#include "utils/Primes.h"
#include <iostream>

using integer_t = utils::Primes::integer_t;
using prime_t = utils::Primes::prime_t;

int main()
{
  integer_t n = 1000000000; // 500000000000;
  utils::Primes generator(n);

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
      ++cnt;
    }
  }
  catch (std::out_of_range const&)
  {
  }

  std::cout << "The " << cnt << "-th prime is " << last_prime << std::endl;
}
