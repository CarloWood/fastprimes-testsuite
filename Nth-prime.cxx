#include "sys.h"
#include "fastprimes/Primes.h"
#include <fstream>
#include <iostream>

using integer_t = fastprimes::Primes::integer_t;
using prime_t = fastprimes::Primes::prime_t;

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
  debug_init_primes();

  integer_t n = 1000000000UL; // 500000000000;
  fastprimes::Primes generator(n);

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

  std::cout << "The " << cnt << "-th prime is " << last_prime << std::endl;
}
