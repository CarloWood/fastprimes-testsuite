#include <array>
#include <vector>
#include <cmath>
#include <cstdint>
#include <iostream>

std::array<unsigned char, 16> inv = { 1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59 };
std::array<unsigned char, 17> first_primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59 };

std::vector<bool> get_primes(unsigned long int limit)
{
  std::vector<bool> primes(limit);
  unsigned long int n;
  bool is_prime = false;

  //First primes
  for (char i = 0; i < first_primes.size(); ++i)
    primes[first_primes[i]] = true;

  for (unsigned long int w = 1; w <= std::ceil(limit / 60) - 1; ++w)
  {
    for (char s = 0; s < inv.size(); ++s)
    {
      n = 60 * w + inv[s];
      is_prime = false;

      unsigned long limit_x;
      unsigned long limit_y;

      // Algorithm 3.1.
      if ((inv[s] % 4 ) == 1)
      {
        limit_x = std::ceil(std::sqrt(n) / 2);
        limit_y = std::ceil(std::sqrt(n));
        for (unsigned long x = 1; x <= limit_x; ++x)
          for (unsigned long y = 1; y <= limit_y; ++y)
            if (4 * x * x + y * y == n)
              is_prime = !is_prime;
      }
      // Algorithm 3.2.
      else if ((inv[s] % 6) == 1)
      {
        limit_x = std::ceil(std::sqrt(n) / 3);
        limit_y = std::ceil(std::sqrt(n));
        for (unsigned long x = 1; x <= limit_x; ++x)
          for (unsigned long y = 1; y <= limit_y; ++y)
            if (3 * x * x + y * y == n)
              is_prime = !is_prime;
      }
      // Algorithm 3.3.
      else if ((inv[s] % 12 ) == 11)
      {
        limit_x = n;
        for (unsigned long x = 2; x <= limit_x; ++x)
          for (unsigned long y = 1; y <= x - 1; ++y)
            if (3 * x * x - y * y == n)
              is_prime = !is_prime;
      }

      // Square Free.
      for (int q = 2; q < std::sqrt(n) && is_prime; ++q)
        if ((primes[q]) && (n % (q * q) == 0))
          is_prime = false;

      if (is_prime)
        primes[n] = true;
    }
  }

  return primes;
}

int main()
{
  auto sieve = get_primes(1000);
  std::vector<uint32_t> primes;
  for (int n = 0; n < sieve.size(); ++n)
    if (sieve[n])
      primes.push_back(n);
  std::cout << "The " << primes.size() << "-th prime is " << primes.back() << std::endl;
}
