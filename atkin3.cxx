#include "sys.h"
#include "fastprimes/Primes.h"
#include <iostream>
#include <cmath>
#include <map>
#include "debug.h"

using prime_t = fastprimes::prime_t;
using integer_t = fastprimes::integer_t;

fastprimes::Primes gen(1000000);
std::vector<prime_t> primes = gen.make_vector();

std::map<prime_t, int> factorize(integer_t n)
{
  std::map<prime_t, int> result;
  for (prime_t p : primes)
  {
    if (n == 1)
      break;
    while (n % p == 0)
    {
      ++result[p];
      n /= p;
    }
  }
  ASSERT(n == 1 || primes.back() >= std::sqrt(n));
  return result;
}

bool is_square_free(integer_t n)
{
  auto pem = factorize(n);
  bool square_free = true;
  for (auto pe : pem)
    if (pe.second > 1)
    {
      square_free = false;
      break;
    }
  return square_free;
}

void print_factorization(std::map<prime_t, int> const& pem)
{
  char const* sep = "";
  for (auto pe : pem)
  {
    std::cout << sep << pe.first;
    if (pe.second > 1)
      std::cout << '^' << pe.second;
    sep = " * ";
  }
}

// Z[i] = { a + bi | a,b \in Z }
//
// Z[i]^* = { 1, -1, i, -i }
//
// I * I --> I, I * N --> N, I * C --> C, I * D --> D
// N * I --> N, N * N --> I, N * C --> D, N * D --> C
// C * I --> C, C * N --> D, C * C --> N, C * D --> I
// D * I --> D, D * N --> C, D * C --> I, D * D --> N
//
// I={   2^k0 (-3)^l0 (-5i)^m0 }
// N={ - 2^k1 (-3)^l1 (-5i)^m1 }
// C={ i 2^k2 (-3)^l2 (-5i)^m2 }
// D={-i 2^k3 (-3)^l3 (-5i)^m3 }
//
//   2 - 3i
//
// S = { 2^k (-3)^l (-5i)^m ... }
//
// Z[i] = { 0, a, -a, bi, -bi, a + bi, a - bi, -a + bi, -a - bi | a,b \in S }
//
int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // 4x^2 + y^2 = n  -->
  //
  // x <= sqrt(n - 1) / 4;
  //
  for (int n = 1; n < 1000000; n += 4)
  {
    if ((n - 1) % 100 == 0)
      std::cout << n << std::endl;
    //if ((n % 12 != 1 && n % 12 != 5))
    //  continue;
    if (!is_square_free(n))
      continue;

    [[maybe_unused]] int solutions = 0;
    int limit_x = std::sqrt((n - 1) / 4);
    for (int x = 1; x <= limit_x; ++x)
    {
      int limit_y = std::sqrt(n - 4 * x * x);
      for (int y = 1; y <= limit_y; ++y)
      {
        if (4 * x * x + y * y == n)
        {
          //std::cout << "4 * " << x << "^2 + " << y << "^2 = " << n << std::endl;
          ++solutions;
        }
      }
    }
    ASSERT((solutions % 2 == 1) == gen.is_prime(n));
  }
}
