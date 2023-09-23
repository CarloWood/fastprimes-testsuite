#include "sys.h"
#include <iostream>
#include <tuple>
#include <array>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cassert>
#include <limits>
#include "debug.h"

using prime_t = uint64_t;
using integer_t = int64_t;
using prime_index_t = int;

int modular_inverse(integer_t n, integer_t m)
{
  integer_t a = n, b = m;
  integer_t x0 = 1, y0 = 0, x1 = 0, y1 = 1;
  while (b > 0)
  {
    integer_t q = a / b;
    std::tie(a, b) = std::make_tuple(b, a % b);
    std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
    std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
  }
  // Ensure the result is positive.
  return (x0 + m) % m;
}

constexpr std::array<prime_t, 7> small_primes = { 2, 3, 5, 7, 11, 13, 17 };

// Returns pₙ# = p₁⋅p₂⋅p₃… pₙ
// Note that p₀ := 1
constexpr int calc_primorial(prime_index_t n)
{
  if (n < 1)
    return 1;
  return small_primes[n - 1] * calc_primorial(n - 1);
}

// Returns Product_{k=1..n} (pₖ-1)
// For example, repeat(4) = (p₁-1)⋅(p₂-1)⋅(p₃-1)⋅(p₄-1) = 1⋅2⋅4⋅6 = 48
constexpr int calc_repeat(prime_index_t n)
{
  if (n < 2)
    return 1;
  return (small_primes[n - 1] - 1) * calc_repeat(n - 1);
}

int const compression = 6;                                      // The number of primes to skip; must be less than 7 or things overflow.
int const compression_primorial = calc_primorial(compression);  // The multiplication of the first `compression` primes.
int const compression_repeat = calc_repeat(compression);        // Same, but subtracting one from each prime first.
int const compression_first_prime = small_primes[compression];  // The first integer that is not divisible by any of the `compression` primes.
int const compression_first_prime_second_row = compression_first_prime + compression_primorial;

constexpr std::array<int, compression_repeat> calc_first_row()
{
  std::array<int, compression_repeat> first_row = {};
  int col = 0;
  int candidate = compression_first_prime;
  do
  {
    bool is_divisible = false;
    for (int i = 1; i < compression; ++i)
      if (candidate % small_primes[i] == 0)
      {
        is_divisible = true;
        break;
      }
    if (!is_divisible)
      first_row[col++] = candidate;
    candidate += 2;
  }
  while (col < compression_repeat);
  return first_row;
}

constexpr std::array<int, compression_repeat> first_row = calc_first_row();

prime_t sieve_index_to_prime(integer_t sieve_index)
{
  integer_t row = sieve_index / compression_repeat;
  int column = sieve_index - row * compression_repeat;
  return row * compression_primorial + first_row[column];
}

constexpr std::array<int, compression_first_prime_second_row> calc_first_row_to_column()
{
  std::array<int, compression_first_prime_second_row> first_row_to_column = {};
  int n = 0;
  for (int col = 0; col <= first_row.size(); ++col)
  {
    while (n < (col == first_row.size() ? compression_first_prime_second_row : first_row[col]))
    {
      first_row_to_column[n] = col - 1;
      ++n;
    }
  }
  return first_row_to_column;
}

constexpr std::array<int, compression_first_prime_second_row> first_row_to_column = calc_first_row_to_column();

integer_t prime_to_sieve_index(uint64_t value)
{
  integer_t row = (value - compression_first_prime) / compression_primorial;
  int column = first_row_to_column[value - row * compression_primorial];
  return row * compression_repeat + column;
}

// Returns an upper bound on the number of primes that are smaller than or equal to n.
// Just a random function that I thought matched quite well (for n > 1000 or so).
// For n larger than 500,000,000 it is off on average 0.0079% (too large, thus).
uint32_t calc_upper_bound_number_of_primes(int n)
{
  double logn = std::log(n);
  assert(logn > 4);
  return std::exp(0.3125 * std::pow(1 / (logn - 4), 1.655) + logn - std::log(logn - 1)) - 4;
}

// Returns all primes less than or equal max_value.
std::vector<prime_t> calculate_primes(uint64_t max_value)
{
  // Pre-allocate the vector that is going to contain the primes.
  uint32_t upper_bound_number_of_primes = calc_upper_bound_number_of_primes(max_value);
  std::vector<prime_t> result;
  result.reserve(upper_bound_number_of_primes);

  // Copy the primes that are compressed away.
  for (int i = 0; i < compression; ++i)
    result.push_back(small_primes[i]);

  // Calculate how many integers are not divisible by the first `compression` number of primes that are less than or equal max_value.
  integer_t sieve_rows = (max_value - compression_first_prime) / compression_primorial + 1;
  // Lets allocate whole rows, even if we only partially need the last one.
  integer_t sieve_size = sieve_rows * compression_repeat;
  // Allocate the sieve and fill it with ones.
  char* sieve = (char*)std::malloc(sieve_size);
  std::memset(sieve, 1, sieve_size);

  integer_t const sqrt_max_value = std::sqrt(max_value);
  integer_t const sqrt_index = prime_to_sieve_index(sqrt_max_value);
  assert(sqrt_max_value >= compression_first_prime_second_row);

  char* next_prime = sieve - 1; // Minus one because we start with an increment.

  for (;;)
  {
    if (*++next_prime)
    {
      prime_t prime = sieve_index_to_prime(next_prime - sieve);

      if (prime >= compression_first_prime_second_row)
      {
        --next_prime;
        break;
      }

      result.push_back(prime);

      integer_t const step = compression_repeat * prime;
      int const offset = compression_first_prime_second_row * prime;
      int compression_primorial_inverse = modular_inverse(compression_primorial, prime);

      for (int col = 0; col < compression_repeat; ++col)
      {
        // The largest value of `offset - first_row[col]` is when `offset` has its largest value
        // and first_row[col] is at its smallest. The latter happens when col = 0 (at which point
        // first_row[col] equals compression_first_prime). The former, `offset` at its largest,
        // happens when `prime` is at its largest, which is `compression_first_prime_second_row`.
        // Note that compression_first_prime_second_row = compression_first_prime + compression_primorial
        //
        // Let P = compression_primorial, F = compression_first_prime.
        // Then compression_first_prime_second_row = P + F.
        //
        // And we can write for the largest values involved:
        //   prime = P + F,
        //   offset = (P + F) * (P + F);
        //   -first_row[col] = -F
        //
        // The largest possible value of compression_primorial_inverse is prime - 1, or P + F - 1.
        //
        // Thus the largest possible value of ((offset - first_row[col]) * compression_primorial_inverse) is (less than or) equal
        //
        //   ((P + F)^2 - F) * (P + F - 1) =
        //   (P^2 + 2 F P + F^2 - F) * (P + F - 1) =
        //   P^3 + 2 F P^2 + F^2 P - F P + F P^2 + 2 F^2 P + F^3 - F^2 - P^2 - 2 F P - F^2 + F =
        //   P^3 + (2F + F - 1) P^2 + (F^2 - F + 2 F^2 - 2F) P + (F^3 - F^2 - F^2 + F) =
        //   P^3 + (3F - 1) P^2 + (3F^2 - 3F) P + (F^3 - 2 F^2 + F) < P^3 + 3 F P^2 + 3 F^2 P + F^3 = (P + F)^3
        //
        // which becomes larger than what fits in an int when P + F > 1290.
        //  compression   P + F
        //  2             2*3 + 5 = 11
        //  3             2*3*5 + 7 = 37
        //  4             2*3*5*7 + 11 = 221
        //  5             2*3*5*7*11 + 13 = 2323
        //
        // This means that if compression_first_prime_second_row is larger than 1290 we need an extra modulo.
        int multiple_index = offset - first_row[col];
        if constexpr (compression_first_prime_second_row > 1290)
          multiple_index %= prime;
        multiple_index *= compression_primorial_inverse;
        multiple_index %= prime;
        for (int i = col + compression_repeat * multiple_index; i < sieve_size; i += step)
          sieve[i] = 0;
      }
    }
  }

  for (;;)
  {
    if (*++next_prime)
    {
      prime_t prime = sieve_index_to_prime(next_prime - sieve);
      assert(prime >= compression_first_prime_second_row);

      if (prime > sqrt_max_value)
        break;

      result.push_back(prime);

      integer_t const step = compression_repeat * prime;
      int compression_primorial_inverse = modular_inverse(compression_primorial, prime);

      for (int col = 0; col < compression_repeat; ++col)
      {
#if 0
        long debug_prime = prime;
        long debug_first_row = first_row[col];
        long debug_inverse = compression_primorial_inverse;
        long debug_multiple_index = prime - first_row[col];
        assert(0 <= debug_multiple_index && debug_multiple_index < prime);
        debug_multiple_index *= debug_inverse;
        assert(debug_multiple_index <= std::numeric_limits<int>::max());
        debug_multiple_index %= debug_prime;
        debug_multiple_index *= compression_repeat;
        debug_multiple_index += col;
        assert(debug_multiple_index <= std::numeric_limits<int>::max());
#endif
        // In this loop prime >= first_row[col] so we can simply subtract first_row from prime.
        int multiple_index = ((prime - first_row[col]) * compression_primorial_inverse) % prime;
        for (int i = col + compression_repeat * multiple_index; i < sieve_size; i += step)
          sieve[i] = 0;
      }
    }
  }

  // Find the last prime that is less than or equal max_value.
  char* last_prime = sieve + sieve_size - 1;
  while (!*last_prime || sieve_index_to_prime(last_prime - sieve) > max_value)
    --last_prime;

  // Add all prime larger than sqrt(max_value) and less than or equal max_value.
  for (char* next_prime = sieve; ; ++next_prime)
    if (*next_prime)
    {
      prime_t prime = sieve_index_to_prime(next_prime - sieve);
      result.push_back(prime);
      if (next_prime == last_prime)
        break;
    }

  std::free(sieve);

  assert(result.size() <= upper_bound_number_of_primes); // Make sure it didn't grow.
  return result;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Dout(dc::notice, "compression = " << compression << "; compression_primorial = " << compression_primorial << "; compression_repeat = " << compression_repeat);

  auto primes = calculate_primes(1000000000);
  std::cout << "The " << primes.size() << "-th prime is " << primes.back() << std::endl;
}
