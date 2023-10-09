#include "sys.h"
#include "utils/macros.h"
#include "threadpool/AIThreadPool.h"
#include "fastprimes/Primes.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <array>
#include <vector>
#include <limits>
#include <iomanip>
#include "debug.h"

// This program calculates the series:
//
// (1 + 1/2) - (1/3 + 1/4 + 1/5) + (1/6 + 1/7 + 1/8 + 1/9 + 1/10) - (1/11 + ... 1/17) + (1/18 + ... 1/29) - ...
//
// That is, the harmonic series, alternating the sign ever next prime number of terms.
//
// The last output of this program (before it coredumped) was:
//
// 2222194682542268053 : 1.150064698632894453398729775362504003041548928880632443254665805
// 2222194692544241084 : 1.1500646941319508200316314461654747013642050885271996578013035
// n * a(455189504) = 2.0487823

using namespace boost::multiprecision;

//=============================================================================
// Factorial
//

std::vector<cpp_int> factorials = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600 };

cpp_int factorial(int n)
{
  if (AI_LIKELY(n < (int)factorials.size()))
    return factorials[n];
  return n * factorial(n - 1);
}

//=============================================================================
// Binomial
//

std::vector<cpp_int> binomials = { 2, 3, 4, 6, 5, 10, 6, 15, 20, 7, 21, 35, 8, 28, 56, 70, 9, 36, 84, 126, 10, 45, 120, 210, 252,
  11, 55, 165, 330, 462, 12, 66, 220, 495, 792, 924, 13, 78, 286, 715, 1287, 1716, 14, 91, 364, 1001, 2002, 3003, 3432, 15, 105,
  455, 1365, 3003, 5005, 6435, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 17, 136, 680, 2380, 6188, 12376, 19448, 24310, 18,
  153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 19, 171, 969, 3876, 11628, 27132, 50388, 75582, 92378, 20, 190, 1140, 4845,
  15504, 38760, 77520, 125970, 167960, 184756, 21, 210, 1330, 5985, 20349, 54264, 116280, 203490, 293930, 352716, 22, 231, 1540,
  7315, 26334, 74613, 170544, 319770, 497420, 646646, 705432, 23, 253, 1771, 8855, 33649, 100947, 245157, 490314, 817190, 1144066,
  1352078, 24, 276, 2024, 10626, 42504, 134596, 346104, 735471, 1307504, 1961256, 2496144, 2704156, 25, 300, 2300, 12650, 53130,
  177100, 480700, 1081575, 2042975, 3268760, 4457400, 5200300, 26, 325, 2600, 14950, 65780, 230230, 657800, 1562275, 3124550,
  5311735, 7726160, 9657700, 10400600, 27, 351, 2925, 17550, 80730, 296010, 888030, 2220075, 4686825, 8436285, 13037895, 17383860,
  20058300, 28, 378, 3276, 20475, 98280, 376740, 1184040, 3108105, 6906900, 13123110, 21474180, 30421755, 37442160, 40116600, 29,
  406, 3654, 23751, 118755, 475020, 1560780, 4292145, 10015005, 20030010, 34597290, 51895935, 67863915, 77558760, 30, 435, 4060,
  27405, 142506, 593775, 2035800, 5852925, 14307150, 30045015, 54627300, 86493225, 119759850, 145422675, 155117520, 31, 465, 4495,
  31465, 169911, 736281, 2629575, 7888725, 20160075, 44352165, 84672315, 141120525, 206253075, 265182525, 300540195, 32, 496, 4960,
  35960, 201376, 906192, 3365856, 10518300, 28048800, 64512240, 129024480, 225792840, 347373600, 471435600, 565722720, 601080390,
  33, 528, 5456, 40920, 237336, 1107568, 4272048, 13884156, 38567100, 92561040, 193536720, 354817320, 573166440, 818809200, 1037158320,
  1166803110, 34, 561, 5984, 46376, 278256, 1344904, 5379616, 18156204, 52451256, 131128140, 286097760, 548354040, 927983760, 1391975640,
  1855967520
};

cpp_int binomial(int n, int k)
{
  ASSERT(n >= 0);

  if (AI_UNLIKELY(2 * k > n))
    k = n - k;

  if (AI_UNLIKELY(k <= 0))
    return k == 0 ? 1 : 0;

  // We're now dealing with just the following table:
  //
  //      k: 0    1    2    3    4    5
  // n: 0                               index-k+1
  //    1
  //    2         2                         0
  //    3         3                         1
  //    4         4    6                    2
  //    5         5   10                    4
  //    6         6   15   20               6
  //    7         7   21   35               9
  //    8         8   28   56   70         12
  //    9         9   36   84  126         16
  //   10        10   45  120  210  252    20 <-- (n^2 - 2n + 1) / 4

  int cache_index = (n * (n - 2) + 1) / 4 + k - 1;

  if (AI_UNLIKELY(cache_index >= (int)binomials.size()))
    binomials.resize(cache_index + 1);

  if (AI_UNLIKELY(binomials[cache_index] == 0))
    binomials[cache_index] = binomial(n - 1, k - 1) + binomial(n - 1, k);

  return binomials[cache_index];
}

//=============================================================================
// Bernoulli numbers
//

std::vector<cpp_rational> Bernoulli_numbers {
   {1, 1}, {1, 6}, {-1, 30}, {1, 42}, {-1, 30}, {5, 66}, {-691, 2730}, {7, 6}, {-3617, 510},
   {43867, 798}, {-174611, 330}, {854513, 138}, {-236364091, 2730}, {8553103, 6}
};

cpp_rational Bplus(int n)
{
  if (AI_UNLIKELY(n == 1))
    return {1, 2};
  if (AI_UNLIKELY(n % 2 == 1))
    return {0, 1};

  int cache_index = n / 2;
  if (AI_UNLIKELY(cache_index >= (int)Bernoulli_numbers.size()))
  {
    // This assumes we call Bplus for incremental values of n/2.
    ASSERT(cache_index == (int)Bernoulli_numbers.size());

    // B_{n} = 1 - Sum_{k=0}^{n-1} Binomial(n, k) B_{k} / (n - k + 1).
    cpp_rational next_B = 1;
    for (int k = 0; k < n; ++k)
    {
      cpp_int denominator = n - k + 1;
      cpp_int numerator = binomial(n, k);
      cpp_rational term = cpp_rational{numerator, denominator} * Bplus(k);
      next_B -= term;
    }

    Bernoulli_numbers.push_back(next_B);
  }

  return Bernoulli_numbers[cache_index];
}

//=============================================================================
// Polygamma
//

constexpr int digits = 64; //256;
using hpfloat_t = number<cpp_dec_float<digits>>;

cpp_int epsilon_inverse = pow(cpp_int{10}, digits + 24);

hpfloat_t digamma(cpp_int n)
{
  // ψ(1 + z) = ln(z) + 1/(2z) - Sum_{j=1}^inf B_{2j} / (2j z^{2j})

  cpp_int z = n - 1;
  cpp_rational s{1, 2 * z};

  cpp_int z2 = z * z;
  cpp_int z2j = 1;
  for (int j2 = 2;; j2 += 2)
  {
    z2j *= z2;
    cpp_rational delta = Bplus(j2) / (j2 * z2j);
    s -= delta;
    if (abs(numerator(delta)) * epsilon_inverse < denominator(delta))
      break;
  }

  return log(static_cast<hpfloat_t>(z)) + static_cast<hpfloat_t>(s);
}

#if 0
// n^m * polygamma(m, n) = (-1)^(m+1) * ( (m-1)! * (1 + m / (2*n)) + sum_{k=1}^inf (2k+m-1)! / (2k)! * B_{2k} / n^(2k) )
//
// Note that (2k+m-1)! / (2k)! = (2k+m-1) * (2k+m-2) * (2k+m-3) * ... * (2k + 1)
// which is for
//     m=1 : 1
//     m=2 : (2k + 1)
//     m=3 : (2k + 1) * (2k + 2)
//     m=4 : (2k + 1) * (2k + 2) * (2k + 3)
// f(m, k) : (2k + 1) * (2k + 2) * (2k + 3) * ... * (2k + m - 1)
//
// Thus f(m, k+1) = (2(k+1) + 1) * (2(k+1) + 2) * (2(k+1) + 3) * ... * (2(k+1) + m - 1) =
//                = (2k + 3)     * (2k + 4)     * (2k + 5)     * ... * (2k     + m + 1) =
//                = f(m, k) * (2k + m) * (2k + m + 1) / ((2k + 1) * (2k + 2))

hpfloat_t polygamma1(int m, cpp_int n)
{
  hpfloat_t n_power = 1.0 / static_cast<hpfloat_t>(n);          // n to the power -1.
  hpfloat_t n_inverse_squared = n_power * n_power;              // n to the power -2.
  hpfloat_t result = 0.5 * n_inverse_squared;
  hpfloat_t eps = epsilon * n_power;

  for (int i = 0;; ++i)
  {
    hpfloat_t delta = static_cast<hpfloat_t>(B2(i)) * n_power;
    result += delta;
    if (i == number_of_B2s - 1 || abs(delta) < eps)
    {
      Dout(dc::notice, "Stop at " << i);
      break;
    }
    n_power *= n_inverse_squared;
  }

  return result;
}
#endif

//=============================================================================
// Partial Harmonic sum
//

// Return Sum_{k=s}^t 1/k
hpfloat_t partial_harmonic_sum(uint64_t s, uint64_t t)
{
  ASSERT(t >= s);
  cpp_rational r{1, s};
  for (uint64_t u = s + 1; u <= t; ++u)
    r += cpp_rational{1, u};
  return static_cast<hpfloat_t>(r);
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  AIThreadPool thread_pool(32);

  constexpr int capacity = 32;
  AIQueueHandle handler = thread_pool.new_queue(capacity);

  using namespace fastprimes;
  Primes generator(10000000000UL, 32, handler);

  hpfloat_t series = 0;
  int n = 1;
  while (n < 137488)
  {
    prime_t p1 = generator.next_prime();
    prime_t p2 = generator.next_prime();

    series += partial_harmonic_sum(n, n + p1 - 1);
    n += p1;
    series -= partial_harmonic_sum(n, n + p2 - 1);
    n += p2;

    std::cout << n << " : " << (n / p2) << " : " << std::setprecision(digits) << series << std::endl;
  }
  cpp_int n2 = n;
  hpfloat_t seriesb;
  int c = 0;
  for (;;)
  {
    prime_t p1 = generator.next_prime();
    prime_t p2 = generator.next_prime();

    series += digamma(n2 + p1) - digamma(n2);
    n2 += p1;
    hpfloat_t an = digamma(n2 + p2) - digamma(n2);

    seriesb = series;
    series -= an;
    n2 += p2;

    if (++c % 512 == 0)
    {
      std::cout << (n2 - p2) << " : " << std::setprecision(digits) << seriesb << std::endl;
      std::cout << n2  << " : " << std::setprecision(digits) << series << std::endl;
      int n = 2 * c;
      std::cout << "n * a(" << n << ") = " << std::setprecision(10) << (n * an) << std::endl;
    }
  }
}
