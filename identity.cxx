#include "sys.h"
#include "utils/macros.h"
#include <boost/multiprecision/cpp_int.hpp>
#include "debug.h"

#define PRINT_TERMS 0

using namespace boost::multiprecision;

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
  if (n < 0)
    return 0;

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

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  int max = 20;
  int nonzero = 0;
  for (int a = 1; a < max; ++a)
  {
    std::vector<cpp_int> sequence;
//    for (int m = 1; m <= max; ++m)
//      for (int n = a; n <= m; ++n)
      for (int n = a; n <= max; ++n)
      {
        int m = 2 * n;
        cpp_int sum_left = 0;
        int sign = 1;
        bool all_terms_zero_left = true;
        char const* sep = "";
        for (int k = 0; k <= n - a; ++k, sign = -sign)
        {
          cpp_int term = sign * binomial(n, k) * binomial(m - k, n);
#if PRINT_TERMS
          std::cout << sep << term;
          sep = " + ";
#endif
          if (term != 0)
            all_terms_zero_left = false;
          sum_left += term;
        }
#if PRINT_TERMS
        std::cout << " = ";
        sep = "";
#endif
        //Dout(dc::notice, "n = " << n << "; m = " << m << "; sum_left = " << sum_left);
        cpp_int sum_right = 0;
        sign = 1;
        bool all_terms_zero_right = true;
        for (int k = 0; k <= n - a; ++k, sign = -sign)
        {
          cpp_int term = sign * binomial(k + a - 1, a - 1) * binomial(m - n + a, k + a);
#if PRINT_TERMS
          std::cout << sep << term;
          sep = " + ";
#endif
          if (term != 0)
            all_terms_zero_right = false;
          sum_right += term;
        }
#if PRINT_TERMS
        std::cout << '\n';
#endif
        ASSERT(sum_left == sum_right);
        if (sum_left != 0)
          ++nonzero;
#if 0
        else
        {
          if (!all_terms_zero_left)
            Dout(dc::notice, "a = " << a << "; n = " << n << "; m = " << m << "; left has non-zero terms!");
          if (!all_terms_zero_right)
            Dout(dc::notice, "a = " << a << "; n = " << n << "; m = " << m << "; right has non-zero terms!");
          if (all_terms_zero_left && all_terms_zero_right)
            Dout(dc::notice, "Trivial zero: a = " << a << "; n = " << n << "; m = " << m);
        }
        if (sum_left != 0 && sum_left != 1)
        {
          Dout(dc::notice, "Sum = " << sum_right << ": a = " << a << "; n = " << n << "; m = " << m);
        }
#endif
        sequence.push_back(sum_right);
      }
    char const* sep = "";
    std::cout << "a=" << a << "; ";
    for (cpp_int e : sequence)
    {
      std::cout << sep << e;
      sep = ",";
    }
    std::cout << std::endl;
  }
}
