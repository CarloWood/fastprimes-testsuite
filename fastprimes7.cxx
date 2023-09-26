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
#include <climits>
#include <new>
#include "debug.h"

using prime_t = uint64_t;
using integer_t = int64_t;
using prime_index_t = int;
using sieve_word_t = uint64_t;

int modular_inverse(integer_t n, integer_t m)
{
  integer_t x0 = 1, x1 = 0;
  integer_t y0 = n, y1 = m;
  while (y1 > 0)
  {
    integer_t q = y0 / y1;
    std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
    std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
  }
  // Ensure the result is positive.
  return (x0 + m) % m;
}

// Called with n < m, n > 0 and m > 0.
int modular_inverse2(integer_t n, integer_t m)
{
  integer_t x0 = 0, x1 = 1;
  integer_t y0 = m, y1 = n;
  do
  {
    integer_t q = y0 / y1;
    std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
    std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
  }
  while (y1 > 0);
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

constexpr int compression = 6;                                      // The number of primes to skip; must be less than 7 or things overflow.
constexpr int compression_primorial = calc_primorial(compression);  // The multiplication of the first `compression` primes.
constexpr int compression_repeat = calc_repeat(compression);        // Same, but subtracting one from each prime first.
constexpr int compression_first_prime = small_primes[compression];  // The first integer that is not divisible by any of the `compression` primes.
constexpr int compression_first_prime_second_row = compression_first_prime + compression_primorial;

constexpr std::array<int, compression_repeat> calc_row0()
{
  std::array<int, compression_repeat> row0 = {};
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
      row0[col++] = candidate;
    candidate += 2;
  }
  while (col < compression_repeat);
  return row0;
}

constexpr std::array<int, compression_repeat> row0 = calc_row0();

prime_t sieve_index_to_prime(int row, int column)
{
  assert(column < row0.size());
  return row * compression_primorial + row0[column];
}

constexpr std::array<int, compression_first_prime_second_row> calc_row0_to_column()
{
  std::array<int, compression_first_prime_second_row> row0_to_column = {};
  int n = 0;
  for (int col = 0; col <= row0.size(); ++col)
  {
    while (n < (col == row0.size() ? compression_first_prime_second_row : row0[col]))
    {
      row0_to_column[n] = col - 1;
      ++n;
    }
  }
  return row0_to_column;
}

constexpr std::array<int, compression_first_prime_second_row> row0_to_column = calc_row0_to_column();

// If, for example, compression_repeat == 13 (which is never the case, but assume it was)
// and sieve_word_t is one byte, we'd have the following situation:
//
//                     111
// Column: 01234567  89012
//        [11111111][11111000] <-- One row of words_per_row (2) sieve_word_t of which the first three bits are unused (repeat = 13).
//         --------       --- <-- unused bits (3)
//         \_ sieve_word_bits = 8

constexpr int sieve_word_bits = sizeof(sieve_word_t) * CHAR_BIT;                                // The number of bits in one sieve_word_t.
constexpr int words_per_row = (compression_repeat + sieve_word_bits - 1) / sieve_word_bits;     // The width of one row in sieve_word_t.
constexpr int unused_bits = words_per_row * sieve_word_bits - compression_repeat;               // Number of unused bits at the end of a row.
// It's too much work to support unused bits (and would needlessly slow down the algorithm).
static_assert(unused_bits == 0);
constexpr sieve_word_t partial_mask = ~sieve_word_t{0} << unused_bits;

// Returns an upper bound on the number of primes that are smaller than or equal to n.
// Just a random function that I thought matched quite well (for n > 1000 or so).
// For n larger than 500,000,000 it is off on average 0.0079% (too large, thus).
uint32_t calc_upper_bound_number_of_primes(integer_t n)
{
  double logn = std::log(n);
  assert(logn > 4);
  return std::exp(0.3125 * std::pow(1 / (logn - 4), 1.655) + logn - std::log(logn - 1)) - 4;
}

#ifdef CWDEBUG
std::vector<prime_t> debug_primes;
void debug_init_primes()
{
  std::ifstream ifs("primes_till_1000000000", std::ios::binary);
  if (!ifs.is_open()) {
    std::cerr << "Failed to open file primes_till_1000000000 for reading." << std::endl;
    assert(false);
    return;
  }

  // Read the size of the vector first.
  size_t size;
  ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

  // Resize the vector and read the data.
  debug_primes.resize(size);
  ifs.read(reinterpret_cast<char*>(debug_primes.data()), size * sizeof(prime_t));

  ifs.close();
}
#endif

#define STORE_PRIME do { \
  result[pi] = prime; \
  ASSERT(prime == debug_primes[pi]); \
  ++pi; \
} while(0)

// Returns all primes less than or equal max_value.
std::vector<prime_t> calculate_primes(uint64_t max_value)
{
  Debug(debug_init_primes());

  // Pre-allocate the vector that is going to contain the primes.
  uint32_t upper_bound_number_of_primes = calc_upper_bound_number_of_primes(max_value);
  std::vector<prime_t> result(upper_bound_number_of_primes);

  // Copy the primes that are compressed away.
  int pi = 0;
  for (; pi < compression;)
  {
    prime_t prime = small_primes[pi];
    STORE_PRIME;
  }

  // Calculate how many integers are not divisible by the first `compression` number of primes that are less than or equal max_value.
  integer_t sieve_rows = (max_value - compression_first_prime) / compression_primorial + 1;

  // Lets allocate whole rows, even if we only partially need the last one.
  integer_t sieve_size = sieve_rows * words_per_row;    // The size of the sieve in words.

  // Allocate the sieve and fill it with ones.
  sieve_word_t* sieve = (sieve_word_t*)std::malloc(sieve_size * sizeof(sieve_word_t));
  if (!sieve)
    throw std::bad_alloc();
  // This assumes that unused_bits is zero.
  std::memset(sieve, 0xff, sieve_size * sizeof(sieve_word_t));
  if constexpr (unused_bits != 0)
  {
    for (sieve_word_t* word = sieve + words_per_row - 1; word < sieve + sieve_size; word += words_per_row)
      *word = partial_mask;
  }

  integer_t const sqrt_max_value = std::sqrt(max_value);
  assert(sqrt_max_value >= compression_first_prime_second_row);

  // sieve is a one dimensional vector, but can best be throught of as a two dimensional
  // table with width compression_repeat.

  // For example, if compression == 3 then compression_repeat == 8, and sieve represents
  // the following integers (all integers not divisible by 2, 3 or 5):
  //
  // row  |col: 0      1      2      3      4      5      6      7
  // -----+--------------------------------------------------------
  //  0 : |  !  7|  @ 11|  # 13|  $ 17|  % 19|  ^ 23|  & 29|  * 31|  <-- row0
  //  1 : |    37|    41|    43|    47|  ! 49|    53|    59|    61|
  //  2 : |    67|    71|    73| @! 77|    79|    83|    89| #! 91|
  //  3 : |    97|   101|   103|   107|   109|   113| $!119|  @121|
  //  4 : |   127|   131| %!133|   137|   139| #@143|   149|   151|
  //  5 : |   157| ^!161|   163|   167|  #169|   173|   179|   181|
  //  6 : | $@187|   191|   193|   197|   199| &!203| %@209|   211|
  //  7 : | *!217| $#221|   223|   227|   229|   233|   239|   241|
  //  8 : | %#247|   251| ^@253|   257|  !259|   263|   269|   271|
  //  9 : |   277|   281|   283|  !287|  $289|   293| ^#299|  !301|
  // 10 : |   307|   311|   313|   317| &@319| %$323|  !329|   331|
  // 11 : |   337| *@341|  !343|   347|   349|   353|   359|  %361|
  // 12 : |   367|  !371|   373| &#377|   379|   383|   389| ^$391|
  // 13 : |   397|   401| *#403|  @407|   409|  !413|   419|   421|
  // 14 : |  !427|   431|   433| ^%437|   439|   443|   449|  @451|
  // 15 : |   457|   461|   463|   467|  !469|  @473|   479|  #481|
  // 16 : |   487|   491| &$493|  !497|   499|   503|   509|  !511|
  // 17 : |  @517|   521|   523| *$527|  ^529|  #533| @!539|   541|
  // 18 : |   547| &%551|  !553|   557|  #559|   563|   569|   571|
  // 19 : |   577|  !581|  @583|   587| *%589|   593|   599|   601|
  // 20 : |   607|  #611|   613|   617|   619|  !623|  $629|   631|
  // 21 : | #!637|   641|   643|   647|  @649|   653|   659|   661|
  // 22 : | &^667|  @671|   673|   677|  !679|   683|  #689|   691|
  // 23 : |  $697|   701|  %703|  !707|   709| *^713|   719|  !721|
  // 24 : |   727|  $731|   733|  @737|   739|   743|  !749|   751|
  // 25 : |   757|   761|  !763|  #767|   769|   773|  %779|  @781|
  // 26 : |   787|  !791|  #793|   797|  $799|  @803|   809|   811|
  // 27 : |  %817|   821|   823|   827|   829| $!833|   839|  &841|
  // 28 : | @!847|  ^851|   853|   857|   859|   863|  @869|  #871|
  // 29 : |   877|   881|   883|   887|  !889|  %893| *&899|  $901|
  // 30 : |   907|   911|  @913|  !917|   919|  #923|   929| %!931|
  // 31 : |   937|   941|  ^943|   947|  #949|   953|  !959|  *961|
  // 32 : | ... etc
  //
  // Here every field contains the integer that is represented by the corresponding bit in sieve.
  // All such integers divisible by 7 are marked with a '!', those that are divisible by 11 are
  // marked with a '@', and so on. Below we'll say "contains a number" while
  // really it's just a single bit in sieve that represents that number.
  //
  // All integers in row 0 have already been calculated during compile time
  // and are stored in row0.
  //
  // In the loop(s) below, we find the next prime (p) and set its corresponding
  // bits to zero. This is done per column; for each column we determine
  // the first row that contains an integer in that column that is divisible
  // by p. For example, for p=13 and col=1, we find that the first number that
  // is divisible by 13 is 221 (= 13*17). While the first number in column 3
  // that is divisible by the prime 11 is 77, row 2.
  //
  // Every subsequent multiple of that prime in the given colum is precisely
  // prime rows lower every time. For example, the vertical distance between
  // every '!' is 7, and between every '$' is 17, etc.
  //
  // Note that the first row does NOT necessarily only contain primes.
  // In this case, the first number that is not a prime is 7^2 = 49 which is
  // larger than 31, the last number in the first row; and thus all numbers in
  // the first row are primes in this case. But that is not the case for
  // compression 4 or larger: for compression 4 the first (top-left) number in
  // the table is 11 and the first non-prime is 11^2 = 121. The first number of
  // row 1 is 11 + (2 * 3 * 5 * 7) = 221. Therefore 121 is part of the first row
  // (and so are 11*13=143, 13*13=169, 11*17=187 and 11*19=209).
  //
  // The first row for a given column that contains a number that is divisible
  // by a prime p is given by: -row0[col] / compression_primorial [MODULO p].
  //
  // For example, we want to know the row in column 2 that contains the first
  // number that is divisible by 19. Note that each time we go one row down,
  // the number is incremented with compression_primorial.
  //
  // row  col=2  MOD 19
  //  0    13      13     = 13 + 0 * 30
  //  1    43       5     = 13 + 1 * 30
  //  2    73      16     = 13 + 2 * 30
  //  3   103       8     = 13 + 3 * 30
  //  4   133       0     = 13 + 4 * 30 --> 4 * 30 = -13 --> 4 = -13 / 30 [MOD 19]
  //
  // The compression_primorial is a constant (only depending on the compression used).
  // But the modulo (p) is not. The row index can be calculated slightly more
  // efficient when p > compression_primorial, which is at least the case for every
  // prime in row 1 or larger (for example, in the case of the above table, with
  // compression = 3, the prime 31 is already larger than the compression_primorial=30,
  // but also 37, the first prime of row 1, is).

  int row = 0;
  int column;
  sieve_word_t column_mask = 1;
  prime_t prime = 0;
  for (;;)
  {
    for (column = 0; column < compression_repeat; ++column)
    {
      int column_word_offset = column / sieve_word_bits;
      sieve_word_t* next_word = sieve + row * words_per_row + column_word_offset;
      sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
      if ((*next_word & column_mask))
      {
        prime = sieve_index_to_prime(row, column);

        if (prime > compression_primorial)
          break;

        STORE_PRIME;

        int const word_step = words_per_row * prime;
        int const offset = compression_first_prime_second_row * prime;
        int compression_primorial_inverse = modular_inverse(compression_primorial, prime);

        for (int col = 0; col < compression_repeat; ++col)
        {
          // The largest value of `offset - row0[col]` is when `offset` has its largest value
          // and row0[col] is at its smallest. The latter happens when col = 0 (at which point
          // row0[col] equals compression_first_prime). The former, `offset` at its largest,
          // happens when `prime` is at its largest, which is `compression_first_prime_second_row`.
          // Note that compression_first_prime_second_row = compression_first_prime + compression_primorial
          //
          // Let P = compression_primorial, F = compression_first_prime.
          // Then compression_first_prime_second_row = P + F.
          //
          // And we can write for the largest values involved:
          //   prime = P + F,
          //   offset = (P + F) * (P + F);
          //   -row0[col] = -F
          //
          // The largest possible value of compression_primorial_inverse is prime - 1, or P + F - 1.
          //
          // Thus the largest possible value of ((offset - row0[col]) * compression_primorial_inverse) is (less than or) equal
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

          // Calculate the first row that has a multiple of this prime in colum `col`.
          int first_row_with_prime_multiple = offset - row0[col];
          if constexpr (compression_first_prime_second_row > 1290)
            first_row_with_prime_multiple %= prime;
          first_row_with_prime_multiple *= compression_primorial_inverse;
          first_row_with_prime_multiple %= prime;
          int col_word_offset = col / sieve_word_bits;
          sieve_word_t col_mask  = sieve_word_t{1} << (col % sieve_word_bits);

          for (int word_index = first_row_with_prime_multiple * words_per_row + col_word_offset;
               word_index < sieve_size;
               word_index += word_step)
          {
            sieve[word_index] &= ~col_mask;

            int debug_row = word_index / words_per_row;
            int debug_col = (word_index % words_per_row) * sieve_word_bits + (col % sieve_word_bits);
            prime_t debug_prime = sieve_index_to_prime(debug_row, debug_col);
            assert(debug_col == col);
            assert(debug_prime % prime == 0);
//            Dout(dc::notice, "Loop1: setting " << debug_prime << " to 0 because it is " << (debug_prime / prime) << " * " << prime);
          }
        }
      }
    }
    if (prime > compression_primorial)
      break;
    // Next row.
    ++row;
    column_mask = 1;
  }

  for (;;)
  {
    for (; column < compression_repeat; ++column)
    {
      int column_word_offset = column / sieve_word_bits;
      sieve_word_t* next_word = sieve + row * words_per_row + column_word_offset;
      sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
      if ((*next_word & column_mask))
      {
        prime = sieve_index_to_prime(row, column);
        assert(prime > compression_primorial);

        STORE_PRIME;

        if (prime > sqrt_max_value)
          break;

        int const word_step = words_per_row * prime;
        int compression_primorial_inverse = modular_inverse2(compression_primorial, prime);

        for (int col = 0; col < compression_repeat; ++col)
        {
#if 0
          long debug_prime = prime;
          long debug_row0 = row0[col];
          long debug_inverse = compression_primorial_inverse;
          long debug_multiple_index = prime - row0[col];
          assert(0 <= debug_multiple_index && debug_multiple_index < prime);
          debug_multiple_index *= debug_inverse;
          assert(debug_multiple_index <= std::numeric_limits<int>::max());
          debug_multiple_index %= debug_prime;
          debug_multiple_index *= compression_repeat;
          debug_multiple_index += col;
          assert(debug_multiple_index <= std::numeric_limits<int>::max());
#endif
          // In this loop prime >= row0[col] so we can simply subtract row0 from prime.
          int first_row_with_prime_multiple = ((prime - row0[col]) * compression_primorial_inverse) % prime;
          int col_word_offset = col / sieve_word_bits;
          sieve_word_t col_mask  = sieve_word_t{1} << (col % sieve_word_bits);

          for (int word_index = first_row_with_prime_multiple * words_per_row + col_word_offset;
               word_index < sieve_size;
               word_index += word_step)
          {
            sieve[word_index] &= ~col_mask;

            int debug_row = word_index / words_per_row;
            int debug_col = (word_index % words_per_row) * sieve_word_bits + (col % sieve_word_bits);
            prime_t debug_prime = sieve_index_to_prime(debug_row, debug_col);
            assert(debug_col == col);
            assert(debug_prime % prime == 0);
//            Dout(dc::notice, "Loop2: setting " << debug_prime << " to 0 because it is " << (debug_prime / prime) << " * " << prime);
          }
        }
      }
    }
    if (prime > sqrt_max_value)
      break;
    // Next row.
    ++row;
    column = 0;
    column_mask = 1;
  }

  // Find the word containing the last prime that is less than or equal max_value.
  int last_word_index = sieve_size;
  while (sieve[--last_word_index] == 0 ||
      sieve_index_to_prime(last_word_index / words_per_row, (last_word_index % words_per_row) * sieve_word_bits) > max_value)
    ;
  int last_row = last_word_index / words_per_row;

  // Add all primes larger than sqrt(max_value) and less than or equal max_value.
  ++column;
  while (row < last_row)
  {
    for (; column < compression_repeat; ++column)
    {
      int column_word_offset = column / sieve_word_bits;
      sieve_word_t* next_word = sieve + row * words_per_row + column_word_offset;
      sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
      if ((*next_word & column_mask))
      {
        prime = sieve_index_to_prime(row, column);
        STORE_PRIME;
      }
    }
    // Next row.
    ++row;
    column = 0;
    column_mask = 1;
  }
  for (; column < compression_repeat; ++column)
  {
    int column_word_offset = column / sieve_word_bits;
    sieve_word_t* next_word = sieve + row * words_per_row + column_word_offset;
    sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
    if ((*next_word & column_mask))
    {
      prime = sieve_index_to_prime(row, column);
      if (prime > max_value)
        break;
      STORE_PRIME;
    }
  }

  std::free(sieve);

  assert(pi <= upper_bound_number_of_primes); // Make sure it didn't overflow.

  result.resize(pi);
  return result;
}

#include <fstream>

void writeVectorToDisk(std::string const& filename, std::vector<prime_t> const& vec)
{
  std::ofstream ofs(filename, std::ios::binary);

  if (!ofs.is_open()) {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  // Write the size of the vector first,
  size_t size = vec.size();
  ofs.write(reinterpret_cast<char const*>(&size), sizeof(size));

  // Then write the vector data.
  ofs.write(reinterpret_cast<char const*>(vec.data()), vec.size() * sizeof(prime_t));

  ofs.close();
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Dout(dc::notice, "compression = " << compression << "; compression_primorial = " << compression_primorial << "; compression_repeat = " << compression_repeat);

  auto primes = calculate_primes(1000000000);
  std::cout << "The " << primes.size() << "-th prime is " << primes.back() << std::endl;
}
