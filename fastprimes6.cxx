#include "sys.h"
#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <cassert>

#include "debug.h"

using prime_t = uint32_t;
using integer_t = uint64_t;
using prime_index_t = int;
using prime_gap_t = prime_t;

#ifdef CWDEBUG
std::array<prime_t, 101> debug_primes = {
  1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541 };

prime_t debug_prime_after(prime_t x)
{
  for (int p : debug_primes)
    if (p >= x)
      return p;
  assert(false);
}
#endif

class Sieve;

class Primes
{
 private:
  std::vector<prime_t> primes_;
  prime_index_t current_sieve_;
  integer_t next_square_;
  Sieve* sieve_;

 public:
  Primes();
  ~Primes();

  prime_t p(prime_index_t n) const
  {
    assert(1 <= n && n < primes_.size());
    return primes_[n];
  }

  prime_t next();

  integer_t primorial(prime_index_t n) const;
  prime_index_t repeat(prime_index_t n) const;
  prime_index_t square(prime_index_t n) const;

#ifdef CWDEBUG
  prime_t debug_p(prime_index_t n) const
  {
    assert(1 <= n);
    if (n < primes_.size())
      return primes_[n];
    if (!(n < debug_primes.size()))
    {
      Dout(dc::warning, "Calling debug_p(" << n << "), returning 0.");
      return 0;
    }
    return debug_primes[n];
  }

  integer_t debug_primorial(prime_index_t n) const;

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, Primes const& primes) { primes.print_on(os); return os; }
#endif
};

// Returns pₙ# = p₁⋅p₂⋅p₃… pₙ
// Note that p₀ := 1
integer_t Primes::primorial(prime_index_t n) const
{
  if (n < 1)
    return 1;
  return p(n) * primorial(n - 1);
}

#ifdef CWDEBUG
integer_t Primes::debug_primorial(prime_index_t n) const
{
  if (n < 1)
    return 1;
  return debug_p(n) * debug_primorial(n - 1);
}
#endif

// Returns Product_{k=1..n} (pₖ-1)
// For example, repeat(4) = (p₁-1)⋅(p₂-1)⋅(p₃-1)⋅(p₄-1) = 1⋅2⋅4⋅6 = 48
prime_index_t Primes::repeat(prime_index_t n) const
{
  if (n < 3)
    return 1;
  return (p(n-1) - 1) * repeat(n - 1);
}

// Returns pₙ²
prime_index_t Primes::square(prime_index_t n) const
{
  return p(n) * p(n);
}

class Sieve
{
 private:
  Primes* primes_;
  Sieve* prev_sieve_;
  prime_index_t begin_;
  prime_t prev_prime_;
  prime_index_t end_;
  prime_index_t n_;
  prime_index_t factor_n_;
  integer_t base_;
  integer_t factor_base_;
  integer_t primorial_;
  integer_t expected_failure_;
#ifdef CWDEBUG
  prime_t debug_lower_bound_;
  prime_t debug_upper_bound_;
  integer_t debug_last_next_;
#endif

 public:
  Sieve(Primes* primes) :
    primes_(primes),
    prev_sieve_(nullptr),
    begin_(1),
    prev_prime_(1),
    end_(2),
    n_(1),
    factor_n_(1),
    base_(0),
    factor_base_(1),
    primorial_(1),
    expected_failure_(0)
#ifdef CWDEBUG
    , debug_lower_bound_(2),
    debug_upper_bound_(5),
    debug_last_next_(0)
#endif
  { }

  ~Sieve()
  {
    delete prev_sieve_;
  }

  // [2]  [3]  [5]  [7]  [11]   [13]
  //
  //  2
  // +1
  //  3  ( 3)
  // ?4 <2*2>
  //      +2
  //       5  ( 5)
  //          (+2)
  //       7  ( 7) ( 7)
  //      ?9 <3*3>
  //           +4  (+4)
  //     ?11   11  (11)  (11)
  //               (+2)  (+2)
  //           13  (13)  (13)   (13)
  //               (+4)  (+4)   (+4)
  //           17  (17)  (17)   (17)
  //               (+2)  (+2)   (+2)
  //           19  (19)  (19)   (19)
  //               (+4)  (+4)   (+4)
  //           23  (23)  (23)   (23)
  //          ?25 <5*5>
  //                +6   (+6)   (+6)
  //          ?29   29   (29)   (29)
  //                +2   (+2)   (+2)
  //          ?31   31   (31)   (31)
  //          ?35 <5*7>
  //                +6   (+6)   (+6)
  //          ?37   37   (37)   (37)
  //                     (+4)   (+4)
  //                41   (41)   (41)
  //                     (+2)   (+2)
  //                43   (43)   (43)
  //                     (+4)   (+4)
  //                47   (47)   (47)
  //               ?49  <7*7>
  //                      +6    (+6)
  //               ?53    53    (53)
  //                      +6    (+6)
  //               ?59    59    (59)
  //                      +2    (+2)
  //               ?61    61    (61)
  //                      +6    (+6)
  //               ?67    67    (67)
  //                      +4    (+4)
  //               ?71    71    (71)
  //                      +2    (+2)
  //               ?73    73    (73)
  //               ?77 <7*11>
  //                      +6    (+6)
  //               ?79    79    (79)
  //                      +4    (+4)
  //               ?83    83    (83)
  //                      +6    (+6)
  //               ?89    89    (89)
  //               ?91 <7*13>
  //                      +8    (+8)
  //               ?97    97    (97)
  //                      +4    (+4)
  //              ?101   101   (101)
  //                      +2    (+2)
  //              ?103   103   (103)
  //                      +4    (+4)
  //              ?107   107   (107)
  //                      +2    (+2)
  //              ?109   109   (109)
  //                      +4    (+4)
  //              ?113   113   (113)
  //              ?119 <7*17>
  //              ?121  ?121 <11*11>
  //                +6    +6    +14
  //              ?127  ?127    127
  //                +4    +4     +4
  //              ?131  ?131    131
  //                +6    +6     +6
  //              ?137  ?137    137
  //                +2    +2     +2
  //              ?139  ?139    139
  //                +4    +4
  //              ?143  ?143 <11*13>
  //                +6    +6    +10
  //              ?149  ?149    149
  //                +2    +2     +2
  //              ?151  ?151    151
  //                +6    +6     +6
  //              ?157  ?157    157
  //                +6    +6     +6
  //              ?163  ?163    163
  //                +4    +4     +4
  //              ?167  ?167    167
  //                +2    +2     +2
  //              ?169  ?169    169
  //                +4    +4     +4
  //              ?173  ?173    173
  //                +6    +6     +6
  //              ?179  ?179    179
  //                +2    +2     +2
  //              ?181  ?181    181
  //                +6    +6
  //              ?187  ?187 <11*17>
  //                +4    +4    +10
  //              ?191  ?191    191
  //                +2    +2     +2
  //              ?193  ?193    193
  //                +4    +4     +4
  //              ?197  ?197    197
  //                +2    +2     +2
  //              ?199  ?199    199
  //               +10   +10
  //              ?209  ?209 <11*19>
  //                +2    +2    +12
  //              ?211  ?211    211
  //               +10   +10    +10
  //              ?221  ?221    221

  Sieve(Sieve* prev_sieve) :
    primes_(prev_sieve->primes()),
    prev_sieve_(prev_sieve),
    begin_(prev_sieve->begin() + 1),
    prev_prime_(primes_->p(begin_ - 1)),
    end_(begin_ + primes_->repeat(begin_)),
    n_(begin_ + 1),
    factor_n_(begin_ + 1),
    base_(0),
    factor_base_(0),
    primorial_(primes_->primorial(begin_ - 1)),
    expected_failure_(prev_prime_ * prev_sieve->next_factor())
#ifdef CWDEBUG
    , debug_lower_bound_(debug_prime_after(primes_->square(begin_ - 1))),
    debug_upper_bound_(primes_->debug_p(begin_ + 1) + primes_->debug_primorial(begin_)),
    debug_last_next_(0)
#endif
  {
    DoutEntering(dc::notice, "Sieve([" << prev_sieve->primes()->debug_p(prev_sieve->begin()) << "]) [" << prev_sieve->primes()->debug_p(begin_) << "]");
    if (n_ == end_)
    {
      n_ = begin_;
      base_ += primorial_;
      expected_failure_ = 0;
    }
    assert(n_ < end_);
    if (factor_n_ == end_)
    {
      factor_n_ = begin_;
      factor_base_ += primorial_;
    }
    assert(factor_n_ < end_);
    Dout(dc::notice, "*this = " << *this);
  }

  Primes* primes() const
  {
    return primes_;
  }

  prime_index_t begin() const
  {
    return begin_;
  }

#ifdef CWDEBUG
  void debug_returning(integer_t next)
  {
    prime_t p = primes_->debug_p(begin_);
    Dout(dc::notice, '[' << p << "] returning: " << next << " (range: [" << debug_lower_bound_ << ", " << debug_upper_bound_ << "])");
    assert(debug_lower_bound_ <= next && next <= debug_upper_bound_);
    if (next == debug_lower_bound_)
      Dout(dc::notice, "This is the first one for [" << p << "].");
    if (next == debug_upper_bound_)
      Dout(dc::notice, "This was the last one for [" << p << "]!");
    debug_last_next_ = next;
  }
#endif

  integer_t next()
  {
    DoutEntering(dc::notice, "Sieve::next() [" << primes_->debug_p(begin_) << "] with expected_failure_ = " << expected_failure_);

#ifdef CWDEBUG
    if (primes_->debug_p(begin_) == 11 && debug_last_next_ == 323)
      Dout(dc::notice, "this goes wrong");
#endif

    if (!expected_failure_)
    {
      Dout(dc::notice, "expected_failure_ == 0; setting next to " << base_ << " + " << primes_->debug_p(n_));
      assert(n_ < end_);
      integer_t next = base_ + primes_->p(n_);
      if (++n_ == end_)
      {
        base_ += primorial_;
        Dout(dc::notice, "Wrapped around; base_ is now " << base_);
        n_ = begin_;
      }
      Debug(debug_returning(next));
      return next;
    }

    integer_t next = prev_sieve_->next();

    if (!base_)
    {
      Dout(dc::notice, "Initial run!");
      base_ = primorial_;
      if (next - primes_->p(begin_) == primorial_)
      {
        Dout(dc::notice, "Finished repeat! Setting expected_failure_ to zero.");
        expected_failure_ = 0;
      }
      Debug(debug_returning(next));
      return next;
    }

    if (next == expected_failure_)
    {
      Dout(dc::notice, "next == expected_failure_ (" << next << "); skipping...");
      next = prev_sieve_->next();
      if (next - primes_->p(begin_) == primorial_)
      {
        Dout(dc::notice, "Finished repeat, setting expected_failure_ to zero.");
        expected_failure_ = 0;
      }
      else
      {
        expected_failure_ = prev_prime_ * prev_sieve_->next_factor();
        Dout(dc::notice, "Set expected_failure_ to " << prev_prime_ << " * " << (expected_failure_ / prev_prime_));
      }
    }
    Debug(debug_returning(next));
    return next;
  }

  integer_t next_factor()
  {
    DoutEntering(dc::notice, "next_factor() [" << primes_->debug_p(begin_) << "]");
    Dout(dc::notice, "factor_base_ = " << factor_base_ << "; factor_n_ = " << factor_n_);
    integer_t next = factor_base_ + primes_->p(factor_n_);
    Dout(dc::notice, "Set next_factor to " << factor_base_ << " + " << primes_->p(factor_n_) << " = " << next);
    if (++factor_n_ == end_)
    {
      factor_base_ += primorial_;
      factor_n_ = begin_;
      Dout(dc::notice, "Wrapped around; factor_base_ is now " << factor_base_ << "; reset factor_n_ to " << begin_);
    }
    return next;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, Sieve const& sieve) { sieve.print_on(os); return os; }
#endif
};

Primes::Primes() : primes_({1, 2}), current_sieve_(1), next_square_(4), sieve_(new Sieve(this))
{
}

Primes::~Primes()
{
  delete sieve_;
}

prime_t Primes::next()
{
  prime_t candidate = sieve_->next();
  if (candidate == next_square_)
  {
    Dout(dc::notice, "Returned candidate is equal to the next square (" << next_square_ << ").");
    sieve_ = new Sieve(sieve_);
    next_square_ = square(sieve_->begin());
    Dout(dc::notice, "New sieve: " << *sieve_ << "; next_square_ = " << next_square_);
    candidate = sieve_->next();
  }

  if (candidate > 2)
  {
    assert(candidate > primes_.back());
    primes_.push_back(candidate);
  }
  return candidate;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Primes primes;

  for (int n = 1; n < 100; ++n)
  {
    prime_t p = primes.next();
    Dout(dc::notice, "NEXT PRIME: " << p);
#ifdef CWDEBUG
    if (p != debug_primes[n])
      DoutFatal(dc::core, "That is NOT the next prime, which is " << debug_primes[n] << "!");
#endif
  }
}

#ifdef CWDEBUG
void Sieve::print_on(std::ostream& os) const
{
  os << '{';
  os << "primes_:" << primes_ <<
      ", prev_sieve_:" << prev_sieve_ <<
      ", begin_:" << begin_ << " (" << primes_->debug_p(begin_) << ")"
      ", prev_prime_:" << prev_prime_ <<
      ", end_:" << end_ << " (" << primes_->debug_p(end_) << ")"
      ", n_:" << n_ <<
      ", factor_n_:" << factor_n_ <<
      ", base_:" << base_ <<
      ", factor_base_:" << factor_base_ <<
      ", primorial_:" << primorial_ <<
      ", expected_failure_:" << expected_failure_ <<
      ", debug_lower_bound_:" << debug_lower_bound_ <<
      ", debug_upper_bound_:" << debug_upper_bound_;
  os << '}';
}

void Primes::print_on(std::ostream& os) const
{
  char const* sep = "";
  os << '{';
  for (prime_t p : primes_)
  {
    os << sep << p;
    sep = ", ";
  }
  os << '}';
}
#endif
