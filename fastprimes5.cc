#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cassert>

using prime_t = uint32_t;
using integer_t = uint64_t;
using prime_index_t = int;
using prime_gap_t = prime_t;

class Sieve;

// Consider integers not divisible by 2, 3, 5 or 7:
//
// SieveWindow starting at 11     Currently ending with 43                                                 11^2
//  |__ start_                          |__ head_                                                           |__ last_
//  v                                   v                                                                   v
// 11 [13  17  19  23  29  31  37  41  43] 47  53  59  61  67  71  73  79  83  89  97 101 103 107 109 113 121
//    { 2   4   2   4   6   2   6   4   2   4   6   6   2   6   4   2   6   4   6   8   4   2   4   2   4   8
//  ^
//  |__ tail_=-1

// SieveWindow starting at 121  Currently ending with 151                                      __ last_
//  |__ start_                    |__ head_                                                   |
//  v                             v                                                           v
// 121 [127 131 137 139 143 149 151]157 163 167 169 173 179 181 187 191 193 197 199 209 211 221
//        6   4   6   2   4   6   2   6   6   4   2   4   6   2   6   4   2   4   2  10   2  10 }
//   ^                                                                                          ^
//   |__ tail_=-1                                                                               |
//                                                                                           starts to repeat here because
//                                                                                           221 == 11 + 2 * 3 * 5 * 7
//
//    223 227 229 233 239 241 247 251 253 257 263 269 271 277 281 283 289 293 299 307 311 313 317 319 323 331 337
//      2   4   2   4   6   2   6   4   2   4   6   6   2   6   4   2   6   4   6   8   4   2   4   2   4   8   6
//      ^
//      |
//    Same gap as 13, etc.

class SieveWindow
{
 private:
  integer_t start_;     // 11 or 11^2, see above.
  integer_t head_;      // Begins at start_, and then becomes equal to the last value returned.
  integer_t last_;      // The last value that should be returned, after which we need to switch to the other window.
  std::vector<prime_gap_t> gaps_;
  int tail_;            // Index into gaps_ of the value that was returned last.
  Sieve* prev_sieve_;
  SieveWindow* other_window_;

 public:
  SieveWindow(integer_t start, integer_t last, Sieve* prev_sieve, SieveWindow* other_window) :
    start_(start), head_(start), last_(last), tail_(-1), prev_sieve_(prev_sieve), other_window_(other_window)
  {
  }

  void destroy_sieve();

  void reset_tail() { tail_ = -1; }
  prime_gap_t read_gap(SieveWindow*& window, bool skip = false);

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, SieveWindow const& sieve_window) { sieve_window.print_on(os); return os; }
};

// Return the "prime factorial" of prime with index pi.
//
// primes must contain { 1, 2, 3, 5, 7, 11, 13, ... }.
// For example, if pi = 3 then returns 5 * 3 * 2.
//
integer_t prime_fac(std::vector<prime_t> const& primes, prime_index_t pi)
{
  integer_t fac = 1;
  for (; pi > 0; --pi)
    fac *= primes[pi];
  return fac;
}

class Sieve
{
 private:
  prime_index_t pi_;                    // For example 5.
  prime_t p_;                           // For example 11 (primes[pi_]).
  SieveWindow window1_;
  SieveWindow window2_;
  SieveWindow* offset_window_;
  SieveWindow* gap_window_;
  prime_t last_prime_;
  integer_t last_offset_;
  integer_t expected_failure_;          // The next expected composite.
  integer_t next_prime_squared_;        // The square of the prime after p_ (for example 169).

  Sieve(Sieve const&) = delete;

 public:
  Sieve(std::vector<prime_t> const& primes) :
    pi_(0),
    p_(1),
    window1_(1, 1, nullptr, &window2_),
    window2_(1, 2, nullptr, &window1_),
    offset_window_(&window1_),
    gap_window_(&window2_),
    last_prime_(1),
    last_offset_(1),
    expected_failure_(4),       // <nothing> --> next_prime_squared_
    next_prime_squared_(4)
  {
    std::cout << "Created " << *this << std::endl;
  }

  // Lets define pₙ to be the n-th prime. p₁=2, p₂=3, p₃=5, p₄=7, ...
  //
  // Let sp(pₙ) = p₍ₙ₋₁₎# = p₁⋅p₂⋅p₃… p₍ₙ₋₁₎
  //
  // Returns sp(pₙ)
  //
  static integer_t sub_primorial(std::vector<prime_t> const& primes, prime_index_t n)
  {
    if (n < 2)
      return 1;
    return sub_primorial(primes, n - 1) * primes[n - 1];        // pₙ = primes[n].
  }

  // For example,
  // p = p₄ = 7
  //
  // window1:
  //                          __ sums up to sp(p) = 30
  //                         /
  //         <---------------------------->
  //        {+4  +2  +4  +2  +4  +6  +2  +6} +4  +2  +4  +2  +4  +6  +2  +6  +4  +2  +4  +2  +4  +6  +2  +6
  //      7  11  13  17  19  23  29  31  37  41  43  47  49  53  59  61  67  71  73  77  79  83  89  91  97
  //      ^                   ^           ^
  //      |                   |           |
  //      p              lined up    min(p + sp(p), p²)
  //                          |
  // window2:                 v
  //                        {+4  +6  +2  +6  +4  +2  +4} +2  +4  +6  +2  +6
  //                     49  53  59  61  67  71  73  77  79  83  89  91  97
  //                      ^               ^
  //                      |               |
  //                      p²        min(p² + p + sp(p) - (p² mod sp(p)), p² + sp(p))
  //
  //      p² + p + sp(p) - (p² mod sp(p)) = 49 + 7 + 30 - (49 mod 30) = 86 - (49 mod 30) = 86 - 19 = 67

  static integer_t window1_last(std::vector<prime_t> const& primes, prime_index_t pi, prime_t p)
  {
    integer_t p2 = p * p;
    return std::min(p + sub_primorial(primes, pi), p2);
  }

  static integer_t window2_last(std::vector<prime_t> const& primes, prime_index_t pi, prime_t p)
  {
    integer_t p2 = p * p;
    integer_t sp = sub_primorial(primes, pi);
    return std::min(p2 + p + sp - (p2 % sp), p2 + sp);
  }

  Sieve(std::vector<prime_t> const& primes, prime_index_t pi, Sieve* prev_sieve) :
    pi_(pi),
    p_(primes[pi]),
    window1_(p_, window1_last(primes, pi, p_), prev_sieve, &window2_),
    window2_(p_ * p_, window2_last(primes, pi, p_), prev_sieve, &window1_),
    offset_window_(&window1_),
    gap_window_(&window2_),
    last_prime_(p_ * p_),
    last_offset_(p_),
    expected_failure_(get_first_expected_failure()),
    next_prime_squared_(primes[pi + 1] * primes[pi + 1])
  {
    std::cout << "Created " << *this << std::endl;
  }

  ~Sieve()
  {
    window1_.destroy_sieve();
  }

  integer_t expected_failure() const
  {
    return expected_failure_;
  }

  integer_t get_first_expected_failure()
  {
    offset_window_->reset_tail();
    last_offset_ += offset_window_->read_gap(offset_window_);
    return p_ * last_offset_;
  }

  void set_next_expected_failure()
  {
    last_offset_ += offset_window_->read_gap(offset_window_);
    expected_failure_ = p_ * last_offset_;
  }

  prime_t next(std::vector<prime_t> const& primes, Sieve*& sieve)
  {
    std::cout << "Entering Sieve::next() [with p_ = " << p_ << "]" << std::endl;
    last_prime_ += gap_window_->read_gap(gap_window_, true);
    std::cout << "Next candidate: " << last_prime_ << " with expected_failure_ = " << expected_failure_ << " and next_prime_squared_ = " << next_prime_squared_ << std::endl;
    if (last_prime_ == next_prime_squared_ || last_prime_ == expected_failure_)
    {
      if (last_prime_ == next_prime_squared_)
      {
        sieve = new Sieve(primes, pi_ + 1, sieve);
        return sieve->next(primes, sieve);
      }
      last_prime_ += gap_window_->read_gap(gap_window_);
      set_next_expected_failure();
      std::cout << "Next candidate: " << last_prime_ << " with expected_failure_ = " << expected_failure_ << " and next_prime_squared_ = " << next_prime_squared_ << std::endl;
      // If expected_failure_ < next_prime_squared_, this can still happen now.
      if (last_prime_ == next_prime_squared_)
      {
        sieve = new Sieve(primes, pi_ + 1, sieve);
        return sieve->next(primes, sieve);
      }
    }
    return last_prime_;
  }

  prime_index_t prime_index() const
  {
    return pi_;
  }

  prime_gap_t next_gap()
  {
    return 1;
  }

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, Sieve const& sieve) { sieve.print_on(os); return os; }
};

void SieveWindow::destroy_sieve()
{
  delete prev_sieve_;
  prev_sieve_ = nullptr;
}

bool divisible_by_any_prime_below(integer_t n, prime_index_t pi)
{
  static std::array<prime_t, 10> primes = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23 };
  assert(pi - 1 < (ssize_t)primes.size());
  for (int i = pi - 1; i > 0; --i)
    if (n % primes[i] == 0)
    {
      std::cout << n << " is divisible by " << primes[i] << "!\n";
      return true;
    }
  return false;
}

//                                   __ tail_ (return +6 upon calling read_gap)
//                                  |            __ gaps_.size()
//                                  |           |
// start_ = 121                     v           v
//                  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21    <-- index into gaps_
//              [  +6  +4  +6  +2  +4  +6  +2 (+6  +6  +4  +2  +4  +6  +2  +6  +4  +2  +4  +2 +10  +2 +10)]  <-- gaps_
//                127 131 137 139 143 149 151 157 163 167 169 173 179 181 187 191 193 197 199 209 211 221 ]
//                                         ^                                                           ^---- last_ = 221
//                                         |__ head_
//
// start_ In the case of window2_ start_ is the last candidate that was generated, and equal to p_ squared.
//        In the case of window1_ it is the factor of the value before the first expected failure,
//        that is a multiple of p_. And since the first expected failure is p_ * 'the next prime after p_',
//        that value is p_^2 and the factor is p_. In other words, in the case of window1_, start_ = p_.
// last_  In the case of window2_ last_ is the last candidate that will be generated after which the gaps
//        start repeating (e.g. after 221 comes 221 + (13 - 11) = 223), where 13 is the first candidate
//        associated with window1_. It is equal to p_ + (2 * 3 * 5 * ... * 'last prime before p_');
//        i.e. 11 + 2 * 3 * 5 * 7 = 11 + 210 = 221.
//        In the case of window1_ it is equal to the start_ of window2_ (p_^2).
// head_  is the factor(window1) / candidate(window2) associated with the last gap in gaps_ (the last gap
//        that was generated). At the moment head_ is returned, tail_ will point to its associated gap.
//        I.e. in the above picture (window2) there are two more tails that can be returned before a new
//        gap has to be generated (+6 and +2).
// tail_  Index into gaps_ of the last gap that was returned (or -1 the next gap to be returned is the first
//        entry in gaps_).

prime_gap_t SieveWindow::read_gap(SieveWindow*& window, bool skip)
{
  if (!prev_sieve_)
  {
    std::cout << "Window " << start_ << " returning +1" << std::endl;
    return 1;
  }

  assert(head_ <= last_);

  if (tail_ >= (ssize_t)gaps_.size())
    std::cout << "Warning: entering window " << start_ << " with tail_ (" << tail_ << ") >= gaps_.size() (" << gaps_.size() << ")." << std::endl;

  if (++tail_ < (ssize_t)gaps_.size())
  {
    std::cout << "Window " << start_ << " returning +" << gaps_[tail_] << " (tail_ = " << tail_ << ")" << std::endl;
    return gaps_[tail_];
  }

  std::cout << "Entered window " << start_ << "; no tail available." << std::endl;
  if (tail_ != gaps_.size())
  {
    std::cout << "SieveWindow is now: " << *this << std::endl;
  }
  assert(tail_ == gaps_.size());

  if (head_ < last_)
  {
    prime_gap_t gap = prev_sieve_->next_gap();
    head_ += gap;
    if (skip && head_ == prev_sieve_->expected_failure())
    {
      prime_gap_t g = prev_sieve_->next_gap();
      gap += g;
      head_ += g;
      prev_sieve_->set_next_expected_failure();
    }
    assert(head_ <= last_);
    gaps_.push_back(gap);
    std::cout << "Window " << start_ << " added +" << gap << " (" << head_ << ") now: " << *this << std::endl;
    assert(!divisible_by_any_prime_below(head_, prev_sieve_->prime_index() + 1));
    assert(tail_ < (ssize_t)gaps_.size());
    return gap;
  }

  // Switch to the beginning of the other window.
  std::cout << "Switching windows from " << start_ << " to " << other_window_->start_ << std::endl;
  window = other_window_;
  window->reset_tail();
  prime_gap_t gap = window->read_gap(window);
  std::cout << "Window " << start_ << " returning +" << gap << " obtained from other window." << std::endl;
  assert(window->tail_ < (ssize_t)window->gaps_.size());
  reset_tail(); // Need for the really small windows.
  if (tail_ >= (ssize_t)gaps_.size())
    std::cout << "Warning: leaving window " << start_ << " with tail_ (" << tail_ << ") >= gaps_.size() (" << gaps_.size() << ")." << std::endl;
  return gap;
}

class PrimeGenerator
{
 private:
  std::vector<prime_t> primes_;
  Sieve* sieve_;

 public:
  PrimeGenerator() : primes_({ 1 }), sieve_(new Sieve(primes_)) { }
  ~PrimeGenerator() { delete sieve_; }

  prime_t next()
  {
    prime_t p = sieve_->next(primes_, sieve_);
    primes_.push_back(p);
    return p;
  }

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, PrimeGenerator const& primes) { primes.print_on(os); return os; }
};

int main()
{
  PrimeGenerator prime_generator;
  std::cout << "gen = " << prime_generator << std::endl;

  for (int i = 0; i < 10; ++i)
  {
    prime_t p = prime_generator.next();
    std::cout << "NEXT PRIME: " << p << std::endl;
    std::cout << "gen = " << prime_generator << std::endl;
  }
}

void SieveWindow::print_on(std::ostream& os) const
{
  os << '{';
  os <<   "start_:" << start_;
  os << ", head_:" << head_;
  os << ", last_:" << last_;
  os << ", tail_:" << tail_;
  os << ", gaps_:{";
  char const* prefix = "";
  for (prime_gap_t gap : gaps_)
  {
    os << prefix << gap;
    prefix = ", ";
  }
  os << "}}";
}

void Sieve::print_on(std::ostream& os) const
{
  os << '{';
  os <<   "pi_:" << pi_;
  os << ", p_:" << p_;
  os << ", window1_:" << window1_;
  os << ", window2_:" << window2_;
  os << ", last_prime_:" << last_prime_;
  os << ", last_offset_:" << last_offset_;
  os << ", expected_failure_:" << expected_failure_;
  os << ", next_prime_squared_:" << next_prime_squared_;
  os << '}';
}

void PrimeGenerator::print_on(std::ostream& os) const
{
  os << '{';
  os <<   "primes_:{";
  char const* prefix = "";
  for (auto prime : primes_)
  {
    os << prefix << prime;
    prefix = ", ";
  }
  os << '}';
  os << ", sieve_:" << *sieve_;
  os << '}';
}
