#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>
#include <deque>

// Prime index.
// Index  Prime
// 0      2
// 1      3
// 2      5
// ...
using prime_index_t = int;

// A prime.
using prime_t = uint32_t;

// A large integer.
using integer_t = uint64_t;     // Must be able to contain the square of the largest prime that we use.

// A distance between two primes.
using prime_gap_t = prime_t;

// List of primes.
using primes_t = std::vector<prime_t>;

// All primes generated so far.
primes_t primes;

int number_of_gaps = 0;

// A class representing a factor for integers that are free
// of prime cofactors up till and including a given prime.
class PrimeBoundCofactorFreeIntegers
{
 private:
  prime_index_t index_;           // The index of the highest prime that is not a factor of these integers.
  prime_t p_;                     // The highest prime that is not a factor of these integers.
  prime_t first_;                 // The first integer in this set (this will be a prime).
  std::vector<prime_gap_t> gaps_; // The gaps between the integers as generated so far.
  integer_t sum_;                 // The sum of all gaps until gaps_ starts to repeat.

  PrimeBoundCofactorFreeIntegers* prev_;        // Points to the PrimeBoundCofactorFreeIntegers of the previous prime p_, or nullptr if p_ == 2.
  integer_t current_value_;                     // The value that will be returned by start(), or the last value that was returned by next().
  int gap_index_;                               // Index into gaps_ of the next gap that should be added to current_value_.
  integer_t gaps_sum_;                          // The sum of all gaps in gaps_ so far.

 public:
  // The default constructor constructs a PrimeBoundCofactorFreeIntegers for the prime 2: all odd integers starting at 3.
  // Aka gaps_ contains { +2 }.
  PrimeBoundCofactorFreeIntegers() :
    index_(-1),                                 // The index of the next (first) prime will be 0.
    p_(0),                                      // Not associated with a prime.
    first_(2),                                  // The first integer larger than 1 is 2.
    gaps_(1, 1),                                // All integers are 1 apart, so the gaps are +1 repeated (or { +1 }).
    sum_(1),                                    // The sum of all gaps after which they will repeat.
    prev_(nullptr),                             // There is no previous PrimeBoundCofactorFreeIntegers.
    current_value_(first_),                     // The first integer to be returned, which is the first prime.
    gap_index_(0),                              // The index of the next gap that should be added.
    gaps_sum_(1)                                // The sum of all gaps in gaps_ so far.
  { }

  // Construct a subsequent PrimeBoundCofactorFreeIntegers.
  PrimeBoundCofactorFreeIntegers(PrimeBoundCofactorFreeIntegers* prev, prime_t first) :
    index_(prev->index() + 1),
    p_(primes[index_]),
    first_(first),
    sum_(prev->sum() * p_),
    prev_(prev),
    current_value_(first_),
    gap_index_(-1),
    gaps_sum_(0)
  { }

  void clear()
  {
    std::cout << "Clearing sieve " << index_ << std::endl;
    std::cout << '{';
    for (auto g : gaps_)
      std::cout << ' ' << g;
    std::cout << " }" << std::endl;
    gaps_.clear();
  }

  integer_t expected_failure() const
  {
    return primes[index_ + 1] * primes[index_ + 1];
  }

  prime_index_t index() const
  {
    return index_;
  }

  integer_t sum() const
  {
    return sum_;
  }

  prime_t start()
  {
    return current_value_;
  }

  integer_t next()
  {
    ++gap_index_;
    if (gaps_sum_ == sum_)
    {
      if (gap_index_ == gaps_.size())
      {
        gap_index_ = 0;
        if (prev_)
        {
          prev_->clear();               // Free up some memory.
          prev_ = nullptr;
        }
      }
    }
    else
    {
      integer_t next_value = prev_->next();
      if (next_value % p_ == 0)
        next_value = prev_->next();
      prime_gap_t next_gap = next_value - current_value_;
      gaps_.push_back(next_gap);
      ++number_of_gaps;
      gaps_sum_ += next_gap;
    }
    current_value_ += gaps_[gap_index_];
    return current_value_;
  }
};

// A list of all sieves, by prime_index_t.
std::deque<PrimeBoundCofactorFreeIntegers> sieves;

int main()
{
  std::cout << "Fast PRIMES v3.0, by Carlo Wood <carlo@alinoe.com>\n" << std::endl;

  // Bootstrap.
  sieves.emplace_back();
  PrimeBoundCofactorFreeIntegers* current_sieve = &sieves.back();
  integer_t expected_failure = 4;

  for (int s = 0; s < 10000; ++s)
  {
    integer_t next_prime = current_sieve->start();
    while (next_prime != expected_failure)
    {
      primes.push_back(next_prime);
      if (primes.size() % 10000 == 0)
        std::cout << "The " << primes.size() << "-th prime is " << *primes.rbegin() << ". s = " << s << "; number_of_gaps = " << number_of_gaps << ".\n";
      next_prime = current_sieve->next();
    }
    sieves.emplace_back(current_sieve, current_sieve->next());
    current_sieve = &sieves.back();
    expected_failure = current_sieve->expected_failure();
  }
}
