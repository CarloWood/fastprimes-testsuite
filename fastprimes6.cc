#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <cassert>

using prime_t = uint32_t;
using integer_t = uint64_t;
using prime_index_t = int;
using prime_gap_t = prime_t;

std::array<prime_t, 23> debug_primes = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79 };

class Primes
{
 private:
  std::vector<prime_t> primes_;

 public:
  Primes() : primes_({1, 2})
  {
  }

  prime_t p(prime_index_t n) const
  {
    assert(1 <= n && n < primes_.size());
    return primes_[n];
  }

  prime_t debug_p(prime_index_t n) const
  {
    assert(1 <= n);
    if (n < primes_.size())
      return primes_[n];
    return debug_primes[n];
  }

  void store(prime_t p)
  {
    assert(p > primes_.back());
    primes_.push_back(p);
  }

  integer_t primorial(prime_index_t n) const;
  prime_index_t repeat(prime_index_t n) const;
  prime_index_t square(prime_index_t n) const;

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, Primes const& primes) { primes.print_on(os); return os; }
};

// Returns pₙ# = p₁⋅p₂⋅p₃… pₙ
// Note that p₀ := 1
integer_t Primes::primorial(prime_index_t n) const
{
  if (n < 1)
    return 1;
  return p(n) * primorial(n - 1);
}

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
  Primes& primes_;
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

 public:
  Sieve(Primes& primes) :
    primes_(primes),
    prev_sieve_(nullptr),
    begin_(1),
    prev_prime_(1),
    end_(2),
    n_(1),
    factor_n_(1),
    base_(0),
    factor_base_(0),
    primorial_(1),
    expected_failure_(0)
  { }

  Sieve(Sieve& prev_sieve) :
    primes_(prev_sieve.primes()),
    prev_sieve_(&prev_sieve),
    begin_(prev_sieve.begin() + 1),
    prev_prime_(prev_sieve.primes().p(begin_ - 1)),
    end_(begin_ + prev_sieve.primes().repeat(begin_)),
    n_(begin_),
    factor_n_(begin_),
    base_(0),
    factor_base_(0),
    primorial_(primes_.primorial(begin_ - 1)),
    expected_failure_(prev_prime_ * prev_sieve_->next_factor())
  { }

  Primes& primes() const
  {
    return primes_;
  }

  prime_index_t begin() const
  {
    return begin_;
  }

  bool done() const
  {
    return expected_failure_ == 0;
  }

  // For example begin_ = 3, end_ = 5
  //
  //  5 <-- p₃ <-- n
  // +2
  //  7
  // +4
  // 11 <-- p₅
  // +2
  // 13
  // +4
  // 17
  // +2
  // 19
  // +4
  // 23
  // +2
  // 25
  //
  integer_t next()
  {
    if (!expected_failure_)
    {
      integer_t next = base_ + primes_.p(n_);
      if (++n_ == end_)
      {
        base_ += primorial_;
        n_ = begin_;
      }
      return next;
    }

#if 1
    static int count = 0;
    ++count;
//    std::cout << "count = " << count << std::endl;
#endif

    bool new_prime = prev_sieve_->done();
    integer_t next = prev_sieve_->next();

    if (!base_)
    {
      base_ = primorial_;
      if (new_prime)
        primes_.store(next);
      return next;
    }

    if (next == expected_failure_)
    {
      new_prime = prev_sieve_->done();
      next = prev_sieve_->next();
      if (next - primes_.p(begin_) == primorial_)
        expected_failure_ = 0;
      else
        expected_failure_ = prev_prime_ * prev_sieve_->next_factor();
    }

    if (new_prime)
      primes_.store(next);
    return next;
  }

  integer_t next_factor()
  {
    integer_t next = factor_base_ + primes_.p(factor_n_);
    if (++factor_n_ == end_)
    {
      factor_base_ += primorial_;
      factor_n_ = begin_;
    }
    return next;
  }

  void print_on(std::ostream& os) const;
  friend std::ostream& operator<<(std::ostream& os, Sieve const& sieve) { sieve.print_on(os); return os; }
};

int main()
{
  {
    Primes primes;
    Sieve s2(primes);
    std::cout << s2 << std::endl;
    for (int i = 0; i < 20; ++i)
      std::cout << s2.next() << ' ';
    std::cout << std::endl;
  }

  Primes primes;
  Sieve s2(primes);
  std::cout << "prime: " << s2.next() << std::endl;
  Sieve s3(s2);
  std::cout << "prime: " << s3.next() << std::endl;
  Sieve s4(s3);
  std::cout << s4 << std::endl;
  for (int i = 0; i < 20; ++i)
    std::cout << s4.next() << ' ';
  std::cout << std::endl;

#if 0
  Sieve s4(s3);
  std::cout << s4 << std::endl;
  for (int i = 0; i < 20; ++i)
    std::cout << s4.next() << ' ';
  std::cout << std::endl;

  Sieve s5(s4);
  std::cout << s5 << std::endl;
  for (int i = 0; i < 20; ++i)
    std::cout << s5.next() << ' ';
  std::cout << std::endl;
#endif
}

void Sieve::print_on(std::ostream& os) const
{
  os << '{';
  os << "primes_:" << primes_ <<
      ", prev_sieve_:" << prev_sieve_ <<
      ", begin_:" << begin_ << " (" << primes_.debug_p(begin_) << ")"
      ", prev_prime_:" << prev_prime_ <<
      ", end_:" << end_ << " (" << primes_.debug_p(end_) << ")"
      ", n_:" << n_ <<
      ", factor_n_:" << factor_n_ <<
      ", base_:" << base_ <<
      ", factor_base_:" << factor_base_ <<
      ", primorial_:" << primorial_ <<
      ", expected_failure_:" << expected_failure_;
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
