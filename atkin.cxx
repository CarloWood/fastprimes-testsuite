#include <bitset>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

void writeVectorToDisk(std::string const& filename, std::vector<uint32_t> const& vec)
{
  std::ofstream ofs(filename, std::ios::binary);

  if (!ofs.is_open())
  {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  // Write the size of the vector first,
  size_t size = vec.size();
  ofs.write(reinterpret_cast<char const*>(&size), sizeof(size));

  // Then write the vector data.
  ofs.write(reinterpret_cast<char const*>(vec.data()), vec.size() * sizeof(uint32_t));

  ofs.close();
}

template<uint64_t lmt>
struct Sieve
{
  std::bitset<lmt / 3>* sieve_;
  mutable uint64_t read_access_{0};
  uint64_t read_write_access_{0};

  Sieve() : sieve_(new std::bitset<lmt / 3>) { }
  ~Sieve() { delete sieve_; }

  void flip(uint64_t n)
  {
    assert(n % 2 != 0 && n % 3 != 0);
    sieve_->flip(n / 3);
    ++read_write_access_;
  }

  void reset(uint64_t n)
  {
    assert(n % 2 != 0 && n % 3 != 0);
    sieve_->reset(n / 3);
    ++read_write_access_;
  }

  bool test(uint64_t n) const
  {
    assert(n % 2 != 0 && n % 3 != 0);
    ++read_access_;
    return sieve_->test(n / 3);
  }

  void print_on(std::ostream& os) const
  {
    os << "read_access_: " << read_access_ << "; read_write_access_: " << read_write_access_;
  }
};

std::vector<uint32_t> SieveOfAtkin()
{
  std::vector<uint32_t> primes;
//  primes.reserve(189961812);

  constexpr uint64_t lmt = 1000; //4000000000;
  Sieve<lmt> sieve;

  for (uint64_t a = 1; a * a < lmt; ++a)
    for (uint64_t b = 1; b * b < lmt; ++b)
    {
      // Main part of Sieve of Atkin.
      uint64_t n = (4 * a * a) + (b * b);
      if (n <= lmt && (n % 12 == 1 || n % 12 == 5))
      {
        if (n == 65)
          std::cout << "WE GET HERE: " << a << ", " << b << std::endl;
        sieve.flip(n);
      }
      n = (3 * a * a) + (b * b);
      if (n <= lmt && n % 12 == 7)
      {
        sieve.flip(n);
      }
      n = (3 * a * a) - (b * b);
      if (a > b && n <= lmt && n % 12 == 11)
      {
        sieve.flip(n);
      }
    }

  uint64_t r = 5, r2 = 25;
  while (r2 < lmt - 2)
  {
    if (sieve.test(r))
    {
      uint64_t i = r2;
      while (i + 4 * r2 < lmt)
      {
        sieve.reset(i);
        i += 4 * r2;
        sieve.reset(i);
        i += 2 * r2;
      }
      if (i < lmt)
        sieve.reset(i);
    }
    ++r;
    r2 += 4 * r;
    ++r;
    if (sieve.test(r))
    {
      uint64_t i = r2;
      while (i + 4 * r2 < lmt)
      {
        sieve.reset(i);
        i += 4 * r2;
        sieve.reset(i);
        i += 2 * r2;
      }
      if (i < lmt)
        sieve.reset(i);
    }
    r += 2;
    r2 += 8 * r;
    r += 2;
  }
  if (r2 < lmt && sieve.test(r))
    sieve.reset(r2);

  primes.push_back(2);
  primes.push_back(3);
  uint64_t x = 5;
  do
  {
    if (sieve.test(x))
      primes.push_back(x);
    x += 2;
    if (sieve.test(x))
      primes.push_back(x);
    x += 4;
  }
  while (x < lmt - 2);
  if (x < lmt && sieve.test(x))
    primes.push_back(x);

  return primes;
}

int main()
{
  auto primes = SieveOfAtkin();

  std::cout << "The " << primes.size() << "-th prime number is " << primes.back() << std::endl;

  for (auto p : primes)
    std::cout << p << std::endl;
//  std::cout << "Writing result to file primes_till_4000000000..." << std::endl;
//  writeVectorToDisk("primes_till_4000000000", primes);
}
