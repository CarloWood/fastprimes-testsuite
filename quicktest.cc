#include <iostream>
#include <vector>
#include <cassert>

std::vector<int> primes = {2, 3, 5, 7};

int prime_to_index(int p)
{
  for (int i = 0; i < primes.size(); ++i)
    if (primes[i] == p)
      return i;
  return -1;
}

int main()
{
  int lk = -1;
  int k = 11;
  int z = 10000000;

  for (int i = 1; (int)primes.size() < z; ++i)
  {
    int lj = i;

    int p1 = primes[i];
    int p2 = primes[i + 1];
    int expected_failure = p1 * p2;
    int failure_index = i + 1;
    int const stop = p2 * p2;

    for (; k <= stop; k += 2)
    {
      bool skip = false;
      for (int prime : primes)
      {
        if (prime >= p1)
          break;
        if (k % prime == 0)
        {
          skip = true;
          break;
        }
      }
      if (skip)
        continue;

      if (lk != -1)
        std::cout << '+' << (k - lk) << ' ';

      bool isComposite = k == expected_failure || k == stop;
      if (!isComposite)
        primes.push_back(k);
      else
        std::cout << '(';

//      std::cout << k;
      lk = k;

      if (isComposite)
        std::cout << ')';

      if (k == expected_failure)
      {
        std::cout << " = " << p1 << " * " << (k / p1);
        int j = prime_to_index(k / p1);
        assert(j == lj + 1);
        lj = j;

        expected_failure = p1 * primes[++failure_index];
      }
      if (k == stop)
        std::cout << " = " << p2 << "^2";

      std::cout << std::endl;
    }
    std::cout << "--------\n";
  }
  std::cout << "The " << z << "-th prime is " << primes[z - 1] << std::endl;
}
