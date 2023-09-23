#include <iostream>
#include <array>
#include <iomanip>
#include <tuple>
#include <cassert>

int modInverse(int n, int N)
{
  int a = n, b = N;
  int x0 = 1, y0 = 0, x1 = 0, y1 = 1;
  while (b > 0)
  {
    int q = a / b;
    std::tie(a, b) = std::make_tuple(b, a % b);
    std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
    std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
  }
  // Ensure the result is positive.
  return (x0 + N) % N;
}

int main()
{
  std::array<int, 8> fr = { 7, 11, 13, 17, 19, 23, 29, 31 };

  int offset = 0;
  int const N = 13;
  int const M = 17;
  int const vo = 2 * 3 * 5;
  int const voi = modInverse(vo, N);
  for (int r = 0; r < vo; ++r)
  {
    std::cout << std::setw(2) << r << " : ";
    for (int n : fr)
    {
      int pr = (-n * voi) % N;
      std::cout << '|';
      if ((offset + n) % N == 0)
      {
        assert((pr - r) % N == 0);
        std::cout << '!';
      }
      else
        std::cout << ' ';
      if ((offset + n) % M == 0)
        std::cout << '#';
      else
        std::cout << ' ';
      std::cout << std::setw(3) << std::setfill(' ') << ((offset + n) /*% N*/);
    }
    offset += vo;
    std::cout << std::endl;
  }
}
