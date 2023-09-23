#include <iostream>
#include <cmath>
#include <iomanip>

long double li(long double x)
{
 static constexpr long double EulerMascheroni = 0.57721566490153286;

 long double logx = std::log(x);

 long double result = EulerMascheroni + std::log(logx);

 // correction sum.
 long double c = 1;
 for (long long n = 1; n < 20; ++n)
 {
   c *= logx / n;
   result += c / n;
 }

 return result;
}

// (((x - 2) / (1 + (x - 2)) - 0.5*(1 -(x-5)/(1+(x-5)))^2.6)

int main()
{
  std::cout << "li(1000000) = " << std::setprecision(18) << li(1000000) << std::endl;
}
