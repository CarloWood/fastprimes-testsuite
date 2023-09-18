#include <iostream>

class Prime {
private:
  long currentPrime;
  unsigned int* p;

public:
  Prime(unsigned int* ptr) : currentPrime(2), p(ptr) { }

  long operator*(void) { return currentPrime; }
  Prime operator++(void)
  {
    unsigned int o;
    while (((o = *p++) & 0x80000000))
      currentPrime += (o | 0x7fffffff);
    currentPrime += o;
    return *this;
  }
  Prime operator++(int)
  {
    Prime p = *this;
    ++(*this);
    return p;
  }
};

using std::cout;

unsigned char control_array[100000];

int main()
{
  cout << "Fast PRIMES v2.0, by Carlo Wood <carlo@alinoe.com>\n" << std::endl;

  unsigned int const arraySize = 100;
  unsigned int offsets[arraySize];
  unsigned int* offsetsEnd = &offsets[arraySize];
  long number = 2;
  offsets[0] = 1;
  offsets[1] = 1;
  offsets[2] = 2;
  unsigned int* loopStart = &offsets[0];
  unsigned int* loopEnd = &offsets[1];
  long sum = 0, loopSum = 1;
  unsigned int offset = 1;
  unsigned int* front = offsets;
  unsigned int* p = loopStart;
  Prime loopSumFactor(&offsets[1]);
  Prime nextPrime(&offsets[1]);
  long prime = *nextPrime++;
  long factor = prime * prime;
  long pfactor = prime;

  for (;;)
  {
    if (number != factor && (number == pfactor || number % pfactor != 0))
    {
      *front++ = offset;
      cout << "\nPrime: " << number << "\t";
      if (control_array[number])
      {
	std::cout << "That is NOT a prime!" << std::endl;
	return 1;
      }
      for (int i = number; i < sizeof(control_array); i += number)
	control_array[i] = 1;
      if (front == offsetsEnd)
        break;
      offset = 0;
    }
    else
    {
      if (number == factor)
      {
	prime = *nextPrime++;
	factor = prime * prime;
      }
      cout << " (hit multiple of " << pfactor;
      pfactor = prime;
      cout << "; new prime = " << pfactor << ")";
    }
    cout << " + " << *p;
    number += *p;
    offset += *p;
    if ((sum += *p++) == loopSum)
    {
      *front = offset;
      loopStart = loopEnd;
      loopEnd = front + 1;
      loopSum *= *loopSumFactor++;
      sum = 0;
      cout << "\n\nArray offsets repeating: [";
      for (unsigned int* q = loopStart; q < loopEnd; ++q)
        cout << " +" << *q;
      cout << " ]" << std::endl;
      p = loopStart;
      continue;
    }
    if (p == loopEnd)
      p = loopStart;
  }
  cout << std::endl;
}
