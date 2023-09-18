#include <iostream>
#include <vector>
#include <list>

using namespace std;

class offset {
private:
  unsigned int M_offset;
public:
  offset(unsigned int __offset) : M_offset(__offset) { }
  friend ostream& operator<<(ostream& os, offset const& __offset) { os << '+' << __offset.M_offset; return os; }
  friend unsigned int operator+(offset const& __offset1, offset const& __offset2) { return __offset1.M_offset + __offset2.M_offset; }
  friend unsigned int operator+(unsigned int __offset1, offset const& __offset2) { return __offset1 + __offset2.M_offset; }
  friend unsigned int operator+(offset const& __offset1, unsigned int __offset2) { return __offset1.M_offset + __offset2; }
  friend unsigned int operator-(offset const& __offset1, offset const& __offset2) { return __offset1.M_offset - __offset2.M_offset; }
  friend unsigned int operator-(unsigned int __offset1, offset const& __offset2) { return __offset1 - __offset2.M_offset; }
  friend unsigned int operator-(offset const& __offset1, unsigned int __offset2) { return __offset1.M_offset - __offset2; }
  friend unsigned int operator*(offset const& __offset1, offset const& __offset2) { return __offset1.M_offset * __offset2.M_offset; }
  friend unsigned int operator*(unsigned int __offset1, offset const& __offset2) { return __offset1 * __offset2.M_offset; }
  friend unsigned int operator*(offset const& __offset1, unsigned int __offset2) { return __offset1.M_offset * __offset2; }
};

class offset_table {
private:
  vector<offset> M_offset_table;
  unsigned int base_prime;
  unsigned int M_offset_table_size;
  // Generation variables:
  offset_table const& previous_offset_table;
  vector<offset>::const_iterator source;
  unsigned int number;
public:
  offset_table(void) :
      base_prime(2),
      M_offset_table_size(1),
      previous_offset_table(*static_cast<offset_table*>(NULL)),
      number(0)
  {
    M_offset_table.push_back(1U);
  }
  offset operator[](unsigned int index)
  {
    if (index >= M_offset_table_size)
      index %= M_offset_table_size;
    while(index >= M_offset_table.size())
    {
      cout << '\n' << base_prime << ": index " << index << " > " << M_offset_table.size() << '\n';
      if (source == previous_offset_table.M_offset_table.end())
        source = previous_offset_table.M_offset_table.begin();
      number = number + *source;
      cout << "Adding: " << *source << "; number now " << number << '\n';
      if (number % previous_offset_table.base_prime != 0)
	M_offset_table.push_back(*source);
      ++source;
    }
    return M_offset_table[index];
  }
  offset_table(offset_table& __previous_offset_table) :
      base_prime(__previous_offset_table.base_prime + __previous_offset_table[1]),
      M_offset_table_size(__previous_offset_table.M_offset_table_size * __previous_offset_table[0]),
      previous_offset_table(__previous_offset_table),
      source(__previous_offset_table.M_offset_table.begin() + 1),
      number(base_prime)
  {
    M_offset_table.push_back(base_prime - 1);
  }
  friend ostream& operator<<(ostream& os, offset_table const& __offset_table)
  {
    os << __offset_table.base_prime << ": {";
    for (vector<offset>::const_iterator i(__offset_table.M_offset_table.begin()); i != __offset_table.M_offset_table.end(); ++i)
      os << ' ' << *i;
    os << " }, size = " << __offset_table.M_offset_table_size;
    if (__offset_table.number > 0)
      cout << "; number = " << __offset_table.number;
    return os;
  }
};

#if 0
class singleton_primes {
private:
  list<offset_table> M_offset_table_list;
};

class primes {
private:
  singleton_primes M_singleton_primes;
};
#endif


int main(void)
{
  offset_table t1;
  cout << t1 << endl;
  offset_table t2(t1);
  cout << t2 << endl;
  offset_table t3(t2);
  cout << t3 << endl;
  offset_table t4(t3);
  cout << t4 << endl;
  offset_table t5(t4);
  cout << t5 << endl;
  cout << '\n';
  cout << t1 << endl;
  cout << t2 << endl;
  cout << t3 << endl;
  cout << t4 << endl;
  cout << t5 << endl;

  return 0;
}
