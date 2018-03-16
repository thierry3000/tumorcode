#include <limits>
#include <array>
#include <iostream>

inline float NANf() { return std::numeric_limits<float>::quiet_NaN(); }

int main()
{
  std::array<float, 42> foo;
  for(int i = 0;i<42;++i)
  {
    foo[i] = NANf();
  }
  for(int i = 0;i<42;++i)
  {
    if( foo[i] > 1.e-3) // compare with nan is always equals false
    {
      std::cout << 2*foo[i] << std::endl;
    }
  }
  std::cout << foo[21] << std::endl;
    
}
