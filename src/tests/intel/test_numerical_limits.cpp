
#include <iostream>
#include <limits>
using namespace std;

int main() 
{
    cout << "Hello, World!" << endl;
    int foo = std::numeric_limits<int>::max();
    cout << foo << endl;
    int bar = foo + 1;
    cout << bar << endl;
    bar = std::numeric_limits<int>::max();
    int foobar = std::min(foo,bar);
    cout << foobar << endl;
    cout << "bla" << endl;
    return 0;
}
