/** compile with
 * gcc -g -rdynamic ./test.c -o test
 */
#include <stdio.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


void handler(int sig) {
  void *array[42];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 42);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(sig);
}

void baz() {
 int *foo = (int*)-1; // make a bad pointer
  printf("%d\n", *foo);       // causes segfault
}

void bar() { baz(); }
void foo() { bar(); }

void numerical_baz()
{
 int number = 42; //use the most famous number
 int number2 = 0;
 number = number/number2;
 std::cout << "Famous number: "<< number << std::endl;
}
void numerical_bar(){numerical_baz();}
void numerical_foo(){numerical_bar();}

int main(int argc, char **argv) {
  signal(SIGSEGV, handler);   // install our handler
  signal(SIGFPE, handler);
  //foo(); // this will call foo, bar, and baz.  baz segfaults.
  numerical_foo(); // this will call numerical_foo, numerical_bar, numerical_baz
}
