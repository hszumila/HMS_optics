#include "myOther.hpp"

#include <cstdio>


// Implementation of reportProgress.

void reportProgressInit() {
  printf("%5.1f%%\b\b\b\b\b\b", 0.0);
  fflush(stdout);
}


void reportProgress(double percent) {
  printf("%5.1f\b\b\b\b\b", percent);
  fflush(stdout);
}


void reportProgress(std::size_t done, std::size_t length) {
  double percent = static_cast<double>(done) / static_cast<double>(length) * 100.0;
  printf("%5.1f\b\b\b\b\b", percent);
  fflush(stdout);
}


void reportProgress(long long int done, long long int length) {
  double percent = static_cast<double>(done) / static_cast<double>(length) * 100.0;
  printf("%5.1f\b\b\b\b\b", percent);
  fflush(stdout);
}


void reportProgressFinish() {
  printf("%5.1f%%\n", 100.0);
}
