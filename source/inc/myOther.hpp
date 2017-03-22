#ifndef myOther_h
#define myOther_h 1

#include <cstddef>



// reportProgress

void reportProgressInit();

void reportProgress(double percent);
void reportProgress(std::size_t done, std::size_t length);
void reportProgress(long long int done, long long int length);

void reportProgressFinish();


#endif  // myOther_h
