#ifndef UTILS_H_
#define UTILS_H_

#include <stdbool.h>

bool isLineEmpty(char* line);
void * safeMalloc(int size);
int findIndex(int * tab, int len, int val);

#endif
