#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

bool isLineEmpty(char* line) {
	int i = 0;
	while(line[i] != '\0') {
		if (!isspace(line[i])) {
			return false;
		}
		i++;
	}
	return true;
}

void * safeMalloc(int size) {
	void * memory = malloc(size);
	if (memory == NULL) {
		fprintf(stderr, "ERROR: malloc returned NULL");
		MPI_Finalize();
		exit(1);
	}
	return memory;
}

int findIndex(int * tab, int len, int val) {
	int left = 1, right = len, pivot;
	while(left < right) {
		pivot = (left + right) / 2;
		if (tab[pivot] < val) {
			left = pivot + 1;
		} else {
			right = pivot;
		}
	}
	if (tab[left] == val) {
		return left;
	}
	return -1;
}

bool arrayContains(int * tab, int len, int val) {
    for (int i = 0; i < len; i++) {
        if (tab[i] == val) {
            return true;
        }
    }
    return false;
}
