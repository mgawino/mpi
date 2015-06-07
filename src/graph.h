#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdbool.h>

typedef struct {
    int ** outEdges;
    int ** outEdgesProc;
    int * outDegrees;
    int ** inEdges;
    int ** inEdgesProc;
    int * inDegrees;
    int * ordering;
    int * parents;
    int * nodesMapping;
    int nodesCount;
} Graph;

Graph* readGraph(FILE* file);
void printGraph(Graph* g);
void freeGraph(Graph *g);

#endif
