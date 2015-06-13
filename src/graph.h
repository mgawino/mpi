#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdbool.h>

typedef struct {
    int ** outEdges;
    int * outDegrees;
    int ** inEdges;
    int * inDegrees;
    int * ordering;
    int * parents;
    int * nodesMapping;
    int nodesCount;
} Graph;

Graph* readGraph(FILE* file);
void freeGraph(Graph *g);

#endif
