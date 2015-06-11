#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include "utils.h"

void printGraph(Graph* g) {
    printf("Graph: %d nodes\n", g->nodesCount);
    for (int i = 1; i <= g->nodesCount; i++) {
		printf("Node %d --> [%d]: ", i, g->outDegrees[i]);
		for (int j = 0; j < g->outDegrees[i]; j++) {
			printf("%d ", g->outEdges[i][j]);
		}
		printf("\n");
		printf("Node %d <-- [%d]: ", i, g->inDegrees[i]);
		for (int j = 0; j < g->inDegrees[i]; j++) {
			printf("%d ", g->inEdges[i][j]);
		}
		printf("\n");
	}
    if (g->parents != NULL) {
    	printf("Graph parents:\n");
    	for (int i = 1; i <= g->nodesCount; i++) {
    		printf("%d ", g->parents[i]);
    	}
    	printf("\n");
    }
    if (g->ordering != NULL) {
		printf("Graph ordering:\n");
		for (int i = 1; i <= g->nodesCount; i++) {
			printf("%d ", g->ordering[i]);
		}
		printf("\n");
    }
    if (g->nodesMapping != NULL) {
		printf("Graph nodesMapping:\n");
		for (int i = 1; i <= g->nodesCount; i++) {
			printf("%d ", g->nodesMapping[i]);
		}
		printf("\n");
	}
}

int dfs(int node, int nodeParent, int nextId, Graph* graph, int* numbering, bool viaReverseEdge) {
    int currentId = nextId;
    numbering[node] = viaReverseEdge ? -currentId : currentId;
    graph->parents[node] = nodeParent;
    for (int i = 0; i < graph->outDegrees[node]; i++) {
        if (numbering[graph->outEdges[node][i]] == 0) {
            currentId = dfs(graph->outEdges[node][i], node, currentId + 1, graph, numbering, false);
        }
    }
    for (int i = 0; i < graph->inDegrees[node]; i++) {
        if (numbering[graph->inEdges[node][i]] == 0) {
            currentId = dfs(graph->inEdges[node][i], node, currentId + 1, graph, numbering, true);
        }
    }
    return currentId;
}

void makeGraphOrdering(Graph * graph) {
	int* numbering = safeMalloc((graph->nodesCount + 1) * sizeof(int));
	memset(numbering, 0, (graph->nodesCount + 1) * sizeof(int));
	dfs(1, -1, 1, graph, numbering, false);
	for (int i = 1; i <= graph->nodesCount; i++) {
		if (numbering[i] > 0) {
			graph->ordering[numbering[i]] = i;
		} else {
			graph->ordering[-numbering[i]] = -i;
	    }
	}
	free(numbering);
}

Graph* readGraph(FILE* file) {
    int neigh, node, neighCount, neighCounter;
    int *edges;
    char line[64];

    Graph* g = safeMalloc(sizeof(Graph));
    fscanf(file, "%d\n", &g->nodesCount);
    g->ordering = safeMalloc((g->nodesCount + 1) * sizeof(int));
    g->parents = safeMalloc((g->nodesCount + 1) * sizeof(int));
    g->inEdges = safeMalloc((g->nodesCount + 1) * sizeof(int*));
    g->outEdges = safeMalloc((g->nodesCount + 1) * sizeof(int*));
    g->outDegrees = safeMalloc((g->nodesCount + 1) * sizeof(int));
    g->inDegrees = safeMalloc((g->nodesCount + 1) * sizeof(int));
    g->nodesMapping = NULL;
    memset(g->outDegrees, 0, (g->nodesCount + 1) * sizeof(int));
    memset(g->inDegrees, 0, (g->nodesCount + 1) * sizeof(int));
    while ((fgets(line, sizeof line, file) != NULL) && !isLineEmpty(line)) {
    	sscanf(line, "%d %d", &node, &neighCount);
        neighCounter = 0;
        g->outDegrees[node] = neighCount;
        edges = safeMalloc(neighCount * sizeof(int));
        while (neighCounter < neighCount) {
        	fgets(line, sizeof line, file);
            sscanf(line, "%d", &neigh);
            edges[neighCounter] = neigh;
            g->inDegrees[neigh]++;
            neighCounter++;
        }
        g->outEdges[node] = edges;
    }
    for (int i = 1; i <= g->nodesCount; i++) {
    	if (g->inDegrees[i] > 0) {
    		g->inEdges[i] = safeMalloc(g->inDegrees[i] * sizeof(int));
    	} else {
    		g->inEdges[i] = NULL;
    	}
    }
    int * counters = safeMalloc((g->nodesCount + 1) * sizeof(int));
    memset(counters, 0, (g->nodesCount + 1) * sizeof(int));
    for (int i = 1; i <= g->nodesCount; i++) {
    	for (int j = 0; j < g->outDegrees[i]; j++) {
    		int target = g->outEdges[i][j];
    		g->inEdges[target][counters[target]++] = i;
    	}
    }
    free(counters);
    makeGraphOrdering(g);
    return g;
}

void freeGraph(Graph *g) {
	for (int i = 1; i <= g->nodesCount; i++) {
		if (g->inDegrees[i] > 0) {
			free(g->inEdges[i]);
		}
		if (g->outDegrees[i] > 0) {
			free(g->outEdges[i]);
		}
	}
	if (g->ordering != NULL) {
		free(g->ordering);
	}
	if (g->parents != NULL) {
		free(g->parents);
	}
	if (g->nodesMapping != NULL) {
		free(g->nodesMapping);
	}
	free(g->inEdges);
	free(g->outEdges);
	free(g->inDegrees);
	free(g->outDegrees);
	free(g);
}
