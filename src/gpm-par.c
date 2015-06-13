#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <assert.h>
#include "graph.h"
#include "utils.h"

#define ROOT 0
#define MAX_MATCH_SIZE 10

#define MPI_PROC_IN_OUT_INFO 70
#define MPI_NODE_IN_OUT_INFO 71
#define MPI_NODE_OUT_NEIGH 72
#define MPI_NODE_IN_INFO 73
#define MPI_NODE_IN_NEIGH 74
#define MPI_END_DISTRIBUTION 75
#define MPI_END_COMPUTATION 76
#define MPI_MATCH_FOUND 77
#define MPI_ASK 78
#define MPI_ASKED_OUT 79
#define MPI_ASKED_IN 80
#define MPI_END 81

typedef struct {
    int matchedNodesCount;
    int matchedNodes[MAX_MATCH_SIZE + 1];
} Match;

typedef struct {
	int inDegree;
	int outDegree;
} NodeInfo;

typedef struct {
	int node;
	int val;
} NodeComp;

typedef struct {
	int nodesWithOutEdges;
	int nodesCount;
	int memory;
} ProcInfo;

void checkArgs(int argc, char* argv[], FILE ** inFile, FILE ** outFile) {
    if (argc != 3) {
    	fprintf(stderr, "ERROR: files not specified\n");
        MPI_Finalize();
        exit(1);
    }
    *inFile = fopen(argv[1], "r");
    if (*inFile == NULL) {
    	fprintf(stderr, "ERROR: file not found %s\n", argv[1]);
    	fclose(*inFile);
    	MPI_Finalize();
    	exit(1);
    }
    *outFile = fopen(argv[2], "w+");
    if (*outFile == NULL) {
		fprintf(stderr, "ERROR: file not found %s\n", argv[2]);
		fclose(*outFile);
		MPI_Finalize();
		exit(1);
    }
}

void writeMatch(int * match, int len, FILE * outFile) {
	for (int i = 1; i < len; i++) {
		fprintf(outFile, "%d ", match[i]);
	}
	fprintf(outFile, "%d\n", match[len]);
}

bool matchContains(int node, Match* match) {
    for (int i = 1; i <= match->matchedNodesCount; i++) {
       if (match->matchedNodes[i] == node) {
           return true;
       }
    }
    return false;
}

bool checkNodeMatches(int graphNode, int patternNode, int * graphOutEdges, int graphOutDegree,
					  int * graphInEdges, int graphInDegree, Graph* pattern, Match* match, int rank) {
    if (matchContains(graphNode, match)) {
        return false;
    }

    for (int i = 0; i < pattern->outDegrees[patternNode]; i++) {
        int patternNeigh = pattern->outEdges[patternNode][i];
        int graphNeigh = match->matchedNodes[patternNeigh];
        if (patternNeigh <= match->matchedNodesCount &&
        	!arrayContains(graphOutEdges, graphOutDegree, graphNeigh)) {
            return false;
        }
    }

    for (int i = 0; i < pattern->inDegrees[patternNode]; i++) {
		int patternNeigh = pattern->inEdges[patternNode][i];
		int graphNeigh = match->matchedNodes[patternNeigh];
		if (patternNeigh <= match->matchedNodesCount &&
			!arrayContains(graphInEdges, graphInDegree, graphNeigh)) {
			return false;
		}
    }
    return true;
}

void trySendEdges(Graph * graph) {
	int send = 0, info[2];
	MPI_Status status;
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ASK, MPI_COMM_WORLD, &send, &status);
	if (send) {
		// graphNode, 1 - all edges, 2 - out edges, 3 - in edges
		MPI_Recv(info, 2, MPI_INT, status.MPI_SOURCE, MPI_ASK, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		int graphNodeIndex = findIndex(graph->nodesMapping, graph->nodesCount, info[0]);
		MPI_Request request = MPI_REQUEST_NULL;
		if (info[1] == 1 || info[1] == 2) {
			MPI_Isend(graph->outEdges[graphNodeIndex], graph->outDegrees[graphNodeIndex],
				MPI_INT, status.MPI_SOURCE, MPI_ASKED_OUT, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);
		}
		if (info[1] == 1 || info[1] == 3) {
			MPI_Isend(graph->inEdges[graphNodeIndex], graph->inDegrees[graphNodeIndex],
				MPI_INT, status.MPI_SOURCE, MPI_ASKED_IN, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);
		}
	}
}

void askForEdges(Graph * graph, int * info, int proc, int ** graphOutEdges, int * graphOutDegree,
				 int ** graphInEdges, int * graphInDegree, int rank) {
    MPI_Request request;
	MPI_Isend(info, 2, MPI_INT, proc, MPI_ASK, MPI_COMM_WORLD, &request);
	int recv = 0;
	do {
		MPI_Test(&request, &recv, MPI_STATUSES_IGNORE);
		trySendEdges(graph);
	} while(!recv);
	int out = 0, in = 0;
	if (info[1] == 1 || info[1] == 2) {
		MPI_Status outStatus;
		do {
			MPI_Iprobe(proc, MPI_ASKED_OUT, MPI_COMM_WORLD, &out, &outStatus);
			trySendEdges(graph);
		} while(!out);
		MPI_Get_count(&outStatus, MPI_INT, graphOutDegree);
		*graphOutEdges = safeMalloc(*graphOutDegree * sizeof(int));
		MPI_Recv(*graphOutEdges, *graphOutDegree, MPI_INT, proc, MPI_ASKED_OUT, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
	if (info[1] == 1 || info[1] == 3) {
		MPI_Status inStatus;
		do {
			MPI_Iprobe(proc, MPI_ASKED_IN, MPI_COMM_WORLD, &in, &inStatus);
			trySendEdges(graph);
		} while(!in);
		MPI_Get_count(&inStatus, MPI_INT, graphInDegree);
		*graphInEdges = safeMalloc(*graphInDegree * sizeof(int));
		MPI_Recv(*graphInEdges, *graphInDegree, MPI_INT, proc, MPI_ASKED_IN, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
}

void exploreMatch(Graph* graph, Graph* pattern, int * nodeProcMap, int rank, Match * match) {
    if (match->matchedNodesCount == pattern->nodesCount) {
    	MPI_Send(match->matchedNodes, match->matchedNodesCount + 1, MPI_INT, ROOT, MPI_MATCH_FOUND, MPI_COMM_WORLD);
    	match->matchedNodes[match->matchedNodesCount] = -1;
        match->matchedNodesCount--;
        return;
    }

    int freeParent = 0;
    int nextPatternNode = pattern->ordering[match->matchedNodesCount + 1];
    int nextPatternNodeParent = pattern->parents[abs(nextPatternNode)];
    int nextGraphNodeParent = match->matchedNodes[nextPatternNodeParent];
    int nextGraphNodeParentIndex = findIndex(graph->nodesMapping, graph->nodesCount, nextGraphNodeParent);
    int * parentEdges;
    int parentEdgesCount;
    if (nextGraphNodeParentIndex == -1) {
    	int info[2];
    	info[0] = nextGraphNodeParent;
    	if (nextPatternNode < 0) {
    		info[1] = 3;
    		askForEdges(graph, info, nodeProcMap[nextGraphNodeParent], NULL, NULL, &parentEdges, &parentEdgesCount, rank);
    	} else {
    		info[1] = 2;
    		askForEdges(graph, info, nodeProcMap[nextGraphNodeParent], &parentEdges, &parentEdgesCount, NULL, NULL, rank);
    	}
    	freeParent = 1;
    } else {
    	parentEdges = nextPatternNode < 0 ? graph->inEdges[nextGraphNodeParentIndex] : graph->outEdges[nextGraphNodeParentIndex];
    	parentEdgesCount = nextPatternNode < 0 ? graph->inDegrees[nextGraphNodeParentIndex] : graph->outDegrees[nextGraphNodeParentIndex];
    }

    int freeEdges;
    int * graphOutEdges, * graphInEdges;
    int graphOutDegree, graphInDegree;
    for (int i = 0; i < parentEdgesCount; i++) {
        int nextGraphNode = parentEdges[i];
        int nextGraphNodeProc = nodeProcMap[nextGraphNode];
        freeEdges = 0;
        if (nextGraphNodeProc == rank) {
        	int graphNodeIndex = findIndex(graph->nodesMapping, graph->nodesCount, nextGraphNode);
        	graphOutEdges = graph->outEdges[graphNodeIndex];
        	graphOutDegree = graph->outDegrees[graphNodeIndex];
        	graphInEdges = graph->inEdges[graphNodeIndex];
        	graphInDegree = graph->inDegrees[graphNodeIndex];
        } else {
        	freeEdges = 1;
        	int info[2];
        	info[0] = nextGraphNode;
        	info[1] = 1; // all edges
        	askForEdges(graph, info, nextGraphNodeProc, &graphOutEdges, &graphOutDegree, &graphInEdges, &graphInDegree, rank);
        }
        if (checkNodeMatches(nextGraphNode, abs(nextPatternNode), graphOutEdges,
        					 graphOutDegree, graphInEdges, graphInDegree, pattern, match, rank)) {
        	match->matchedNodes[++match->matchedNodesCount] = nextGraphNode;
            exploreMatch(graph, pattern, nodeProcMap, rank, match);
        }
        if (freeEdges) {
        	free(graphOutEdges);
        	free(graphInEdges);
        }
    }
    if (freeParent) {
    	free(parentEdges);
    }
    match->matchedNodes[match->matchedNodesCount] = -1;
    match->matchedNodesCount--;
}

void findMatches(Graph *graph, Graph* pattern, int * nodeProcMap, int rank) {
	for (int i = 1; i <= graph->nodesCount; i++) {
		Match match;
		memset(match.matchedNodes, -1, (MAX_MATCH_SIZE + 1) * sizeof(int));
		match.matchedNodesCount = 1;
		match.matchedNodes[1] = graph->nodesMapping[i];
		exploreMatch(graph, pattern, nodeProcMap, rank, &match);
	}
	int end = 0;
	MPI_Request request;
	MPI_Isend(NULL, 0, MPI_INT, ROOT, MPI_END_COMPUTATION, MPI_COMM_WORLD, &request);
	MPI_Request_free(&request);
	while (!end) {
		trySendEdges(graph);
		MPI_Iprobe(ROOT, MPI_END, MPI_COMM_WORLD, &end, MPI_STATUSES_IGNORE);
	}
	MPI_Recv(NULL, 0, MPI_INT, ROOT, MPI_END, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}


// ---------------------------- distribiute -----------------------------------------

NodeInfo * preprocessGraph(FILE * file, int * nodesCount) {
	char line[64];
	int node, neighCount, neighbour;
	fscanf(file, "%d\n", nodesCount);
	NodeInfo * nodeInfo = safeMalloc(((*nodesCount) + 1) * sizeof(NodeInfo));
	for (int i = 1 ; i <= *nodesCount; i++) {
		nodeInfo[i].inDegree = 0;
		nodeInfo[i].outDegree = 0;
	}
	while ((fgets(line, sizeof line, file) != NULL) && !isLineEmpty(line)) {
		sscanf(line, "%d %d", &node, &neighCount);
		nodeInfo[node].outDegree = neighCount;
		for (int i = 0; i < neighCount; i++) {
			fgets(line, sizeof line, file);
			sscanf(line, "%d", &neighbour);
			nodeInfo[neighbour].inDegree++;
		}
	}
	return nodeInfo;
}

int compare (const void * a, const void * b) {
	NodeComp * node1 = (NodeComp *) a;
	NodeComp * node2 = (NodeComp *) b;
	return node1->val - node2->val;
}

int intCompare(const void *a, const void *b) {
    int *val1 = (int *)a;
    int *val2 = (int *)b;
    return (*val1) - (*val2);
}

ProcInfo * assignNodesToProcesses(int nodesCount, int procNum, NodeInfo * nodeInfo, int ** nodeProcMap) {
	ProcInfo * procInfo = (ProcInfo *) safeMalloc(procNum * sizeof(ProcInfo));
	for (int i = 1; i < procNum; i++) {
		procInfo[i].memory = 500000000;
		procInfo[i].nodesWithOutEdges = 0;
		procInfo[i].nodesCount = 0;
	}
	*nodeProcMap = malloc((nodesCount + 1) * sizeof(int));
	memset(*nodeProcMap, -1, (nodesCount + 1) * sizeof(int));
	NodeComp * nodeComp = (NodeComp *) safeMalloc(nodesCount * sizeof(NodeComp));
	for (int i = 1; i <= nodesCount; i++) {
		nodeComp[i-1].val = nodeInfo[i].inDegree + nodeInfo[i].outDegree;
		nodeComp[i-1].node = i;
	}
	qsort(nodeComp, nodesCount, sizeof(NodeComp), compare);
	int proc = 1;
	int direction = 1;
	int i = 1;
	while (i <= nodesCount) {
		int node = nodeComp[i-1].node;
		if (nodeComp[i-1].val * 8 + 8 < procInfo[proc].memory) {
			(*nodeProcMap)[node] = proc;
			procInfo[proc].memory -= (nodeComp[i-1].val * 8 + 8);
			if (nodeInfo[node].outDegree > 0) {
				procInfo[proc].nodesWithOutEdges++;
			}
			procInfo[proc].nodesCount++;
			i++;
		}
		proc += direction;
		if (proc == procNum) {
			direction = -1;
			proc--;
		} else if (proc == 0) {
			direction = 1;
			proc++;
		}
	}
	free(nodeComp);
	return procInfo;
}

void distribiutePattern(Graph * pattern) {
	int nodesCount = pattern->nodesCount;
	MPI_Bcast(&nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->outDegrees, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->inDegrees, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->parents, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->ordering, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	for (int i = 1; i <= nodesCount; i++) {
		int outDegree = pattern->outDegrees[i];
		int inDegree = pattern->inDegrees[i];
		if (outDegree > 0) {
			MPI_Bcast(pattern->outEdges[i], outDegree, MPI_INT, ROOT, MPI_COMM_WORLD);
		}
		if (inDegree > 0) {
			MPI_Bcast(pattern->inEdges[i], inDegree, MPI_INT, ROOT, MPI_COMM_WORLD);
		}
	}
}

Graph * gatherPattern() {
	Graph * pattern = safeMalloc(sizeof(Graph));
	MPI_Bcast(&pattern->nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	pattern->inEdges = safeMalloc((pattern->nodesCount + 1) * sizeof(int*));
	pattern->outEdges = safeMalloc((pattern->nodesCount + 1) * sizeof(int*));
	pattern->inDegrees = safeMalloc((pattern->nodesCount + 1) * sizeof(int));
	pattern->outDegrees = safeMalloc((pattern->nodesCount + 1) * sizeof(int));
	pattern->parents = safeMalloc((pattern->nodesCount + 1) * sizeof(int));
	pattern->ordering = safeMalloc((pattern->nodesCount + 1) * sizeof(int));
	pattern->nodesMapping = NULL;
	MPI_Bcast(pattern->outDegrees, pattern->nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->inDegrees, pattern->nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->parents, pattern->nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	MPI_Bcast(pattern->ordering, pattern->nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	for (int i = 1; i <= pattern->nodesCount; i++) {
		int outDegree = pattern->outDegrees[i];
		int inDegree = pattern->inDegrees[i];
		if (outDegree > 0) {
			pattern->outEdges[i] = safeMalloc(outDegree * sizeof(int));
			MPI_Bcast(pattern->outEdges[i], outDegree, MPI_INT, ROOT, MPI_COMM_WORLD);
		} else {
			pattern->outEdges[i] = NULL;
		}
		if (inDegree > 0) {
			pattern->inEdges[i] = safeMalloc(inDegree * sizeof(int));
			MPI_Bcast(pattern->inEdges[i], inDegree, MPI_INT, ROOT, MPI_COMM_WORLD);
		} else {
			pattern->inEdges[i] = NULL;
		}
	}
	return pattern;
}

void distribiuteGraph(FILE * file, int *nodeProcMap, NodeInfo * nodeInfo, ProcInfo * procInfo, int procNum) {
	int nodesCount, receivier;
	int * neighbours, * neighboursProc;
	char line[64];
	// Sending info about nodes to all processes
	int procInOutInfo[2];
	for (int i = 1; i < procNum; i++) {
		procInOutInfo[0] = procInfo[i].nodesWithOutEdges;
		procInOutInfo[1] = procInfo[i].nodesCount;
		MPI_Send(procInOutInfo, 2, MPI_INT, i, MPI_PROC_IN_OUT_INFO, MPI_COMM_WORLD);
	}
	int nodeInOutInfo[3];
	fscanf(file, "%d\n", &nodesCount);
	while ((fgets(line, sizeof line, file) != NULL) && !isLineEmpty(line)) {
			sscanf(line, "%d %d", &nodeInOutInfo[0], &nodeInOutInfo[1]);
			nodeInOutInfo[2] = nodeInfo[nodeInOutInfo[0]].inDegree;
			receivier = nodeProcMap[nodeInOutInfo[0]];
			// Sending: node, outDegree, inDegree
			MPI_Send(nodeInOutInfo, 3, MPI_INT, receivier, MPI_NODE_IN_OUT_INFO, MPI_COMM_WORLD);
			neighbours = safeMalloc(nodeInOutInfo[1] * sizeof(int));
			neighboursProc = safeMalloc(nodeInOutInfo[1] * sizeof(int));
			for (int i = 0; i < nodeInOutInfo[1]; i++) {
				fgets(line, sizeof line, file);
				sscanf(line, "%d", &neighbours[i]);
				neighboursProc[i] = nodeProcMap[neighbours[i]];
			}
			MPI_Send(neighbours, nodeInOutInfo[1], MPI_INT, receivier, MPI_NODE_OUT_NEIGH, MPI_COMM_WORLD);
	}
	// Sending info of nodes without output edges
	int nodeInInfo[2];
	for (int i = 1; i <= nodesCount; i++) {
		if (nodeInfo[i].outDegree == 0) {
			nodeInInfo[0] = i;
			nodeInInfo[1] = nodeInfo[i].inDegree;
			// Sending: node, inDegree
			MPI_Send(nodeInInfo, 2, MPI_INT, nodeProcMap[i], MPI_NODE_IN_INFO, MPI_COMM_WORLD);
		}
	}
}

Graph * gatherGraph(int * inNodes, int * outNodes, int * nodeProcMap, int rank) {
	Graph * g = safeMalloc(sizeof(Graph));
	*inNodes = 0;
	*outNodes = 0;
	int procInOutInfo[2];
	// Receiving: nodesWithOutEdges, nodesCount
	MPI_Recv(procInOutInfo, 2, MPI_INT, ROOT, MPI_PROC_IN_OUT_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	int nodesWithOutEdges = procInOutInfo[0];
	int nodesCount = procInOutInfo[1];

	g->nodesCount = nodesCount;
	g->nodesMapping = safeMalloc((nodesCount + 1) * sizeof(int));
	int k = 1, l = 1;
	while (k <= nodesCount) {
		if (nodeProcMap[l] == rank) {
			g->nodesMapping[k] = l;
			k++;
		}
		l++;
	}
	g->nodesMapping[0] = -1;
	qsort(g->nodesMapping, g->nodesCount + 1, sizeof(int), intCompare);

	g->outDegrees = safeMalloc((nodesCount + 1) * sizeof(int));
	g->outEdges = safeMalloc((nodesCount + 1) * sizeof(int*));

	g->inDegrees = safeMalloc((nodesCount + 1) * sizeof(int));
	g->inEdges = safeMalloc((nodesCount + 1) * sizeof(int*));

	// Receiving: node, outDegree, inDegree
	int nodeInOutInfo[3];
	for (int i = 1; i <= nodesWithOutEdges; i++) {
		MPI_Recv(&nodeInOutInfo, 3, MPI_INT, ROOT, MPI_NODE_IN_OUT_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int node = nodeInOutInfo[0];
		int outDegree = nodeInOutInfo[1];
		int inDegree = nodeInOutInfo[2];

		*outNodes += outDegree;
		*inNodes += inDegree;
		int nodeIndex = findIndex(g->nodesMapping, g->nodesCount, node);
		g->outDegrees[nodeIndex] = outDegree;
		g->outEdges[nodeIndex] = safeMalloc(outDegree * sizeof(int));

		g->inDegrees[nodeIndex] = inDegree;
		if (inDegree > 0) {
			g->inEdges[nodeIndex] = safeMalloc(inDegree * sizeof(int));
		} else {
			g->inEdges[nodeIndex] = NULL;
		}
		MPI_Recv(g->outEdges[nodeIndex], nodeInOutInfo[1], MPI_INT, ROOT, MPI_NODE_OUT_NEIGH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// Receiving info about nodes without output edges
	int nodeInInfo[2];
	for (int i = 0; i < (nodesCount - nodesWithOutEdges); i++) {
		// Receive: node, inDegree
		MPI_Recv(&nodeInInfo, 2, MPI_INT, ROOT, MPI_NODE_IN_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int node = nodeInInfo[0];
		int inDegree = nodeInInfo[1];
		int nodeIndex = findIndex(g->nodesMapping, g->nodesCount, node);
		g->outDegrees[nodeIndex] = 0;
		g->outEdges[nodeIndex] = NULL;
		if (inDegree > 0) {
			g->inDegrees[nodeIndex] = inDegree;
			g->inEdges[nodeIndex] = safeMalloc(inDegree * sizeof(int));
		} else {
			g->inDegrees[nodeIndex] = 0;
			g->inEdges[nodeIndex] = NULL;
		}
		*inNodes += inDegree;
	}
	return g;
}

void transpose(Graph * g, int * nodeProcMap, int inNodes, int outNodes, int rank) {
	int * nodeInInfo = safeMalloc(outNodes * 2 * sizeof(int));
	int msgCounter = 0;
	MPI_Request request;
	for (int i = 1; i <= g->nodesCount; i++) {
		int node = g->nodesMapping[i];
		for (int j = 0; j < g->outDegrees[i]; j++) {
			nodeInInfo[msgCounter] = g->outEdges[i][j];
			nodeInInfo[msgCounter + 1] = node;
			// Sending outNode, node, rank
			MPI_Isend(nodeInInfo + msgCounter, 2, MPI_INT, nodeProcMap[g->outEdges[i][j]], MPI_NODE_IN_NEIGH, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);
			msgCounter += 2;
		}
	}
	int * counters = safeMalloc((g->nodesCount + 1) * sizeof(int));
	memset(counters, 0, (g->nodesCount + 1) * sizeof(int));
	// Receiving input edges
	int result[2];
	for (int i = 0; i < inNodes; i++) {
		MPI_Recv(&result, 2, MPI_INT, MPI_ANY_SOURCE, MPI_NODE_IN_NEIGH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int source = findIndex(g->nodesMapping, g->nodesCount, result[0]);
		g->inEdges[source][counters[source]++] = result[1];
	}
	free(nodeInInfo);
	free(counters);
}

void distribiute(int argc, char *argv[], int rank, int procNum, Graph ** graph, Graph ** pattern, int ** nodeProcMap, int * patternSize, FILE **outFile) {
	int nodesCount;
	if (rank == ROOT) {
		double startDistTime = MPI_Wtime();
		FILE *inFile;
		checkArgs(argc, argv, &inFile, outFile);
		NodeInfo * nodeInfo = preprocessGraph(inFile, &nodesCount);
		ProcInfo * procInfo = assignNodesToProcesses(nodesCount, procNum, nodeInfo, nodeProcMap);
		// Brodcast node --> process mapping
		MPI_Bcast(&nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(*nodeProcMap, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		// Distribute pattern
		*pattern = readGraph(inFile);
		*patternSize = (*pattern)->nodesCount;
		distribiutePattern(*pattern);
		freeGraph(*pattern);
		// Distribiute graph
		fseek(inFile, 0, SEEK_SET);
		distribiuteGraph(inFile, *nodeProcMap, nodeInfo, procInfo, procNum);
		// Free arrays
		free(procInfo);
		free(nodeInfo);
		free(*nodeProcMap);
		fclose(inFile);

		for (int i = 1; i <= procNum - 1; i++) {
			MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_END_DISTRIBUTION, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
		printf("Distribution time[s]: %d\n", (int)(MPI_Wtime() - startDistTime));
	} else {
		// Receive node --> process mapping
		MPI_Bcast(&nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		*nodeProcMap = safeMalloc((nodesCount + 1) * sizeof(int));
		MPI_Bcast(*nodeProcMap, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		// Receive pattern
		*pattern = gatherPattern();
		// Receive graph
		int inNodes, outNodes;
		*graph = gatherGraph(&inNodes, &outNodes, *nodeProcMap, rank);
		transpose(*graph, *nodeProcMap, inNodes, outNodes, rank);
		// Receive end of distribution
		MPI_Request request;
		MPI_Isend(NULL, 0, MPI_INT, ROOT, MPI_END_DISTRIBUTION, MPI_COMM_WORLD, &request);
		MPI_Request_free(&request);
	}
}

void matchNodes(int rank, int procNum, Graph * graph, Graph * pattern, int * nodeProcMap, int patternSize, FILE * outFile) {
	if (rank == ROOT) {
		double startComputationTime = MPI_Wtime();
		int flag;
		int * match = malloc((patternSize + 1) * sizeof(int));
		int endProcesses = 0;
		while(endProcesses != procNum - 1) {
			flag = 0;
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_MATCH_FOUND, MPI_COMM_WORLD, &flag, MPI_STATUSES_IGNORE);
			if (flag) {
				MPI_Recv(match, patternSize + 1, MPI_INT, MPI_ANY_SOURCE, MPI_MATCH_FOUND, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				writeMatch(match, patternSize, outFile);
				continue;
			}
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_END_COMPUTATION, MPI_COMM_WORLD, &flag, MPI_STATUSES_IGNORE);
			if (flag) {
				MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_END_COMPUTATION, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				endProcesses += 1;
			}
		}
		fclose(outFile);
		for (int i = 1; i < procNum; i++) {
			MPI_Send(NULL, 0, MPI_INT, i, MPI_END, MPI_COMM_WORLD);
		}
		printf("Computations time[s]: %d\n", (int)(MPI_Wtime() - startComputationTime));
	} else {
		findMatches(graph, pattern, nodeProcMap, rank);
		freeGraph(graph);
		free(pattern);
		free(nodeProcMap);
	}
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, procNum;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Graph * pattern, * graph;
    int * nodeProcMap;
    int patternSize;
    FILE * outFile;
    distribiute(argc, argv, rank, procNum, &graph, &pattern, &nodeProcMap, &patternSize, &outFile);
    matchNodes(rank, procNum, graph, pattern, nodeProcMap, patternSize, outFile);
    MPI_Finalize();
    return 0;
}
