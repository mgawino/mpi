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

#define MPI_NODE_IN_INFO 74
#define MPI_NODE_IN_NEIGH 75

#define MPI_END_DISTRIBUTION 100
#define MPI_MATCH 101
#define MPI_LAST_SEND_COMPLETE 102
#define MPI_END_COMPUTATION 103
#define MPI_MATCH_FOUND 104

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
	int sumDegrees;
} NodeComp;

typedef struct {
	int nodesWithOutEdges;
	int nodesCount;
	int memory;
} ProcInfo;

// ------------------ utils ----------------------------

FILE * checkArgs(int argc, char* argv[]) {
    if (argc != 3) {
    	fprintf(stderr, "ERROR: files not specified\n");
        MPI_Finalize();
        exit(1);
    }
    FILE* file;
    file = fopen(argv[1], "r");
    if (file == NULL) {
    	fprintf(stderr, "ERROR: file not found %s\n", argv[1]);
    	fclose(file);
    	MPI_Finalize();
    	exit(1);
    }
    return file;
}

// ---------------------- matching ------------------------------------

void printArray(int * tab, int len) {
	for (int i = 1; i <= len; i++) {
		printf("%d ", tab[i]);
	}
	printf("\n");
}

bool matchContains(int node, Match* match) {
    for (int i = 1; i <= match->matchedNodesCount; i++) {
       if (match->matchedNodes[i] == node) {
           return true;
       }
    }
    return false;
}

bool contains(int * nodes, int nodesCount, int target) {
    for (int i = 0; i < nodesCount; i++) {
        if (nodes[i] == target) {
            return true;
        }
    }
    return false;
}

bool checkNodeMatches(int graphNode, int patternNode, Graph* graph, Graph* pattern, Match* match) {
    if (matchContains(graphNode, match)) {
        return false;
    }

    int graphNodeIndex = findIndex(graph->nodesMapping, graph->nodesCount, graphNode);
    for (int i = 0; i < pattern->outDegrees[patternNode]; i++) {
        int patternNeigh = pattern->outEdges[patternNode][i];
        int graphNeigh = match->matchedNodes[patternNeigh];
        if (patternNeigh <= match->matchedNodesCount &&
        	!contains(graph->outEdges[graphNodeIndex], graph->outDegrees[graphNodeIndex], graphNeigh)) {
            return false;
        }
    }

    for (int i = 0; i < pattern->inDegrees[patternNode]; i++) {
		int patternNeigh = pattern->inEdges[patternNode][i];
		int graphNeigh = match->matchedNodes[patternNeigh];
		if (patternNeigh <= match->matchedNodesCount &&
			!contains(graph->inEdges[graphNodeIndex], graph->inDegrees[graphNodeIndex], graphNeigh)) {
			return false;
		}
    }
    return true;
}

void exploreReceivedMatch(Graph * graph, Graph * pattern, int * nodeProcMap, int rank, Match * match, MPI_Request * lastSendRequest);

void tryReceive(Graph * graph, Graph * pattern, int * nodeProcMap, int rank, MPI_Request * lastSendRequest) {
	int count, flag;
	MPI_Status recvStatus;
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_MATCH, MPI_COMM_WORLD, &flag, &recvStatus);
	if (flag) {
		MPI_Get_count(&recvStatus, MPI_INT, &count);
		Match match;
		match.matchedNodesCount = count - 1;
		MPI_Recv(&match.matchedNodes, count, MPI_INT, MPI_ANY_SOURCE, MPI_MATCH, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		exploreReceivedMatch(graph, pattern, nodeProcMap, rank, &match, lastSendRequest);
	}
}

void exploreMatch(Graph* graph, Graph* pattern, int * nodeProcMap, int rank, Match * match, MPI_Request * lastSendRequest) {
    if (match->matchedNodesCount == pattern->nodesCount) {
    	printf("MATCH\n");
    	MPI_Send(match->matchedNodes, match->matchedNodesCount + 1, MPI_INT, ROOT, MPI_MATCH_FOUND, MPI_COMM_WORLD);
        match->matchedNodesCount--;
        return;
    }

    int nextPatternNode = pattern->ordering[match->matchedNodesCount + 1];
    int nextPatternNodeParent = pattern->parents[abs(nextPatternNode)];
    int nextGraphNodeParent = match->matchedNodes[nextPatternNodeParent];
    int nextGraphNodeParentIndex = findIndex(graph->nodesMapping, graph->nodesCount, nextGraphNodeParent);
    int * parentEdges = nextPatternNode < 0 ? graph->inEdges[nextGraphNodeParentIndex] : graph->outEdges[nextGraphNodeParentIndex];
    int parentEdgesCount = nextPatternNode < 0 ? graph->inDegrees[nextGraphNodeParentIndex] : graph->outDegrees[nextGraphNodeParentIndex];

    int canSend = 1;
    MPI_Request sendRequest;
    int sendBuffer[MAX_MATCH_SIZE + 1];
    for (int i = 0; i < parentEdgesCount; i++) {
    	tryReceive(graph, pattern, nodeProcMap, rank, lastSendRequest);
        int nextGraphNode = parentEdges[i];
        int nextGraphNodeProc = nodeProcMap[nextGraphNode];
        if (nextGraphNodeProc == rank) {
        	if (checkNodeMatches(nextGraphNode, abs(nextPatternNode), graph, pattern, match)) {
        		match->matchedNodes[++match->matchedNodesCount] = nextGraphNode;
        		exploreMatch(graph, pattern, nodeProcMap, rank, match, lastSendRequest);
        	}
        } else {
        	while(!canSend) {
        		MPI_Test(&sendRequest, &canSend, MPI_STATUSES_IGNORE);
        		tryReceive(graph, pattern, nodeProcMap, rank, lastSendRequest);
        	}
        	match->matchedNodes[0] = nextGraphNode;
        	memcpy(sendBuffer, match->matchedNodes, (match->matchedNodesCount + 1) * sizeof(int));
        	MPI_Isend(sendBuffer, match->matchedNodesCount + 1, MPI_INT, nextGraphNodeProc, MPI_MATCH, MPI_COMM_WORLD, &sendRequest);
        	*lastSendRequest = sendRequest;
        	canSend = 0;
        }
    }
    match->matchedNodesCount--;
}

void exploreReceivedMatch(Graph * graph, Graph * pattern, int * nodeProcMap, int rank, Match * match, MPI_Request * lastSendRequest) {
	int nextGraphNode = match->matchedNodes[0];
	int nextPatternNode = pattern->ordering[match->matchedNodesCount + 1];
	if (checkNodeMatches(nextGraphNode, abs(nextPatternNode), graph, pattern, match)) {
		match->matchedNodes[++match->matchedNodesCount] = nextGraphNode;
		exploreMatch(graph, pattern, nodeProcMap, rank, match, lastSendRequest);
	}
}

void findMatches(Graph *graph, Graph* pattern, int * nodeProcMap, int rank) {
	MPI_Request lastSendRequest = MPI_REQUEST_NULL;
	for (int i = 1; i <= graph->nodesCount; i++) {
		Match match;
		memset(match.matchedNodes, -1, (MAX_MATCH_SIZE + 1) * sizeof(int));
		match.matchedNodesCount = 1;
		match.matchedNodes[1] = graph->nodesMapping[i];
		exploreMatch(graph, pattern, nodeProcMap, rank, &match, &lastSendRequest);
	}
	int compleated = lastSendRequest == MPI_REQUEST_NULL;
	MPI_Request request, dummyRequest;
	while(!compleated) {
		MPI_Test(&lastSendRequest, &compleated, MPI_STATUSES_IGNORE);
		tryReceive(graph, pattern, nodeProcMap, rank, &dummyRequest);
	}
	MPI_Isend(NULL, 0, MPI_INT, ROOT, MPI_LAST_SEND_COMPLETE, MPI_COMM_WORLD, &request);
	MPI_Request_free(&request);
	int end = 0;
	while(!end) {
		MPI_Iprobe(ROOT, MPI_END_COMPUTATION, MPI_COMM_WORLD, &end, MPI_STATUSES_IGNORE);
		tryReceive(graph, pattern, nodeProcMap, rank, &dummyRequest);
	}
	MPI_Recv(NULL, 0, MPI_INT, ROOT, MPI_END_COMPUTATION, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
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
	if(node1->sumDegrees < node2->sumDegrees) {
		return -1;
	} else if(node1->sumDegrees == node2->sumDegrees) {
		return 0;
	}
	return 1;
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
		nodeComp[i-1].sumDegrees = nodeInfo[i].inDegree + nodeInfo[i].outDegree;
		nodeComp[i-1].node = i;
	}
	qsort(nodeComp, nodesCount, sizeof(NodeComp), compare);
	int proc = 1;
	int direction = 1;
	int i = 1;
	while (i <= nodesCount) {
		int node = nodeComp[i-1].node;
		if (nodeComp[i-1].sumDegrees * 8 + 8 < procInfo[proc].memory) {
			(*nodeProcMap)[node] = proc;
			procInfo[proc].memory -= (nodeComp[i-1].sumDegrees * 8 + 8);
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

Graph * gatherGraph(int * inNodes, int * outNodes) {
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

	g->outDegrees = safeMalloc((nodesCount + 1) * sizeof(int));
	g->outEdges = safeMalloc((nodesCount + 1) * sizeof(int*));

	g->inDegrees = safeMalloc((nodesCount + 1) * sizeof(int));
	g->inEdges = safeMalloc((nodesCount + 1) * sizeof(int*));

	int nodeIndex;
	// Receiving: node, outDegree, inDegree
	int nodeInOutInfo[3];
	for (nodeIndex = 1; nodeIndex <= nodesWithOutEdges; nodeIndex++) {
		MPI_Recv(&nodeInOutInfo, 3, MPI_INT, ROOT, MPI_NODE_IN_OUT_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int node = nodeInOutInfo[0];
		int outDegree = nodeInOutInfo[1];
		int inDegree = nodeInOutInfo[2];

		*outNodes += outDegree;
		*inNodes += inDegree;

		g->nodesMapping[nodeIndex] = node;
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

		g->outDegrees[nodeIndex] = 0;
		g->outEdges[nodeIndex] = NULL;
		g->nodesMapping[nodeIndex] = node;
		if (inDegree > 0) {
			g->inDegrees[nodeIndex] = inDegree;
			g->inEdges[nodeIndex] = safeMalloc(inDegree * sizeof(int));
		} else {
			g->inDegrees[nodeIndex] = 0;
			g->inEdges[nodeIndex] = NULL;
		}
		*inNodes += inDegree;
		nodeIndex++;
	}
	return g;
}

void transpose(Graph * g, int * nodeProcMap, int inNodes, int outNodes, int rank) {
	int * counters = safeMalloc((g->nodesCount + 1) * sizeof(int));
	int * nodeInInfo = safeMalloc(outNodes * 3 * sizeof(int));
	int msgCounter = 0;
	MPI_Request request;
	memset(counters, 0, (g->nodesCount + 1) * sizeof(int));
	for (int i = 1; i <= g->nodesCount; i++) {
		int node = g->nodesMapping[i];
		for (int j = 0; j < g->outDegrees[i]; j++) {
			nodeInInfo[msgCounter] = g->outEdges[i][j];
			nodeInInfo[msgCounter + 1] = node;
			nodeInInfo[msgCounter + 2] = rank;
			// Sending outNode, node, rank
			MPI_Isend(nodeInInfo + msgCounter, 3, MPI_INT, nodeProcMap[g->outEdges[i][j]], MPI_NODE_IN_NEIGH, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);
			msgCounter += 3;
		}
	}

	// Receiving input edges
	int result[3];
	for (int i = 0; i < inNodes; i++) {
		MPI_Recv(&result, 3, MPI_INT, MPI_ANY_SOURCE, MPI_NODE_IN_NEIGH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int source = findIndex(g->nodesMapping, g->nodesCount, result[0]);
		g->inEdges[source][counters[source]++] = result[1];
	}
	free(nodeInInfo);
	free(counters);
}

void distribiute(int argc, char *argv[], int rank, int procNum, Graph ** graph, Graph ** pattern, int ** nodeProcMap, int * patternSize) {
	int nodesCount;
	if (rank == ROOT) {
		double startTime = MPI_Wtime();
		FILE * file = checkArgs(argc, argv);
		NodeInfo * nodeInfo = preprocessGraph(file, &nodesCount);
		ProcInfo * procInfo = assignNodesToProcesses(nodesCount, procNum, nodeInfo, nodeProcMap);
		// Brodcast node --> process mapping
		printf("Node assign time: %f\n", MPI_Wtime() - startTime);
		MPI_Bcast(&nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(*nodeProcMap, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		printf("Proc map brodcast time: %f\n", MPI_Wtime() - startTime);
		// Distribute pattern
		*pattern = readGraph(file);
		*patternSize = (*pattern)->nodesCount;
		distribiutePattern(*pattern);
		freeGraph(*pattern);
		// Distribiute graph
		fseek(file, 0, SEEK_SET);
		distribiuteGraph(file, *nodeProcMap, nodeInfo, procInfo, procNum);
		// Free arrays
		free(procInfo);
		free(nodeInfo);
		free(*nodeProcMap);
		fclose(file);

		// Check if all processes ends - TO REMOVE!!!
		printf("Graph sent time: %f\n", MPI_Wtime() - startTime);
		MPI_Gather(NULL, 0, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
		printf("Distribution time: %f\n", MPI_Wtime() - startTime);
	} else {
		// Receive node --> process mapping
		MPI_Bcast(&nodesCount, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		*nodeProcMap = safeMalloc((nodesCount + 1) * sizeof(int));
		MPI_Bcast(*nodeProcMap, nodesCount + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		// Receive pattern
		*pattern = gatherPattern();
		// Receive graph
		int inNodes, outNodes;
		*graph = gatherGraph(&inNodes, &outNodes);
		transpose(*graph, *nodeProcMap, inNodes, outNodes, rank);
		// Receive end of distribution
		MPI_Request request;
		MPI_Isend(NULL, 0, MPI_INT, ROOT, MPI_END_DISTRIBUTION, MPI_COMM_WORLD, &request);
		MPI_Request_free(&request);
	}
}

void matchNodes(int rank, int procNum, Graph * graph, Graph * pattern, int * nodeProcMap, int patternSize) {
	if (rank == ROOT) {
		int lastSendsCompleated = 0;
		int flag;
		int * match = malloc((patternSize + 1) * sizeof(int));
		while(lastSendsCompleated != procNum - 1) {
			flag = 0;
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_MATCH_FOUND, MPI_COMM_WORLD, &flag, MPI_STATUSES_IGNORE);
			if (flag) {
				MPI_Recv(match, patternSize + 1, MPI_INT, MPI_ANY_SOURCE, MPI_MATCH_FOUND, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				printArray(match, patternSize);
				continue;
			}
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_LAST_SEND_COMPLETE, MPI_COMM_WORLD, &flag, MPI_STATUSES_IGNORE);
			if (flag) {
				MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_LAST_SEND_COMPLETE, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				lastSendsCompleated += 1;
			}
		}
		MPI_Request request;
		for (int i = 1; i < procNum; i++) {
			MPI_Isend(NULL, 0, MPI_INT, i, MPI_END_COMPUTATION, MPI_COMM_WORLD, &request);
			MPI_Request_free(&request);
		}
	} else {
		findMatches(graph, pattern, nodeProcMap, rank);
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
    distribiute(argc, argv, rank, procNum, &graph, &pattern, &nodeProcMap, &patternSize);
    matchNodes(rank, procNum, graph, pattern, nodeProcMap, patternSize);
    MPI_Finalize();
    return 0;
}
