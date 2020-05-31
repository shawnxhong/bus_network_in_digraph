#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define MAX 46
#define INFINITE 1000
 
typedef struct Edge{
	struct Edge *next;//next edge with the same source
	int length; // the length of this edge
	int destID; // to whom it is pointing 
}Edge; 

// vertex is a bus Stop
typedef struct Vertex
{
	char *name; // the name of this stop
	struct Edge *firstEdge; // the first edge this vertex point out 
}Vertex;

//
typedef struct Graph
{
	struct Vertex *array; // array storing all vertices in the graph
	int vexNum; // number of vertices
	int edgeNum; // number of edges
}Graph;

typedef struct MinHeapNode
{
	int v; // vertex ID
	int distance; // the distance between this node and the source node
}MinHeapNode;

typedef struct MinHeap
{
	struct MinHeapNode **array; // pointer to the array storing all heap nodes
	int size; // size of the heap
	int capacity; // capacity of the heap
	int *pos; // array that store the position of vertex v in this heap
	
}MinHeap;

// create a new vertex
Vertex *newVertex(char *stopName)
{
	Vertex *new;
	new = malloc(sizeof(Vertex));
	assert(new != NULL);
	new->name = strdup(stopName);
	new->firstEdge = NULL;
	return new;
}

// create a new edge
Edge *newEdge(int dest, int distance)
{
	Edge *new;
	new = malloc(sizeof(Edge));
	assert(new != NULL);
	new->destID = dest;
	new->length = distance;
	new->next = NULL;
	return new;
}

// create a new graph
Graph *newGraph(int VNums){
	Graph *G;
	G = malloc(sizeof(Graph));
	assert(G != NULL);
	G->vexNum = VNums;
	G->edgeNum = 0;	
	G->array = (Vertex *)malloc(VNums * sizeof(Vertex));
	int i;
	for (i = 0;i < VNums;++i)
	{
		G->array[i].firstEdge= NULL;
	}
	return G;
}

MinHeapNode *newMinHeapNode(int v, int dist)
{
	MinHeapNode *new;
	new = malloc(sizeof(MinHeapNode));
	assert(new != NULL);
	new->v = v;
	new->distance = dist;
	return new;
}

MinHeap *newMinHeap(int c)
{
	MinHeap *H;
	H = malloc(sizeof(MinHeap));
	assert(H != NULL);
	H->size = 0;
	H->capacity = c;
	H->array = (MinHeapNode **)malloc(c * sizeof(MinHeapNode*));
	H->pos = (int *)malloc(c * sizeof(int));
	return H;
}

void swapMinHeapNode(MinHeapNode **a, MinHeapNode **b)
{
	MinHeapNode *temp = *a;
	*a = *b;
	*b = temp;
}

// A standard function to heapify at given idx 
// This function also updates position of nodes when they are swapped. 
// Position is needed for decreaseKey() 
// O(logN)
void minHeapify(MinHeap* minHeap, int idx) 
{ 
    int smallest, left, right; 
    smallest = idx; 
    left = 2 * idx + 1;
    right = 2 * idx + 2; 

    if (left < minHeap->size && minHeap->array[left]->distance < minHeap->array[smallest]->distance )
      smallest = left;
  
    if (right < minHeap->size && minHeap->array[right]->distance < minHeap->array[smallest]->distance ) 
      smallest = right; 
  
    if (smallest != idx) 
    { 
        // The nodes to be swapped in min heap 
        MinHeapNode *smallestNode = minHeap->array[smallest]; 
        MinHeapNode *idxNode = minHeap->array[idx]; 
  
        // Swap positions 
        minHeap->pos[smallestNode->v] = idx; 
        minHeap->pos[idxNode->v] = smallest; 
  
        // Swap nodes 
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]); 
  
        minHeapify(minHeap, smallest); //heapify downwards recursively 
    }
}
  
// Standard function to extract minimum node from heap 
// O(logN)
MinHeapNode* extractMin(MinHeap* minHeap) 
{ 
    if (minHeap->size == 0) 
        return NULL; 
  
    // Store the root node 
    MinHeapNode* root = minHeap->array[0]; 
  
    // Replace root node with last node 
    MinHeapNode* lastNode = minHeap->array[minHeap->size - 1]; 
    minHeap->array[0] = lastNode; 
  
    // Update position of last node 
    minHeap->pos[root->v] = minHeap->size-1;  //pos�������vertex�ڶ��������λ�� 
    minHeap->pos[lastNode->v] = 0;
  
    // Reduce heap size and heapify root 
    --minHeap->size; 
    minHeapify(minHeap, 0);
  
    return root; 
} 
  
// Function to decreasy dist value of a given vertex v. This function 
// uses pos[] of min heap to get the current index of node in min heap 
// O(1)
void decreaseKey(MinHeap* minHeap, int v, int dist) 
{ 
    // Get the index of v in  heap array 
    int i = minHeap->pos[v]; 
  
    // Get the node and update its dist value 
    minHeap->array[i]->distance = dist; 
  
    // Travel up while the complete tree is not hepified. 
    // This is a O(Logn) loop 
    while (i && minHeap->array[i]->distance < minHeap->array[(i - 1) / 2]->distance) 
    { 
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i-1)/2; 
        minHeap->pos[minHeap->array[(i-1)/2]->v] = i; 
        swapMinHeapNode(&minHeap->array[i],  &minHeap->array[(i - 1) / 2]); 
  
        // move to parent index 
        i = (i - 1) / 2; 
    }
} 
  
// A utility function to check if a given vertex 
// 'v' is in min heap or not 
// O(1)
int isInMinHeap(struct MinHeap *minHeap, int v) 
{ 
   if (minHeap->pos[v] < minHeap->size) 
    	return 1; 
   return 0; 
}

//get the number of bus stop from the input file 
//O (n) 
int busStopsNum(const char *busStops)
{
	FILE *fp = fopen(busStops, "r");
	char buffer[48];

	int num = 0;
	while (fgets(buffer, 48, fp) !=NULL)
	{
		num++;
	}
	fclose(fp);
	return num;
}

// time complexity O (n) 
// given the name of the input file, the number of vertices and an empty array
// return an array of vertices
Vertex *verticeArray(const char *busStops, int num, Vertex *array[]){
	FILE *fp = fopen(busStops, "r");
	char buffer[48];
	
	char stopName[20];
	int key, ret;
	
	while (fgets(buffer, 48, fp) !=NULL)
	{
		ret = sscanf(buffer, "%d:%[^\n]", &key, stopName); //read characters until \n is found
		Vertex *new = newVertex(stopName);
		array[key] = new;
	}
	fclose(fp);
	
	return array[0];
}

//given a file, store all vertex information into a graph without the edges
//O(n)
Graph *createGraph(const char *busStops)
{
	// get the size of the graph first 
	int StopNum;
	StopNum = busStopsNum(busStops);
	
	//make an array of all the vertices
	Vertex *Vertices[StopNum];
	Vertices[0] = verticeArray(busStops, StopNum, Vertices);
	//create an empty graph
	Graph *G = newGraph(StopNum);	
	G->vexNum = StopNum;
	//put the vertice array into the graph
	int i = 0;
	while (i < StopNum)
	{
		G->array[i].name = Vertices[i]->name;
		i++;
	}
	/*i = 0;
	while (i<StopNum)
	{
		printf("%d: %s\n", i, G->array[i].name);
		i++;
	}*/ 
	return G;
}

//insert edges into the graph
//O(m), m is the number of edges
void *InsertEdge(Graph *G, const char *Distance)
{
	FILE *fp = fopen(Distance, "r");
	int head, tail, distance;
	while (fscanf(fp, "%d-%d:%d", &head, &tail, &distance)!=EOF)
	{	
		Edge *new = newEdge(tail, distance);
		new->next = G->array[head].firstEdge;
		G->array[head].firstEdge = new;
		G->edgeNum += 1;
		//printf("%d-%d:%d\n", head, tail, distance);
	}
	
	fclose(fp);
}

void printGraph(Graph *G)
{
	int i;
	for (i=0; i < G->vexNum; i++)
	{
		Edge *currE = G->array[i].firstEdge;
		printf("\n Adjacency list of vertex %d\n edge ", i);
		while (currE != NULL)
		{
			printf("-> %d", currE->destID);
			currE = currE->next;
		}
	}
}

//generated a file to store edges information in reversed orientation
char *reverseDistanceFile(const char *Distance)
{	
	char *filename = "DistanceReverse.txt";
	FILE *fp1 = fopen(Distance, "r");
	FILE *fp2 = fopen(filename, "a");
	int head, tail, distance;
	while (fscanf(fp1, "%d-%d:%d", &head, &tail, &distance)!=EOF)
	{	
		fprintf(fp2, "%d-%d:%d\n", tail, head, distance);
	}
	
	fclose(fp1);
	fclose(fp2);
	return filename;
}

// an array to store boolean value of whether a vertex has been visited
int visited[MAX];

//O(n)
void DFS(Graph *G, int v)
{
	//printf("visiting v %d: %s\n", v, G->array[v].name);
	visited[v] = 1;
	Edge *currE = G->array[v].firstEdge;

	while (currE!=NULL) // dfs for unvisited adjacent vertices of v
	{
		if( visited[currE->destID] == 0)
			DFS(G, currE->destID);
		currE = currE->next;
	}
}

//do a DFS traverse starting from vertex v; 
// O(n)
void DFSTraverse(Graph *G, int v)
{
	int n = G->vexNum;
	int i;
	for (i=0; i<n; i++)
	{
		visited[i] = 0;
	} //mark allthe vertices as unvisited
	
	/*for (v=0; v<n; v++)
	{
		if (visited[v] == 0)
			DFS(G,v);
	}*/
	DFS(G,v); 
}

//given those input files, decide whether the graph is strongly connected or not
// worest case time complexity O(m+n):
// create the graph O(n)
// insert all the edges O(m)
// DFS travesal twice O(n)
// altogether O(m+n)
int StrongConnectivity(const char *busStops, const char * BusNames, const char *BusRoutes, const char *Distance)
{	
	int allVisited = 1;
	
	Graph *G;
	G = createGraph(busStops);
	int n = G->vexNum; 
	InsertEdge(G, Distance);
	
	char *reverseDistance;
	reverseDistance = reverseDistanceFile(Distance);

	Graph *reverseG;
	reverseG = createGraph(busStops);
	InsertEdge(reverseG, reverseDistance);
	
	//printGraph(G);
	//printGraph(reverseG);
	
	//DFS from vertex 0
	DFSTraverse(G, 0);
	int i;
	for( i = 0; i < n; i++){
		//printf("check v %d, visited?: %d\n",i, visited[i]);
		if(visited[i] == 0){
			allVisited = 0;
			break;
		}
	}
	
	//DFS the reversed graph from vertex 0
	DFSTraverse(reverseG, 0);
	for( i = 0; i < n; i++){
		//printf("check v %d, visited?: %d\n",i, visited[i]);
		if(visited[i] == 0){
			allVisited = 0;
			break;
		}
	}
	
	if (allVisited == 0)
		return 0;
	else
		return 1;
}

// print those maximal strongly connected components separately
// creat the graph with the edges O(m+n)
// DFS travers all the SCC, altogether O(n)
// So O(m+n)
void maximalStronglyComponents(const char *busStops, const char *BusNames, const char *BusRoutes, const char *Distance)
{
	Graph *G;
	G = createGraph(busStops);
	int n = G->vexNum; 
	InsertEdge(G, Distance);
	
	int SCCNum = 1;
	int currentVertex = 0;
	
	while (currentVertex < n){
		DFSTraverse(G, currentVertex);
		printf("\nStrongly connected component %d:\n",SCCNum);
		
		while (visited[currentVertex] !=0 )
		{
			printf("bus stop %d: %s\n", currentVertex, G->array[currentVertex].name);
			currentVertex++;
			if (currentVertex > n -1){
				break;
			}
		}
		SCCNum++;
	}
}


// given the name of a bus stop and the graph, return the serial number of that stop
// worest case O(n)
int searchStop(const char *StopName, Graph *G)
{
	int v;
	int n = G->vexNum;
	
	for(v=0; v<n; v++)
	{
		if (strcmp(StopName, G->array[v].name) == 0)
			return v;
	}
}

// do a BFS traversal starting from vertex v;
// O(n)
void BFSTraverse(Graph *G, int v)
{
	int i, n;
	n = G->vexNum;
	for(i=0; i < n; i++)
	{
		visited[i] = 0;
	}
	
	int queue[n];
	int Front = 0;
	int Rear = 0;
	int currV = v;
	printf("bus stop %d: %s\n",currV, G->array[currV].name);
	visited[currV] = 1; 
	queue[Rear] = currV;
	Rear++;
	
	while (Front != Rear)
	{
		currV = queue[Front];
		Front++;
		
		Edge *currE = G->array[currV].firstEdge;
		while (currE !=NULL)
		{
			if (visited[currE->destID] == 0){
				printf("bus stop %d: %s\n",currE->destID, G->array[currE->destID].name);
				visited[currE->destID] = 1;
				queue[Rear] = currE->destID;
				Rear++;
			}
			currE = currE->next;
		}
	}
}
 
// print all the reachable stops from the source stop
// worst case O (n + m) 
void reachableStops(const char *sourceStop, const char *busStops, const char *BusNames, const char *BusRoutes, const char *Distance)
{
	Graph *G;
	G = createGraph(busStops);
	int n = G->vexNum; 
	InsertEdge(G, Distance);
	
	int v;
	v = searchStop(sourceStop, G);
	
	printf("\nBus stops that are reachable from %s:\n", sourceStop);
	BFSTraverse(G, v);
}



// Dijstra algrithom using a priority queue : O((m+n)logn)
void TravelRoute(const char *sourceStop, const char *destinationStop, const char *busStops, const char *BusNames, const char *BusRoutes, const char *Distance)
{	
	// read in the bus route and generate an 2D array storing which bus to take between two adjacent stops
	/*int Choice[46][46];
	FILE *fp1 = fopen(BusRoutes, "r");
	char buffer[256];
	char route[256];
	int bus, stop, ret;
	
	while (fgets(buffer, 256, fp1) !=NULL)
	{
		ret = sscanf(buffer, "%d:%s", &bus, route); //read characters until \n is not found
		//printf("%d: %s\n", bus, route);
		char *token = strtok(route, ",#");
		int Q[2] = {-1,-1};
		while (token != NULL){
			int x = atoi(token);
			token = strtok(NULL, ",#");
			
			if (Q[0] == -1 && Q[1] == -1)
			{
				Q[0] = x;
			}
			else if (Q[0] != -1 && Q[1] == -1)
			{
				Q[1] = x;
				Choice[Q[0]][Q[1]] = bus;
			}
			else if (Q[0] != -1 && Q[1] != -1)
			{
				Q[0] = Q[1];
				Q[1] = x;
				Choice[Q[0]][Q[1]] = bus;
			}
		}
	}
	
	fclose(fp1);*/
	
	// build the Graph
	Graph *G;
	G = createGraph(busStops);
	int n = G->vexNum; 
	InsertEdge(G, Distance);
	
	//match the stop ID and its name
	int sourceV = searchStop(sourceStop, G);
	int destV = searchStop(destinationStop, G);
		
	int dist[n]; // the shortest distance from this stop to the source stop
	int prev[n]; // the previous stop ID on the shortest path from this stop to the source stop
	MinHeap *Queue = newMinHeap(n);
	int i;
	for(i=0; i<n; i++)
	{
		dist[i] = INFINITE;
		Queue->array[i] = newMinHeapNode(i, dist[i]);
		Queue->pos[i] = i;
		prev[i] = -1;
	}
	
	Queue->pos[sourceV] = sourceV;
	dist[sourceV] = 0;
	decreaseKey(Queue, sourceV, dist[sourceV]); 
	
	Queue->size = n;
	while( Queue->size > 0)
	{
		MinHeapNode *minNode = extractMin(Queue);
		int u = minNode->v;
		Edge *currE = G->array[u].firstEdge;
		int v;
		while(currE !=NULL)
		{
			v = currE->destID;
			if (isInMinHeap(Queue, v) && dist[u] != INFINITE && currE->destID + dist[u] < dist[v])
			{
				dist[v] = dist[u] + currE->length;
				prev[v] = u;
				decreaseKey(Queue, v, dist[v]);
			}
			currE = currE->next;
		}
	}
	
	/*for(i=0; i<n; i++)
	{
		printf("vertex: %d, dist[%d]= %d, prev[%d] = %d\n", i,i,dist[i],i,prev[i]);
	}*/
	
	if (dist[destV] == INFINITE)
	{
		printf("\nNo route exists from %s to %s\n", sourceStop, destinationStop);
	}
	else
	{	
		int currStop = destV;
		int prevStop;
		
		int stack[MAX];
		for(i=0; i < MAX; i++)
		{
			stack[i] = -1;
		}
		
		stack[0] = destV;
		i = 1;
		
		while (prev[currStop] != -1)
		{
			stack[i] = prev[currStop];
			i++;
			currStop = prev[currStop];
		}
		int stackSize = i;
		
		int path[stackSize];
		int k=0;
		for ( i = stackSize-1 ; i >= 0; i--)
		{
			path[k] = stack[i];
			k++;
		}
		
		printf("\nshortest path from %s to %s: ", sourceStop, destinationStop);
		for(k=0; k<stackSize; k++)
		{
			printf("%d, ", path[k]);
		}
		
		/*for(q = 0; q < i; q++)
		{
			printf("from %d to %d, bus: %d\n",path[q], path[q+1], Choice[path[q]][path[q+1]]); //bus name
		}*/
	}
}


int main(){
	
	int c;
	c = StrongConnectivity("busStops.txt", "busNames.txt", "busRoute.txt", "Distance.txt");
	
	printf("connectivity: %d\n", c);
	
	maximalStronglyComponents("busStops.txt", "busNames.txt", "busRoute.txt", "Distance.txt");
	
	reachableStops("General Holmes", "busStops.txt", "busNames.txt", "busRoute.txt", "Distance.txt");
	
	TravelRoute("Bathurst Street","Broadway","busStops.txt", "busNames.txt", "busRoute.txt", "Distance.txt");
	
	
	return 0;
}
