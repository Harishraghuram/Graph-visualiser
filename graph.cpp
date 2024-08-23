#include <graphviz/gvc.h>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <stack>
#include <sstream>
#include <algorithm>

struct Edge {
    int src, dest;
    double cost;
};

// Function to add a description node to the graph
void addGraphDescription(Agraph_t* g, int numNodes, int numEdges, int sourceNode, int destinationNode, double minCost, const std::vector<int>& topoOrder, const std::vector<int>& articulationPoints, const std::vector<int>& path) {
    std::stringstream description;
    description << "Number of nodes: " << numNodes << "\\n";
    description << "Number of edges: " << numEdges << "\\n";
    description << "Minimum cost from Node " << (sourceNode + 1) << " to Node " << (destinationNode + 1) << ": " << minCost << "\\n";
    
    description << "Path: ";
    for (size_t i = 0; i < path.size(); ++i) {
        description << path[i] << " ";
    }
    description << "\\n";
    
    description << "Topological Order: ";
    for (size_t i = 0; i < topoOrder.size(); ++i) {
        description << (topoOrder[i] + 1) << " ";
    }
    
    description << "\\nArticulation Points: ";
    for (size_t i = 0; i < articulationPoints.size(); ++i) {
        description << (articulationPoints[i] + 1) << " ";
    }
    
    Agnode_t *descNode = agnode(g, const_cast<char*>("Description"), 1);
    agsafeset(descNode, const_cast<char*>("label"), const_cast<char*>(description.str().c_str()), const_cast<char*>(""));
    agsafeset(descNode, const_cast<char*>("shape"), const_cast<char*>("note"), const_cast<char*>(""));
    agsafeset(descNode, const_cast<char*>("color"), const_cast<char*>("black"), const_cast<char*>(""));
    agsafeset(descNode, const_cast<char*>("style"), const_cast<char*>("dashed"), const_cast<char*>(""));
}

// Helper function for articulation point detection
void articulationPointsUtil(int u, int parent, std::vector<int>& disc, std::vector<int>& low, std::vector<bool>& visited, std::vector<bool>& ap, const std::vector<std::vector<int> >& adj) {
    static int time = 0;
    int children = 0;
    visited[u] = true;
    disc[u] = low[u] = ++time;
    
    for (size_t i = 0; i < adj[u].size(); ++i) {
        int v = adj[u][i];
        if (!visited[v]) {
            ++children;
            articulationPointsUtil(v, u, disc, low, visited, ap, adj);
            low[u] = std::min(low[u], low[v]);
            if (parent != -1 && low[v] >= disc[u]) {
                ap[u] = true;
            }
        } else if (v != parent) {
            low[u] = std::min(low[u], disc[v]);
        }
    }
    
    if (parent == -1 && children > 1) {
        ap[u] = true;
    }
}

// Function to find articulation points
std::vector<int> findArticulationPoints(int numNodes, const std::vector<Edge>& edges) {
    std::vector<std::vector<int> > adj(numNodes);
    for (size_t i = 0; i < edges.size(); ++i) {
        const Edge& edge = edges[i];
        adj[edge.src].push_back(edge.dest);
        adj[edge.dest].push_back(edge.src); // For undirected graph
    }
    
    std::vector<int> disc(numNodes, -1);
    std::vector<int> low(numNodes, -1);
    std::vector<bool> visited(numNodes, false);
    std::vector<bool> ap(numNodes, false);
    
    for (int i = 0; i < numNodes; ++i) {
        if (!visited[i]) {
            articulationPointsUtil(i, -1, disc, low, visited, ap, adj);
        }
    }
    
    std::vector<int> articulationPoints;
    for (int i = 0; i < numNodes; ++i) {
        if (ap[i]) {
            articulationPoints.push_back(i);
        }
    }
    
    return articulationPoints;
}

// Depth-First Search helper function for topological sorting
void topologicalSortUtil(int v, std::vector<bool>& visited, std::stack<int>& Stack, const std::vector<std::vector<int> >& adj) {
    visited[v] = true;
    
    for (size_t i = 0; i < adj[v].size(); ++i) {
        int i_node = adj[v][i];
        if (!visited[i_node]) {
            topologicalSortUtil(i_node, visited, Stack, adj);
        }
    }
    
    Stack.push(v);
}

// Function to perform topological sorting
std::vector<int> topologicalSort(int numNodes, const std::vector<Edge>& edges) {
    std::vector<std::vector<int> > adj(numNodes);
    for (size_t i = 0; i < edges.size(); ++i) {
        const Edge& edge = edges[i];
        adj[edge.src].push_back(edge.dest);
    }
    
    std::stack<int> Stack;
    std::vector<bool> visited(numNodes, false);
    
    for (int i = 0; i < numNodes; ++i) {
        if (!visited[i]) {
            topologicalSortUtil(i, visited, Stack, adj);
        }
    }
    
    std::vector<int> topoOrder;
    while (!Stack.empty()) {
        topoOrder.push_back(Stack.top());
        Stack.pop();
    }
    
    return topoOrder;
}

int main() {
    // Initialize a Graphviz context
    GVC_t *gvc = gvContext();

    // Create a new directed graph
    Agraph_t *g = agopen(const_cast<char*>("g"), Agdirected, 0);

    int numNodes, numEdges;
    std::cout << "Enter the number of nodes: ";
    std::cin >> numNodes;
    std::cout << "Enter the number of edges: ";
    std::cin >> numEdges;

    std::vector<Agnode_t*> nodes;
    std::vector<Edge> edges;

    // Create nodes
    for (int i = 0; i < numNodes; ++i) {
        std::string nodeName = "Node " + std::to_string(i + 1);
        Agnode_t *node = agnode(g, const_cast<char*>(nodeName.c_str()), 1);
        nodes.push_back(node);
    }

    // Create edges
    for (int i = 0; i < numEdges; ++i) {
        int src, dest;
        double cost;
        std::cout << "Enter the source node (1 to " << numNodes << ") for edge " << (i + 1) << ": ";
        std::cin >> src;
        std::cout << "Enter the destination node (1 to " << numNodes << ") for edge " << (i + 1) << ": ";
        std::cin >> dest;
        std::cout << "Enter the cost for edge " << (i + 1) << ": ";
        std::cin >> cost;

        if (src < 1 || src > numNodes || dest < 1 || dest > numNodes) {
            std::cout << "Invalid node number. Please enter values between 1 and " << numNodes << "." << std::endl;
            --i; // Decrement i to retry this edge
            continue;
        }

        Edge edge;
        edge.src = src - 1;
        edge.dest = dest - 1;
        edge.cost = cost;
        edges.push_back(edge);

        Agedge_t *graphEdge = agedge(g, nodes[src - 1], nodes[dest - 1], 0, 1);

        // Set an optional attribute for the edge
        std::string edgeLabel = "Cost: " + std::to_string(cost);
        agsafeset(graphEdge, const_cast<char*>("label"), const_cast<char*>(edgeLabel.c_str()), const_cast<char*>(""));
    }

    // Bellman-Ford algorithm to find the minimum cost
    std::vector<double> distances(numNodes, std::numeric_limits<double>::infinity());
    std::vector<int> predecessors(numNodes, -1);

    int sourceNode, destinationNode;
    std::cout << "Enter the source node (1 to " << numNodes << ") for cost calculation: ";
    std::cin >> sourceNode;
    std::cout << "Enter the destination node (1 to " << numNodes << ") for cost calculation: ";
    std::cin >> destinationNode;

    sourceNode--;
    destinationNode--;

    distances[sourceNode] = 0;

    // Relax edges
    for (int i = 1; i <= numNodes - 1; ++i) {
        for (size_t j = 0; j < edges.size(); ++j) {
            const Edge& edge = edges[j];
            if (distances[edge.src] + edge.cost < distances[edge.dest]) {
                distances[edge.dest] = distances[edge.src] + edge.cost;
                predecessors[edge.dest] = edge.src;
            }
        }
    }

    // Check for negative-weight cycles
    for (size_t j = 0; j < edges.size(); ++j) {
        const Edge& edge = edges[j];
        if (distances[edge.src] + edge.cost < distances[edge.dest]) {
            std::cerr << "Graph contains a negative-weight cycle" << std::endl;
            return 1;
        }
    }

    // Print minimum cost
    double minCost = distances[destinationNode];
    std::vector<int> path;
    if (minCost == std::numeric_limits<double>::infinity()) {
        std::cout << "No path from Node " << (sourceNode + 1) << " to Node " << (destinationNode + 1) << std::endl;
    } else {
        std::cout << "Minimum cost from Node " << (sourceNode + 1) << " to Node " << (destinationNode + 1) << " is " << minCost << std::endl;
        std::cout << "Path: ";
        std::stack<int> pathStack;
        int v = destinationNode;
        while (v != -1) {
            pathStack.push(v + 1);
            v = predecessors[v];
        }
        while (!pathStack.empty()) {
            path.push_back(pathStack.top());
            pathStack.pop();
            std::cout << path.back() << " ";
        }
        std::cout << std::endl;

        // Highlight the shortest path in the graph
        for (size_t i = 0; i < path.size() - 1; ++i) {
            Agedge_t *edge = agedge(g, nodes[path[i] - 1], nodes[path[i + 1] - 1], 0, 0);
            if (edge) {
                agsafeset(edge, const_cast<char*>("color"), const_cast<char*>("red"), const_cast<char*>(""));
                agsafeset(edge, const_cast<char*>("penwidth"), const_cast<char*>("3"), const_cast<char*>(""));
            }
        }
        for (int i = 0; i < numNodes; ++i) {
            if (i == sourceNode || i == destinationNode) {
                agsafeset(nodes[i], const_cast<char*>("color"), const_cast<char*>("blue"), const_cast<char*>(""));
                agsafeset(nodes[i], const_cast<char*>("style"), const_cast<char*>("bold"), const_cast<char*>(""));
            }
        }
    }

    // Perform topological sorting
    std::vector<int> topoOrder = topologicalSort(numNodes, edges);
    
    // Find articulation points
    std::vector<int> articulationPoints = findArticulationPoints(numNodes, edges);

    // Print details to the terminal
    std::cout << "Number of nodes: " << numNodes << std::endl;
    std::cout << "Number of edges: " << numEdges << std::endl;
    std::cout << "Topological Order: ";
    for (size_t i = 0; i < topoOrder.size(); ++i) {
        std::cout << (topoOrder[i] + 1) << " ";
    }
    std::cout << std::endl;

    std::cout << "Articulation Points: ";
    for (size_t i = 0; i < articulationPoints.size(); ++i) {
        std::cout << (articulationPoints[i] + 1) << " ";
    }
    std::cout << std::endl;

    // Add description to the graph
    addGraphDescription(g, numNodes, numEdges, sourceNode, destinationNode, minCost, topoOrder, articulationPoints, path);

    // Render the graph to a file
    gvLayout(gvc, g, "dot");
    gvRenderFilename(gvc, g, "png", "graph.png");

    // Cleanup
    gvFreeLayout(gvc, g);
    agclose(g);
    gvFreeContext(gvc);

    std::cout << "Graph has been rendered to graph.png" << std::endl;

    return 0;
}




















