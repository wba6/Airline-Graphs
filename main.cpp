#include <climits>
#include <iostream>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

// Structure to represent an edge in the graph.
struct Edge {
    int u, v;       // u and v are indices into the 'cities' vector (0-indexed)
    int distance;   // distance in miles
    double cost;    // price of the route
};

// Union-Find (Disjoint Set) class to support Kruskal’s algorithm.
class UnionFind {
public:
    std::vector<int> parent;
    std::vector<int> rank;

    UnionFind(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    // Find with path compression.
    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    // Union by rank. Returns false if x and y are already in the same set.
    bool unionSets(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY)
            return false;
        if (rank[rootX] < rank[rootY])
            std::swap(rootX, rootY);
        parent[rootY] = rootX;
        if (rank[rootX] == rank[rootY])
            rank[rootX]++;
        return true;
    }
};

// Graph class: reads the file, stores the cities and edges, prints routes, and computes the MST.
class Graph {
public:
    int numCities;
    std::vector<std::string> cities;
    std::vector<Edge> edges;

    Graph() : numCities(0) {}

    // Reads the graph data from a file.
    bool readFromFile(const std::string &filename) {
        std::ifstream infile(filename);
        if (!infile) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return false;
        }
        std::string line;
        // First line: number of cities.
        getline(infile, line);
        if (line.empty())
            return false;
        numCities = stoi(line);

        // Read city names.
        for (int i = 0; i < numCities; i++) {
            getline(infile, line);
            // Remove trailing carriage return if present
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            cities.push_back(line);
        }


        // Read each route (edge). Format: city1 city2 distance price.
        while (getline(infile, line)) {
            if (line.empty())
                continue;
            std::istringstream iss(line);
            int u, v, dist;
            double price;
            if (!(iss >> u >> v >> dist >> price))
                break;
            Edge e;
            // Adjust indices from 1-indexed to 0-indexed.
            e.u = u - 1;
            e.v = v - 1;
            e.distance = dist;
            e.cost = price;
            edges.push_back(e);
        }
        infile.close();
        return true;
    }

    // Prints all direct routes in a well-formatted manner.
    void printGraph() {
        // Print summary
        std::cout << "There are " << numCities << " cities and "
         << edges.size() << " direct connections." << std::endl << std::endl;

        // For each city, list all routes (edges) that include that city
        for (int i = 0; i < numCities; i++) {
            std::cout << "(" << cities[i] << ")" << std::endl;
            for (auto &edge : edges) {
                if (edge.u == i) {
                    // edge.u is the current city; edge.v is the "other" city
                    std::cout << "... " << cities[i] << "-" << cities[edge.v] << ", "
                     << edge.distance << " miles, $"
                     << std::fixed << std::setprecision(2) << edge.cost << std::endl;
                } else if (edge.v == i) {
                    // edge.v is the current city; edge.u is the "other" city
                    std::cout << "... " << cities[i] << "-" << cities[edge.u] << ", "
                     << edge.distance << " miles, $"
                     << std::fixed << std::setprecision(2) << edge.cost << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    // Computes the MST edges using Kruskal's algorithm (based on distance).
    std::vector<Edge> computeMSTEdges() {
        std::vector<Edge> mst;
        // Make a copy of edges and sort them by distance.
        std::vector<Edge> sortedEdges = edges;
        sort(sortedEdges.begin(), sortedEdges.end(), [](const Edge &a, const Edge &b) {
            return a.distance < b.distance;
        });
        UnionFind uf(numCities);
        for (auto &edge : sortedEdges) {
            if (uf.unionSets(edge.u, edge.v)) {
                mst.push_back(edge);
            }
        }
        return mst;
    }

    std::vector<Edge> computeShortestEdges(std::string startCity, std::string endCity) {
        std::vector<int> dist{INT_MAX};
        std::vector<int> parent{-1};

        // get the integer representations of both cities
        auto startIndex = find(cities.begin(), cities.end(), startCity);
        if (startIndex == cities.end()) {
            std::cerr << "Unable to locate city: " << startCity << std::endl; 
        }

        auto endIndex = find(cities.begin(), cities.end(), endCity);
        if (endIndex == cities.end()) {
            std::cerr << "Unable to locate city: " << endCity << std::endl;
        }

        dist[std::distance(std::begin(cities), startIndex)] = 0;

    }

};

int main(int argc, char* argv[]) {
    std::string filename;
    // If no filename is given as a command-line argument, prompt the user.
    if (argc < 2) {
        std::cout << "Enter the route file name: ";
        std::cin >> filename;
    } else {
        filename = argv[1];
    }

    Graph graph;
    if (!graph.readFromFile(filename)) {
        std::cerr << "Failed to read graph from file." << std::endl;
        return 1;
    }

    // Part I, Query 1: Show the entire list of direct routes.
    graph.printGraph();

    // Part I, Query 2: Compute the MST (minimum spanning tree) based on distances.
    std::vector<Edge> mstEdges = graph.computeMSTEdges();

    // Group MST edges by connected component.
    UnionFind uf(graph.numCities);
    for (auto &edge : mstEdges) {
        uf.unionSets(edge.u, edge.v);
    }
    std::map<int, std::vector<Edge>> components;
    for (auto &edge : mstEdges) {
        int rep = uf.find(edge.u);
        components[rep].push_back(edge);
    }

    // Display the MST(s).
    if (components.size() == 1) {
        std::cout << "Minimum Spanning Tree (based on distances):" << std::endl;
        int totalDistance = 0;
        for (auto &edge : mstEdges) {
            std::cout << graph.cities[edge.u] << " <-> " << graph.cities[edge.v]
                 << ", Distance: " << edge.distance << " miles"
                 << ", Price: $" << edge.cost << std::endl;
            totalDistance += edge.distance;
        }
        std::cout << "Total distance of MST: " << totalDistance << " miles" << std::endl;
    } else {
        int compNum = 1;
        for (auto &comp : components) {
            std::cout << "Minimum Spanning Tree for Component " << compNum << ":" << std::endl;
            int totalDistance = 0;
            for (auto &edge : comp.second) {
                std::cout << graph.cities[edge.u] << " <-> " << graph.cities[edge.v]
                     << " | Distance: " << edge.distance << " miles"
                     << " | Price: $" << edge.cost << std::endl;
                totalDistance += edge.distance;
            }
            std::cout << "Total distance: " << totalDistance << " miles" << std::endl << std::endl;
            compNum++;
        }
    }

    return 0;
}

