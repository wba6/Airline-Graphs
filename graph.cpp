#include "graph.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <map>

/**
* @brief Constructor that initializes n disjoint sets.
* @param n Number of sets (typically, the number of cities).
*/
UnionFind::UnionFind(int n) {
    parent.resize(n);
    rank.resize(n, 0);
    for (int i = 0; i < n; ++i)
        parent[i] = i;
}


/**
 * @brief Finds the representative of the set that contains x.
 * @param x The element to find.
 * @return The representative (root) of the set.
 */
int UnionFind::find(int x) {
    if (parent[x] != x)
        parent[x] = find(parent[x]);
    return parent[x];
}

/**
 * @brief Unites the sets that contain x and y.
 *
 * Uses union by rank.
 *
 * @param x First element.
 * @param y Second element.
 * @return True if the sets were merged; false if they were already in the same set.
 */
bool UnionFind::unionSets(int x, int y) {
    int rootX = find(x), rootY = find(y);
    if (rootX == rootY)
        return false;
    if (rank[rootX] < rank[rootY])
        std::swap(rootX, rootY);
    parent[rootY] = rootX;
    if (rank[rootX] == rank[rootY])
        rank[rootX]++;
    return true;
}

/**
 * @brief Reads the graph data from a file.
 *
 * File format:
 *   - First line: number of cities (N)
 *   - Next N lines: each city name
 *   - Remaining lines: routes in the format:
 *         city1 city2 distance cost
 *     (City numbers are 1-indexed in the file.)
 *
 * @param filename The path to the route file.
 * @return True if reading was successful; false otherwise.
 */
bool Graph::readFromFile(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;
    // Read number of cities.
    if (!std::getline(infile, line)) {
        std::cerr << "File is empty or invalid." << std::endl;
        return false;
    }
    numCities = std::stoi(line);

    // Read city names.
    cities.resize(numCities);
    for (int i = 0; i < numCities; i++) {
        if (!std::getline(infile, line)) {
            std::cerr << "Not enough city names provided." << std::endl;
            return false;
        }
        // Remove any trailing carriage return (Windows line endings).
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        cities[i] = line;
        cityToIndex[line] = i;
    }

    // Initialize adjacency list.
    adjList.assign(numCities, std::vector<EdgeInfo>());

    // Read route data.
    while (std::getline(infile, line)) {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        int c1, c2, dist;
        double price;
        if (!(iss >> c1 >> c2 >> dist >> price)) {
            // Skip malformed lines.
            continue;
        }
        // Adjust from 1-indexed to 0-indexed.
        int u = c1 - 1;
        int v = c2 - 1;

        // Store edge (only once) for MST.
        edges.push_back({u, v, dist, price});

        // For undirected graph, add edge to both u and v lists.
        adjList[u].push_back({v, dist, price});
        adjList[v].push_back({u, dist, price});
    }

    infile.close();
    return true;
}

/**
 * @brief Displays all direct routes in the graph.
 *
 * Each city is printed along with its direct routes.
 */
void Graph::printGraph() const {
    std::cout << "----- Direct Routes -----\n";
    std::cout << "Total Cities: " << numCities << "\n";
    std::cout << "Total Direct Routes: " << edges.size() << "\n\n";
    for (int i = 0; i < numCities; ++i) {
        std::cout << "(" << cities[i] << ")\n";
        // Display all routes for this city.
        for (const auto &edge : adjList[i]) {
            std::cout << "  -> " << cities[i] << " - " << cities[edge.neighbor]
                      << " | Distance: " << edge.distance << " miles"
                      << " | Price: $" << std::fixed << std::setprecision(2) << edge.cost << "\n";
        }
        std::cout << "\n";
    }
}

/**
 * @brief Computes the Minimum Spanning Tree (MST) using Kruskal's algorithm.
 *
 * @return A vector of edges that are in the MST.
 */
std::vector<Edge> Graph::computeMSTEdges() const {
    std::vector<Edge> mst;
    // Copy edges for sorting.
    std::vector<Edge> sortedEdges = edges;
    std::sort(sortedEdges.begin(), sortedEdges.end(), [](const Edge &a, const Edge &b) {
        return a.distance < b.distance;
    });

    UnionFind uf(numCities);
    for (const auto &edge : sortedEdges) {
        if (uf.unionSets(edge.u, edge.v)) {
            mst.push_back(edge);
        }
    }
    return mst;
}

/**
 * @brief Displays the Minimum Spanning Tree (MST) of the graph. Uses Kruskal's
 * 
 * If the graph is disconnected, the MSTs for each connected component are displayed.
 */
void Graph::printMST() const {
    std::vector<Edge> mstEdges = computeMSTEdges();
    // Group MST edges by connected component using UnionFind.
    UnionFind uf(numCities);
    for (const auto &edge : mstEdges) {
        uf.unionSets(edge.u, edge.v);
    }
    std::map<int, std::vector<Edge>> components;
    for (const auto &edge : mstEdges) {
        int comp = uf.find(edge.u);
        components[comp].push_back(edge);
    }

    std::cout << "----- Minimum Spanning Tree(s) (by distance) -----\n";
    if (components.size() == 1) {
        int totalDistance = 0;
        for (const auto &edge : mstEdges) {
            std::cout << "  " << cities[edge.u] << " <-> " << cities[edge.v]
                      << " | Distance: " << edge.distance << " miles"
                      << " | Price: $" << std::fixed << std::setprecision(2) << edge.cost << "\n";
            totalDistance += edge.distance;
        }
        std::cout << "Total MST Distance: " << totalDistance << " miles\n\n";
    } else {
        int compNum = 1;
        for (const auto &entry : components) {
            std::cout << "Component " << compNum << ":\n";
            int totalDistance = 0;
            for (const auto &edge : entry.second) {
                std::cout << "  " << cities[edge.u] << " <-> " << cities[edge.v]
                          << " | Distance: " << edge.distance << " miles"
                          << " | Price: $" << std::fixed << std::setprecision(2) << edge.cost << "\n";
                totalDistance += edge.distance;
            }
            std::cout << "Total Distance: " << totalDistance << " miles\n\n";
            compNum++;
        }
    }
}

/**
 * @brief Prints the shortest path based on total distance between two cities.
 *
 * This function uses the unified Dijkstra's algorithm with edge distance as the weight function.
 *
 * @param srcIndex The index of the source city.
 * @param destIndex The index of the destination city.
 */
void Graph::printShortestPathByDistance(int srcIndex, int destIndex) {
    double totalDistance;
    std::vector<int> parentDistance = unifiedDijkstra(srcIndex, destIndex,
        [](const EdgeInfo &edge) -> double { return edge.distance; },
        totalDistance);
    std::cout << "\n----- Shortest Path by Distance -----\n";
    if (totalDistance == std::numeric_limits<double>::infinity()) {
        std::cout << "No path exists between " << cities.at(srcIndex) << " and " << cities.at(destIndex)<< ".\n";
    } else {
        std::vector<int> path = reconstructPath(parentDistance, destIndex);
        for (std::size_t i = 0; i < path.size(); ++i) {
            std::cout << cities[path[i]];
            if (i < path.size() - 1) {
                double legDistance = 0;
                for (const auto &edge : adjList[path[i]]) {
                    if (edge.neighbor == path[i+1]) {
                        legDistance = edge.distance;
                        break;
                    }
                }
                std::cout << " (" << legDistance << " miles) -> ";
            }
        }

        std::cout << "\nTotal Distance: " << totalDistance << " miles\n";
    }
}

/**
 * @brief Prints the shortest path based on total price between two cities.
 *
 * This function uses the unified Dijkstra's algorithm with edge cost as the weight function.
 *
 * @param srcIndex The index of the source city.
 * @param destIndex The index of the destination city.
 */
void Graph::printShortestPathByPrice(int srcIndex, int destIndex) {
    double totalPrice;
    std::vector<int> parentPrice = unifiedDijkstra(srcIndex, destIndex,
        [](const EdgeInfo &edge) -> double { return edge.cost; },
        totalPrice);
    std::cout << "\n----- Shortest Path by Price -----\n";
    if (totalPrice == std::numeric_limits<double>::infinity()) {
        std::cout << "No path exists between " << cities.at(srcIndex) << " and " << cities.at(destIndex)<< ".\n";
    } else {
        std::vector<int> path = reconstructPath(parentPrice, destIndex);
        for (std::size_t i = 0; i < path.size(); ++i) {
            std::cout << cities[path[i]];
            if (i < path.size() - 1) {
                double legPrice = 0;
                for (const auto &edge : adjList[path[i]]) {
                    if (edge.neighbor == path[i+1]) {
                        legPrice = edge.cost;
                        break;
                    }
                }
                std::cout << " ($" << std::fixed << std::setprecision(2) << legPrice << ") -> ";
            }
        }

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\nTotal Price: $" << totalPrice << "\n";
    }
}

/**
 * @brief Prints the shortest path based on the fewest number of hops (edges) between two cities.
 *
 * This function uses a breadth-first search (BFS) to determine the path with the fewest hops
 * from the source city to the destination city.
 *
 * @param srcIndex The index of the source city.
 * @param destIndex The index of the destination city.
 */
void Graph::printShortestPathByJumps(int srcIndex, int destIndex) {
    // BFS for fewest hops
    std::vector<int> parentHops = shortestPathBFS(srcIndex, destIndex);

    // Reconstruct path using the same method as with Dijkstra
    std::vector<int> bfsPath = reconstructPath(parentHops, destIndex);

    std::cout << "\n----- Shortest Path by Number of Hops (BFS) -----\n";
    if (bfsPath.empty() || bfsPath[0] != srcIndex) {
        std::cout << "No path exists between " << cities.at(srcIndex) << " and " << cities.at(destIndex)<< ".\n";
    } else {
        // Print path
        for (std::size_t i = 0; i < bfsPath.size(); ++i) {
            std::cout << cities[bfsPath[i]];
            if (i < bfsPath.size() - 1)
                std::cout << " -> ";
        }
        // The number of hops is path.size() - 1
        int fewestHops = static_cast<int>(bfsPath.size()) - 1;
        std::cout << "\nTotal Hops: " << fewestHops << "\n";
    }
}

/**
 * @brief Unified Dijkstra algorithm for shortest path queries.
 *
 * This function uses a lambda function (weightFunc) to determine the weight
 * of an edge. It can be used for different metrics:
 *  - For distance queries: weightFunc(edge) = edge.distance.
 *  - For price queries:    weightFunc(edge) = edge.cost.
 *  - For hops queries:     weightFunc(edge) = 1.
 *
 * @param start Index of the starting city.
 * @param goal Index of the destination city.
 * @param weightFunc A lambda that takes an EdgeInfo and returns its weight.
 * @param outTotalCost Returns the total cumulative weight for the computed path.
 * @return A vector representing the parent of each node in the shortest path tree.
 */
std::vector<int> Graph::unifiedDijkstra(int start, int goal, std::function<double(const EdgeInfo&)> weightFunc, double &outTotalCost) {
    std::vector<double> cost(numCities, std::numeric_limits<double>::infinity());
    std::vector<int> parent(numCities, -1);
    cost[start] = 0.0;

    // Priority queue: (cumulative cost, city index)
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    pq.push({0.0, start});

    while (!pq.empty()) {
        auto [currCost, u] = pq.top();
        pq.pop();

        // Skip stale entries.
        if (currCost > cost[u])
            continue;
        if (u == goal)
            break;  // Destination reached.

        // Relax adjacent edges.
        for (const auto &edge : adjList[u]) {
            int v = edge.neighbor;
            double newCost = cost[u] + weightFunc(edge);
            if (newCost < cost[v]) {
                cost[v] = newCost;
                parent[v] = u;
                pq.push({newCost, v});
            }
        }
    }

    outTotalCost = cost[goal];
    return parent;
}

/**
 * @brief Reconstructs the path from the source to destination.
 *
 * Given the parent vector from the shortest path algorithm, this function
 * reconstructs the path from the source city to the destination city.
 *
 * @param parent The parent vector.
 * @param dest The destination city index.
 * @return A vector containing the indices of the cities along the path.
 */
std::vector<int> Graph::reconstructPath(const std::vector<int> &parent, int dest) const {
    std::vector<int> path;
    for (int cur = dest; cur != -1; cur = parent[cur])
        path.push_back(cur);
    std::reverse(path.begin(), path.end());
    return path;
}

 /**
 * @brief Finds and prints all trips starting from the given city within the budget.
 *
 * @param startCity The starting city's name.
 * @param budget The maximum total cost allowed for a trip.
 */
void Graph::findAllTripsFrom(const std::string &startCity, double budget) {
    if (cityToIndex.find(startCity) == cityToIndex.end()) {
        std::cerr << "City " << startCity << " not found.\n";
        return;
    }
    int startIdx = cityToIndex[startCity];
    std::vector<bool> visited(numCities, false);
    std::vector<int> path;
    visited[startIdx] = true;
    path.push_back(startIdx);

    std::cout << "Trips starting from (" << startCity << "):\n";  
    findTrips(startIdx, 0.0, budget, visited, path);
}

 /**
 * @brief Recursively finds and prints all trips (paths) starting from the current city 
 *        that have a total cost within the given budget.
 *
 * @param current The index of the current city.
 * @param currentCost The cumulative cost from the starting city along the current path.
 * @param budget The maximum allowed cost for a trip.
 * @param visited A vector indicating which cities have already been visited.
 * @param path The current trip path (sequence of city indices).
 */
void Graph::findTrips(int current, double currentCost, double budget,
               std::vector<bool> &visited, std::vector<int> &path) const {
    // If path has more than one city and the cost is within the budget, print the trip.
    if (path.size() > 1 && currentCost <= budget) {
        std::cout << "\t";
        for (size_t i = 0; i < path.size(); ++i) {
            std::cout << cities[path[i]];
            if (i < path.size() - 1) {
                double legCost = 0, miles = 0;
                for (const auto &edge : adjList[path[i]]) {
                    if (edge.neighbor == path[i+1]) {
                        legCost = edge.cost;
                        miles = edge.distance;
                        break;
                    }
                }
                std::cout << " ($" << std::fixed << std::setprecision(2) << legCost << ", " << miles << " miles) -> ";
            }
        }
        std::cout << " | Total Cost: $" << std::fixed << std::setprecision(2) << currentCost << "\n";
    }

    // Explore each neighboring city.
    for (const auto &edge : adjList[current]) {
        int neighbor = edge.neighbor;
        double newCost = currentCost + edge.cost;
        // Prune if neighbor is already visited or if adding this edge exceeds the budget.
        if (!visited[neighbor] && newCost <= budget) {
            visited[neighbor] = true;
            path.push_back(neighbor);
            findTrips(neighbor, newCost, budget, visited, path);
            path.pop_back();
            visited[neighbor] = false;
        }
    }
}

/*
 * @breif finds all trips based on a given budget
 *
 * @param budget The macimum allowed cost
 */
void Graph::printAllTripsFromBuget(const double budget) {
    std::cout << "All trips under " << budget << std::endl;

    // Find trips from all cities
    for(const auto& city: cities) {
       findAllTripsFrom(city, budget); 
    }    

}

/**
 * @brief Finds the fewest-hops path from start to goal using BFS.
 *
 * @param start The index of the starting city.
 * @param goal  The index of the destination city.
 * @return A parent vector for path reconstruction.
 */
std::vector<int> Graph::shortestPathBFS(int start, int goal) {
    // Create a visited array and a parent array
    std::vector<bool> visited(numCities, false);
    std::vector<int> parent(numCities, -1);

    // BFS setup
    std::queue<int> q;
    visited[start] = true;
    q.push(start);

    // BFS loop
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == goal) {
            // reached the goal stop BFS
            break;
        }

        // Explore neighbors
        for (auto &edge : adjList[u]) {
            int v = edge.neighbor;
            if (!visited[v]) {
                visited[v] = true;
                parent[v] = u;
                q.push(v);
            }
        }
    }

    return parent;
}


