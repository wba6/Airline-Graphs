#include "graph.hpp"
#include <limits>
#include <iomanip>



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
    // 1. Read number of cities.
    if (!std::getline(infile, line)) {
        std::cerr << "File is empty or invalid." << std::endl;
        return false;
    }
    numCities = std::stoi(line);

    // 2. Read city names.
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

    // 3. Initialize adjacency list.
    adjList.assign(numCities, std::vector<EdgeInfo>());

    // 4. Read route data.
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
 * @brief Displays the Minimum Spanning Tree (MST) of the graph.
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

