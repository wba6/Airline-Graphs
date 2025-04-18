#include <functional>
#include <unordered_map>

/**
 * @brief Structure representing a route for MST computation.
 */
struct Edge {
    int u;          ///< Index of the first city (0-indexed)
    int v;          ///< Index of the second city (0-indexed)
    int distance;   ///< Distance in miles
    double cost;    ///< Price of the route
};

/**
 * @brief Structure representing an edge in the adjacency list.
 *
 * Used for shortest-path queries.
 */
struct EdgeInfo {
    int neighbor;   ///< Index of the neighboring city
    int distance;   ///< Distance in miles
    double cost;    ///< Price of the route
};

/**
 * @brief Disjoint-set (Union-Find) data structure for Kruskal's algorithm.
 */
class UnionFind {
public:
    std::vector<int> parent; ///< Parent array for each set
    std::vector<int> rank;   ///< Rank array for union by rank

    /**
     * @brief Constructor that initializes n disjoint sets.
     * @param n Number of sets (typically, the number of cities).
     */
    UnionFind(int n);

    /**
     * @brief Finds the representative of the set that contains x.
     * @param x The element to find.
     * @return The representative (root) of the set.
     */
    int find(int x);

    /**
     * @brief Unites the sets that contain x and y.
     *
     * Uses union by rank.
     *
     * @param x First element.
     * @param y Second element.
     * @return True if the sets were merged; false if they were already in the same set.
     */
    bool unionSets(int x, int y);
};

/**
 * @brief Graph class encapsulating cities and routes.
 *
 * This class reads graph data from a file and provides methods for:
 *  - Displaying direct routes.
 *  - Computing and printing the Minimum Spanning Tree (MST).
 *  - Computing shortest paths using a unified Dijkstra algorithm.
 */
class Graph {
public:
    // Number of cities
    int numCities;                                      
    // City names
    std::vector<std::string> cities;
    //Map from city name to its index
    std::unordered_map<std::string, int> cityToIndex;
    // Adjacency list
    std::vector<std::vector<EdgeInfo>> adjList;
    // Edge list for MST
    std::vector<Edge> edges;

    /**
     * @brief Default constructor.
     */
    Graph() : numCities(0) {}

    /**
     * @brief Reads the graph data from a file.
     *
     * File format:
     *   - First line: number of cities (N)
     *   - Next N lines: each city name
     *   - Remaining lines: routes in the format:
     *         city1 city2 distance cost
     *
     * @param filename The path to the route file.
     * @return True if reading was successful; false otherwise.
     */
    bool readFromFile(const std::string &filename);

    /**
     * @brief Displays all direct routes in the graph.
     *
     * Each city is printed along with its direct routes.
     */
    void printGraph() const;

    /**
     * @brief Computes the Minimum Spanning Tree (MST) using Kruskal's algorithm.
     *
     * @return A vector of edges that are in the MST.
     */
    std::vector<Edge> computeMSTEdges() const;

    /**
     * @brief Displays the Minimum Spanning Tree (MST) of the graph.
     *
     * If the graph is disconnected, the MSTs for each connected component are displayed.
     */
    void printMST() const;

    /*
     * @breif finds all trips based on a given budget and prints them
     *
     * @param budget The macimum allowed cost
     */
    void printAllTripsFromBuget(const double budget);

    /**
     * @brief Prints the shortest path based on total distance between two cities.
     *
     * This function uses the unified Dijkstra's algorithm with edge distance as the weight function.
     *
     * @param srcIndex The index of the source city.
     * @param destIndex The index of the destination city.
     */
    void printShortestPathByDistance(int srcIndex, int destIndex);

    /**
     * @brief Prints the shortest path based on total price between two cities.
     *
     * This function uses the unified Dijkstra's algorithm with edge cost as the weight function.
     *
     * @param srcIndex The index of the source city.
     * @param destIndex The index of the destination city.
     */
    void printShortestPathByPrice(int srcIndex, int destIndex);

    /**
     * @brief Prints the shortest path based on the fewest number of hops (edges) between two cities.
     *
     * This function uses a breadth-first search (BFS) to determine the path with the fewest hops
     * from the source city to the destination city.
     *
     * @param srcIndex The index of the source city.
     * @param destIndex The index of the destination city.
     */
    void printShortestPathByJumps(int srcIndex, int destIndex);

private:
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
    std::vector<int> unifiedDijkstra(int start, int goal, std::function<double(const EdgeInfo&)> weightFunc, double &outTotalCost);

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
    std::vector<int> reconstructPath(const std::vector<int> &parent, int dest) const;
    
    /**
    * @brief Finds and prints all trips starting from the given city within the budget.
    *
    * @param startCity The starting city's name.
    * @param buget The maximum toatal cost allowed for a trip.
    */
    void findAllTripsFrom(const std::string& startCity, double buget);  

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
    void findTrips(int current, double currentCost, double budget,
               std::vector<bool> &visited, std::vector<int> &path) const; 

    /**
     * @brief Finds the fewest-hops path from start to goal using BFS.
     *
     * @param start The index of the starting city.
     * @param goal  The index of the destination city.
     * @return A parent vector for path reconstruction.
     */
    std::vector<int> shortestPathBFS(int start, int goal);
};


