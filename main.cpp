/**
 * @file airline_project.cpp
 * @brief Implementation for Airline Graph Project (Parts 1 and 2).
 *
 * This program reads an airline route file containing the number of cities,
 * the list of city names, and the routes between cities (each with distance
 * and price). It then displays:
 *  - The direct routes (Part 1 Query 1)
 *  - The Minimum Spanning Tree (MST) (Part 1 Query 2)
 *  - The shortest path between two cities based on three metrics:
 *      * Total distance (miles)
 *      * Total price ($)
 *      * Number of hops
 *
 * Dijkstra’s algorithm is used for all three metrics by parameterizing the
 * weight with a lambda function.
 */

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <iomanip>
#include "graph.hpp"

/**
 * @brief Main entry point.
 *
 * This function:
 *  - Prompts the user for the route file name.
 *  - Reads the graph data.
 *  - Displays the direct routes and the MST.
 *  - Prompts for the source and destination cities.
 *  - Computes and displays the shortest path based on distance, price, and hops.
 *
 * @return int Exit status.
 */
int main() {
    Graph graph;
    std::string filename;

    std::cout << "Enter route file name (e.g., airline1.txt): ";
    std::cin >> filename;

    if (!graph.readFromFile(filename)) {
        std::cerr << "Error reading graph data from file. Exiting.\n";
        return 1;
    }

    // Display direct routes and MST.
    std::cout << "\n";
    graph.printGraph();
    graph.printMST();

    // Ask for the source and destination cities.
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear newline
    std::string srcCity, dstCity;
    std::cout << "Enter source city: ";
    std::getline(std::cin, srcCity);
    std::cout << "Enter destination city: ";
    std::getline(std::cin, dstCity);

    // Validate cities.
    if (graph.cityToIndex.find(srcCity) == graph.cityToIndex.end() ||
        graph.cityToIndex.find(dstCity) == graph.cityToIndex.end()) {
        std::cerr << "Error: One or both cities not found in the graph.\n";
        return 1;
    }
    int srcIdx = graph.cityToIndex[srcCity];
    int dstIdx = graph.cityToIndex[dstCity];

    // --- Compute and display shortest path by distance ---
    double totalDistance;
    std::vector<int> parentDistance = graph.unifiedDijkstra(srcIdx, dstIdx,
        [](const EdgeInfo &edge) -> double { return edge.distance; },
        totalDistance);
    std::cout << "\n----- Shortest Path by Distance -----\n";
    if (totalDistance == std::numeric_limits<double>::infinity()) {
        std::cout << "No path exists between " << srcCity << " and " << dstCity << ".\n";
    } else {
        std::vector<int> path = graph.reconstructPath(parentDistance, dstIdx);
        for (std::size_t i = 0; i < path.size(); ++i) {
            std::cout << graph.cities[path[i]];
            if (i < path.size() - 1)
                std::cout << " -> ";
        }
        std::cout << "\nTotal Distance: " << totalDistance << " miles\n";
    }

    // --- Compute and display shortest path by price ---
    double totalPrice;
    std::vector<int> parentPrice = graph.unifiedDijkstra(srcIdx, dstIdx,
        [](const EdgeInfo &edge) -> double { return edge.cost; },
        totalPrice);
    std::cout << "\n----- Shortest Path by Price -----\n";
    if (totalPrice == std::numeric_limits<double>::infinity()) {
        std::cout << "No path exists between " << srcCity << " and " << dstCity << ".\n";
    } else {
        std::vector<int> path = graph.reconstructPath(parentPrice, dstIdx);
        for (std::size_t i = 0; i < path.size(); ++i) {
            std::cout << graph.cities[path[i]];
            if (i < path.size() - 1)
                std::cout << " -> ";
        }
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\nTotal Price: $" << totalPrice << "\n";
    }

    // --- Compute and display shortest path by number of hops ---
    double totalHops;
    std::vector<int> parentHops = graph.unifiedDijkstra(srcIdx, dstIdx,
        [](const EdgeInfo &edge) -> double { return 1.0; },
        totalHops);
    std::cout << "\n----- Shortest Path by Number of Hops -----\n";
    if (totalHops == std::numeric_limits<double>::infinity()) {
        std::cout << "No path exists between " << srcCity << " and " << dstCity << ".\n";
    } else {
        std::vector<int> path = graph.reconstructPath(parentHops, dstIdx);
        for (std::size_t i = 0; i < path.size(); ++i) {
            std::cout << graph.cities[path[i]];
            if (i < path.size() - 1)
                std::cout << " -> ";
        }
        std::cout << "\nTotal Hops: " << totalHops << "\n";
    }

    return 0;
}

