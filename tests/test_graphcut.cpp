#include "gpf/graphcut.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    std::vector<double> source_cap = {10.0, 5.0, 15.0};
    std::vector<double> sink_cap = {8.0, 10.0, 7.0};
    
    gpf::GraphCut gc(source_cap, sink_cap);
    
    gc.add_edge(0, 1, 10.0, 10.0);
    gc.add_edge(1, 2, 5.0, 5.0);
    gc.add_edge(0, 2, 15.0, 15.0);
    
    double flow = gc.max_flow();
    
    assert(flow > 0.0);
    std::cout << "Max flow: " << flow << "\n";
    std::cout << "Graph cut tests passed!\n";
    return 0;
}
