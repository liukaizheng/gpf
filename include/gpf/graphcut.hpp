#pragma once

#include <vector>
#include <limits>

namespace gpf {

namespace graphcut_detail {
constexpr size_t INVALID_IND_GC = std::numeric_limits<size_t>::max();
constexpr size_t ORPHAN = INVALID_IND_GC - 1;
constexpr size_t TERMINAL = INVALID_IND_GC - 2;
}

struct Arc {
    size_t head;
    size_t next;
    size_t sister;
    double r_cap;

    Arc() : head(graphcut_detail::INVALID_IND_GC), next(graphcut_detail::INVALID_IND_GC), sister(graphcut_detail::INVALID_IND_GC), r_cap(0.0) {}
    
    Arc(size_t h, size_t n, size_t s, double rc) 
        : head(h), next(n), sister(s), r_cap(rc) {}
};

class GraphCut {
public:
    GraphCut(const std::vector<double>& source_cap, const std::vector<double>& sink_cap);

    void add_edge(size_t i, size_t j, double cap, double rev_cap);
    double max_flow();

    std::vector<bool> is_sink;
    double flow;

private:
    void set_active(size_t i);
    size_t next_active();
    void set_orphan_front(size_t i);
    void set_orphan_rear(size_t i);
    void augment(size_t middle_arc_id);
    void process_source_orphan(size_t i);
    void process_sink_orphan(size_t i);

    size_t time;
    size_t queue_first[2];
    size_t queue_last[2];
    size_t orphan_first;
    size_t orphan_last;
    std::vector<double> tr_cap;
    std::vector<Arc> arcs;
    std::vector<size_t> dist;
    std::vector<size_t> first_arc;
    std::vector<size_t> next;
    std::vector<size_t> next_orphan;
    std::vector<size_t> parent;
    std::vector<size_t> ts;
};

template<typename T>
class ArcBuilder {
public:
    std::vector<std::tuple<size_t, size_t, T>> arcs;

    ArcBuilder(const std::vector<T>& source_caps, const std::vector<T>& sink_caps) {
        size_t sink = source_caps.size() + 1;
        for (size_t i = 0; i < source_caps.size(); ++i) {
            arcs.push_back({0, i + 1, source_caps[i]});
        }
        for (size_t i = 0; i < sink_caps.size(); ++i) {
            arcs.push_back({i + 1, sink, sink_caps[i]});
        }
    }

    void add_arc(size_t from, size_t to, T cap, bool rev = false) {
        from += 1;
        to += 1;
        arcs.push_back({from, to, cap});
        if (rev) {
            arcs.push_back({to, from, cap});
        }
    }
};

} // namespace gpf
