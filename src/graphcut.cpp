#include "gpf/graphcut.hpp"
#include <algorithm>

namespace gpf {

GraphCut::GraphCut(const std::vector<double>& source_cap, const std::vector<double>& sink_cap) 
    : time(0), orphan_first(graphcut_detail::INVALID_IND_GC), orphan_last(graphcut_detail::INVALID_IND_GC), flow(0.0) {
    
    size_t n_nodes = source_cap.size();
    queue_first[0] = graphcut_detail::INVALID_IND_GC;
    queue_first[1] = graphcut_detail::INVALID_IND_GC;
    queue_last[0] = graphcut_detail::INVALID_IND_GC;
    queue_last[1] = graphcut_detail::INVALID_IND_GC;
    
    dist.resize(n_nodes, 1);
    is_sink.resize(n_nodes, false);
    next.resize(n_nodes, graphcut_detail::INVALID_IND_GC);
    next_orphan.resize(n_nodes, graphcut_detail::INVALID_IND_GC);
    parent.resize(n_nodes, graphcut_detail::TERMINAL);
    ts.resize(n_nodes, 0);
    first_arc.resize(n_nodes, graphcut_detail::INVALID_IND_GC);
    
    tr_cap.resize(n_nodes);
    for (size_t i = 0; i < n_nodes; ++i) {
        flow += std::min(source_cap[i], sink_cap[i]);
        tr_cap[i] = source_cap[i] - sink_cap[i];
    }
    
    for (size_t i = 0; i < n_nodes; ++i) {
        is_sink[i] = tr_cap[i] < 0.0;
        if (tr_cap[i] == 0.0) {
            parent[i] = graphcut_detail::INVALID_IND_GC;
        } else {
            set_active(i);
        }
    }
}

void GraphCut::set_active(size_t i) {
    if (next[i] != graphcut_detail::INVALID_IND_GC) {
        return;
    }
    if (queue_last[1] == graphcut_detail::INVALID_IND_GC) {
        queue_first[1] = i;
    } else {
        next[queue_last[1]] = i;
    }
    queue_last[1] = i;
    next[i] = i;
}

size_t GraphCut::next_active() {
    while (true) {
        size_t i = queue_first[0];
        if (i == graphcut_detail::INVALID_IND_GC) {
            i = queue_first[1];
            queue_first[0] = queue_first[1];
            queue_last[0] = queue_last[1];
            queue_first[1] = graphcut_detail::INVALID_IND_GC;
            queue_last[1] = graphcut_detail::INVALID_IND_GC;
            if (i == graphcut_detail::INVALID_IND_GC) {
                return graphcut_detail::INVALID_IND_GC;
            }
        }
        if (next[i] == i) {
            queue_first[0] = graphcut_detail::INVALID_IND_GC;
            queue_last[0] = graphcut_detail::INVALID_IND_GC;
        } else {
            queue_first[0] = next[i];
        }
        next[i] = graphcut_detail::INVALID_IND_GC;
        if (parent[i] != graphcut_detail::INVALID_IND_GC) {
            return i;
        }
    }
}

void GraphCut::set_orphan_front(size_t i) {
    parent[i] = graphcut_detail::ORPHAN;
    next_orphan[i] = orphan_first;
    orphan_first = i;
}

void GraphCut::set_orphan_rear(size_t i) {
    parent[i] = graphcut_detail::ORPHAN;
    if (orphan_last != graphcut_detail::INVALID_IND_GC) {
        next_orphan[orphan_last] = i;
    } else {
        orphan_first = i;
    }
    orphan_last = i;
    next_orphan[i] = graphcut_detail::INVALID_IND_GC;
}

void GraphCut::add_edge(size_t i, size_t j, double cap, double rev_cap) {
    size_t arc_id = arcs.size();
    size_t sister_arc_id = arc_id + 1;
    arcs.push_back(Arc(j, first_arc[i], sister_arc_id, cap));
    first_arc[i] = arc_id;
    arcs.push_back(Arc(i, first_arc[j], arc_id, rev_cap));
    first_arc[j] = sister_arc_id;
}

void GraphCut::process_source_orphan(size_t i) {
    size_t a0_min = graphcut_detail::INVALID_IND_GC;
    size_t d_min = graphcut_detail::INVALID_IND_GC;
    size_t a0 = first_arc[i];
    
    while (a0 != graphcut_detail::INVALID_IND_GC) {
        if (arcs[arcs[a0].sister].r_cap == 0.0) {
            a0 = arcs[a0].next;
            continue;
        }
        size_t j = arcs[a0].head;
        if (is_sink[j] || parent[j] == graphcut_detail::INVALID_IND_GC) {
            a0 = arcs[a0].next;
            continue;
        }
        
        size_t d = 0;
        while (true) {
            if (ts[j] == time) {
                d += dist[j];
                break;
            }
            size_t a = parent[j];
            d++;
            if (a == graphcut_detail::TERMINAL) {
                ts[j] = time;
                dist[j] = 1;
                break;
            }
            if (a == graphcut_detail::ORPHAN) {
                d = graphcut_detail::INVALID_IND_GC;
                break;
            }
            j = arcs[a].head;
        }
        
        if (d != graphcut_detail::INVALID_IND_GC) {
            if (d < d_min) {
                a0_min = a0;
                d_min = d;
            }
            j = arcs[a0].head;
            while (true) {
                if (ts[j] == time) {
                    break;
                }
                ts[j] = time;
                dist[j] = d;
                d--;
                j = arcs[parent[j]].head;
            }
        }
        a0 = arcs[a0].next;
    }
    
    parent[i] = a0_min;
    if (parent[i] != graphcut_detail::INVALID_IND_GC) {
        ts[i] = time;
        dist[i] = d_min + 1;
    } else {
        a0 = first_arc[i];
        while (a0 != graphcut_detail::INVALID_IND_GC) {
            size_t j = arcs[a0].head;
            size_t a = parent[j];
            if (!is_sink[j] && a != graphcut_detail::INVALID_IND_GC) {
                if (arcs[arcs[a0].sister].r_cap != 0.0) {
                    set_active(j);
                }
                if (a != graphcut_detail::TERMINAL && a != graphcut_detail::ORPHAN && arcs[a].head == i) {
                    set_orphan_rear(j);
                }
            }
            a0 = arcs[a0].next;
        }
    }
}

void GraphCut::process_sink_orphan(size_t i) {
    size_t a0_min = graphcut_detail::INVALID_IND_GC;
    size_t d_min = graphcut_detail::INVALID_IND_GC;
    size_t a0 = first_arc[i];
    
    while (a0 != graphcut_detail::INVALID_IND_GC) {
        if (arcs[a0].r_cap == 0.0) {
            a0 = arcs[a0].next;
            continue;
        }
        size_t j = arcs[a0].head;
        if (!is_sink[j] || parent[j] == graphcut_detail::INVALID_IND_GC) {
            a0 = arcs[a0].next;
            continue;
        }
        
        size_t d = 0;
        while (true) {
            if (ts[j] == time) {
                d += dist[j];
                break;
            }
            size_t a = parent[j];
            d++;
            if (a == graphcut_detail::TERMINAL) {
                ts[j] = time;
                dist[j] = 1;
                break;
            }
            if (a == graphcut_detail::ORPHAN) {
                d = graphcut_detail::INVALID_IND_GC;
                break;
            }
            j = arcs[a].head;
        }
        
        if (d != graphcut_detail::INVALID_IND_GC) {
            if (d < d_min) {
                a0_min = a0;
                d_min = d;
            }
            j = arcs[a0].head;
            while (true) {
                if (ts[j] == time) {
                    break;
                }
                ts[j] = time;
                dist[j] = d;
                d--;
                j = arcs[parent[j]].head;
            }
        }
        a0 = arcs[a0].next;
    }
    
    parent[i] = a0_min;
    if (parent[i] != graphcut_detail::INVALID_IND_GC) {
        ts[i] = time;
        dist[i] = d_min + 1;
    } else {
        a0 = first_arc[i];
        while (a0 != graphcut_detail::INVALID_IND_GC) {
            size_t j = arcs[a0].head;
            size_t a = parent[j];
            if (is_sink[j] && a != graphcut_detail::INVALID_IND_GC) {
                if (arcs[a0].r_cap != 0.0) {
                    set_active(j);
                }
                if (a != graphcut_detail::TERMINAL && a != graphcut_detail::ORPHAN && arcs[a].head == i) {
                    set_orphan_rear(j);
                }
            }
            a0 = arcs[a0].next;
        }
    }
}

void GraphCut::augment(size_t middle_arc_id) {
    double bottle_neck = arcs[middle_arc_id].r_cap;
    size_t middle_arc_sister_id = arcs[middle_arc_id].sister;
    size_t i = arcs[middle_arc_sister_id].head;
    size_t aid;
    
    while (true) {
        aid = parent[i];
        if (aid == graphcut_detail::TERMINAL) {
            break;
        }
        const Arc& sister = arcs[arcs[aid].sister];
        if (bottle_neck > sister.r_cap) {
            bottle_neck = sister.r_cap;
        }
        i = arcs[aid].head;
    }
    
    if (bottle_neck > tr_cap[i]) {
        bottle_neck = tr_cap[i];
    }
    
    i = arcs[middle_arc_id].head;
    while (true) {
        aid = parent[i];
        if (aid == graphcut_detail::TERMINAL) {
            break;
        }
        const Arc& arc = arcs[aid];
        if (bottle_neck > arc.r_cap) {
            bottle_neck = arc.r_cap;
        }
        i = arc.head;
    }
    
    if (bottle_neck > -tr_cap[i]) {
        bottle_neck = -tr_cap[i];
    }
    
    arcs[middle_arc_sister_id].r_cap += bottle_neck;
    arcs[middle_arc_id].r_cap -= bottle_neck;
    
    i = arcs[middle_arc_sister_id].head;
    while (true) {
        aid = parent[i];
        if (aid == graphcut_detail::TERMINAL) {
            break;
        }
        
        arcs[aid].r_cap += bottle_neck;
        size_t sister_id = arcs[aid].sister;
        arcs[sister_id].r_cap -= bottle_neck;
        if (arcs[sister_id].r_cap == 0.0) {
            set_orphan_front(i);
        }
        i = arcs[aid].head;
    }
    
    tr_cap[i] -= bottle_neck;
    if (tr_cap[i] == 0.0) {
        set_orphan_front(i);
    }
    
    i = arcs[middle_arc_id].head;
    while (true) {
        aid = parent[i];
        if (aid == graphcut_detail::TERMINAL) {
            break;
        }
        size_t sister_id = arcs[aid].sister;
        arcs[sister_id].r_cap += bottle_neck;
        arcs[aid].r_cap -= bottle_neck;
        if (arcs[aid].r_cap == 0.0) {
            set_orphan_front(i);
        }
        i = arcs[aid].head;
    }
    tr_cap[i] += bottle_neck;
    if (tr_cap[i] == 0.0) {
        set_orphan_front(i);
    }
    flow += bottle_neck;
}

double GraphCut::max_flow() {
    size_t current_node = graphcut_detail::INVALID_IND_GC;
    
    while (true) {
        size_t i = current_node;
        if (i != graphcut_detail::INVALID_IND_GC) {
            next[i] = graphcut_detail::INVALID_IND_GC;
            if (parent[i] == graphcut_detail::INVALID_IND_GC) {
                i = graphcut_detail::INVALID_IND_GC;
            }
        }
        if (i == graphcut_detail::INVALID_IND_GC) {
            i = next_active();
            if (i == graphcut_detail::INVALID_IND_GC) {
                break;
            }
        }
        
        size_t aid;
        if (!is_sink[i]) {
            aid = first_arc[i];
            while (aid != graphcut_detail::INVALID_IND_GC) {
                if (arcs[aid].r_cap != 0.0) {
                    size_t j = arcs[aid].head;
                    if (parent[j] == graphcut_detail::INVALID_IND_GC) {
                        is_sink[j] = false;
                        parent[j] = arcs[aid].sister;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                        set_active(j);
                    } else if (is_sink[j]) {
                        break;
                    } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                        parent[j] = arcs[aid].sister;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                    }
                }
                aid = arcs[aid].next;
            }
        } else {
            aid = first_arc[i];
            while (aid != graphcut_detail::INVALID_IND_GC) {
                size_t sister_id = arcs[aid].sister;
                if (arcs[sister_id].r_cap != 0.0) {
                    size_t j = arcs[aid].head;
                    if (parent[j] == graphcut_detail::INVALID_IND_GC) {
                        is_sink[j] = true;
                        parent[j] = sister_id;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                        set_active(j);
                    } else if (!is_sink[j]) {
                        aid = sister_id;
                        break;
                    } else if (ts[j] <= ts[i] && dist[j] > dist[i]) {
                        parent[j] = sister_id;
                        ts[j] = ts[i];
                        dist[j] = dist[i] + 1;
                    }
                }
                aid = arcs[aid].next;
            }
        }
        
        time++;
        if (aid != graphcut_detail::INVALID_IND_GC) {
            next[i] = i;
            current_node = i;
            augment(aid);
            
            i = orphan_first;
            while (i != graphcut_detail::INVALID_IND_GC) {
                size_t next_node = next_orphan[i];
                next_orphan[i] = graphcut_detail::INVALID_IND_GC;
                i = orphan_first;
                while (i != graphcut_detail::INVALID_IND_GC) {
                    orphan_first = next_orphan[i];
                    if (orphan_first == graphcut_detail::INVALID_IND_GC) {
                        orphan_last = graphcut_detail::INVALID_IND_GC;
                    }
                    if (is_sink[i]) {
                        process_sink_orphan(i);
                    } else {
                        process_source_orphan(i);
                    }
                    i = orphan_first;
                }
                orphan_first = next_node;
                i = orphan_first;
            }
        } else {
            current_node = graphcut_detail::INVALID_IND_GC;
        }
    }
    return flow;
}

} // namespace gpf
