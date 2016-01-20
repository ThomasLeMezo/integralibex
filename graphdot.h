#ifndef GRAPHDOT_H
#define GRAPHDOT_H

#include "graph.h"
#include <iostream>
#include <fstream>

class GraphDot
{
public:
    GraphDot(Graph *g);

    void write(string namefile);
    void parse_node_definition(Pave *p, ofstream &file, int node_id);
    void parse_node_link(Pave *p, ofstream &file);
    void parse_graph(ofstream &file);

private:
    Graph *m_graph;
};

#endif // GRAPHDOT_H
