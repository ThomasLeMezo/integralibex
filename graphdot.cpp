#include "graphdot.h"
#include <iostream>
#include <fstream>

GraphDot::GraphDot(Graph *g){
    m_graph = g;
}

void GraphDot::write(string namefile){
    ofstream file;
    file.open(namefile);
    parse_graph(file);
    file.close();
}

void GraphDot::parse_graph(ofstream &file){

    file << "digraph " << "Graph" << m_graph->get_graph_id() << " {" << endl;

    for(int i=0; i<m_graph->get_node_list().size(); i++){
        parse_node_definition(m_graph->get_node_list()[i], file, i);
    }

    for(int i=0; i<m_graph->get_node_list().size(); i++){
        parse_node_link(m_graph->get_node_list()[i], file);
    }

    file << "}" << endl;
}

void GraphDot::parse_node_definition(Pave *p, ofstream &file, int node_id){
    file << "subgraph { " << endl;
    file << "\"" << p << "\"" << " [label=\"" << node_id << " " << p->get_position() << "\"];" << endl;
    for(int face=0; face<4; face++){
        file << "\"" << p->get_border(face) << "\"" << " [label=\"face " << face << "\", shape=box];" << endl;
    }
    for(int face = 0; face < 4; face++){
        file << "\"" << p << "\"" << " -> \"" << p->get_border(face) << "\" [color=blue];" << endl;
    }
    file << "}" << endl;
}

void GraphDot::parse_node_link(Pave *p, ofstream &file){
    for(int face = 0; face < 4; face++){
        for(int i=0; i<p->get_border(face)->get_inclusions().size(); i++){
            file << "\"" << p->get_border(face) << "\"" << " -> \"" << p->get_border(face)->get_inclusion(i).get_border() << "\" ;" << endl;
        }
    }
}



