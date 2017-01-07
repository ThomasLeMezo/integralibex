#include "graphdot.h"
#include <iostream>
#include <fstream>

GraphDot::GraphDot(Graph *g){
    m_graph = g;
}

void GraphDot::write(std::string namefile){
    ofstream file;
    file.open(namefile);
    parse_graph(file);
    file.close();
}

void GraphDot::parse_graph(ofstream &file){

    file << "digraph " << "Graph" << m_graph->get_graph_id() << " {" << endl;

    for(int i=m_graph->get_node_list().size()-1; i>=0; i--){
        parse_node_definition(m_graph->get_node_list()[i], file, i);
    }

    for(int i=0; i<m_graph->get_node_list().size(); i++){
        parse_node_link(m_graph->get_node_list()[i], file);
    }


    file << "}" << endl;
}

void GraphDot::parse_node_definition(Pave *p, ofstream &file, int node_id){
    file << "\"" << p << "\"" << "[shape=record,label=\"<"<<
            p->get_border(3) << "> face 3 \\n "<< p->get_border(3) << "|{ <" <<
            p->get_border(2) << "> face 2 \\n "<< p->get_border(2) << "|{"<<
            node_id << "\\n" << p->get_position() <<"\\n "<< p << "}| <" <<
            p->get_border(0) << "> face 0 \\n "<< p->get_border(0) << "}| <" <<
            p->get_border(1) << "> face 1\\n "<< p->get_border(1) << "\"];" << endl;

}

void GraphDot::parse_node_link(Pave *p, ofstream &file){
    for(int face = 0; face < 2*p->get_dim(); face++){
        for(int i=0; i<p->get_border(face)->get_inclusions().size(); i++){
            file << "\"" << p << "\"" << ":\"" << p->get_border(face) << "\"" << " -> \"" << p->get_border(face)->get_inclusion(i)->get_border()->get_pave() << "\"" << ":\"" << p->get_border(face)->get_inclusion(i)->get_border() << "\" ;" << endl;
        }

        for(int i=0; i<p->get_border(face)->get_inclusions_receving().size(); i++){
            file << "\"" << p << "\"" << ":\"" << p->get_border(face) << "\""
                 << " -> "
                 << "\"" << p->get_border(face)->get_inclusions_receving()[i]->get_owner()->get_pave() << "\"" << ":\"" << p->get_border(face)->get_inclusions_receving()[i]->get_owner() << "\""
                 <<" [color=blue, style=dotted];" << endl;
        }
    }
}



