#include "graph.h"
#include "vibes.h"

using namespace std;
using namespace ibex;

Graph::Graph(const IntervalVector &box, ibex::Function *f, Utils *utils, int graph_id=0){
    Pave *p = new Pave(box, f);
    m_node_list.push_back(p);
    m_graph_id = graph_id;
    m_drawing_cpt = 0;
    m_utils = utils;
}

Graph::Graph(Graph* g, int graph_id){
    for(auto &node:g->get_node_list()){
        Pave *p = new Pave(node);
        node->set_copy_node(p);
        m_node_list.push_back(p);
    }

    for(int i=0; i<m_node_list.size(); i++){
        Pave* pave_root = g->get_node_list()[i];
        Pave* pave_copy = m_node_list[i];

        for(int face = 0; face<4; face++){
            for(int j=0; j<pave_root->get_border(face)->get_brothers().size(); j++){
                pave_copy->get_border(face)->set_brother(pave_root->get_border(face)->get_brothers()[j]->get_pave()->get_copy_node()->get_border((face+2)%4), j);
            }
        }
    }

    m_graph_id = graph_id;
    m_drawing_cpt = 0;
    m_utils = g->get_utils();
}

Graph::Graph(Graph* g, Pave* activated_node, int graph_id) : Graph(g, graph_id){

    for(auto &node:m_node_list){
        node->set_empty();
    }

    Pave* copy_node = activated_node->get_copy_node();
    copy_node->set_full();
    *(copy_node) &= *(activated_node);
    for(int face=0; face<4; face++){
        vector<Pave*> brothers_pave = copy_node->get_brothers(face);
        for(auto &pave : brothers_pave){
            if(!pave->is_in_queue()){
                m_node_queue.push_back(pave);
                pave->set_in_queue(true);
            }
        }
    }
}

Graph::~Graph(){
    for(auto &node:m_node_list){
        delete(node);
    }
    for(auto &node:m_node_empty_list){
        delete(node);
    }
}

void Graph::sivia(double epsilon_theta, int iterations_max, bool backward, bool bisect_empty){
    int iterations = 0;
    vector<Pave *> tmp_pave_list(m_node_list);
    m_node_list.clear();

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<iterations_max){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        double diam = tmp->get_theta_diam();

        if(diam < epsilon_theta || tmp->is_empty() && bisect_empty){// || (not_full_test && tmp->is_full() && diam < M_PI)){
            m_node_list.push_back(tmp);
            iterations++;
        }
        else{
            tmp->bisect(tmp_pave_list);
        }
    }

    for(int i=0; i<tmp_pave_list.size(); i++){
        if(backward){
            tmp_pave_list[i]->set_full();
            m_node_queue.push_back(tmp_pave_list[i]);
        }
        m_node_list.push_back(tmp_pave_list[i]);
    }
}

int Graph::process(int max_iterations, bool backward){
    int iterations = 0;
    while(m_node_queue.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = m_node_queue.front();
        m_node_queue.erase(m_node_queue.begin());
        pave->set_in_queue(false);

        bool change = m_utils->CtcContinuity(pave, backward);
        if(change){
            m_utils->CtcPaveConsistency(pave, backward);

            // Warn scheduler to process new pave
            for(int face=0; face<4; face++){
                vector<Pave*> brothers_pave = pave->get_brothers(face);
                for(int i=0; i<brothers_pave.size(); i++){
                    if(brothers_pave[i]->is_in_queue() == false){
                        m_node_queue.push_back(brothers_pave[i]);
                        brothers_pave[i]->set_in_queue(true);
                    }
                }
            }
        }

        //        if(iterations%100 == 0){
        //            cout << iterations << endl;
        //        }
    }

    m_node_queue.clear();
    return iterations;
}

void Graph::set_full(){
    for(auto &node : m_node_list){
        node->set_full();
    }
}

Pave* Graph::get_pave(double x, double y){
    IntervalVector position(2);
    position[0] = Interval(x);
    position[1] = Interval(y);

    for(auto &node:m_node_list){
        if(!(position & node->get_position()).is_empty()){
            return node;
        }
    }
    return NULL;
}

std::vector<Pave*> Graph::get_pave(const ibex::IntervalVector &box){
    std::vector<Pave*> node_list_inter;
    for(auto &node:m_node_list){
        if(!(box & node->get_position()).is_empty()){
            node_list_inter.push_back(node);
        }
    }
    return node_list_inter;
}

void Graph::set_active_pave(const IntervalVector &box){
    vector<Pave*> pave_activated = get_pave(box);

    for(auto &pave : pave_activated){
        pave->set_full();
        for(int face=0; face<4; face++){
            vector<Pave*> pave_brother_list = pave->get_brothers(face);
            for(auto &pave_brother : pave_brother_list){
                if(!pave_brother->is_in_queue()){
                    m_node_queue.push_back(pave_brother);
                    pave_brother->set_in_queue(true);
                }
            }
        }
    }
}

std::vector<Pave*> Graph::get_node_list(){
    return m_node_list;
}

Pave* Graph::get_node_const(int i) const{
    return m_node_list[i];
}

void Graph::draw(int size, bool filled){

    stringstream ss;
    ss << "integralIbex" << m_graph_id<< "-" << m_drawing_cpt;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    for(auto &node:m_node_empty_list){
        node->draw(filled, "gray[]");
    }

    for(auto &node:m_node_list){
        node->draw(filled);
    }
    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    m_drawing_cpt++;
}

void Graph::print_pave_info(double x, double y, string color){

    Pave* p = get_pave(x, y);
    if(p==NULL){
        cout << "PAVE NOT FOUND" << endl;
        return;
    }

    cout << "BOX = " << p->get_position() << endl;
    cout << p << endl;
    for(int i= 0; i<p->get_borders().size(); i++){
        cout << "border " << i << '\t' << p->get_border(i)->get_position() << '\t' << p->get_border(i)->get_segment_in() << p->get_border(i)->get_segment_out()<< endl;
    }
    cout << "theta " << p->get_theta(0) << " " << p->get_theta(1) << endl;

    for(int i=0; i<p->get_borders().size(); i++){
        if(p->get_border(i)->get_brothers().size()!=0){
            for(int j = 0; j<p->get_border(i)->get_brothers().size(); j++){
                cout << "border=" << i << " (" << p->get_border(i) << ") brother=" << j << " " << p->get_border(i)->get_brother(j) << " pave(" << p->get_border(i)->get_brother(j)->get_pave() << ")" << endl;
            }
        }
        else{
            cout << "border=" << i << " (" << p->get_border(i) << ")" << endl;
        }

    }

    double r=0.5*min(p->get_position()[0].diam(), p->get_position()[1].diam())/2.0;

    vibes::drawCircle(p->get_position()[0].mid(), p->get_position()[1].mid(), r, color);

    cout << endl;
}

Utils* Graph::get_utils(){
    return m_utils;
}

void Graph::set_y_symetry(){
    vector<Pave*> pave_list_symetry;

    IntervalVector box_symetry(2);
    box_symetry[0] = Interval::ZERO;
    box_symetry[1] = Interval::ALL_REALS;

    for(auto &pave:m_node_list){
        if(!(pave->get_position() & box_symetry).is_empty()){
            pave_list_symetry.push_back(pave);
        }
    }

    for(auto &pave:pave_list_symetry){
        for(auto &pave_brother:pave_list_symetry){
            if(!(pave->get_position()[1] & -pave_brother->get_position()[1]).is_empty() && pave_brother!=pave){
                if(((pave->get_border(3)->get_segment_in() | -pave_brother->get_border(3)->get_segment_out()) & pave->get_border(3)->get_segment_full())
                        != pave->get_border(3)->get_segment_in()){
                    if(!pave->is_in_queue()){
                        m_node_queue.push_back(pave);
                        pave->set_in_queue(true);
                    }
                }
                pave->get_borders_symetry()[3].set_segment_in(-pave_brother->get_border(3)->get_segment_out(), false);
//                pave->m_borders_symetry[3].set_segment_out(-pave_brother->m_borders[3].segment_in(), false);
            }
        }
    }

}

int Graph::size() const{
    return m_node_list.size();
}

void Graph::remove_empty_node(){
    for(int i=0; i<m_node_list.size(); i++){
        if(m_node_list[i]->is_empty()){
            m_node_list[i]->remove_from_brothers();
            m_node_empty_list.push_back(m_node_list[i]);
            m_node_list.erase(m_node_list.begin() + i);
            i--;
        }
    }
}

bool Graph::is_empty(){
    bool empty = true;
    for(auto &node:m_node_list){
        node->reset_full_empty();
        if(!node->is_empty())
            empty = false;
    }
    return empty;
}

Pave* Graph::get_semi_full_node(){
    for(auto &node:m_node_list){
        if(!node->is_empty() && !node->is_full()){
            return node;
        }
    }

    // Case all full or empty
    for(auto &node:m_node_list){
        if(!node->is_empty()){
            return node;
        }
    }

    return NULL;
}

bool Graph::diff(const Graph &g){
    bool change = false;
    if(this->size() == g.size()){
        for(int i=0; i<g.size(); i++){
            if(m_node_list[i]->diff(*(g.get_node_const(i)))){
                m_node_queue.push_back(m_node_list[i]);
                change = true;
            }
        }
    }
    return change;
}

bool Graph::inter(const Graph &g){
    bool change = false;
    if(this->size() == g.size()){
        for(int i=0; i<g.size(); i++){
            if(m_node_list[i]->inter(*(g.get_node_const(i)))){
                m_node_queue.push_back(m_node_list[i]);
                change = true;
            }
        }
    }
    return change;
}

void Graph::set_empty(){
    for(auto &pave : m_node_list){
        pave->set_empty();
    }
}
