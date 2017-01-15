#include "graph.h"
#include "vibes.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace ibex;

Graph::Graph(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, Utils *utils, int graph_id, bool diseable_singleton):
    m_search_box(2)
{
    Pave *p = new Pave(box, f_list, diseable_singleton);
    m_count_alive = 1;
    m_search_box = box;
    m_node_list.push_back(p);
    m_graph_id = graph_id;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
    m_compute_inner = false;
    m_inner_mode = false;
    m_positive_invariant = false;
}

Graph::Graph(Utils *utils, int graph_id=0):
    m_search_box(2)
{
    m_count_alive = 0;
    m_graph_id = graph_id;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
    m_compute_inner = false;
    m_inner_mode = false;
    m_positive_invariant = false;
}

Graph::Graph(Graph* g, int graph_id):
    m_search_box(2)
{
#pragma omp parallel for schedule(static) ordered
    for(std::vector<Pave *>::iterator pave = g->get_node_list().begin(); pave < g->get_node_list().end(); ++pave){
        Pave *p_copy = new Pave(*pave);
        (*pave)->set_copy_node(p_copy);
#pragma omp ordered
        m_node_list.push_back(p_copy);
    }

    for(Pave *node_border:g->get_border_list()){
        Pave *p = new Pave(node_border);
        node_border->set_copy_node(p);
        m_node_border_list.push_back(p);
    }

    for(Pave *pave : g->get_node_queue()){
        add_to_all_queue(pave->get_copy_node());
    }

    // For node list
#pragma omp parallel for schedule(dynamic)
    for(int node_nb=0; node_nb<m_node_list.size(); ++node_nb){
        Pave* pave_root = g->get_node_list()[node_nb];
        Pave* pave_copy = pave_root->get_copy_node();

        for(int face = 0; face<4; ++face){
            for(int j=0; j<pave_root->get_border(face)->get_inclusions().size(); ++j){
                Inclusion *i = new Inclusion(pave_root->get_border(face)->get_inclusion(j));
                i->set_border(pave_root->get_border(face)->get_inclusion(j)->get_border()->get_pave()->get_copy_node()->get_border(pave_root->get_border(face)->get_inclusion(j)->get_brother_face()));
                i->set_owner(pave_root->get_border(face)->get_inclusion(j)->get_owner()->get_pave()->get_copy_node()->get_border(pave_root->get_border(face)->get_inclusion(j)->get_owner()->get_face()));

                pave_copy->get_border(face)->add_inclusion(i);
            }
        }
    }

    // For border node list
    for(int node_nb=0; node_nb<m_node_border_list.size(); node_nb++){
        Pave* pave_root = g->get_border_list()[node_nb];
        Pave* pave_copy = m_node_border_list[node_nb];

        for(int face = 0; face<4; face++){
            for(int j=0; j<pave_root->get_border(face)->get_inclusions().size(); j++){
                Inclusion *i = new Inclusion(pave_root->get_border(face)->get_inclusion(j));
                i->set_border(pave_root->get_border(face)->get_inclusion(j)->get_border()->get_pave()->get_copy_node()->get_border(pave_root->get_border(face)->get_inclusion(j)->get_brother_face()));
                i->set_owner(pave_root->get_border(face)->get_inclusion(j)->get_owner()->get_pave()->get_copy_node()->get_border(pave_root->get_border(face)->get_inclusion(j)->get_owner()->get_face()));

                pave_copy->get_border(face)->add_inclusion(i);
            }
        }
    }

    m_graph_id = graph_id;
    m_utils = g->get_utils();
    m_search_box = g->get_search_box();
    m_count_alive = g->get_alive_node();

    debug_marker1 = false;
    debug_marker2 = false;

    m_compute_inner = g->get_compute_inner();
    m_inner_mode = g->get_inner_mode();
    m_positive_invariant = g->get_positive_invariant();
    m_pos_attractor_list = g->get_pos_attractor_list();

    m_inside_curve_list = g->get_inside_curve_list();
}

Graph::Graph(Graph* g, Pave* activated_node, int graph_id) : Graph(g, graph_id){
    //    cout << "COPY GRAPH size = " << size() << endl;
    for(Pave *node:m_node_list){
        node->set_empty();
    }

    Pave* copy_node = activated_node->get_copy_node();
    copy_node->set_full();
    *(copy_node) &= *(activated_node);

    vector<Pave*> brothers_pave = copy_node->get_all_brothers();
    for(Pave *p:brothers_pave){
        if(!p->is_in_queue()){
            add_to_all_queue(p);
        }
    }
}

Graph::~Graph(){
    for(Pave *node:m_node_list){
        delete(node);
    }
    for(Pave *node:m_node_border_list){
        delete(node);
    }
}

void Graph::clear_node_queue(){
    if(!m_compute_inner)
        clear_node_queue_outer();
    else{
        if(m_inner_mode)
            clear_node_queue_inner();
        else
            clear_node_queue_outer();
    }
    for(Pave *node:m_node_list){
        node->set_in_queue(false);
    }
}

void Graph::clear_node_queue_all(){
    clear_node_queue_inner();
    clear_node_queue_outer();

    for(Pave *node:m_node_list){
        node->set_in_queue(false);
    }
}

void Graph::clear_node_queue_inner(){
    for(Pave *node:m_node_list){
        node->set_in_queue_inner(false);
    }
    m_node_queue_inner.clear();
}

void Graph::clear_node_queue_outer(){
#pragma omp parallel for
    for(auto node = m_node_list.begin(); node < m_node_list.end(); ++node){
        (*node)->set_in_queue_outer(false);
    }
    m_node_queue_outer.clear();
}

void Graph::sivia(int nb_node, GRAPH_BW_FW_DIRECTION direction, bool do_not_bisect_empty, bool do_not_bisect_full, bool apply_heuristic){
    //    if(nb_node<=m_count_alive)
    //        return;
    cout << "--> node_list.size() = " << m_node_list.size() << " m_count_alive = " << m_count_alive << endl;
    vector<Pave *> tmp_pave_list(m_node_list);
    m_node_list.clear();
    m_node_list.reserve(nb_node);

    while(tmp_pave_list.size()!=0 && m_count_alive<nb_node){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        if(do_not_bisect_empty || do_not_bisect_full)
            tmp->reset_full_empty();

        if(!tmp->is_active() || tmp->is_removed_pave_union()
                || ((do_not_bisect_empty && tmp->is_empty_inter()) || (do_not_bisect_full && tmp->is_full_inter()))
                || tmp->get_theta_diam_max()<0.0){
            m_node_list.push_back(tmp);
            if(!tmp->get_zone_propagation())
                compute_propagation_zone(tmp);
        }
        else{
            tmp->bisect(tmp_pave_list, direction, apply_heuristic);
            delete(tmp);
            m_count_alive++;
        }
    }

    for(Pave *p:tmp_pave_list){
        if(direction && !p->is_removed_pave()){
            add_to_queue_outer(p);
            if(m_compute_inner){
                add_to_queue_inner(p);
            }
        }
        m_node_list.push_back(p);
    }

    if(m_compute_inner){
#pragma omp parallel for
        for(vector<Pave *>::iterator p=tmp_pave_list.begin(); p<tmp_pave_list.end(); ++p){
            compute_propagation_zone((*p));
        }
    }

    cout << "--> sivia (" << m_node_list.size() << ") outer(" << m_node_queue_outer.size() << ") inner(" <<  m_node_queue_inner.size() << ")" << endl;
}

int Graph::process(int max_iterations, GRAPH_BW_FW_DIRECTION direction, bool union_functions){
    int iterations_final=0;
    int iterations_thread = 0;
    omp_lock_t queue_lock;
    omp_init_lock(&queue_lock);

#pragma omp parallel for private(iterations_thread) schedule(dynamic)
    for(int iterations=0; iterations<max_iterations; iterations++){
        ++iterations_thread;
        if(iterations<max_iterations)
        {
            Pave *pave;

            omp_set_lock(&queue_lock);
            if(!is_empty_node_queue())
            {
                if(omp_get_thread_num()%2==0)
                {
                    pave = get_node_queue_access().front();
                    pop_front_queue();
                }
                else
                {
                    pave = get_node_queue_access().back();
                    pop_back_queue();
                }
            }
            omp_unset_lock(&queue_lock);

            if(pave != NULL)
            {
                pave->lock_pave();
                pave->lock_pave_queue();
                pave->set_in_queue(false);
                pave->unlock_pave_queue();

                /// ******* PROCESS CONTINUITY *******
                bool change;
                pave->lock_border();
                change = m_utils->CtcContinuity(pave, direction);

                pave->increment_cpt_continuity();
                /// TODO improve external inifinity border consistency
                if(pave->is_active()
                        && !pave->is_removed_pave() && !pave->is_removed_pave_outer()
                        && (change || pave->get_first_process()))
                {


                    /// ******* PROCESS CONSISTENCY *******
                    std::vector<bool> change_tab;
                    for(int i=0; i<4; i++)
                        change_tab.push_back(false);


                    if(!pave->is_external_border()) // Temporary (until Ctc ok)
                        m_utils->CtcConsistency(pave, direction, change_tab, union_functions);
                    pave->unlock_border();

                    pave->set_first_process(false);
                    pave->increment_cpt_consistency();

                    /// ******* PUSH BACK NEW PAVES *******
                    // Warn scheduler to process new pave
                    for(int face=0; face<4; face++){
                        if(change_tab[face] || pave->is_external_border()){ // Temporary
                            for(Pave *p:pave->get_brothers(face)){
                                omp_set_lock(&queue_lock);
                                p->lock_pave_queue();
                                if(!p->is_in_queue())
                                    add_to_queue(p);
                                p->unlock_pave_queue();
                                omp_unset_lock(&queue_lock);
                            }
                        }
                    }
                }
                else{
                    pave->unlock_border();
                }
                pave->unlock_pave();
            }
        }

        omp_set_lock(&queue_lock);
        if(is_empty_node_queue()){
            if(iterations<max_iterations)
                iterations_final+=iterations_thread;
            iterations=max_iterations;
        }
        omp_unset_lock(&queue_lock);
    }

    clear_node_queue();
    if(get_inner_mode())
        cout << "--> processing inner (" << iterations_final << ")" << endl;
    else
        cout << "--> processing outer (" << iterations_final << ")" << endl;
    //    myfile.close();
    return 0;
}

void Graph::set_full(){
#pragma omp parallel for
    for(vector<Pave *>::iterator node = m_node_list.begin(); node<m_node_list.end(); ++node){
        if((*node)->is_active() && !(*node)->is_removed_pave())
            (*node)->set_full();
    }
}

void Graph::add_all_to_queue(){
    m_node_queue_outer.clear();
    if(m_compute_inner)
        m_node_queue_inner.clear();

    for(Pave *node: m_node_list){
        m_node_queue_outer.push_back(node);
        if(m_compute_inner)
            m_node_queue_inner.push_back(node);
    }
}

Pave* Graph::get_pave(double x, double y) const{
    IntervalVector position(2);
    position[0] = Interval(x);
    position[1] = Interval(y);

    for(Pave *node:m_node_list){
        if(!(position & node->get_position()).is_empty()){
            return node;
        }
    }
    for(Pave *node:m_node_border_list){
        if(!(position & node->get_position()).is_empty()){
            return node;
        }
    }
    return NULL;
}

const std::vector<Pave *> Graph::get_pave(const ibex::IntervalVector &box) const{
    std::vector<Pave*> node_list_inter;
    for(Pave *node:m_node_list){
        if(!(box & node->get_position()).is_empty()){
            node_list_inter.push_back(node);
        }
    }
    return node_list_inter;
}

void Graph::initialize_queues_with_initial_condition(const ibex::IntervalVector &box){
    std::vector<ibex::IntervalVector> box_list;
    box_list.push_back(box);
    initialize_queues_with_initial_condition(box_list);
}

void Graph::initialize_queues_with_initial_condition(const std::vector<ibex::IntervalVector> &box_list){
    clear_node_queue_inner();
    clear_node_queue_outer();

#pragma omp parallel for
    for(vector<Pave*>::iterator pave=m_node_list.begin(); pave<m_node_list.end(); ++pave){
        (*pave)->reset_full_empty();
    }

    //    std::for_each(m_node_list.begin(), m_node_list.end(), std::bind(&Pave::reset_full_empty, std::placeholders::_1));

#pragma omp parallel for
    for(std::vector<Pave*>::iterator pave = m_node_list.begin(); pave<m_node_list.end(); ++pave){
        //    for(Pave *pave:m_node_list){
        if((*pave)->is_active() /*&& !pave->is_removed_pave_outer()*/){
            for(IntervalVector box:box_list){
                if(!(box & (*pave)->get_position()).is_empty()){

                    // Outer
                    (*pave)->set_full_outer();

                    // Inner
                    if((*pave)->get_position().is_strict_interior_subset(box)){
                        if(!(*pave)->is_empty_inner()){
                            (*pave)->set_empty_inner_in(); // Do not set removed pave inner !!! => bc inner out is not empty
                            // pave->set_removed_pave_inner(true);
                            add_to_queue_inner((*pave));
                        }
                        (*pave)->set_bassin(true);
                    }
                }
            }
        }
    }

#pragma omp parallel for
    for(std::vector<Pave*>::iterator pave = m_node_list.begin(); pave<m_node_list.end(); ++pave){
        // Inner pave
        //        if(pave->is_removed_pave_inner() && !pave->is_removed_pave_outer()){
        //            for(int face=0; face<4; face++){
        //                vector<Pave*> pave_brother_list = pave->get_brothers(face);
        //                for(Pave *pave_brother:pave_brother_list){
        //                    if(!pave->is_removed_pave_outer() && !pave->is_removed_pave_inner())
        //                        add_to_queue_inner(pave_brother);
        //                }
        //            }
        //        }
        if(!(*pave)->is_removed_pave_outer())
            add_to_queue_inner((*pave));

        // Outer pave
        if((*pave)->is_full_outer()){
            for(int face=0; face<4; face++){
                vector<Pave*> pave_brother_list = (*pave)->get_brothers(face);
                for(Pave *pave_brother:pave_brother_list){
                    if(!pave_brother->is_full_outer() && !(*pave)->is_removed_pave_outer())
                        add_to_queue_outer(pave_brother);
                }
            }
        }
    }

    // Temp
    //    for(Pave *pave : m_node_border_list){
    //        pave->set_full_outer();
    //        pave->set_full_inner();
    //        add_to_queue_outer(pave);
    //    }
}

void Graph::set_active_pave(const IntervalVector &box){
    vector<Pave*> pave_activated = get_pave(box);

    for(Pave *pave:pave_activated){
        pave->set_full_outer();
        for(int face=0; face<4; face++){
            vector<Pave*> pave_brother_list = pave->get_brothers(face);
            for(Pave *pave_brother : pave_brother_list){
                if(!pave_brother->is_in_queue_outer()){
                    add_to_queue_outer(pave_brother);
                }
            }
        }
    }
}

std::vector<Pave *> &Graph::get_node_list() {
    return m_node_list;
}

std::vector<Pave *>& Graph::get_border_list(){
    return m_node_border_list;
}

void Graph::push_back(Pave* p){
    m_node_list.push_back(p);
    if(p->is_active())
        m_count_alive++;
}

void Graph::push_back_external_border(Pave *p){
    m_node_border_list.push_back(p);
}

const std::list<Pave *>& Graph::get_node_queue() const {
    if(!m_compute_inner)
        return m_node_queue_outer;
    else{
        if(m_inner_mode)
            return m_node_queue_inner;
        else
            return m_node_queue_outer;
    }
}

std::list<Pave *>& Graph::get_node_queue_access(){
    if(!m_compute_inner)
        return m_node_queue_outer;
    else{
        if(m_inner_mode)
            return m_node_queue_inner;
        else
            return m_node_queue_outer;
    }
}

Pave* Graph::get_node_const(int i) const{
    return m_node_list[i];
}

Pave& Graph::operator[](int id){
    return *(m_node_list[id]);
}

void Graph::draw(int size, bool filled, string comment, bool inner_only, int position, bool pos_invariant){
    // Magenta = #FF00FF
    // Gray light =  #D3D3D3
    // Blue = #4C4CFF

    stringstream ss;
    ss << "cycleSOLVER" << m_graph_id << " " << comment;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",position,"y",position,"width",size,"height",size));

    for(Pave *node:m_node_list){
        node->draw(filled, inner_only);
    }

    for(Pave *node:m_node_border_list){
        node->draw(filled);
    }

    // Pos Invariant
    if(pos_invariant){
        if(m_pos_attractor_list.size()!=0){
            for(vector<vector< vector<IntervalVector>>> &attractor:m_pos_attractor_list){
                for(vector< vector<IntervalVector>> &segment_list:attractor){
                    for(vector<IntervalVector> &segment:segment_list){
                        vector<double> x, y;
                        for(IntervalVector &pt:segment){
                            x.push_back(pt[0].mid());
                            y.push_back(pt[1].mid());
                        }
                        vibes::drawLine(x, y, "g",vibesParams("LineWidth",10.0));
                    }
                }
            }
        }
    }

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));

    vibes::axisLimits(get_search_box()[0].lb()-1.0,get_search_box()[0].ub()+1.0, get_search_box()[1].lb()-1.0,get_search_box()[1].ub()+1.0);
}

void Graph::drawInner(bool filled){
    for(Pave *node:m_node_list){
        node->draw(filled, true);
    }
}

void Graph::print_pave_info(double x, double y, string color) const{

    Pave* p = get_pave(x, y);
    if(p==NULL){
        cout << "PAVE NOT FOUND" << endl;
        return;
    }

    cout << "*******************" << endl;
    cout << "BOX = " << p->get_position() << endl;
    cout << p << endl;
    cout << "active = " << p->is_active() << endl;
    p->print_theta_list();
    cout << "inner mode = " << p->get_inner_mode() << endl;
    p->reset_full_empty();
    cout << "is_full_inner = " << p->is_full_inner() << " is_empty_inner = " << p->is_empty_inner() << endl;
    cout << "is_full_outer = " << p->is_full_outer() << " is_empty_outer = " << p->is_empty_outer() << endl;
    cout << "nb\t" << "in_inner\t" << "out_inner\t" << "in_outer\t" << "out_outer\t" << endl;
    cout << "cpt_continuity_inner = " << p->get_cpt_continuity_inner() << " cpt_continuity_outer = " << p->get_cpt_continuity_outer() << endl;
    cout << "cpt_consistency_inner = " << p->get_cpt_consistency_inner() << " cpt_consistency_outer = " << p->get_cpt_consistency_outer() << endl;
    cout << "id, \tin_inner, out_inner, in_outer, out_outer" << endl;
    for(int i= 0; i<p->get_borders().size(); i++){
        cout << i << '\t'
                //<< p->get_border(i)->get_position() << "       \t"
             << p->get_border(i)->get_segment_in_inner() << '\t'
             << p->get_border(i)->get_segment_out_inner() << '\t'
             << p->get_border(i)->get_segment_in_outer() << '\t'
             << p->get_border(i)->get_segment_out_outer() << '\t'
             << endl;
    }


    for(int i=0; i<4; i++){
        if(p->get_border(i)->get_inclusions().size()!=0){
            for(int j = 0; j<p->get_border(i)->get_inclusions().size(); j++){
                cout << "border=" << i << " (" << p->get_border(i)
                     << ") brother=" << j << " " << p->get_border(i)->get_inclusion(j)->get_border()
                     << " pave(" << p->get_border(i)->get_inclusion(j)->get_border()->get_pave() << ")"
                     << endl;
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

void Graph::print() const{
    cout << "********" << endl;
    cout << "GRAPH id= " << m_graph_id << endl;
    for(Pave *node:m_node_list){
        node->print();
    }
}

Utils* Graph::get_utils(){
    return m_utils;
}

int Graph::size() const{
    return m_node_list.size();
}

void Graph::mark_empty_node(){
#pragma omp parallel for
    for (vector<Pave*>::iterator pave = m_node_list.begin(); pave < m_node_list.end(); ++pave){
        if((*pave)->is_active()){

            // Removing pave that satisfies a special condition (given by a function < 0)
            bool removed_inside_curve = false;
            if(m_inside_curve_list.size()>0 && (*pave)->is_full_outer()){
                for(ibex::Function *f_curve:m_inside_curve_list){
                    Interval result = f_curve->eval((*pave)->get_position());
                    if(result.is_subset(Interval::NEG_REALS)){
                        //                        removed_inside_curve = true;
                        (*pave)->set_empty_inner();
                        (*pave)->set_bassin(true);
                        (*pave)->set_external_border(true);
                    }
                }
            }

            // Analyze if a (*pave) is empty and should be removed
            bool test_two_pi = (*pave)->is_theta_more_than_two_pi();
            if(!test_two_pi || removed_inside_curve){
                (*pave)->reset_full_empty();
                bool empty_outer = false;
                bool empty_inner = false;

                // Outer
                if((*pave)->is_removed_pave_outer() || (*pave)->is_empty_outer()){
                    (*pave)->set_removed_pave_outer(true);
                    empty_outer = true;
                }

                // Inner
                if(m_compute_inner && ((*pave)->is_removed_pave_inner() || (*pave)->is_empty_inner())){
                    (*pave)->set_removed_pave_inner(true);
                    empty_inner = true;
                }

                if(empty_outer || (m_compute_inner && empty_inner)){
                    (*pave)->set_active(false);
                    m_count_alive--;
                }
            }
            else if(test_two_pi){
                (*pave)->set_full();
            }
        }
    }
}

bool Graph::is_empty(){
    if(get_alive_node()<=0)
        return true;
    for(Pave *node:m_node_list){
        node->reset_full_empty();
        if(!node->is_empty_inter())
            return false;
    }
    m_count_alive = 0;
    return true;
}

Pave* Graph::get_semi_full_node(){
    int end = m_node_list.size();
    int start = ceil(end/1.2);

    for(int i=start; i!=(start-1)%end; i = (i+1)%end){
        Pave *node = m_node_list[i];
        //    for(Pave *node:m_node_list){
        if(!node->is_removed_pave() && node->is_border() && node->get_theta_diam()<M_PI/2.0){
            cout << "FIRST RETURN, START = " << start << " i= "<< i << endl;
            return node;
        }
    }

    //    for(Pave *node:m_node_list){
    for(int i=start; i!=(start-1)%end; i = (i+1)%end){
        Pave *node = m_node_list[i];
        if(!node->is_removed_pave() && !node->is_empty() && !node->is_full()){
            cout << "2ND RETURN" << endl;
            return node;
        }
    }

    // Case all full or empty
    cout << "WARNING - get_semi_full_node : ALL FULL/EMPTY NODE" << endl;
    for(Pave *node:m_node_list){
        if(!node->is_removed_pave() && !node->is_empty()){
            return node;
        }
    }

    return NULL;
}

void Graph::diff(const Graph &g){
    if(this->size() == g.size()){
#pragma omp parallel for
        for(int i=0; i<g.size(); ++i){
            m_node_list[i]->diff(*(g.get_node_const(i)));
        }
    }
}

void Graph::inter(const Graph &g, bool with_bwd){
    if(this->size() == g.size()){
#pragma omp parallel for
        for(int i=0; i<g.size(); ++i){
            m_node_list[i]->inter(*(g.get_node_const(i)), with_bwd);
        }
    }
}

void Graph::set_empty(){
    for(Pave *node : m_node_list){
        if(node->is_active() && !node->is_removed_pave())
            node->set_empty();
    }
    for(Pave *node : m_node_border_list){
        node->set_empty();
    }
}

void Graph::set_empty_outer_full_inner(){
#pragma omp parallel for
    for(vector<Pave*>::iterator pave = m_node_list.begin(); pave < m_node_list.end(); ++pave){
        if(!(*pave)->is_removed_pave_outer())
            (*pave)->set_empty_outer();
        if(!(*pave)->is_removed_pave_inner()){
            (*pave)->set_full_inner();
        }
        else{
            (*pave)->set_full_outer();
        }

        if((*pave)->is_active())
            (*pave)->set_first_process_all(true);
    }

    for(Pave *pave : m_node_border_list){
        pave->set_empty_outer();
        pave->set_full_inner();
        pave->set_first_process_all(true);
    }
}

void Graph::set_symetry(Function* f, int face_in, int face_out){
    Inclusion *i = new Inclusion(m_node_list[0]->get_border(face_in), f, face_in);
    Pave *p = m_node_list[0];
    // ToDo : not optimal (remove all other inclusions but should only remove ones seleced by user)
    p->get_border(face_out)->remove_inclusion(-1);
    p->get_border(face_out)->add_inclusion(i);
}

//void Graph::set_all_first_process(){
//    for(Pave *pave:m_node_list){
//        pave->set_first_process(true);
//    }
//}

int Graph::get_graph_id(){
    return m_graph_id;
}

IntervalVector Graph::get_bounding_box() const{
    IntervalVector boundingBox(2);
    boundingBox[0] = Interval::EMPTY_SET;
    boundingBox[1] = Interval::EMPTY_SET;

    for(Pave *p:m_node_list){
        boundingBox |= p->get_position();
    }

    return boundingBox;
}

void Graph::build_graph(){

    // ToDo : To Improve !!!
    for(int i=0; i<m_node_list.size()-1; i++){
        Pave *p1 = m_node_list[i];
        for(int j=i+1; j<m_node_list.size(); j++){
            Pave *p2 = m_node_list[j];

            IntervalVector inter = p1->get_position() & p2->get_position();
            //            cout << "== p1=" << p1->get_position() << '\t' << "p2=" << p2->get_position() << '\t' << "result="<< inter << '\t' << "test=" << !inter.is_empty() << endl;
            if(!inter.is_empty()){
                // Find common border

                for(int face = 0; face < 4; face++){
                    Border *b1 = p1->get_border(face);
                    Border *b2 = p2->get_border((face+2)%4);

                    IntervalVector inter_b = b1->get_position() & b2->get_position();

                    if(!inter_b.is_empty() ){
                        //                        cout << "b1 = " << b1->get_position() << '\t' << "b2 = " << b2->get_position() << '\t' << "inter = " << inter_b << '\t' << "test = " << inter_b.is_empty() << endl;

                        Inclusion *inclusion_to_b1 = new Inclusion(b1, b1->get_face());
                        if(!b2->add_inclusion(inclusion_to_b1))
                            delete(inclusion_to_b1);

                        Inclusion *inclusion_to_b2 = new Inclusion(b2, b2->get_face());
                        if(!b1->add_inclusion(inclusion_to_b2))
                            delete(inclusion_to_b2);
                    }
                }
            }
        }
    }

    for(int i = 0; i<m_node_list.size(); i++){
        if(m_node_list[i]->is_external_border()){
            m_node_border_list.push_back(m_node_list[i]);
            m_node_list.erase(m_node_list.begin()+i);
            i--;
        }
    }
}

void Graph::update_queue(bool border_condition, bool empty_condition){
    clear_node_queue();
    omp_lock_t add_queue;
    omp_init_lock(&add_queue);
#pragma omp parallel for
    for(vector<Pave *>::iterator p=m_node_list.begin(); p<m_node_list.end(); ++p){
        (*p)->set_first_process(true);
        (*p)->set_in_queue(false);
        (*p)->reset_full_empty();

        if((*p)->is_active() && (
                    (!(*p)->is_full() && !(*p)->is_empty())
                    || (border_condition && (*p)->is_border())
                    || (empty_condition && (*p)->is_near_empty()))){
            omp_set_lock(&add_queue);
            add_to_queue((*p));
            omp_unset_lock(&add_queue);
        }
    }
    if(get_inner_mode())
        cout << "--> queue inner (" << m_node_queue_inner.size() << ")" << endl;
    else
        cout << "--> queue outer (" << m_node_queue_outer.size() << ")" << endl;
}

int Graph::get_f_size() const{
    if(m_node_list.size()>0){
        return m_node_list[0]->get_f_list().size();
    }
    else{
        return -1;
    }
}

void Graph::set_active_f(int id){
#pragma omp parallel for
    for(vector<Pave *>::iterator p=m_node_list.begin(); p<m_node_list.end(); ++p){
        (*p)->set_active_function(id);
    }
}

bool Graph::identify_attractor(){
    /// ToDo : correction of angle condition : should be real angle instead of M_PI !!!
    /// Get the angle of the transversal
    if(m_node_list.size()==0)
        return true;
    reset_marker_attractor();
    bool found_all_attractor = true;

    for(Pave *p:m_node_list){

        // Study the subgraph starting with p
        if(p->is_active() && !p->is_marked_attractor() && !p->is_removed_pave()){
            vector<Pave*> attractor_list;
            attractor_list.push_back(p);
            get_recursive_attractor(p, attractor_list);

            // Test if the subgraph is an attractor
            bool test_attractor = true;

            IntervalVector bounding_box(2, Interval::EMPTY_SET);
            for(Pave *p:attractor_list){
                if(!p->is_positive_invariant()){
                    test_attractor = false;
                }
                bounding_box |= p->get_bounding_pave();
            }

            vector<IntervalVector> list_search_box = m_utils->get_segment_from_box(m_search_box);
            for(IntervalVector &iv:list_search_box){
                if(!((bounding_box & iv).is_empty()))
                    test_attractor = false;
            }

            if(test_attractor){
                cout << "FIND AN ATTRACTOR" << endl;
                for(Pave *p:attractor_list){
                    p->set_active(false);
                }
            }
            else{
                found_all_attractor = false;
            }
        }
    }
    return found_all_attractor;
}

void Graph::get_recursive_attractor(Pave* p, vector<Pave*> &list){
    for(int face=0; face<4; face++){
        vector<Pave*> brothers_face = p->get_brothers(face);
        for(Pave *p_brother:brothers_face){
            if(!p_brother->is_removed_pave() && !p_brother->is_marked_attractor()){
                list.push_back(p_brother);
                p_brother->set_marker_attractor(true);
                get_recursive_attractor(p_brother, list);
            }
        }
    }
}

void Graph::reset_marker_attractor(){
    for(Pave *p:m_node_list){
        p->set_marker_attractor(false);
    }
}

void Graph::reset_marker(vector<Pave*> list){
    for(Pave *p:list){
        p->set_marker(false);
    }
}

void Graph::set_marker(vector<Pave*> list, bool val){
    for(Pave *p:list){
        p->set_marker(val);
    }
}

IntervalVector Graph::get_search_box() const{
    return m_search_box;
}

void Graph::complementaire(){
    for(Pave *p:m_node_list){
        p->complementaire();
    }
}

void Graph::set_external_boundary(bool in, bool out){
    if(m_node_border_list.size()==0){
        cout << "ERROR : border_list.size()==0" << endl;
        exit(1);
    }
    cout << m_node_border_list.size() << " boundary boxes" << endl;
    for(Pave *p:m_node_border_list){
        p->set_segment(in, out);
    }
}

void Graph::set_all_active(){
    for(Pave *p:m_node_list){
        p->set_active(true);
        p->set_removed_pave(false);
    }
}

int Graph::get_alive_node() const{
    return m_count_alive;
}

bool Graph::is_no_active_function(){
    for(Pave *p:m_node_list){
        if(p->get_active_function()!=-1){
            return false;
        }
    }
    return true;
}

void Graph::add_to_all_queue(Pave *p){
    add_to_queue_outer(p);
    add_to_queue_inner(p);
}

bool Graph::get_compute_inner(){
    return m_compute_inner;
}

bool Graph::get_inner_mode(){
    return m_inner_mode;
}

void Graph::set_compute_inner(bool val){
    m_compute_inner = val;
#pragma omp parallel for
    for(vector<Pave*>::iterator pave = m_node_list.begin(); pave < m_node_list.end(); ++pave){
        (*pave)->set_compute_inner(val);
    }
}

void Graph::set_inner_mode(bool val){
    if(!m_compute_inner)
        set_compute_inner(true);
    m_inner_mode = val;
    for(Pave *p:m_node_list){
        p->set_inner_mode(val);
        p->set_first_process(true);
    }
    for(Pave *p:m_node_border_list){
        p->set_inner_mode(val);
    }
}

void Graph::add_to_queue_inner(Pave *p){
    if(!p->is_in_queue_inner()){
        m_node_queue_inner.push_back(p);
        p->set_in_queue_inner(true);
    }
}

void Graph::add_to_queue_outer(Pave *p){
    if(!p->is_in_queue_outer()){
        m_node_queue_outer.push_back(p);
        p->set_in_queue_outer(true);
    }
}

bool Graph::is_empty_node_queue(){
    if(!m_compute_inner)
        return (m_node_queue_outer.size() == 0);
    else{
        if(m_inner_mode)
            return (m_node_queue_inner.size() == 0);
        else
            return (m_node_queue_outer.size() == 0);
    }
}

void Graph::add_to_queue(Pave *p){
    if(!m_compute_inner)
        add_to_queue_outer(p);
    else{
        if(m_inner_mode)
            add_to_queue_inner(p);
        else
            add_to_queue_outer(p);
    }
}

void Graph::pop_front_queue(){
    if(!m_compute_inner){
        m_node_queue_outer.pop_front();
    }
    else{
        if(m_inner_mode){
            m_node_queue_inner.pop_front();
        }
        else{
            m_node_queue_outer.pop_front();
        }
    }
}

void Graph::pop_back_queue(){
    if(!m_compute_inner){
        m_node_queue_outer.pop_back();
    }
    else{
        if(m_inner_mode){
            m_node_queue_inner.pop_back();
        }
        else{
            m_node_queue_outer.pop_back();
        }
    }
}

void Graph::copy_to_inner(){
    for(Pave *p:m_node_list)
        p->copy_to_inner();
}

void Graph::compute_propagation_zone(Pave *p, bool compute_anyway){
    if(!compute_anyway && (p->get_f_list().size()==1 || p->get_zone_propagation()))
        return;
    bool bwd_function_save = p->get_backward_function();
    for(int fwd=0; fwd<2; fwd++){
        Pave *p_copy = new Pave(p);
        p_copy->set_inner_mode(false);
        p_copy->set_compute_inner(false);

        p_copy->set_backward_function(fwd);
        p->set_backward_function(fwd);
        for(int id_function = 0; id_function<p_copy->get_f_list().size(); id_function++){
            p_copy->set_active_function(id_function);
            p_copy->set_full();
            vector<bool> change_tab;
            for(int i=0; i<4; i++)
                change_tab.push_back(false);
            m_utils->CtcPaveBackward(p_copy, true, change_tab);

            for(int face = 0; face<4; face++){
                if(!p_copy->get_border(face)->get_segment_in().is_empty())
                    p->get_border(face)->push_back_zone_function_in(true);
                else
                    p->get_border(face)->push_back_zone_function_in(false);

                if(!p_copy->get_border(face)->get_segment_out().is_empty())
                    p->get_border(face)->push_back_zone_function_out(true);
                else
                    p->get_border(face)->push_back_zone_function_out(false);
            }
        }
    }
    p->set_zone_propagation(true);
    p->set_backward_function(bwd_function_save);
}

void Graph::compute_all_propagation_zone(bool compute_anyway){
    for(Pave *p:m_node_list)
        compute_propagation_zone(p, compute_anyway);
}

void Graph::reset_computation_zone(){
    for(Pave *p:m_node_list)
        p->reset_computation_zone();
}

void Graph::set_backward_function(bool val){
    for(Pave *p:m_node_list)
        p->set_backward_function(val);
}

double Graph::get_area_outer(){
    double area = 0.0;
    for(Pave *p:m_node_list){
        if(!p->is_removed_pave_outer()){
            area += p->get_area_outer();
        }
    }
    return area;
}

void Graph::get_recursive_zone(Pave* p, vector<Pave*> &list){
    for(int face=0; face<4; face++){
        vector<Pave*> brothers_face = p->get_brothers(face);
        for(Pave *p_brother:brothers_face){
            if(!p_brother->is_removed_pave() && !p_brother->is_marked()){
                list.push_back(p_brother);
                p_brother->set_marker(true);
                get_recursive_zone(p_brother, list);
            }
        }
    }
}

void Graph::get_recursive_contour(Pave* p, vector<Pave*> &list){
    for(int face=0; face<4; face++){

        for(Inclusion *i:p->get_border(face)->get_inclusions()){
            Pave *p_brother = i->get_border()->get_pave();
            // Brother pave not full nor empty
            if(!p_brother->is_removed_pave() && !p_brother->is_marked()
                    && !p_brother->is_empty() && !p_brother->is_full()){

                // Link between paves
                Interval seg_inter = p->get_border(face)->get_segment_in_union_out()
                        & i->get_border()->get_segment_in_union_out();
                if(seg_inter!=Interval::EMPTY_SET){
                    list.push_back(p_brother);
                    p_brother->set_marker(true);
                    get_recursive_contour(p_brother, list);
                }
            }
        }
    }
}

vector<vector<Pave*>> Graph::get_contour_nodes(){
    vector<vector<Pave*>> contours_list;

    if(m_node_list.size()==0)
        return contours_list;
    reset_marker(m_node_list);

    for(Pave *p:m_node_list){
        // Study the subgraph starting with p
        if(p->is_active() && !p->is_marked() && !p->is_removed_pave_outer()){
            vector<Pave*> zone_list;
            zone_list.push_back(p);
            get_recursive_zone(p, zone_list);

            // Isolate the boundaries of the list
            // Find the first none full node of the list
            reset_marker(zone_list);
            for(Pave *p_zone:zone_list){
                if(!p_zone->is_marked() && !p_zone->is_full() && !p_zone->is_empty()){
                    vector<Pave *> contour_list;
                    contour_list.push_back(p_zone);
                    get_recursive_contour(p_zone, contour_list);
                    contours_list.push_back(contour_list);
                }
            }
            set_marker(zone_list, true);
        }
    }
    return contours_list;
}

vector<double> Graph::get_perimeters(){
    vector<double> perimeters;
    vector<vector<Pave*>> contours_list = get_contour_nodes();
    for(vector<Pave*> contour:contours_list){
        double perimeter = 0.0;
        for(Pave *p:contour){
            perimeter += p->get_perimeter();
        }
        perimeters.push_back(perimeter);
    }
    return perimeters;
}

bool Graph::is_positive_invariant(){
    bool one_pave_not_empty = false;
    for(Pave *p:m_node_list){
        if(!p->is_positive_invariant())
            return false;
        if(!p->is_empty())
            one_pave_not_empty = true;
    }
    return one_pave_not_empty && true;
}

void Graph::push_back_pos_attractor(){
    vector<vector< vector<IntervalVector>>> segment_list;
    for(Pave *p:m_node_list){
        if(p->get_segment_list().size()!=0)
            segment_list.push_back(p->get_segment_list());
    }
    m_pos_attractor_list.push_back(segment_list);
}

bool Graph::get_positive_invariant() const{
    return m_positive_invariant;
}

void Graph::set_positive_invariant(bool val){
    m_positive_invariant = val;
}

void Graph::reset_pave_segment_list(){
#pragma omp parallel for
    for(vector<Pave *>::iterator p = m_node_list.begin(); p<m_node_list.end(); ++p){
        (*p)->reset_segment_list();
    }
}

void Graph::reset_full_empty(){
    for(Pave *p:m_node_list){
        p->reset_full_empty();
    }
}

std::vector<std::vector<std::vector< std::vector<ibex::IntervalVector>>>> Graph::get_pos_attractor_list() const{
    return m_pos_attractor_list;
}

bool Graph::is_sufficiently_discretized(){
    double diam_x = m_search_box[0].diam();
    double diam_y = m_search_box[1].diam();

    for(Pave *p:m_node_list){
        if(p->get_position()[0].diam() == diam_x || p->get_position()[1].diam() == diam_y)
            return false;
    }
    return true;
}

void Graph::reset_queues(){
    for(Pave *p:m_node_queue_inner){
        p->set_in_queue_inner(false);
    }
    for(Pave *p:m_node_queue_outer){
        p->set_in_queue_outer(false);
    }
    m_node_queue_inner.clear();
    m_node_queue_outer.clear();
}

std::list<ibex::Function*>  Graph::get_inside_curve_list() const{
    return m_inside_curve_list;
}

void Graph::push_back_inside_curve(ibex::Function* curve){
    m_inside_curve_list.push_back(curve);
}
