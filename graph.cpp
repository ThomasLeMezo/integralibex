#include "graph.h"
#include "vibes.h"

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
    m_drawing_cpt = 0;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
    m_compute_inner = false;
    m_inner_mode = false;
}

Graph::Graph(Utils *utils, int graph_id=0):
    m_search_box(2)
{
    m_graph_id = graph_id;
    m_drawing_cpt = 0;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
    m_compute_inner = false;
    m_inner_mode = false;
}

Graph::Graph(Graph* g, int graph_id):
    m_search_box(2)
{
    for(Pave *node:g->get_node_list()){
        Pave *p = new Pave(node);
        node->set_copy_node(p);
        m_node_list.push_back(p);
    }
    for(Pave *node:g->get_node_queue()){
        add_to_all_queue(node->get_copy_node());
    }

    for(int i=0; i<m_node_list.size(); i++){
        Pave* pave_root = g->get_node_list()[i];
        Pave* pave_copy = m_node_list[i];

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
    m_drawing_cpt = 0;
    m_utils = g->get_utils();
    m_search_box = g->get_search_box();
    m_count_alive = g->get_alive_node();

    debug_marker1 = false;
    debug_marker2 = false;

    m_compute_inner = g->get_compute_inner();
    m_inner_mode = g->get_inner_mode();
}

Graph::Graph(Graph* g, Pave* activated_node, int graph_id) : Graph(g, graph_id){
    cout << "COPY GRAPH size = " << size() << endl;
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
    for(Pave *node:m_node_list){
        node->set_in_queue_outer(false);
    }
    m_node_queue_outer.clear();
}

void Graph::sivia(int nb_node, bool backward, bool do_not_bisect_empty, bool do_not_bisect_full){
    //    if(nb_node<=m_count_alive)
    //        return;
    vector<Pave *> tmp_pave_list(m_node_list);
    m_node_list.clear();
    m_node_list.reserve(nb_node);

    while(((int)tmp_pave_list.size()!=0) & (m_count_alive<nb_node)){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        if(do_not_bisect_empty || do_not_bisect_full)
            tmp->reset_full_empty();
        //        if(m_utils->m_imageIntegral_activated)
        //            tmp->set_inner(m_utils->m_imageIntegral->testBox(tmp->get_position()));

        if(!tmp->is_active() || tmp->is_removed_pave_union()
                || ((do_not_bisect_empty && tmp->is_empty_inter()) || (do_not_bisect_full && tmp->is_full_inter()))
                || tmp->get_theta_diam_max()<0.0){
            m_node_list.push_back(tmp);
            if(!tmp->get_zone_propagation())
                compute_propagation_zone(tmp);
        }
        else{
            tmp->bisect(tmp_pave_list, backward);
            delete(tmp);
            m_count_alive++;
        }
    }

    for(Pave *p:tmp_pave_list){
        if(backward && !p->is_removed_pave()){
            add_to_queue_outer(p);
            if(m_compute_inner){
                add_to_queue_inner(p);
            }
        }
        if(!p->is_removed_pave() && m_compute_inner)
            compute_propagation_zone(p);
        m_node_list.push_back(p);
    }
    cout << "SIVIA outer(" << m_node_queue_outer.size() << ") inner(" <<  m_node_queue_inner.size() << ")" << endl;
}

int Graph::process(int max_iterations, bool backward, bool union_functions){
    int iterations = 0;
    while(!is_empty_node_queue() & iterations < max_iterations){
        iterations++;
        Pave *pave = get_node_queue_access().front();
        pop_front_queue();
        pave->set_in_queue(false);

        /// ******* PROCESS CONTINUITY *******
        bool change = m_utils->CtcContinuity(pave, backward);
        if(pave->is_active() && !pave->is_removed_pave() && (change || pave->get_first_process())){

            /// ******* PROCESS CONSISTENCY *******
            std::vector<bool> change_tab;
            for(int i=0; i<4; i++)
                change_tab.push_back(false);
            m_utils->CtcConsistency(pave, backward, change_tab, union_functions);

            /// ******* PUSH BACK NEW PAVES *******
            // Warn scheduler to process new pave
            if(backward && !pave->get_zone_propagation())
                compute_propagation_zone(pave);
            for(int face=0; face<4; face++){
                if(change_tab[face]){
                    for(Pave *p:pave->get_brothers(face)){
                        if(p->is_in_queue() == false){
                            add_to_queue(p);
                        }
                    }
                }
            }

            pave->set_first_process_false();
        }
    }

    clear_node_queue();
    if(get_inner_mode())
        cout << "--> processing inner (" << iterations << ")" << endl;
    else
        cout << "--> processing outer (" << iterations << ")" << endl;
    return iterations;
}

void Graph::set_full(){
    for(Pave *node : m_node_list){
        if(node->is_active() && !node->is_removed_pave())
            node->set_full();
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

void Graph::set_active_outer_inner(const ibex::IntervalVector &box){
    std::vector<ibex::IntervalVector> box_list;
    box_list.push_back(box);
    set_active_outer_inner(box_list);
}

void Graph::set_active_outer_inner(const std::vector<ibex::IntervalVector> &box_list){
    for(Pave *pave:m_node_list){
        for(IntervalVector box:box_list){
            if(!(box & pave->get_position()).is_empty()){
                // Outer
                pave->set_full_outer();

                // Inner
                bool inner=false;
                if(pave->get_position().is_strict_interior_subset(box)){
                    pave->set_empty_inner_in();
                    pave->set_removed_pave_inner(true);
                    pave->set_active(false);
                    inner = true;
                    m_count_alive--;
                }

                compute_propagation_zone(pave, true);

                // Update queue
                for(int face=0; face<4; face++){
                    vector<Pave*> pave_brother_list = pave->get_brothers(face);
                    for(Pave *pave_brother:pave_brother_list){
                        if(!pave_brother->is_in_queue_outer()){
                            add_to_queue_outer(pave_brother);
                        }
                        if(inner && !pave_brother->is_in_queue_inner()
                                && !pave->get_border(face)->get_zone_function_in(pave->get_active_function())){
                            add_to_queue_inner(pave_brother);
                        }
                    }
                }
            }
        }
    }
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

std::list<Pave *> &Graph::get_node_queue_access(){
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

void Graph::draw(int size, bool filled, string comment, bool inner_only){

    // Magenta = #FF00FF
    // Gray light =  #D3D3D3
    // Blue = #4C4CFF

    stringstream ss;
    ss << "integralIbex" << m_graph_id<< "-" << m_drawing_cpt << " " << comment;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    for(Pave *node:m_node_list){
        node->draw(filled, inner_only);
    }

    for(Pave *node:m_node_border_list){
        node->draw(filled);
    }

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    //    m_drawing_cpt++;
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
    p->print_theta_list();
    cout << "inner mode = " << p->get_inner_mode() << endl;
    p->reset_full_empty();
    cout << "is_full_inner = " << p->is_full_inner() << " is_empty_inner = " << p->is_empty_inner() << endl;
    cout << "is_full_outer = " << p->is_full_outer() << " is_empty_outer = " << p->is_empty_outer() << endl;
    cout << "nb\t" << "in_inner\t" << "out_inner\t" << "in_outer\t" << "out_outer\t" << endl;
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
    for(Pave *pave:m_node_list){
        if(!pave->is_removed_pave() && pave->is_active()){
            pave->reset_full_empty();
            bool empty = false;
            if(pave->is_empty_outer()){
                pave->set_removed_pave_outer(true);
                empty = true;
            }
            if(m_compute_inner && pave->is_empty_inner()){
                pave->set_removed_pave_inner(true);
                empty = true;
            }
            if(empty){
                pave->set_active(false);
                m_count_alive--;
            }
        }
    }
}

bool Graph::is_empty(){
    if(get_alive_node()==0)
        return true;
    bool empty = true;
    for(Pave *node:m_node_list){
        node->reset_full_empty();
        if(!node->is_empty_inter())
            empty = false;
    }
    return empty;
}

Pave* Graph::get_semi_full_node(){
    for(Pave *node:m_node_list){
        if(!node->is_removed_pave() && node->is_border() && node->get_theta_diam()<M_PI/2.0){
            return node;
        }
    }

    for(Pave *node:m_node_list){
        if(!node->is_removed_pave() && !node->is_empty() && !node->is_full()){
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
        for(int i=0; i<g.size(); i++){
            m_node_list[i]->diff(*(g.get_node_const(i)));
        }
    }
}

void Graph::inter(const Graph &g, bool with_bwd){
    if(this->size() == g.size()){
        for(int i=0; i<g.size(); i++){
            m_node_list[i]->inter(*(g.get_node_const(i)), with_bwd);
        }
    }
}

void Graph::set_empty(){
    for(Pave *pave : m_node_list){
        if(pave->is_active())
            pave->set_empty();
    }
}

void Graph::set_empty_outer_full_inner(){
    for(Pave *pave : m_node_list){
        if(pave->is_active()){
            pave->set_empty_outer();
            pave->set_full_inner();
        }
    }
}

void Graph::set_symetry(Function* f, int face_in, int face_out){
    Inclusion *i = new Inclusion(m_node_list[0]->get_border(face_in), f, face_in);
    m_node_list[0]->get_border(face_out)->add_inclusion(i);
}

void Graph::set_all_first_process(){
    for(Pave *pave:m_node_list){
        pave->set_first_process_true();
    }
}

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
                        b2->add_inclusion(inclusion_to_b1);

                        Inclusion *inclusion_to_b2 = new Inclusion(b2, b2->get_face());
                        b1->add_inclusion(inclusion_to_b2);
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
    for(Pave *p:m_node_list){
        p->set_first_process_true();
        p->set_in_queue(false);
        p->reset_full_empty();

        if(p->is_active() && (
                    (!p->is_full() && !p->is_empty())
                    || (border_condition && p->is_border())
                    || (empty_condition && p->is_near_empty()))){
            add_to_queue(p);
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
    for(Pave *p:m_node_list){
        p->set_active_function(id);
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
                if(p->is_border() && p->get_theta_diam()>=M_PI){
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

IntervalVector Graph::get_search_box() const{
    return m_search_box;
}

void Graph::complementaire(){
    for(Pave *p:m_node_list){
        p->complementaire();
    }
}

void Graph::set_external_boundary(bool in, bool out){
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

int Graph::get_alive_node(){
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
    for(Pave *p:m_node_list){
        p->set_compute_inner(val);
    }
}

void Graph::set_inner_mode(bool val){
    if(!m_compute_inner)
        set_compute_inner(true);
    m_inner_mode = val;
    for(Pave *p:m_node_list){
        p->set_inner_mode(val);
    }
    for(Pave *p:m_node_border_list){
        p->set_inner_mode(val);
    }
}

void Graph::add_to_queue_inner(Pave *p){
    m_node_queue_inner.push_back(p);
    p->set_in_queue_inner(true);
}

void Graph::add_to_queue_outer(Pave *p){
    m_node_queue_outer.push_back(p);
    p->set_in_queue_outer(true);
}

bool Graph::is_empty_node_queue(){
    if(!m_compute_inner)
        return m_node_queue_outer.size() == 0;
    else{
        if(m_inner_mode)
            return m_node_queue_inner.size() == 0;
        else
            return m_node_queue_outer.size() == 0;
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
    if(!m_compute_inner)
        m_node_queue_outer.pop_front();
    else{
        if(m_inner_mode)
            m_node_queue_inner.pop_front();
        else
            m_node_queue_outer.pop_front();
    }
}

void Graph::copy_to_inner(){
    for(Pave *p:m_node_list)
        p->copy_to_inner();
}

void Graph::compute_propagation_zone(Pave *p, bool compute_anyway){
    if(!compute_anyway && p->get_f_list().size()==1 && !p->get_zone_propagation())
        return;
    for(int fwd=0; fwd<2; fwd++){
        Pave *p_copy = new Pave(p);
        p_copy->set_inner_mode(false);
        p_copy->set_compute_inner(false);

        p_copy->set_backward_function(fwd);
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
