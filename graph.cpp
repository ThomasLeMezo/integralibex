#include "graph.h"
#include "vibes.h"

using namespace std;
using namespace ibex;

Graph::Graph(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, Utils *utils, const IntervalVector &u, int graph_id, bool diseable_singleton):
    m_search_box(2)
{
    Pave *p = new Pave(box, f_list, u, diseable_singleton);
    m_count_alive = 1;
    m_search_box = box;
    m_node_list.push_back(p);
    m_graph_id = graph_id;
    m_drawing_cpt = 0;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
}

Graph::Graph(Utils *utils, int graph_id=0):
    m_search_box(2)
{
    m_graph_id = graph_id;
    m_drawing_cpt = 0;
    m_utils = utils;

    debug_marker1 = false;
    debug_marker2 = false;
}

Graph::Graph(Graph* g, int graph_id):
    m_search_box(2)
{
    for(auto &node:g->get_node_list()){
        Pave *p = new Pave(node);
        node->set_copy_node(p);
        m_node_list.push_back(p);
    }
    for(auto &node:g->get_node_queue()){
        m_node_queue.push_back(node->get_copy_node());
        node->get_copy_node()->set_in_queue(true);
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
}

Graph::Graph(Graph* g, Pave* activated_node, int graph_id) : Graph(g, graph_id){
    cout << "COPY GRAPH size = " << size() << endl;
    for(auto &node:m_node_list){
        node->set_empty();
    }

    Pave* copy_node = activated_node->get_copy_node();
    copy_node->set_full();
    *(copy_node) &= *(activated_node);

    vector<Pave*> brothers_pave = copy_node->get_all_brothers();
    for(auto &p:brothers_pave){
        if(!p->is_in_queue()){
            m_node_queue.push_back(p);
            p->set_in_queue(true);
        }
    }
}

Graph::~Graph(){
    for(auto &node:m_node_list){
        delete(node);
    }
    for(auto &node:m_node_border_list){
        delete(node);
    }
}

void Graph::clear_node_queue(){
    for(auto &node:m_node_list){
        node->set_in_queue(false);
    }
    m_node_queue.clear();
}

void Graph::sivia(int nb_node, bool backward, bool do_not_bisect_empty, bool do_not_bisect_full, double theta_limit){
    int iterations = 0;
    m_count_alive = 0;
    vector<Pave *> tmp_pave_list(m_node_list);
    m_node_list.clear();
    m_node_list.reserve(nb_node);

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<nb_node){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        if(do_not_bisect_empty || do_not_bisect_full)
            tmp->reset_full_empty();
        if(m_utils->m_imageIntegral_activated)
            tmp->set_inner(m_utils->m_imageIntegral->testBox(tmp->get_position()));

        if(!tmp->is_active() || tmp->is_removed_pave()
                || ((do_not_bisect_empty && tmp->is_empty()) || (do_not_bisect_full && tmp->is_full()))
                || tmp->get_theta_diam()<theta_limit){
            m_node_list.push_back(tmp);
            if(!tmp->is_removed_pave()){
                iterations++;
                m_count_alive++;
            }
        }
        else{
            tmp->bisect(tmp_pave_list, backward);
            delete(tmp);
        }
    }

    for(int i=0; i<tmp_pave_list.size(); i++){
        if(backward && !tmp_pave_list[i]->is_removed_pave()){
            tmp_pave_list[i]->set_in_queue(true);
            m_node_queue.push_back(tmp_pave_list[i]);
        }
        m_node_list.push_back(tmp_pave_list[i]);
        m_count_alive++;
    }
    cout << "SIVIA node_queue.size() = " << m_node_queue.size() << endl;
}

int Graph::process(int max_iterations, bool backward, bool enable_function_iteration, bool inner){
    int iterations = 0;
    while(!m_node_queue.empty() & iterations < max_iterations){
        iterations++;
        Pave *pave = m_node_queue.front();
        m_node_queue.erase(m_node_queue.begin());
        pave->set_in_queue(false);

        IntervalVector test(2);
        test[0] = Interval(0.4);
        test[1] = Interval(-1.7);
//        if(debug_marker2 && !(test & pave->get_position()).is_empty()){
//            draw(1024, "debug");
//            print_pave_info(test[0].mid(), test[1].mid(), "b[b]");
//            cout << "DEBUG" << endl;
//        }

        /// ******* PROCESS CONTINUITY *******
        bool change = m_utils->CtcContinuity(pave, backward);
        if(pave->is_active() && !pave->is_removed_pave() && (change || pave->get_first_process())){

//            if(debug_marker2 && !(test & pave->get_position()).is_empty()){
//                draw(1024, "debug");
//                print_pave_info(test[0].mid(), test[1].mid(), "b[b]");
//                cout << "DEBUG" << endl;
//            }

            /// ******* PROCESS CONSISTENCY *******
            std::vector<bool> change_tab;
            for(int i=0; i<4; i++)
                change_tab.push_back(false);
            m_utils->CtcPaveConsistency(pave, backward, change_tab, enable_function_iteration, inner);

            /// ******* PUSH BACK NEW PAVES *******
            // Warn scheduler to process new pave
            for(int face=0; face<4; face++){
                if(change_tab[face]){
                    vector<Pave*> brothers_pave = pave->get_brothers(face);
                    for(int i=0; i<brothers_pave.size(); i++){
                        if(brothers_pave[i]->is_in_queue() == false){
                            m_node_queue.push_back(brothers_pave[i]);
                            brothers_pave[i]->set_in_queue(true);
                        }
                    }
                }
            }

            pave->set_first_process_false();
        }

//        if(debug_marker1 && iterations%100==0){
//            draw(4048, true);
//            cout << "SAVE image " << iterations << endl;
//            stringstream ss; ss << "/home/lemezoth/Images/VIBES/iteration_" << iterations << ".png";
//            vibes::saveImage(ss.str());
//        }
//        if(debug_marker2 && !(test & pave->get_position()).is_empty()){
//            draw(1024, "debug");
//            print_pave_info(test[0].mid(), test[1].mid(), "b[b]");
//            cout << "DEBUG" << endl;
//        }
    }

    m_node_queue.clear();
    return iterations;
}

void Graph::set_full(){
    for(auto &node : m_node_list){
        if(node->is_active() && !node->is_removed_pave())
            node->set_full();
    }
}

void Graph::add_all_to_queue(){
    m_node_queue.clear();
    for(auto &node: m_node_list){
        m_node_queue.push_back(node);
    }
}

Pave* Graph::get_pave(double x, double y) const{
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

const std::vector<Pave *> Graph::get_pave(const ibex::IntervalVector &box) const{
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

const std::vector<Pave *> Graph::get_node_queue() const {
    return m_node_queue;
}

std::vector<Pave *> &Graph::get_node_queue_access(){
    return m_node_queue;
}

Pave* Graph::get_node_const(int i) const{
    return m_node_list[i];
}

Pave& Graph::operator[](int id){
    return *(m_node_list[id]);
}

void Graph::draw(int size, bool filled, string comment){

    // Magenta = #FF00FF
    // Gray light =  #D3D3D3
    // Blue = #4C4CFF

    stringstream ss;
    ss << "integralIbex" << m_graph_id<< "-" << m_drawing_cpt << " " << comment;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    for(auto &node:m_node_list){
        //        if(node->is_active()){
        if(node->is_near_inner())
            node->draw(filled, "#D3D3D3[#FF00FF]"); // magenta
        else
            node->draw(filled, "#D3D3D3[#4C4CFF]"); // blue
        //        }
        //        else{
        //            node->draw(filled, "#D3D3D3[blue]");
        //        }
    }

    for(auto &node:m_node_border_list){
        node->draw(filled, "gray[gray]");
    }

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
//    m_drawing_cpt++;
}

void Graph::drawInner(bool filled){
    for(auto &node:m_node_list){
        node->draw(filled, "[]", true);
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
    cout << "nb\t" << "position\t" << "in\t" << "out\t" << "contaminated_in\t" << "contaminated_out\t" << "continuity_in\t" << "continuity_out" << endl;
    for(int i= 0; i<p->get_borders().size(); i++){
        cout << i << '\t'
             //<< p->get_border(i)->get_position() << "       \t"
             << p->get_border(i)->get_segment_in() << "       \t"
             << p->get_border(i)->get_segment_out() << '\t'
             << p->get_border(i)->get_continuity_in() << '\t'
             << p->get_border(i)->get_continuity_out() << '\t'
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
    for(auto &node:m_node_list){
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
    for(auto &pave:m_node_list){
        if(!pave->is_removed_pave() && pave->is_active()){
            pave->reset_full_empty();
            if(pave->is_empty()){
                if(pave->is_near_inner()){
                    pave->set_inner(true);
                }
                pave->set_removed_pave(true);
                pave->set_active(false);
                m_count_alive--;
            }
        }
    }
    propagate_inner();
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
        if(!node->is_removed_pave() && node->is_border() && node->get_theta_diam()<M_PI/2.0){
            return node;
        }
    }

    for(auto &node:m_node_list){
        if(!node->is_removed_pave() && !node->is_empty() && !node->is_full()){
            return node;
        }
    }

    // Case all full or empty
    for(auto &node:m_node_list){
        if(!node->is_removed_pave() && !node->is_empty()){
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
        if(pave->is_active())
            pave->set_empty();
    }
}

void Graph::set_symetry(Function* f, int face_in, int face_out){
    Inclusion *i = new Inclusion(m_node_list[0]->get_border(face_in), f, face_in);
    m_node_list[0]->get_border(face_out)->add_inclusion(i);
}

void Graph::set_all_first_process(){
    for(auto& node:m_node_list){
        node->set_first_process_true();
    }
}

int Graph::get_graph_id(){
    return m_graph_id;
}

IntervalVector Graph::get_bounding_box() const{
    IntervalVector boundingBox(2);
    boundingBox[0] = Interval::EMPTY_SET;
    boundingBox[1] = Interval::EMPTY_SET;

    for(auto &p:m_node_list){
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

    m_node_queue.clear();
    for(auto &p:m_node_list){
        p->set_first_process_true();
        p->set_in_queue(false);
        p->reset_full_empty();

        if((!p->is_full() && !p->is_empty()) || (border_condition && p->is_border()) || (empty_condition && p->is_near_empty())){
            p->set_in_queue(true);
            m_node_queue.push_back(p);
        }
    }
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
    for(auto &p:m_node_list){
        p->set_active_function(id);
    }
}

bool Graph::identify_attractor(){
    if(m_node_list.size()==0)
        return true;
    reset_marker_attractor();
    bool found_all_attractor = true;

    for(auto &p:m_node_list){

        // Study the subgraph starting with p
        if(p->is_active() && !p->is_marked_attractor() && !p->is_removed_pave()){
            vector<Pave*> attractor_list;
            attractor_list.push_back(p);
            get_recursive_attractor(p, attractor_list);

            // Test if the subgraph is an attractor
            bool test_attractor = true;

            IntervalVector bounding_box(2, Interval::EMPTY_SET);
            for(auto &p:attractor_list){
                if(p->is_border() && p->get_theta_diam()>=M_PI){
                    test_attractor = false;
                }
                bounding_box |= p->get_bounding_pave();
            }

            vector<IntervalVector> list_search_box = m_utils->get_segment_from_box(m_search_box);
            for(auto &b:list_search_box){
                if(!((bounding_box & b).is_empty()))
                    test_attractor = false;
            }

            if(test_attractor){
                cout << "FIND AN ATTRACTOR" << endl;
                for(auto &p:attractor_list){
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
        for(auto &p_brother:brothers_face){
            if(!p_brother->is_removed_pave() && !p_brother->is_marked_attractor()){
                list.push_back(p_brother);
                p_brother->set_marker_attractor(true);
                get_recursive_attractor(p_brother, list);
            }
        }
    }
}

void Graph::reset_marker_attractor(){
    for(auto &p:m_node_list){
        p->set_marker_attractor(false);
    }
}

IntervalVector Graph::get_search_box() const{
    return m_search_box;
}

void Graph::complementaire(){
    for(auto &p:m_node_list){
        p->complementaire();
    }
}

void Graph::set_external_boundary(bool in, bool out){
    for(auto &p:m_node_border_list){
        p->set_segment(in, out);
    }
}

void Graph::set_all_active(){
    for(auto &p:m_node_list){
        p->set_active(true);
        p->set_removed_pave(false);
    }
}

void Graph::mark_full_pave_as_inner(){
    for(auto &p:m_node_list){
        if(p->is_full()){
            p->set_inner(true);
        }
    }
}

void Graph::propagate_inner(){
    for(auto &p:m_node_list){
        if(p->is_inner()){
            recursive_mark_inner(p);
        }
    }
}

void Graph::recursive_mark_inner(Pave* p){
    vector<Pave *> brothers = p->get_all_brothers();
    for(auto &bro:brothers){
        if(!bro->is_inner() && bro->is_empty()){
            bro->set_inner(true);
            recursive_mark_inner(bro);
        }
    }
}

int Graph::get_alive_node(){
    return m_count_alive;
}

bool Graph::is_no_active_function(){
    for(auto &p:m_node_list){
        if(p->get_active_function()!=-1){
            return false;
        }
    }
    return true;
}

bool Graph::set_all_inner(bool inner){
    for(auto &p:m_node_list){
        p->set_inner(inner);
    }
}
