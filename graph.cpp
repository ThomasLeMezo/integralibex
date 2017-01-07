#include "graph.h"
#include "vibes.h"

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkAppendPolyData.h>

using namespace std;
using namespace ibex;

using namespace Parma_Polyhedra_Library::IO_Operators;

Graph::Graph(const IntervalVector &box, ibex::Function *f, const ibex::IntervalVector &u, int graph_id){
    Pave *p = new Pave(box, f, u);
    m_node_list.push_back(p);
    m_graph_id = graph_id;
    m_drawing_cpt = 0;
}

Graph::Graph(Graph* g, int graph_id){
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

        for(int face = 0; face<pave_root->get_size(); face++){
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
}

Graph::Graph(Graph* g, Pave* activated_node, int graph_id) : Graph(g, graph_id){

    for(auto &node:m_node_list){
        node->set_empty();
    }

    Pave* copy_node = activated_node->get_copy_node();
    copy_node->set_full();
    *(copy_node) &= *(activated_node);
    for(int face=0; face<copy_node->get_size(); face++){
        vector<Pave*> brothers_pave = copy_node->get_brothers(face);
        for(int i=0; i<brothers_pave.size(); i++){
            if(!brothers_pave[i]->is_in_queue()){
                m_node_queue.push_back(brothers_pave[i]);
                brothers_pave[i]->set_in_queue(true);
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

void Graph::clear_node_queue(){
    for(auto &node:m_node_list){
        node->set_in_queue(false);
    }
    m_node_queue.clear();
}

void Graph::sivia(int nb_node, bool backward, bool do_not_bisect_empty, bool do_not_bisect_full){
    int iterations = 0;
    vector<Pave *> tmp_pave_list(m_node_list);
    m_node_list.clear();
    m_node_list.reserve(nb_node);

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<nb_node){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        if(do_not_bisect_empty)
            tmp->reset_full_empty();
        if((do_not_bisect_empty && tmp->is_empty()) || (do_not_bisect_full && tmp->is_full())){// || (not_full_test && tmp->is_full() && diam < M_PI)){
            m_node_list.push_back(tmp);
            //            iterations++;
        }
        else{
            tmp->bisect(tmp_pave_list);
            delete(tmp);
        }
    }

    for(int i=0; i<tmp_pave_list.size(); i++){
        if(backward){
            m_node_queue.push_back(tmp_pave_list[i]);
            tmp_pave_list[i]->set_in_queue(true);
        }
        m_node_list.push_back(tmp_pave_list[i]);
    }
}

int Graph::process(int max_iterations, bool backward, bool inner){
    int iterations = 0;

    while(!m_node_queue.empty() & iterations < max_iterations){
        iterations++;
        Pave *pave = m_node_queue.front();
        m_node_queue.erase(m_node_queue.begin());
        pave->set_in_queue(false);

        bool change = CtcContinuity(pave, backward);
        if(change || pave->get_first_process()){
            CtcPaveConsistency(pave, backward, inner);

            // Warn scheduler to process new pave
            for(int face=0; face<pave->get_size(); face++){
                vector<Pave*> brothers_pave = pave->get_brothers(face);
                for(int i=0; i<brothers_pave.size(); i++){
                    if(brothers_pave[i]->is_in_queue() == false){
                        m_node_queue.push_back(brothers_pave[i]);
                        brothers_pave[i]->set_in_queue(true);
                    }
                }
            }
            pave->set_first_process_false();
        }

//        if(iterations%100==0){
//            cout << '\r' << "ITERATIONS = " << iterations << " / " << max_iterations
//                 << " " << "node queue size = " << m_node_queue.size() << "/" << m_node_list.size() << flush;
//        }
    }

    cout << '\r' << "ITERATIONS = " << iterations << " / " << max_iterations << endl;
    m_node_queue.clear();
    return iterations;
}

void Graph::set_full(){
    for(auto &node : m_node_list){
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
    position[0] = ibex::Interval(x);
    position[1] = ibex::Interval(y);

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
        for(int face=0; face<pave->get_size(); face++){
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

const std::vector<Pave *> &Graph::get_node_list() const {
    return m_node_list;
}

const std::vector<Pave *> Graph::get_node_queue() const {
    return m_node_queue;
}

Pave* Graph::get_node_const(int i) const{
    return m_node_list[i];
}

Pave& Graph::operator[](int id){
    return *(m_node_list[id]);
}

void Graph::draw_vtk(string filename){
    cout << "************ DRAWING ************" << endl;
    bool polygon = false;
    bool box_empty = false;

    if(m_node_list.size()>0)
        polygon = true;
    if(m_node_empty_list.size()>0)
        box_empty = true;

    vtkSmartPointer<vtkAppendPolyData> polyData_polygon = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> polyData_box_active = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> polyData_box_empty = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> polyData_vector_field = vtkSmartPointer<vtkAppendPolyData>::New();

    if(polygon){
        cout << "m_node_list.size() = " << m_node_list.size() << endl;
        for(auto &node:m_node_list){
            node->reset_full_empty();
            if(!node->is_empty()){
                node->draw_vtk(polyData_polygon);
            }
        }
        polyData_polygon->Update();

        for(auto &node:m_node_list){
            node->draw_box(polyData_box_active);
        }
        polyData_box_active->Update();

        for(auto &node:m_node_list){
            node->draw_vector_field(polyData_vector_field);
        }
        polyData_vector_field->Update();
    }

    if(box_empty){
        cout << "m_node_empty_list.size() = " << m_node_empty_list.size() << endl;
        for(auto &node:m_node_empty_list){
            node->draw_box(polyData_box_empty);
        }
        polyData_box_empty->Update();
    }


    vtkSmartPointer<vtkXMLPolyDataWriter> outputWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    string file = filename + to_string(m_graph_id) + "-" + to_string(m_drawing_cpt);

    string filePolygon = file + "polygon.vtp";
    string fileBox_active = file + "box_active.vtp";
    string fileBox_empty = file + "box_empty.vtp";
    string fileVectorField = file + "vector_field.vtp";

    if(polygon){
        outputWriter->SetFileName(filePolygon.c_str());
        outputWriter->SetInputData(polyData_polygon->GetOutput());
        outputWriter->Write();

        outputWriter->SetFileName(fileBox_active.c_str());
        outputWriter->SetInputData(polyData_box_active->GetOutput());
        outputWriter->Write();

        outputWriter->SetFileName(fileVectorField.c_str());
        outputWriter->SetInputData(polyData_vector_field->GetOutput());
        outputWriter->Write();
    }

    if(box_empty){
        outputWriter->SetFileName(fileBox_empty.c_str());
        outputWriter->SetInputData(polyData_box_empty->GetOutput());
        outputWriter->Write();
    }
}

int Graph::size() const{
    return m_node_list.size();
}

void Graph::remove_empty_node(){
    for(int i=0; i<m_node_list.size(); i++){
        m_node_list[i]->reset_full_empty();
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

void Graph::set_symetry(Function* f, int axis_in, int side_in, int axis_out, int side_out){
    Inclusion *i = new Inclusion(m_node_list[0]->get_border(2*axis_in + side_in), f, axis_in, side_in);
    m_node_list[0]->get_border(2*axis_out+side_out)->add_inclusion(i);
}

void Graph::set_all_first_process(){
    for(auto& node:m_node_list){
        node->set_first_process_true();
    }
}

int Graph::get_graph_id(){
    return m_graph_id;
}


void Graph::print_pave_info(double x, double y) const{
    Pave *p = get_pave(x, y);
    if(p!=NULL){
        cout << "*****" << endl;
        cout << "pave = " << p->get_position() << endl;
        cout << "volume_in = " << p->get_volume_in().generators() << endl;
        cout << "volume_out = " << p->get_volume_out().generators() << endl;

        cout << "theta = " << p->get_theta() << endl;
    }
    else{
        cout << "pave not found" << endl;
    }
}

void Graph::reset_full_empty(){
    for(auto &p:m_node_list){
        p->reset_full_empty();
    }
}
