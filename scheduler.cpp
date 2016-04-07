#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"
#include "omp.h"

#include "graphdot.h"
#include "imageintegral.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(const IntervalVector &box, ibex::Function *f, const IntervalVector &u){
    m_graph_list.push_back(new Graph(box, f, &m_utils, u, 0));
}

Scheduler::~Scheduler(){
    for(auto &graph:m_graph_list){
        delete(graph);
    }
}

Scheduler::Scheduler(const IntervalVector &box, const IntervalVector &box_remove, ibex::Function *f, const IntervalVector &u){
    Graph *g = new Graph(&m_utils, 0);
    m_graph_list.push_back(g);

    // Creates paves
    IntervalVector* box_diff;
    int nb_box = box.diff(box_remove, box_diff);

    for(int i=0; i<nb_box; i++){
        Pave* p = new Pave(box_diff[i], f, u);
        g->get_node_list().push_back(p);
    }

    // Rebuilt continuity inside graph
    g->build_graph();

    // Diseable continuity on bounding box
    for(int i=0; i<nb_box; i++){
        for(auto &b:g->get_node_list()[i]->get_borders()){
            IntervalVector test_inter = b->get_position() & box;
            if(!(test_inter[0].is_degenerated() && test_inter[1].is_degenerated())){
                b->set_continuity(false);
            }
        }
    }
}

void Scheduler::set_symetry(Function *f, int face_in, int face_out){
    m_graph_list[0]->set_symetry(f, face_in, face_out);
}

void Scheduler::cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe, bool inner){
    vector<IntervalVector> initial_boxes;
    initial_boxes.push_back(initial_boxe);
    cameleon_propagation(iterations_max, process_iterations_max, initial_boxes, inner);
}

void Scheduler::cameleon_propagation(int iterations_max, int process_iterations_max, vector<IntervalVector> &initial_boxes, bool inner){
    if(m_graph_list.size()!=1 && m_graph_list[0]->size() !=1)
        return;
    int iterations = 0;

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        m_graph_list[0]->sivia(4,false, false, false); // Start with 4 boxes
        m_graph_list[0]->set_empty();
        for(auto &initial_box:initial_boxes)
            m_graph_list[0]->set_active_pave(initial_box);
        m_graph_list[0]->process(process_iterations_max, false, false);
        m_graph_list[0]->remove_empty_node();
        iterations++;
    }
    int nb_graph = 0;

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        m_graph_list[0]->remove_empty_node();
        m_graph_list[0]->sivia(2*m_graph_list[0]->size(), false, true, true);

        m_graph_list[0]->set_empty();
        for(auto &initial_box:initial_boxes)
            m_graph_list[0]->set_active_pave(initial_box);

        // Process the forward with the subpaving
        cout << "GRAPH No "<< nb_graph << " (" << m_graph_list[0]->size() << ")" << endl;
        m_graph_list[0]->process(process_iterations_max, false, false);

        // Remove empty pave & Test if the graph is empty

        if(inner && iterations == iterations_max -1){
            if(m_graph_inner_list.size()==1){
                delete(m_graph_inner_list[0]);
                m_graph_inner_list.clear();
            }
            Graph *graph_inner = new Graph(m_graph_list[0]);
            graph_inner->set_empty();
            for(auto &initial_box:initial_boxes)
                graph_inner->set_active_pave(initial_box);
            graph_inner->process(process_iterations_max, false, true);
            m_graph_inner_list.push_back(graph_inner);
        }

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool inner, bool do_not_bisect_inside){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;

    int iterations = 0;
    m_graph_list[0]->set_full();

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        m_graph_list[0]->sivia(4,true, false, false); // Start with 4 boxes
        m_graph_list[0]->process(process_iterations_max, true, false); // ? Usefull ??? ToDo
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        //        if(m_utils.m_imageIntegral_activated){
        //            m_utils.m_imageIntegral->set_box(m_graph_list[0]->get_bounding_box());
        //        }

        for(int nb_graph=0; nb_graph<m_graph_list.size(); nb_graph++){

            if(m_graph_list[nb_graph]->size()==0 || m_graph_list.size()==0)
                break;
            m_graph_list[nb_graph]->clear_node_queue();
            m_graph_list[nb_graph]->sivia(2*m_graph_list[nb_graph]->size(), true, false, do_not_bisect_inside);

            // Process the backward with the subpaving
            cout << "GRAPH No "<< nb_graph << " (" << m_graph_list[nb_graph]->size() << ")" << endl;

            int graph_list_process_cpt = m_graph_list[nb_graph]->process(process_iterations_max, true, false);
            cout << "--> processing outer = " << graph_list_process_cpt << endl;

            if(inner && iterations>0){
                cout << "COMPUTE INNER" << endl;
                if(m_graph_inner_list.size()==m_graph_list.size()){
                    delete(m_graph_inner_list[nb_graph]);
                    m_graph_inner_list.erase(m_graph_inner_list.begin() + nb_graph);
                }
                Graph* graph_inner;
                graph_inner = new Graph(m_graph_list[nb_graph]);
                graph_inner->add_all_to_queue();
                graph_inner->set_all_first_process();
                int graph_inner_process_cpt = graph_inner->process(process_iterations_max, true, true);
                cout << "--> processing inner = " << graph_inner_process_cpt << endl;
                m_graph_inner_list.insert(m_graph_inner_list.begin()+nb_graph, graph_inner);

                //                if(!graph_inner->is_empty()){
                //                    cout << "GRAPH inner not empty, iteration = " << iterations << endl;
                //                    break;
                //                }
            }

            // Remove empty pave
            m_graph_list[nb_graph]->remove_empty_node();

            // Test if the graph is empty
            if(m_graph_list[nb_graph]->is_empty()){
                if(m_graph_inner_list.size()==m_graph_list.size()){
                    delete(m_graph_inner_list[nb_graph]);
                    m_graph_inner_list.erase(m_graph_inner_list.begin() + nb_graph);
                }

                delete(m_graph_list[nb_graph]);
                m_graph_list.erase(m_graph_list.begin()+nb_graph);
                cout << " REMOVE EMPTY GRAPH" << endl;
                if(nb_graph!=0)
                    nb_graph--;
                break;
            }


            // ***************************************************
            //              REMOVE INSIDE PROCEDURE
            // Copy graph & propagate one Pave + intersect with cycle
            // Find the first non-full & non-empty Pave
            if(remove_inside && m_graph_list.size() < graph_max){
                Pave *pave_start = m_graph_list[nb_graph]->get_semi_full_node(); // Find a pave semi full
                if(pave_start == NULL)
                    break;

                Graph* graph_propagation = new Graph(m_graph_list[nb_graph], pave_start); // copy graph with 1 activated node (pave_start)
                //                Graph* graph_diff = new Graph(m_graph_list[nb_graph], m_graph_list.size());

                graph_propagation->process(process_iterations_max, false, false); // process forward

                m_graph_list[nb_graph]->inter(*graph_propagation); // intersect the graph with the propagation graph


                //                graph_diff->diff(*m_graph_list[nb_graph]);

                //                if(!graph_diff->is_empty() && false){ // If there is an inside, add to graph_list
                //                    cout << " REMOVE INSIDE" << endl;
                //                    m_graph_list.push_back(graph_diff);
                //                    graph_diff->clear_node_queue();
                //                    cout << " ADD NEW GRAPH No " << m_graph_list.size()-1 << endl;
                //                    //m_graph_list.back()->print();
                //                }
                //                else{
                //                    delete(graph_diff);
                //                }
                m_graph_list[nb_graph]->remove_empty_node();

                delete(graph_propagation);

            }
        }
        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Scheduler::draw(int size, bool filled){
    for(int i=0; i<m_graph_list.size(); i++){
        m_graph_list[i]->draw(size, filled);

        if(m_graph_inner_list.size()==m_graph_list.size())
            m_graph_inner_list[i]->drawInner(filled);
    }
}

Graph* Scheduler::get_graph_list(int i){
    return m_graph_list[i];
}

void Scheduler::print_pave_info(int graph, double x, double y, string color){
    m_graph_list[graph]->print_pave_info(x, y, color);
}

void Scheduler::set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY){
    m_utils.m_imageIntegral = new imageIntegral(range, f, t_range, nbBisectionT, nbBisectionXY);
    m_utils.m_imageIntegral_activated = true;
}
