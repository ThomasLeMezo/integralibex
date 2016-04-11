#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"
#include "omp.h"

#include "graphdot.h"
#include "imageintegral.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, const IntervalVector &u, bool diseable_singleton){
    m_graph_list.push_back(new Graph(box, f_list, &m_utils, u, 0, diseable_singleton));
}

Scheduler::~Scheduler(){
    for(auto &graph:m_graph_list){
        delete(graph);
    }
}

Scheduler::Scheduler(const IntervalVector &box, const vector<IntervalVector> &remove_boxes, const std::vector<ibex::Function *> &f_list, const IntervalVector &u, bool diseable_singleton){
    Graph *g = new Graph(&m_utils, 0);
    m_graph_list.push_back(g);

    // Creates paves
    vector<IntervalVector> list_boxes;
    list_boxes.push_back(box);

    for(auto &box_remove:remove_boxes){
        vector<IntervalVector> list_boxes_tmp;
        for(auto &b:list_boxes){
            IntervalVector* box_result;
            int nb_boxes = b.diff(box_remove, box_result);

            for(int i=0; i<nb_boxes; i++)
                list_boxes_tmp.push_back(box_result[i]);
        }
        list_boxes.swap(list_boxes_tmp);
        list_boxes_tmp.clear();
    }

    for(auto &b:list_boxes){
        Pave* p = new Pave(b, f_list, u, diseable_singleton);
        g->get_node_list().push_back(p);
    }
    // Diseable continuity on the bounding box
    define_continuity(g, box, false, false);

//    box[0] = Interval(-1.0, 13.0);
//    box[1] = Interval(-16, 16);
//    vector<IntervalVector> list_border;
//    double size_border = 0.1;
//    IntervalVector test(2);
//    test[0] = Interval(box[0].lb(), box[0].ub());
//    test[1] = Interval(box[1].lb()-size_border, box[1].lb());
//    list_border.push_back(test);
//    test[0] = Interval(box[0].ub(), box[0].ub()+size_border);
//    test[1] = Interval(box[1].lb(), box[1].ub());
//    list_border.push_back(test);
//    test[0] = Interval(box[0].lb(), box[0].ub());
//    test[1] = Interval(box[1].ub(), box[1].ub()+size_border);
//    list_border.push_back(test);
//    test[0] = Interval(box[0].lb()-size_border, box[0].lb());
//    test[1] = Interval(box[1].lb(), box[1].ub());
//    list_border.push_back(test);

//    for(auto &b:list_border){
//        Pave* p = new Pave(b, f_list, u, diseable_singleton, false);
////        p->set_full_out();
//        p->set_full_in();
////        p->set_continuity_out(false);
//        p->set_continuity_in(false);
//        g->get_node_list().push_back(p);
//    }

    for(auto &b:remove_boxes){
        Pave* p = new Pave(b, f_list, u, diseable_singleton, false);
        p->set_full_out();
        p->set_continuity_out(false);
        g->get_node_list().push_back(p);
    }

    // Rebuilt continuity inside graph
    g->build_graph();
}

void Scheduler::define_continuity(Graph *g, const IntervalVector &box, bool continuity_in, bool continuity_out){
    IntervalVector box0(2), box1(2), box2(2), box3(2);
    box0[0] = box[0];                   box0[1] = Interval(box[1].lb());
    box1[0] = Interval(box[0].ub());    box1[1] = box[1];
    box2[0] = box[0];                   box2[1] = Interval(box[1].ub());
    box3[0] = Interval(box[0].lb());    box3[1] = box[1];

    for(int i=0; i<g->get_node_list().size(); i++){
        for(auto &b:g->get_node_list()[i]->get_borders()){
            IntervalVector inter0 = b->get_position() & box0;
            IntervalVector inter1 = b->get_position() & box1;
            IntervalVector inter2 = b->get_position() & box2;
            IntervalVector inter3 = b->get_position() & box3;

            bool test0 = !inter0.is_empty() && (!inter0[0].is_degenerated() || !inter0[1].is_degenerated());
            bool test1 = !inter1.is_empty() && (!inter1[0].is_degenerated() || !inter1[1].is_degenerated());
            bool test2 = !inter2.is_empty() && (!inter2[0].is_degenerated() || !inter2[1].is_degenerated());
            bool test3 = !inter3.is_empty() && (!inter3[0].is_degenerated() || !inter3[1].is_degenerated());

            if(test0 || test1 || test2 || test3){
                b->set_continuity_in(continuity_in);
                b->set_continuity_out(continuity_out);
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
