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

Scheduler::Scheduler(const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list, const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out){
    Graph *g = new Graph(&m_utils, 0);
    m_graph_list.push_back(g);

    if(bassin_boxes.size()>0){
        vector<IntervalVector> list_boxes;
        list_boxes.push_back(box);

        // Build diff boxes
        for(auto &box_bassin:bassin_boxes){
            vector<IntervalVector> list_boxes_tmp;
            for(auto &b:list_boxes){
                std::vector<IntervalVector> box_result = m_utils.diff(b, box_bassin);

                for(int i=0; i<box_result.size(); i++)
                    list_boxes_tmp.push_back(box_result[i]);
            }
            list_boxes.swap(list_boxes_tmp);
            list_boxes_tmp.clear();
        }

        // Build Paves and push back them
        for(auto &b:list_boxes){
            Pave* p = new Pave(b, f_list, u, diseable_singleton);
            g->get_node_list().push_back(p);
        }

        /// ****** ADD BASSIN BOXES *******

        for(auto &b:bassin_boxes){
            Pave* p = new Pave(b, f_list, u, diseable_singleton, false);
            p->set_full_out(); // WARNING : Requiered when initial box is too large, and some trajectories can leave !!
            g->get_node_list().push_back(p);
        }
    }

    /// ****** CREATE BORDER EXTRA BOXES *******
    vector<IntervalVector> list_border = m_utils.get_segment_from_box(box, 0.1);

    for(auto &b:list_border){
        Pave* p = new Pave(b, f_list, u, diseable_singleton, false);
        if(border_out)
            p->set_full_out();
        if(border_in)
            p->set_full_in();
        g->get_node_list().push_back(p);
    }

    // ****** REBUILD GRAPH *******
    g->build_graph(); // connect boxes according to continuity
}

Scheduler::Scheduler(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out):
    Scheduler(box, f_list, u, diseable_singleton)
{
    Graph *g = m_graph_list[0];

    /// ****** CREATE BORDER EXTRA BOXES *******
    vector<IntervalVector> list_border = m_utils.get_segment_from_box(box, 0.1);

    for(auto &b:list_border){
        Pave* p = new Pave(b, f_list, u, diseable_singleton, false);
        p->set_external_border(true);
        if(border_out)
            p->set_full_out();
        if(border_in)
            p->set_full_in();
        g->get_node_list().push_back(p);
    }

    // ****** REBUILD GRAPH *******
    g->build_graph(); // connect boxes according to continuity and remove external boundary boxes

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
        m_graph_list[0]->process(process_iterations_max, false);
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
        m_graph_list[0]->process(process_iterations_max, false);

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::compute_attractor(int iterations_max, int process_iterations_max){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;

    int iterations = 0;
    m_graph_list[0]->set_full();

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        m_graph_list[0]->sivia(4,true, false, false, false); // Start with 4 boxes
        m_graph_list[0]->process(process_iterations_max, true, false); // ? Usefull ??? ToDo
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        if(m_graph_list[0]->size()==0 || m_graph_list.size()==0)
            break;
        m_graph_list[0]->clear_node_queue();
        m_graph_list[0]->sivia(2*m_graph_list[0]->size(), true, false, false, false);

        for(int nb_f=0; nb_f<m_graph_list[0]->get_f_size(); nb_f++){
            m_graph_list[0]->set_active_f(nb_f);
            if(nb_f>0)
                m_graph_list[0]->update_queue();

            const clock_t sivia_time = clock();
            cout << "--> time (sivia) = " << float( sivia_time - begin_time ) /  CLOCKS_PER_SEC << endl;

            // Process the backward with the subpaving
            cout << "GRAPH No "<< 0 << " (" << m_graph_list[0]->size() << ")" << endl;
            int graph_list_process_cpt = m_graph_list[0]->process(process_iterations_max, true, false);

            cout << "--> processing outer = " << graph_list_process_cpt << endl;
            cout << "--> time (processing) = " << float( clock() - sivia_time ) /  CLOCKS_PER_SEC << endl;

            // Remove empty pave
            m_graph_list[0]->remove_empty_node();

            // Test if the graph is empty
            if(m_graph_list[0]->is_empty()){
                delete(m_graph_list[0]);
                m_graph_list.erase(m_graph_list.begin());
                cout << "REMOVE EMPTY GRAPH" << endl;
                break;
            }
        }
        if(m_graph_list.size()!=0){
            if(m_graph_list[0]->identify_attractor()){
                cout << "NO MORE ATTRACTOR TO FIND" << endl;
                break;
            }
        }

        cout << "--> time (total) = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;

    }
}

void Scheduler::cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside, bool near_bassin){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;

    int iterations = 0;
    m_graph_list[0]->set_full();

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        m_graph_list[0]->sivia(4,true, false, false, near_bassin); // Start with 4 boxes
        m_graph_list[0]->process(process_iterations_max, true); // ? Usefull ??? ToDo
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        for(int nb_graph=0; nb_graph<m_graph_list.size(); nb_graph++){

            if(m_graph_list[nb_graph]->size()==0 || m_graph_list.size()==0)
                break;
            m_graph_list[nb_graph]->clear_node_queue();
            m_graph_list[nb_graph]->sivia(2*m_graph_list[nb_graph]->size(), true, false, do_not_bisect_inside, near_bassin);
            const clock_t sivia_time = clock();
            cout << "--> time (sivia) = " << float( sivia_time - begin_time ) /  CLOCKS_PER_SEC << endl;

            if(near_bassin){
                m_graph_list[nb_graph]->update_queue();
            }
            // Process the backward with the subpaving
            cout << "GRAPH No "<< nb_graph << " (" << m_graph_list[nb_graph]->size() << ")" << endl;
            int graph_list_process_cpt = m_graph_list[nb_graph]->process(process_iterations_max, true);

            cout << "--> processing outer = " << graph_list_process_cpt << endl;
            cout << "--> time (processing) = " << float( clock() - sivia_time ) /  CLOCKS_PER_SEC << endl;

            // Remove empty pave
            m_graph_list[nb_graph]->remove_empty_node();

            // Test if the graph is empty
            if(m_graph_list[nb_graph]->is_empty()){
                m_graph_list.erase(m_graph_list.begin()+nb_graph);
                delete(m_graph_list[nb_graph]);
                cout << " REMOVE EMPTY GRAPH" << endl;
                if(nb_graph!=0)
                    nb_graph--;
                break;
            }


            // ***************************************************
            //              REMOVE INSIDE PROCEDURE
            // Copy graph & propagate one Pave + intersect with cycle
            // Find the first non-full & non-empty Pave (i.e. on the border)
            /// TODO : improve the algorithm
            if(remove_inside && m_graph_list.size() < graph_max){
                Pave *pave_start = m_graph_list[nb_graph]->get_semi_full_node(); // Find a pave semi full (border pave)
                if(pave_start == NULL)
                    break;

                Graph* graph_propagation = new Graph(m_graph_list[nb_graph], pave_start); // copy graph with 1 activated node (pave_start)
                //                Graph* graph_diff = new Graph(m_graph_list[nb_graph], m_graph_list.size());

                graph_propagation->process(process_iterations_max, false); // process forward

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
        cout << "--> time (total) = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Scheduler::draw(int size, bool filled, string comment){
    for(int i=0; i<m_graph_list.size(); i++){
        m_graph_list[i]->draw(size, filled, comment);
    }
}

Graph* Scheduler::get_graph_list(int i){
    return m_graph_list[i];
}

void Scheduler::print_pave_info(int graph, double x, double y, string color){
    if(m_graph_list.size()>graph){
        m_graph_list[graph]->print_pave_info(x, y, color);
    }
    else{
        cout << "GRAPH NOT FOUND" << endl;
    }
}

void Scheduler::set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY){
    m_utils.m_imageIntegral = new imageIntegral(range, f, t_range, nbBisectionT, nbBisectionXY);
    m_utils.m_imageIntegral_activated = true;
}
