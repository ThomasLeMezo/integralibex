#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(const IntervalVector &box, ibex::Function *f){
    m_graph_list.push_back(new Graph(box, f, &m_utils, 0));
}

Scheduler::~Scheduler(){
    for(auto &graph:m_graph_list){
        delete(graph);
    }
}

void Scheduler::cameleon_propagation(int iterations_max, int process_iterations_max, IntervalVector &initial_box, int max_symetry){
    if(m_graph_list.size()!=1 && m_graph_list[0]->size() !=1)
        return;
    int iterations = 0;

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        m_graph_list[0]->sivia(0.0,4,false, false); // Start with 4 boxes
        m_graph_list[0]->process(process_iterations_max, false);
        iterations++;
    }
    int nb_graph = 0;

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        m_graph_list[0]->sivia(0.0, 2*m_graph_list[0]->size(), false, true);

        m_graph_list[0]->set_empty();
        m_graph_list[0]->set_active_pave(initial_box);

        // Process the forward with the subpaving
        cout << " GRAPH No "<< nb_graph << " (" << m_graph_list[0]->size() << ")" << endl;
        int symetry = 0;

        m_graph_list[0]->process(process_iterations_max, false);
//        while(!=0 && symetry < max_symetry){
//            m_graph_list[0]->set_y_symetry();
//            symetry++;
//        }

        // Test if the graph is empty
        if(m_graph_list[0]->is_empty()){
            cout << "GRAPH EMPTY" << endl;
            break;
        }

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;

    int iterations = 0;
    m_graph_list[0]->set_full();

    if(iterations < iterations_max && this->m_graph_list[0]->size()<4){
        m_graph_list[0]->sivia(0.0,4,false, false); // Start with 4 boxes
        m_graph_list[0]->process(process_iterations_max, false);
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        for(int nb_graph=0; nb_graph<m_graph_list.size(); nb_graph++){

            if(m_graph_list[nb_graph]->size()==0 || m_graph_list.size()==0)
                break;

            m_graph_list[nb_graph]->sivia(0.0, 2*m_graph_list[nb_graph]->size(), true, true);

            // Process the backward with the subpaving
            cout << " GRAPH No "<< nb_graph << " (" << m_graph_list[nb_graph]->size() << ")" << endl;
            m_graph_list[nb_graph]->process(process_iterations_max, true);

            // Remove empty pave
            m_graph_list[nb_graph]->remove_empty_node();

            // Test if the graph is empty
            if(m_graph_list[nb_graph]->is_empty()){
                m_graph_list.erase(m_graph_list.begin()+nb_graph);
                if(nb_graph!=0)
                    nb_graph--;
                break;
            }

            // ***************************************************
            // Copy graph & propagate one Pave + intersect with cycle
            // Find the first non-full & non-empty Pave
            if(remove_inside && m_graph_list.size() < graph_max){
                cout << "REMOVE INSIDE" << endl;
                Pave *pave_start = m_graph_list[nb_graph]->get_semi_full_node(); // Find a pave semi full
                if(pave_start == NULL)
                    break;

                Graph* graph_propagation = new Graph(m_graph_list[nb_graph], pave_start); // copy graph with 1 activated node (pave_start)
                Graph* graph_diff = new Graph(m_graph_list[nb_graph], m_graph_list.size());

                graph_propagation->process(process_iterations_max, false); // process forward

                m_graph_list[nb_graph]->inter(*graph_propagation); // intersect the graph with the propagation graph

                graph_diff->diff(*m_graph_list[nb_graph]);

                if(!graph_diff->is_empty()){ // If there is an inside, add to graph_list
                    m_graph_list.push_back(graph_diff);
                }
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
    for(auto &graph:m_graph_list){
        graph->draw(size, filled);
    }
}
