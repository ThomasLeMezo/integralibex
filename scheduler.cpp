#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"
#include "omp.h"

#include "graphdot.h"
#include "imageintegral.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton){
    m_graph_id = -1;
    m_graph_list.push_back(new Graph(box, f_list, &m_utils, get_graph_id(), diseable_singleton));
}

Scheduler::~Scheduler(){
    for(Graph *graph:m_graph_list){
        delete(graph);
    }
}

Scheduler::Scheduler(const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_in, bool border_out){
    m_graph_id = -1;
    Graph *g = new Graph(&m_utils, get_graph_id());
    m_graph_list.push_back(g);

    vector<IntervalVector> list_boxes;
    list_boxes.push_back(box);

    // Build diff boxes
    for(IntervalVector box_bassin:bassin_boxes){
        vector<IntervalVector> list_boxes_tmp;
        for(IntervalVector b:list_boxes){
            std::vector<IntervalVector> box_result = m_utils.diff(b, box_bassin);

            for(int i=0; i<box_result.size(); i++)
                list_boxes_tmp.push_back(box_result[i]);
        }
        list_boxes.swap(list_boxes_tmp);
        list_boxes_tmp.clear();
    }

    // Build Paves and push back them
    for(IntervalVector b:list_boxes){
        Pave* p = new Pave(b, f_list, diseable_singleton);
        p->set_full_all();
        g->push_back(p);
    }

    /// ****** ADD BASSIN BOXES *******

    for(IntervalVector b:bassin_boxes){
        Pave* p = new Pave(b, f_list, diseable_singleton, false);
        p->set_external_border(true);
        p->set_inner_mode(false);
        p->set_full();

        p->set_inner_mode(true);
        p->set_full_out(); // WARNING : Requiered when initial box is too large, and some trajectories can leave !!
        //        p->set_continuity_in(false);
        //        p->set_continuity_out(false);

        p->set_bassin(true);

        g->push_back(p); // Inactive box
    }

    /// ****** CREATE BORDER EXTRA BOXES *******
    vector<IntervalVector> list_border = m_utils.get_segment_from_box(box, 0.1);

    for(IntervalVector b:list_border){
        Pave* p = new Pave(b, f_list, diseable_singleton, false);
        p->set_external_border(true);
        p->set_inner_mode(true); p->set_full();

        if(border_out){
            p->set_inner_mode(false); p->set_full_out();
        }
        if(border_in){
            p->set_inner_mode(false); p->set_full_in();
        }
        g->push_back(p); // Box is put into external box list when building the graph
    }

    // ****** REBUILD GRAPH *******
    g->build_graph(); // connect boxes according to continuity
    g->set_inner_mode(true);
}

Scheduler::Scheduler(const IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_inner_in, bool border_inner_out, bool border_outer_in, bool border_outer_out):
    Scheduler(box, f_list, diseable_singleton)
{
    Graph *g = m_graph_list[0];

    /// ****** CREATE BORDER EXTRA BOXES *******
    vector<IntervalVector> list_border = m_utils.get_segment_from_box(box, 0.1);

    for(IntervalVector b:list_border){
        Pave* p = new Pave(b, f_list, diseable_singleton, false);
        p->set_external_border(true);
        if(border_inner_in)
            p->set_full_inner_in();
        if(border_inner_out)
            p->set_full_inner_out();
        if(border_outer_in)
            p->set_full_outer_in();
        if(border_outer_out)
            p->set_full_outer_out();
        g->push_back(p); // Box is put into external box list when building the graph
    }

    // ****** REBUILD GRAPH *******
    g->build_graph(); // connect boxes according to continuity and remove external boundary boxes

}

void Scheduler::set_symetry(Function *f, int face_in, int face_out){
    Graph *g = m_graph_list[0];
    g->set_symetry(f, face_in, face_out);
}

void Scheduler::push_back_inside_curve(ibex::Function *curve){
    Graph *g = m_graph_list[0];
    g->push_back_inside_curve(curve);
}

void Scheduler::cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe){
    vector<IntervalVector> initial_boxes;
    initial_boxes.push_back(initial_boxe);
    cameleon_propagation(iterations_max, process_iterations_max, initial_boxes);
}

void Scheduler::cameleon_propagation(int iterations_max, int process_iterations_max, const vector<IntervalVector> &initial_boxes){

    Graph *graph = m_graph_list[0];
    if(m_graph_list.size()!=1 && graph->size() !=1)
        return;
    int iterations = 0;

    if(iterations < iterations_max && graph->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph->sivia(4,false, false, false); // Start with 4 boxes
        graph->set_empty();
        for(IntervalVector initial_box:initial_boxes)
            graph->set_active_pave(initial_box);
        graph->process(process_iterations_max, false);
        graph->mark_empty_node();
        iterations++;
    }
    int nb_graph = 0;

    while(iterations < iterations_max){

        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph->mark_empty_node();
        graph->sivia(2*graph->get_alive_node(), false, false, false);

        graph->set_empty();
        for(IntervalVector initial_box:initial_boxes)
            graph->set_active_pave(initial_box);

        // Process the forward with the subpaving
        cout << "GRAPH No "<< nb_graph << " (" << graph->size() << ")" << endl;
        graph->process(process_iterations_max, false);

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe){
    vector<IntervalVector> initial_boxes;
    initial_boxes.push_back(initial_boxe);
    cameleon_propagation_with_inner(iterations_max, process_iterations_max, initial_boxes);
}

void Scheduler::cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, const vector<IntervalVector> &initial_boxes){
    Graph *graph = m_graph_list[0];
    if(m_graph_list.size()!=1 && graph->size() !=1)
        return;
    int iterations = 0;

    graph->set_compute_inner(true);

    graph->set_inner_mode(false);
    graph->set_external_boundary(false, false);

    graph->set_inner_mode(true);
    graph->set_external_boundary(true, true);
    graph->set_empty_outer_full_inner();

    if(iterations < iterations_max && graph->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        graph->sivia(4, false, false, false); // Start with 4 boxes
        graph->set_empty_outer_full_inner();
        graph->set_active_outer_inner(initial_boxes);

        // Inner
        graph->set_inner_mode(true);
        graph->set_backward_function(true);
        graph->process(process_iterations_max, true);

        // Outer
        graph->set_inner_mode(false);
        graph->set_backward_function(false);
        graph->process(process_iterations_max, false);

        graph->mark_empty_node();
        iterations++;
    }
    int nb_graph = 0;

    while(iterations < iterations_max){

        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph->mark_empty_node();
        //        if(iterations>12 && iterations <16)
        //            graph->sivia(2*graph->get_alive_node(), false, false, false, true);
        //        else
        graph->sivia(2*graph->get_alive_node(), false, false, false);
        graph->set_empty_outer_full_inner();
        graph->set_active_outer_inner(initial_boxes);

        // Process the forward with the subpaving
        cout << "GRAPH No "<< nb_graph << " (" << graph->size() << ")" << endl;
        // Inner
        graph->set_inner_mode(true);
        graph->set_backward_function(true);
        graph->process(process_iterations_max, true);

        // Outer
        graph->set_inner_mode(false);
        graph->set_backward_function(false);
        graph->process(process_iterations_max, false);

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::compute_attractor(int iterations_max, int process_iterations_max){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;

    int iterations = 0;
    Graph *graph = m_graph_list[0];
    ///////////////////////////// INNER MODE //////////////////////////////
    graph->set_full();

    if(iterations < iterations_max && graph->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph->sivia(4,true, false, false); // Start with 4 boxes
        graph->process(process_iterations_max, true); // ? Usefull ??? ToDo
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        if(graph->get_alive_node()==0 || m_graph_list.size()==0)
            break;
        graph->clear_node_queue();
        graph->sivia(2*graph->get_alive_node(), true, false, false);

        for(int nb_f=0; nb_f<graph->get_f_size(); nb_f++){
            graph->set_active_f(nb_f);
            if(nb_f>0)
                graph->update_queue();

            const clock_t sivia_time = clock();
            cout << "--> time (sivia) = " << float( sivia_time - begin_time ) /  CLOCKS_PER_SEC << endl;

            // Process the backward with the subpaving
            cout << "GRAPH No " << 0 << " (" << graph->size() << ")" << endl;
            graph->process(process_iterations_max, true);
            cout << "--> time (processing) = " << float( clock() - sivia_time ) /  CLOCKS_PER_SEC << endl;

            // Remove empty pave
            graph->mark_empty_node();

        }

        if(graph->identify_attractor()){
            cout << "NO MORE ATTRACTOR TO FIND" << endl;
            graph->push_back_pos_attractor();
            break;
        }

        cout << "--> time (total) = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;

    }
}

void Scheduler::cameleon_viability(int iterations_max, int process_iterations_max, bool border_condition){
    cout << endl << endl;
    cout << "************ COMPUTE KERNEL ************" << endl;
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;
    int iterations = 0;

    Graph *graph = m_graph_list[0];
    graph->set_active_f(-1);
    graph->compute_all_propagation_zone();
    //    graph->debug_marker2 = true;

    if(iterations < iterations_max && graph->size()>4){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        cout << "GRAPH No 0 (" << graph->size() << ")" << endl;
        graph->mark_empty_node();

        graph->set_inner_mode(true);
        graph->update_queue(border_condition, true);
        graph->process(process_iterations_max, true);

        graph->set_inner_mode(false);
        graph->update_queue(true, true);
        graph->process(process_iterations_max, true);

        cout << "--> time (processing) = " << float( clock() - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        if(graph->get_alive_node()==0)
            break;

        graph->sivia(2*graph->get_alive_node()+graph->size(), true, false, false);
        const clock_t sivia_time = clock();
        cout << "--> time (sivia) = " << float( sivia_time - begin_time ) /  CLOCKS_PER_SEC << endl;

        //        if(iterations == 0)
        //            graph->debug_marker2 = true;

        // Process the backward with the subpaving
        cout << "GRAPH No 0 (" << graph->size() << ")" << endl;
        graph->mark_empty_node();

        graph->set_inner_mode(true);
        graph->update_queue(border_condition, true);
        graph->process(process_iterations_max, true);

        graph->set_inner_mode(false);
        graph->update_queue(true, true);
        graph->process(process_iterations_max, true, true);

        cout << "--> time (processing) = " << float( clock() - sivia_time ) /  CLOCKS_PER_SEC << endl;

        // Remove empty pave
        graph->mark_empty_node();

        // Test if the graph is empty
        if(graph->is_empty()){
            cout << "EMPTY GRAPH" << endl;
            break;
        }

        cout << "--> time (total) = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

void Scheduler::cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside, bool stop_first_pos_invariant){
    if(this->m_graph_list.size()<1 && this->m_graph_list[0]->size() <1)
        return;
    int iterations = 0;
    Graph *graph_initial = m_graph_list[0];
    graph_initial->get_node_list()[0]->set_full();

    if(iterations < iterations_max && graph_initial->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph_initial->sivia(4,true, false, false); // Start with 4 boxes
        while(!graph_initial->is_sufficiently_discretized()){
            graph_initial->reset_queues();
            graph_initial->sivia(2*graph_initial->get_alive_node(),true, false, false);
        }
        graph_initial->process(process_iterations_max, true); // ? Usefull ??? ToDo
        iterations++;
    }

    while(iterations < iterations_max){
        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        for(int nb_graph=0; nb_graph<m_graph_list.size(); nb_graph++){
            Graph *graph = m_graph_list[nb_graph];
            cout << "GRAPH No "<< graph->get_graph_id() << " (" << graph->size() << ")" << endl;

            if(graph->get_alive_node()==0 || m_graph_list.size()==0)
                break;
            if(!(stop_first_pos_invariant && graph->get_positive_invariant())){
                graph->clear_node_queue();
                graph->sivia(2*graph->get_alive_node(), true, false, do_not_bisect_inside);
                const clock_t sivia_time = clock();
                cout << "--> time (sivia) = " << float( sivia_time - begin_time ) /  CLOCKS_PER_SEC << endl;

                // Process the backward with the subpaving
                int graph_list_process_cpt = graph->process(process_iterations_max, true);

                cout << "--> processing outer = " << graph_list_process_cpt << endl;
                cout << "--> time (processing) = " << float( clock() - sivia_time ) /  CLOCKS_PER_SEC << endl;

                // Remove empty pave
                graph->mark_empty_node();

                if(graph->is_empty() && m_graph_list.size()>1){
                    m_graph_list.erase(m_graph_list.begin()+nb_graph);
                    cout << "--> remove empty graph" << endl;
                }
                else{
                    // Test if positive invariant
                    if(graph->is_positive_invariant()){
                        cout << "--> graph IS positive invariant" << endl;
                        graph->push_back_pos_attractor();
                        graph->set_positive_invariant(true);
                        emit publishLog("Invariant found at step " + QString::number(iterations+1));
                    }
                    else{
                        cout << "--> graph IS NOT positive invariant" << endl;
                    }

                    // ***************************************************
                    //              REMOVE INSIDE PROCEDURE
                    // Copy graph & propagate one Pave + intersect with cycle
                    // Find the first non-full & non-empty Pave (i.e. on the border)
                    /// TODO : improve the algorithm
                    if(remove_inside && m_graph_list.size() < graph_max && !(stop_first_pos_invariant && graph->get_positive_invariant())){
                        Pave *pave_start = graph->get_semi_full_node(); // Find a pave semi full (border pave)
                        if(pave_start == NULL){
                            cout << " REMOVE INSIDE - NO START PAVE FOUND" << endl;
                            break;
                        }
                        else{
                            cout << "--> pave selected = " << pave_start->get_position() << endl;
                        }

                        Graph* graph_propagation = new Graph(graph, pave_start); // copy graph with 1 activated node (pave_start)
                        Graph* graph_diff = new Graph(graph, get_graph_id());

                        graph_propagation->process(process_iterations_max, false); // process forward
                        graph->inter(*graph_propagation); // intersect the graph with the propagation graph
                        graph->reset_pave_segment_list();
                        graph_diff->diff(*graph);
                        graph_diff->mark_empty_node();

                        if(!graph_diff->is_empty()){ // If there is an inside, add to graph_list
                            cout << "--> sucess to separate cycle" << endl;
                            graph_diff->set_positive_invariant(false);
                            graph_diff->reset_pave_segment_list();
                            graph->set_positive_invariant(false);
                            m_graph_list.push_back(graph_diff);
                            graph_diff->clear_node_queue();
                            cout << " ADD NEW GRAPH No " << graph_diff->get_graph_id() << endl;
                            //m_graph_list.back()->print();
                        }
                        else{
                            cout << "--> not possible to separate cycle" << endl;
                            delete(graph_diff);
                        }
                        graph->mark_empty_node();

                        delete(graph_propagation);

                    }
                }
            }
        }
        cout << "--> time (total) = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        emit iteration_status(iterations+1, iterations_max);
        iterations++;
    }

    if(m_graph_list.size()>1)
        emit publishLog(QString::number(m_graph_list.size()) + " possible cycles were found");
//    else
        //        emit publishLog(QString::number(m_graph_list.size()) + " possible cycle was found");
}

void Scheduler::find_path(int iterations_max, int process_iterations_max, const ibex::IntervalVector &boxA, const ibex::IntervalVector &boxB){
    Graph *graph = m_graph_list[0];
    if(m_graph_list.size()!=1 && graph->size() !=1)
        return;
    int iterations = 0;
    graph->set_compute_inner(true);

    graph->set_inner_mode(false);
    graph->set_external_boundary(false, false);

    graph->set_inner_mode(true);
    graph->set_external_boundary(true, true);

    if(iterations < iterations_max && graph->size()<4){
        cout << "************ ITERATION = " << iterations << " ************" << endl;

        graph->sivia(4, false, false, false); // Start with 4 boxes
        graph->set_empty_outer_full_inner();
        graph->clear_node_queue();
        Graph *graph_backward = new Graph(graph);

        //// Outer Graph
        graph->set_active_outer_inner(boxA);
        // INNER
        graph->set_inner_mode(true);
        graph->set_backward_function(true);
        graph->process(process_iterations_max, true);
        // OUTER
        graph->set_inner_mode(false);
        graph->set_backward_function(false);
        graph->process(process_iterations_max, false, true);

        /// Backward graph
        graph_backward->set_active_outer_inner(boxB);
        //        // INNER
        graph_backward->set_inner_mode(true);
        graph_backward->set_backward_function(false); // Invert bc of bwd
        graph_backward->process(process_iterations_max, true);
        // OUTER
        graph_backward->set_inner_mode(false);
        graph_backward->set_backward_function(true); // Invert bc of bwd
        graph_backward->process(process_iterations_max, false, true);

        // Intersect graph
        graph->inter(graph_backward, true);
        delete(graph_backward);
        graph->mark_empty_node();
        iterations++;
    }

    while(iterations < iterations_max){

        const clock_t begin_time = clock();
        cout << "************ ITERATION = " << iterations << " ************" << endl;
        graph->mark_empty_node();
        if(graph->is_empty()){
            cout << "THERE IS NO PATH TO LINK THE TWO BOXES" << endl;
            cout << "size = " << graph->size() << endl;
            break;
        }
        graph->sivia(2*graph->get_alive_node(), false, false, false);

        graph->set_empty_outer_full_inner();
        graph->clear_node_queue();
        Graph *graph_backward = new Graph(graph);

        /// Forward graph
        cout << "GRAPH FWD (" << graph->size() << ")" << endl;
        graph->set_active_outer_inner(boxA);
        // INNER
        graph->set_inner_mode(true);
        graph->set_backward_function(true);
        graph->process(process_iterations_max, true);
        // OUTER
        graph->set_inner_mode(false);
        graph->set_backward_function(false);
        graph->process(process_iterations_max, false, true);

        /// Backward graph
        //        cout << "GRAPH BWD (" << graph_backward->size() << ")" << endl;
        graph_backward->set_active_outer_inner(boxB);
        //        // INNER
        graph_backward->set_inner_mode(true);
        graph_backward->set_backward_function(false); // Invert bc of bwd
        graph_backward->process(process_iterations_max, true);
        // OUTER
        graph_backward->set_inner_mode(false);
        graph_backward->set_backward_function(true); // Invert bc of bwd
        graph_backward->process(process_iterations_max, false, true);

        //        graph->draw(1024, true, "fwd");
        //        graph_backward->draw(1024, true, "bwd");
        //        cin.ignore();

        graph->inter(graph_backward, true);
        delete(graph_backward);

        //        graph->draw(1024, true, "inter");

        cout << "--> graph_time = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
        iterations++;
    }
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Scheduler::draw(int size, bool filled, string comment, bool positive_invariant){
    int position = 0;
    for(Graph *g:m_graph_list){
        g->draw(size, filled, comment, false, position, positive_invariant);
        position += 100;
    }
}

Graph* Scheduler::get_graph_list(int i){
    return m_graph_list[i];
}

void Scheduler::print_pave_info(int graph, double x, double y, string color){
    if(m_graph_list.size()>graph){
        Graph *g = m_graph_list[graph];
        g->print_pave_info(x, y, color);
    }
    else{
        cout << "GRAPH NOT FOUND" << endl;
    }
}

//void Scheduler::set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY){
//    m_utils.m_imageIntegral = new imageIntegral(range, f, t_range, nbBisectionT, nbBisectionXY);
//    m_utils.m_imageIntegral_activated = true;
//}

void Scheduler::attractor_to_kernel(){
    cout << endl << endl;
    cout << "*******************" << endl;
    cout << "ATTRACTOR TO KERNEL" << endl;

    // reset removed and active pave
    Graph *graph = m_graph_list[0];
    graph->set_all_active();

    /// OUTER
    cout << "-> configure outer" << endl;
    graph->set_inner_mode(false);
    graph->copy_to_inner();
    graph->set_full();
    graph->set_external_boundary(false, true);

    /// INNER
    cout << "-> configure inner" << endl;
    graph->set_inner_mode(true);
    graph->complementaire();
    graph->set_external_boundary(true, true);

    cout << "-> compute propagation zone" << endl;
    graph->compute_all_propagation_zone();

    graph->mark_empty_node();
}

int Scheduler::get_graph_id(){
    m_graph_id++;
    return m_graph_id;
}

void Scheduler::set_inner_mode(bool val){
    for(Graph *g:m_graph_list)
        g->set_inner_mode(val);
}

void Scheduler::set_external_boundary(bool in, bool out){
    for(Graph *g:m_graph_list)
        g->set_external_boundary(in, out);
}
