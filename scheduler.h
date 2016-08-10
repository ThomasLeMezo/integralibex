#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "graph.h"
#include "ibex.h"
#include "utils.h"

class Pave;
class Border;
class Scheduler
{
/***************** Functions ******************/
public:
    Scheduler(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton);
    Scheduler(const ibex::IntervalVector &box, const vector<ibex::IntervalVector> &remove_boxes, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_in=true, bool border_out=true);
    Scheduler(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_inner_in, bool border_inner_out, bool border_outer_in, bool border_outer_out);
    ~Scheduler();

    void cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside=false);
    void cameleon_propagation(int iterations_max, int process_iterations_max, const vector<ibex::IntervalVector> &initial_boxes);
    void cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe);
    void cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, const vector<ibex::IntervalVector> &initial_boxes);
    void cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe);
    void compute_attractor(int iterations_max, int process_iterations_max);
    void cameleon_viability(int iterations_max, int process_iterations_max, bool border_condition=false);
    void find_path(int iterations_max, int process_iterations_max, const ibex::IntervalVector &boxA, const ibex::IntervalVector &boxB);

    void attractor_to_kernel();

    void set_symetry(ibex::Function *f, int face_in, int face_out);
//    void set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY);

    // ******** Drawing functions ********
    void draw(int size, bool filled, string comment="");
    void print_pave_info(int graph, double x, double y, string color="black[g]");

    // Getter
    Graph* get_graph_list(int i);


/***************** Variables ******************/
private:
    std::vector<Graph*> m_graph_list;

    Utils m_utils;
};

#endif // SCHEDULER_H
