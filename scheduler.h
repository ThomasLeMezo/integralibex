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
    Scheduler(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, const ibex::IntervalVector &u, bool diseable_singleton=false);
    Scheduler(const ibex::IntervalVector &box, const vector<ibex::IntervalVector> &remove_boxes, const std::vector<ibex::Function *> &f_list, const ibex::IntervalVector &u, bool diseable_singleton=false);
    ~Scheduler();

    void cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool inner, bool do_not_bisect_inside=false);
    void cameleon_propagation(int iterations_max, int process_iterations_max, vector<ibex::IntervalVector> &initial_boxes, bool inner=false);
    void cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe, bool inner=false);

    void set_symetry(ibex::Function *f, int face_in, int face_out);
    void set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY);

    void define_continuity(Graph *g, const ibex::IntervalVector &test_box, bool continuity_in, bool continuity_out);
    // ******** Drawing functions ********
    void draw(int size, bool filled);
    void print_pave_info(int graph, double x, double y, string color);

    // Getter
    Graph* get_graph_list(int i);


/***************** Variables ******************/
private:
    std::vector<Graph*> m_graph_list;
    std::vector<Graph*> m_graph_inner_list;

    Utils m_utils;
};

#endif // SCHEDULER_H
