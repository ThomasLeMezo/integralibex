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
    Scheduler(const ibex::IntervalVector &box, ibex::Function *f);
    ~Scheduler();

    void cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside);
    void cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_box, int max_symetry);

    void graph_symetry(vector<Pave *> &pave_list, vector<Pave *> &pave_queue);

    // ******** Drawing functions ********
    void draw(int size, bool filled);


/***************** Variables ******************/
public:
    std::vector<Graph*> m_graph_list;

    Utils m_utils;
};

#endif // SCHEDULER_H
