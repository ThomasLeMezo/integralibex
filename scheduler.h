#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "pave.h"
#include "ibex.h"
#include "utils.h"

class Pave;
class Border;
class Scheduler
{
/***************** Functions ******************/
public:
    Scheduler();
    ~Scheduler();

    void process(std::vector<Pave*> &pave_queue, int max_iterations, bool backward);
    void SIVIA(std::vector<Pave*> &pave_list, std::vector<Pave*> &pave_queue, double epsilon_theta, int iterations_max, bool backward);
    void process_SIVIA_cycle(int iterations_max, int pave_max, int process_iterations_max);

    void set_full(std::vector<Pave *> &pave_list);
    void set_initial_pave(const ibex::IntervalVector &box, ibex::Function *f);
    Pave* get_pave(std::vector<Pave*> &pave_list, double x, double y);
    void activate_pave(std::vector<Pave*> &pave_list, std::vector<Pave *> &pave_queue, double x, double y);

    vector<Pave *> copy_graph(int graph, bool empty, Pave *pave_keep);

    // ******** Drawing functions ********
    void draw(int size, bool filled);
    void draw(vector<Pave *> pave_list, int size, bool filled);

    // ******** Utils functions ********
    void print_pave_info(int graph, double x, double y, string color="g[]");

/***************** Variables ******************/
public:
    std::vector<std::vector<Pave*> > m_global_pave_list;
    std::vector<std::vector<Pave*> > m_global_pave_list_empty;
    std::vector<std::vector<Pave*> > m_global_pave_queue;

    int m_draw_nb;
    Utils m_utils;
};

#endif // SCHEDULER_H
