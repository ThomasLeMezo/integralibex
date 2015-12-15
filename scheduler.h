#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "pave.h"
#include "ibex.h"
#include "border.h"
#include "utils.h"

class Pave;
class Border;
class Scheduler
{
public:
    Scheduler();
    ~Scheduler(){}

    void add_to_queue(Pave* pave, bool forward);
    void draw(int size, bool filled);
    void process(int max_iterations);
    void process_graph(int iterations_max, int pave_max);
    void process_backward(int max_iterations);

    void SIVIA(double epsilon_theta, int iterations_max);
    void add_segment(int id_box);
    void add_segment(double x, double y);
    void set_initial_pave(const ibex::IntervalVector &box);

    Pave* get_pave(double x, double y);
    void print_pave_info(double x, double y, string color);

    void set_full_continuity();
    void set_empty_borders();

    Utils utils;

public:
    std::vector<Pave*> pave_list;
    std::vector<Pave*> pave_queue_forward,pave_queue_backward;

    int draw_nb;
};

#endif // SCHEDULER_H
