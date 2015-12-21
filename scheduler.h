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
    ~Scheduler(){}

    void process(int max_iterations);
    void SIVIA(double epsilon_theta, int iterations_max, bool not_full_test);
    void process_SIVIA_cycle(int iterations_max, int pave_max, int process_iterations_max);

    void set_full_continuity();
    void set_initial_pave(const ibex::IntervalVector &box, ibex::Function *f);
    Pave* get_pave(double x, double y);

    // ******** Drawing functions ********
    void draw(int size, bool filled);

    // ******** Utils functions ********
    void print_pave_info(double x, double y, string color);



/***************** Variables ******************/
public:
    std::vector<Pave*> pave_list;
    std::vector<Pave*> pave_list_empty;
    std::vector<Pave*> pave_queue;

    int draw_nb;
    Utils utils;
};

#endif // SCHEDULER_H
