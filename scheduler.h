#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "pave.h"
#include "ibex.h"
#include "border.h"

class Pave;
class Border;
class Scheduler
{
public:
    Scheduler();
    ~Scheduler(){}

    void add_to_queue(Pave* pave);
    void draw();
    void process(int max_iterations);
    void SIVIA(double epsilon_theta, int iterations_max);
    void add_segment();
    void set_initial_pave(const ibex::IntervalVector &box);

public:
    std::vector<Pave*> pave_list;
    std::vector<Pave*> pave_queue;
};

#endif // SCHEDULER_H
