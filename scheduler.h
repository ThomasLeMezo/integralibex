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
    void add_segment(int id_box);
    void add_segment(double x, double y);
    void set_initial_pave(const ibex::IntervalVector &box);
    std::vector<ibex::Interval> rotate(ibex::Interval theta, ibex::Interval x, ibex::Interval y);

    Pave* get_pave(double x, double y);
    void print_pave_info(double x, double y);


    // Contractor
    void CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dy);
    void CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy);

public:
    std::vector<Pave*> pave_list;
    std::vector<Pave*> pave_queue;

    ibex::CtcPolar contract_polar;
};

#endif // SCHEDULER_H
