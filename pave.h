#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>
#include <scheduler.h>

class Border;
class Scheduler;
class Pave
{
/***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &box, Scheduler *scheduler);
    ~Pave(){}

    void draw() const;
    void process();
    void bisect(vector<Pave *> &result);
    void computePropagation(ibex::Interval seg_in, int face);
    void push_queue(Border &b);
    void warn_scheduler();
    void activate_pave();

/***************** Variables ******************/

public:
    Scheduler *scheduler;

    vector<ibex::Interval> theta;
    ibex::Interval speed;

    ibex::IntervalVector box;

    ibex::Interval table_rotation[4] = {-ibex::Interval::HALF_PI, ibex::Interval::PI, ibex::Interval::HALF_PI, ibex::Interval(0.0)};

    std::vector<Border> queue;
    std::vector<Border> borders;
};

#endif // PAVE_H
