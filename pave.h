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

    void process();
    void bisect(vector<Pave *> &result);

    void add_new_segment(Border &b);
    void warn_scheduler();
    void activate_pave();

    void set_theta(ibex::Interval theta);

    void draw(bool filled);
    void draw_borders(bool filled);

    void compute_successors();

    void add_precursors(Pave* p);
    void add_successors(Pave* p);
    void clear_graph();

    bool test_cycle(Pave *p_test, int depth, int depth_max);

/***************** Variables ******************/

public:
    Scheduler *scheduler;

    ibex::Interval theta[2];
    ibex::Interval speed;

    ibex::IntervalVector box;

    std::vector<Border> queue;
    std::vector<Border> borders;

    std::vector<Pave*> precursors;
    std::vector<Pave*> successors;

    bool visited_node;
};

#endif // PAVE_H
