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
    Pave(const ibex::IntervalVector &box, ibex::Function f);
    ~Pave(){}

    bool copy_segment(Pave *p);
    bool equal_segment(Pave *p);

    void process_forward_old();
    void process_backward_old();

    void bisect(vector<Pave *> &result);

    void activate_pave();
    void set_full_continuity();

    void set_theta(ibex::Interval theta);
    ibex::IntervalVector get_border_position(int face);
    double get_theta_diam();

    // Drawing functions
    void draw(bool filled, string color="b[]");
    void draw_borders(bool filled);

    void compute_flow();

    bool get_brother_empty(int level=1);
    void remove_from_brothers();
    void remove_brothers(Pave* p, int face);
    bool all_brothers_full(int level=1);

    bool is_empty();
    bool is_full();
    void set_empty(bool val);
    bool set_full(bool val);

    bool netwon_test();

    vector<Pave*> get_brothers(int face);
    std::vector<ibex::Interval> rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y);

/***************** Variables ******************/

public:
    ibex::Interval theta[2];
    ibex::IntervalVector box;
    std::vector<Border> borders;

private:
    bool empty;
    bool full;
};

#endif // PAVE_H
