#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>

class Border;
class Pave
{
/***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &box);
    ~Pave(){}

    void draw();
    void process();
    void cut(std::vector<Pave> *result);
    void project(std::vector<ibex::Interval> &seg_out, ibex::Interval seg_in, ibex::Interval theta, ibex::Interval c0, ibex::Interval c1);
    void computePropagation(ibex::Interval seg_in, int face);

/***************** Variables ******************/

public:
    ibex::Interval theta;
    ibex::Interval speed;

    ibex::IntervalVector box;

    ibex::Interval table_rotation[4] = {ibex::Interval(0.0), ibex::Interval::PI/2.0, ibex::Interval::PI, 3.0*ibex::Interval::PI/2.0};

    std::vector<Border> queue;
    std::vector<Border> borders;
};

#endif // PAVE_H
