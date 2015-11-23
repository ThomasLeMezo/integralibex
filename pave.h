#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>

class Pave
{
/***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &b);
    ~Pave(){}

    void draw();
    void process();
    void cut(std::vector<Pave> *result);

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
