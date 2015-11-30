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
    void bisect(std::vector<Pave> *result);
    void computePropagation(ibex::Interval seg_in, int face);

/***************** Variables ******************/

public:
    ibex::Interval theta;
    ibex::Interval speed;

    ibex::IntervalVector box;

    ibex::Interval table_rotation[4] = {-ibex::Interval::HALF_PI, ibex::Interval::PI, ibex::Interval::HALF_PI, ibex::Interval(0.0)};

    std::vector<Border> queue;
    std::vector<Border> borders;
};

#endif // PAVE_H
