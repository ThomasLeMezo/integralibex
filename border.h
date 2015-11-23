#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>

class Pave;
class Border
{
public:
    Border(const ibex::IntervalVector& p, int face, Pave *pave);

// State Variable
public:
    std::vector<ibex::Interval> segments; // List of impacted interval of the segment
    int face;
    std::vector<Border> borthers;
    double offset;
    ibex::IntervalVector position;

    Pave *pave;
};

#endif // BORDER_H
