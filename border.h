#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <pave.h>

class Pave;

class Border
{
public:
    Border(const ibex::IntervalVector& position, int face, Pave *pave);
    Border(const ibex::Interval &segment, int face);

    void draw();
    void add_segement(ibex::Interval seg);
    void publish_to_borthers(ibex::Interval seg);

// State Variable
public:
    std::vector<ibex::Interval> segments;   // List of impacted interval of the segment
    int face;                               // Number of the face (0=bottom, 1=right, ...)
    std::vector<Border*> borthers;          // Pointer to brothers Borders
    ibex::IntervalVector position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *pave;                             // Pointer to its container
};

#endif // BORDER_H
