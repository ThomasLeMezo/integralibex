#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <pave.h>

class Pave;
class Scheduler;
class Border
{
public:
    Border(const ibex::IntervalVector& position, const int face, Pave *pave);
    Border(const ibex::Interval &segment, const int face);

    void draw() const;
    ibex::Interval add_segment(ibex::Interval seg);
    void publish_to_borthers(ibex::Interval seg);
    void add_brothers(std::vector<Border *> brother_list);
    void update_brothers(Border* border_pave1, Border* border_pave2);

// State Variable
public:
    std::vector<ibex::Interval> segments;   // List of impacted interval of the segment
    int face;                               // Number of the face (0=bottom, 1=right, ...)
    std::vector<Border*> brothers;          // Pointer to brothers Borders
    ibex::IntervalVector position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *pave;                             // Pointer to its container
};

#endif // BORDER_H
