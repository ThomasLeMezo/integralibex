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
    Border(const ibex::IntervalVector& position, const int face, const ibex::Interval &segment);
    ~Border(){}

    void draw() const;
    void get_points(std::vector<double> &x, std::vector<double> &y);

    std::vector<ibex::Interval> add_segment(const ibex::Interval &seg);
    bool plug_segment(ibex::Interval &input, ibex::Interval &position, bool modify);

    void publish_to_borthers(ibex::Interval seg, bool forward);
    void add_brothers(std::vector<Border *> brother_list);
    void update_brothers(Border* border_pave1, Border* border_pave2);

    void set_full();
    bool is_empty();
    bool is_full();

// State Variable
public:
    ibex::Interval segment;   // List of impacted interval of the segment
    // We suppose that there is only ONE segment and not a list (otherwise, SIVIA requiered)

    int face;                               // Number of the face (0=bottom, 1=right, ...)
    std::vector<Border*> brothers;          // Pointer to brothers Borders
    ibex::IntervalVector position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *pave;                             // Pointer to its container

    bool flow_in;
    bool flow_out[4];

private:
    bool empty;
    bool full;
};

#endif // BORDER_H
