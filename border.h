#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <pave.h>

class Pave;
class Border
{
/***************** Functions ******************/
public:
    Border(const ibex::IntervalVector& position, const int face, Pave *pave);
    Border(const ibex::IntervalVector& position, const int face, const ibex::Interval &segment);
    ~Border(){}

    // ******** Drawing functions ********
    void draw() const;

    // ******** Graph building ********
    void add_brothers(std::vector<Border *> brother_list);
    void update_brothers(Border* border_pave1, Border* border_pave2);

    // ******** Border Properties ********
    // Setters
    void set_full();
    bool set_full_continuity();

    // Getters
    void get_points(std::vector<double> &x, std::vector<double> &y);

    // Tests
    bool is_empty();
    bool is_full();

/***************** Variables ******************/
public:
    ibex::Interval segment;
    ibex::Interval segment_full;

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
