#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>

class Border;
class Pave
{

/***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &box, ibex::Function *f);
    ~Pave(){}

    bool copy_segment(Pave *p);
    bool equal_segment(Pave *p);

    // ******** Drawing functions ********
    void draw(bool filled, string color="b[]");
    void draw_borders(bool filled);

    // ******** Graph building ********
    void bisect(vector<Pave *> &result);
    void remove_from_brothers();
    void remove_brothers(Pave* p, int face);

    // ******** Pave Properties ********
    // Tests
    bool is_empty();
    bool is_full();
    bool is_one_brother_empty(int level=1);
    bool is_all_brothers_full(int level=1);

    // Setter
    void set_empty(bool val);
    bool set_full(bool val);
    void set_full();
    bool set_full_continuity();
    void set_theta(ibex::Interval theta);
    void set_same_properties(Pave *p);

    // Getters
    ibex::IntervalVector get_border_position(int face);
    double get_theta_diam();
    vector<Pave*> get_brothers(int face);

    // ******** Utils functions ********
    std::vector<ibex::Interval> rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y);


/***************** Variables ******************/
public:
    ibex::Interval theta[2];
    ibex::IntervalVector box;
    std::vector<Border> borders;

    ibex::Function *f;

private:
    bool empty;
    bool full;
};

#endif // PAVE_H
