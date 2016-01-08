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
    Pave(const Pave *p);
    ~Pave();

    Pave& operator&=(const Pave &p);
    bool inter(const Pave &p);
    void diff(const Pave &p);

    // ******** Drawing functions ********
    void draw(bool filled, string color="black[]");
    void draw_borders(bool filled);

    // ******** Graph building ********
    void bisect(vector<Pave *> &result);
    void remove_from_brothers();
    void remove_brothers(Pave* p, int face);

    // ******** Pave Properties ********
    // Tests
    bool is_empty();
    bool is_full();

    // Setter
    void set_full();
    void set_empty();
    void set_theta(ibex::Interval theta);
    void reset_full_empty();

    // Getters
    ibex::IntervalVector get_border_position(int face);
    double get_theta_diam();
    vector<Pave*> get_brothers(int face);

    // ******** Utils functions ********
    std::vector<ibex::Interval> rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y);

    // ******** Tarjan functions ********
//    void tarjan_compute_successors();
//    void strongconnect(int &index, std::vector<Pave *> *S, std::vector<std::vector<Pave *> > *SCC);

/***************** Variables ******************/
public:
    ibex::Interval m_theta[2];
    ibex::IntervalVector m_box;
    std::vector<Border> m_borders;
    std::vector<Border> m_borders_symetry;
    bool m_in_queue;

    ibex::Function *m_f;

    Pave* m_copy_node;

private:
    bool m_empty;
    bool m_full;
};

#endif // PAVE_H
