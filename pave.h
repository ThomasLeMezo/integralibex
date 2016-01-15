#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>

class Border;
class Pave
{

    /***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &position, ibex::Function *f, ibex::Interval u=ibex::Interval::EMPTY_SET);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    bool                        inter(const Pave &p);
    bool                        diff(const Pave &p);

    // ******** Drawing functions ********
    void                        draw(bool filled, string color="black[]", bool inner=false);
    void                        draw_borders(bool filled, string color_polygon="y[y]");
    void                        print();

    // ******** Graph building ********
    void                        bisect(vector<Pave *> &result);
    void                        remove_from_brothers();
    void                        remove_brothers(Pave* p, int face);

    // ******** Pave Properties ********
    // Tests
    bool                        is_empty();
    bool                        is_full();
    bool                        is_in_queue() const;

    // Setter
    void                        set_full();
    void                        set_empty();
    void                        set_theta(ibex::Interval theta);
    void                        set_in_queue(bool flag);
    void                        set_copy_node(Pave *p);

    void                        reset_full_empty();

    // Getters
    const ibex::IntervalVector&         get_border_position(int face) const;
    double                              get_theta_diam();
    const std::vector<Pave*>            get_brothers(int face);
    const ibex::Interval&               get_theta(int i) const;
    const std::vector<ibex::Interval>   get_theta() const;
    const ibex::Interval&               get_u() const;
    const ibex::IntervalVector&         get_position() const;

    const std::vector<Border>&          get_borders();
    Border*                             get_border(int face);
    const Border*                       get_border_const(int face) const;

    Pave*                               get_copy_node();
    ibex::Function*                     get_f() const;

    Border& operator[](int face);

    /***************** Variables ******************/
private:
    std::vector<ibex::Interval> m_theta;
    ibex::Interval              m_u;
    ibex::IntervalVector        m_position;
    std::vector<Border>         m_borders;

    ibex::Function              *m_f;
    Pave*                       m_copy_node;

private:
    bool                        m_empty;
    bool                        m_full;

    bool                        m_in_queue;
};

#endif // PAVE_H
