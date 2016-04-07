#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>

class Border;
class Pave
{

    /***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &position, ibex::Function *f, const ibex::IntervalVector &u);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    bool                        inter(const Pave &p);
    bool                        diff(const Pave &p);

    // ******** Drawing functions ********
    void                        draw(bool filled, string color="black[]", bool borders_only=false, bool cmd_u=false);
    void                        draw_borders(bool filled, string color_polygon="g[g]");
    void                        draw_position();
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
    void                        set_first_process_true();
    void                        set_first_process_false();
    void                        set_inner(bool inner);

    void                        reset_full_empty();

    // Getters
    double                              get_theta_diam();
    const std::vector<Pave*>            get_brothers(int face);
    const ibex::Interval&               get_theta(int i) const;
    const std::vector<ibex::Interval>   get_theta() const;
    const ibex::Interval&               get_u(int i) const;
    const std::vector<ibex::Interval>   get_u() const;
    const ibex::IntervalVector&         get_u_iv() const;
    const ibex::IntervalVector&         get_position() const;

    const std::vector<Border *> &       get_borders();
    Border*                             get_border(int face);
    const Border*                       get_border_const(int face) const;

    Pave*                               get_copy_node();
    ibex::Function*                     get_f() const;

    bool                                get_first_process() const;
    bool                                get_inner() const;

    Border* operator[](int face);

    /***************** Variables ******************/
private:
    std::vector<ibex::Interval> m_theta;
    std::vector<ibex::Interval> m_u;
    ibex::IntervalVector        m_u_iv;
    ibex::IntervalVector        m_position;
    std::vector<Border*>        m_borders;

    ibex::Function              *m_f;
    Pave*                       m_copy_node;

private:
    bool                        m_empty;
    bool                        m_full;

    bool                        m_in_queue;

    bool                        m_first_process;
    bool                        m_inner;
};

#endif // PAVE_H
