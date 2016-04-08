#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>
#include <string>

class Border;
class Pave
{

    /***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &position, const std::vector<ibex::Function *> &f_list, const ibex::IntervalVector &u);
    Pave(const ibex::IntervalVector &position, ibex::Function* f, const ibex::IntervalVector &u);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    Border*                     operator[](int face);
    bool                        inter(const Pave &p);
    bool                        diff(const Pave &p);

    // ******** Drawing functions ********
    void                        draw(bool filled, std::string color="black[]", bool borders_only=false, bool cmd_u=false);
    void                        draw_borders(bool filled, std::string color_polygon="g[g]");
    void                        draw_position();
    void                        print();

    // ******** Graph building ********
    void                        bisect(std::vector<Pave *> &result);
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

    void                        set_active_function(int id);

    void                        reset_full_empty();

    // Getters
    double                              get_theta_diam();
    const std::vector<Pave*>            get_brothers(int face);
    const ibex::Interval&               get_theta(int i) const;
    const std::vector<ibex::Interval>   get_theta() const;
    std::vector<std::vector<ibex::Interval>>  get_theta_list() const;

    const ibex::Interval&               get_u(int i) const;
    const std::vector<ibex::Interval>   get_u() const;    
    const ibex::IntervalVector&         get_u_iv() const;
    const ibex::IntervalVector&         get_position() const;

    const std::vector<Border *> &       get_borders();
    Border*                             get_border(int face);
    const Border*                       get_border_const(int face) const;

    Pave*                               get_copy_node();
    ibex::Function*                     get_f() const;
    std::vector<ibex::Function *>       get_f_list() const;
    int                                 get_active_function() const;

    bool                                get_first_process() const;
    bool                                get_inner() const;

    // Other functions
    const std::vector<ibex::Interval> compute_theta(ibex::Function *f);



    /***************** Variables ******************/
private:
//    std::vector<ibex::Interval> m_theta;
    std::vector< std::vector<ibex::Interval>> m_theta_list;

    std::vector<ibex::Interval> m_u;
    ibex::IntervalVector        m_u_iv;
    ibex::IntervalVector        m_position;
    std::vector<Border*>        m_borders;

    std::vector<ibex::Function*> m_f_list;
    int                         m_active_function;
    Pave*                       m_copy_node;

private:
    bool                        m_empty;
    bool                        m_full;

    bool                        m_in_queue;

    bool                        m_first_process;
    bool                        m_inner;
};

#endif // PAVE_H
