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
    Pave(const ibex::IntervalVector &position, const std::vector<ibex::Function *> &f_list, const ibex::IntervalVector &u, bool diseable_singeleton=false, bool active=true);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    Pave&                       operator|=(const Pave &p);
    Border*                     operator[](int face);
    bool                        inter(const Pave &p);
    bool                        diff(const Pave &p);

    void                        combine(std::vector<Pave *> &pave_list);
    void                        combine(Pave &p);
    void                        complementaire();
    ibex::IntervalVector        bounding_pave_in() const;
    ibex::IntervalVector        bounding_pave_out() const;
    void                        intersect_face(const ibex::IntervalVector &box_in, const ibex::IntervalVector &box_out);

    // ******** Drawing functions ********
    void                        draw(bool filled, std::string color="black[]", bool borders_only=false) const;
    void                        draw_borders(bool filled, std::string color_polygon="g[g]") const;
    void                        draw_test(int size, std::string comment) const;
    void                        draw_theta() const;
    void                        print();
    void                        print_theta_list();

    // ******** Graph building ********
    void                        bisect(std::vector<Pave *> &result, bool backward=false);
    void                        remove_from_brothers();
    void                        remove_brothers(Pave* p, int face);

    // ******** Pave Properties ********
    // Tests
    bool                        is_empty();
    bool                        is_full();
    bool                        is_full_geometricaly() const;
    bool                        is_fully_full();
    bool                        is_in_queue() const;
    bool                        is_active() const;
    bool                        is_border() const;
    bool                        is_test(int face) const;
    bool                        is_marked_attractor() const;
    bool                        is_external_border() const;
    bool                        is_removed_pave() const;
    bool                        is_inner() const;
    bool                        is_near_inner();

    // Setter
    void                        set_full();
    void                        set_full_in();
    void                        set_full_out();
    void                        set_empty();
    void                        set_segment(bool in, bool out);
    void                        set_active(bool val);
    void                        set_theta(std::vector<ibex::Interval> theta_list);
    void                        set_theta(ibex::Interval theta);
    void                        set_in_queue(bool flag);
    void                        set_copy_node(Pave *p);
    void                        set_first_process_true();
    void                        set_first_process_false();
    void                        set_inner(bool inner);

    void                        set_continuity_in(bool enable);
    void                        set_continuity_out(bool enable);

    void                        set_active_function(int id);

    void                        set_marker_attractor(bool val);
    void                        set_external_border(bool val);
    void                        set_removed_pave(bool val);

    void                        reset_full_empty();

    // Getters
    double                              get_theta_diam(int active_function=-1) const;
    double                              get_theta_diam_min();
    const std::vector<Pave*>            get_brothers(int face);
    const std::vector<Pave *>           get_all_brothers();
    const ibex::Interval&               get_theta(int i) const;
    const std::vector<ibex::Interval>   get_theta() const;
    std::vector<std::vector<ibex::Interval>>  get_theta_list() const;
    const std::vector<ibex::Interval>   get_all_theta() const;

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

    bool                                get_diseable_singelton() const;

    ibex::IntervalVector                get_bounding_pave() const;

    // Other functions
    const std::vector<ibex::Interval>   compute_theta(ibex::Function *f);



    /***************** Variables ******************/
private:
//    std::vector<ibex::Interval> m_theta;
    std::vector< std::vector<ibex::Interval>> m_theta_list;
    std::vector< ibex::Interval> m_theta;

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
    bool                        m_fully_full;

    bool                        m_in_queue;

    bool                        m_first_process;
    bool                        m_inner;

    bool                        m_active;
    bool                        m_diseable_singeleton;
    bool                        m_marker_attractor;
    bool                        m_external_border;
    bool                        m_removed_pave;
};

#endif // PAVE_H
