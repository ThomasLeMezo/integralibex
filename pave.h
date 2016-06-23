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
    Pave(const ibex::IntervalVector &position, const std::vector<ibex::Function *> &f_list, bool diseable_singeleton=false, bool active=true);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    Pave&                       operator|=(const Pave &p);
    Border*                     operator[](int face);
    bool                        inter(const Pave &p, bool with_bwd=false);
    bool                        diff(const Pave &p);
    bool                        inter_inner(const std::vector<Pave*> pave_list);

    void                        complementaire();
    void                        copy_to_inner();

    // ******** Drawing functions ********
    void                        draw(bool filled, bool inner_only=false);
    void                        draw_borders(bool filled, std::string color_polygon="g[g]", bool complementary=false) const;
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
    bool                        is_empty_inter();
    bool                        is_full();
    bool                        is_full_inter();
    bool                        is_full_geometricaly() const;
    bool                        is_fully_full();
    bool                        is_in_queue() const;
    bool                        is_in_queue_outer() const;
    bool                        is_in_queue_inner() const;
    bool                        is_active() const;
    bool                        is_border() const;
    bool                        is_marked_attractor() const;
    bool                        is_external_border() const;
    bool                        is_removed_pave() const;
    bool                        is_removed_pave_inner() const;
    bool                        is_removed_pave_outer() const;
    bool                        is_removed_pave_union() const;
    bool                        is_near_empty();
    bool                        is_empty_inner();
    bool                        is_empty_outer();
    bool                        is_full_inner();
    bool                        is_full_outer();

    // Setter
    void                        set_full();
    void                        set_full_all();
    void                        set_full_in();
    void                        set_full_out();
    void                        set_empty();
    void                        set_empty_outer();
    void                        set_full_inner();
    void                        set_empty_inner();
    void                        set_empty_inner_in();
    void                        set_empty_inner_out();
    void                        set_full_outer();
    void                        set_segment(bool in, bool out);
    void                        set_active(bool val);
    void                        set_theta(std::vector<ibex::Interval> theta_list);
    void                        set_theta(ibex::Interval theta);
    void                        set_in_queue(bool flag);
    void                        set_in_queue_inner(bool flag);
    void                        set_in_queue_outer(bool flag);
    void                        set_copy_node(Pave *p);
    void                        set_first_process_true();
    void                        set_first_process_false();

    void                        set_continuity_in(bool enable);
    void                        set_continuity_out(bool enable);

    void                        set_active_function(int id);

    void                        set_marker_attractor(bool val);
    void                        set_external_border(bool val);
    void                        set_removed_pave(bool val);
    void                        set_removed_pave_inner(bool val);
    void                        set_removed_pave_outer(bool val);

    void                        set_compute_inner(bool val);
    void                        set_inner_mode(bool val);

    void                        set_zone_propagation(bool val);

    void                        set_backward_function(bool val);

    void                        reset_full_empty();

    // Getters
    double                              get_theta_diam(int active_function=-1) const;
    double                              get_theta_diam_min();
    double                              get_theta_diam_max() const;
    const std::vector<Pave*>            get_brothers(int face);
    const std::vector<Pave *>           get_all_brothers();
    const ibex::Interval&               get_theta(int i) const;
    const std::vector<ibex::Interval>   get_theta() const;

    std::vector<std::vector<ibex::Interval>>    get_theta_list() const;
    std::vector<ibex::Interval>                 get_theta_list(int function_id) const;
    const std::vector<ibex::Interval>           get_all_theta(bool all=false) const;

    std::vector<std::vector<ibex::Interval>>    get_theta_list_bwd() const;
    std::vector<std::vector<ibex::Interval>>    get_theta_list_fwd() const;
    std::vector<ibex::Interval>                 get_all_theta_bwd() const;
    std::vector<ibex::Interval>                 get_all_theta_fwd() const;

    const ibex::IntervalVector&         get_position() const;

    std::vector<Border *> &             get_borders();
    const std::vector<Border *>         get_borders_const() const;
    Border*                             get_border(int face);
    const Border*                       get_border_const(int face) const;

    Pave*                               get_copy_node();
    ibex::Function*                     get_f() const;
    std::vector<ibex::Function *>       get_f_list() const;
    int                                 get_active_function() const;

    bool                                get_first_process() const;

    bool                                get_diseable_singelton() const;

    ibex::IntervalVector                get_bounding_pave() const;
    bool                                get_compute_inner() const;
    bool                                get_inner_mode() const;

    bool                                get_zone_propagation() const;
    void                                reset_computation_zone();
    bool                                get_backward_function() const;

    // Other functions
    const std::vector<ibex::Interval>   compute_theta(ibex::Function *f, bool backward_function=false);



    /***************** Variables ******************/
private:
//    std::vector<ibex::Interval> m_theta;
    std::vector< std::vector<ibex::Interval>>   m_theta_list;
    std::vector< ibex::Interval>                m_theta;

    std::vector< std::vector<ibex::Interval>>   m_theta_list_bwd;
    std::vector< ibex::Interval>                m_theta_bwd;

    bool                        m_backward_function;


    ibex::IntervalVector        m_position;
    std::vector<Border*>        m_borders;

    std::vector<ibex::Function*> m_f_list;
    int                         m_active_function;
    Pave*                       m_copy_node;

private:
    bool                        m_empty_outer;
    bool                        m_empty_inner;
    bool                        m_full_outer;
    bool                        m_full_inner;

    bool                        m_in_queue_inner;
    bool                        m_in_queue_outer;

    bool                        m_first_process;

    bool                        m_active;
    bool                        m_diseable_singeleton;
    bool                        m_marker_attractor;
    bool                        m_external_border;
    bool                        m_removed_pave_inner;
    bool                        m_removed_pave_outer;

    bool                        m_compute_inner;
    bool                        m_inner_mode;

    bool                        m_zone_propagation;

};

#endif // PAVE_H
