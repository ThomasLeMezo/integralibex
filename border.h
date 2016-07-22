#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <pave.h>
#include <inclusion.h>

class Pave;
class Inclusion;
class Border
{
/***************** Functions ******************/
public:
    Border(const ibex::IntervalVector& position, const int face, Pave *pave);
    Border(const Border *border);
    ~Border();

    // ******** Drawing functions ********
    void                            draw(bool same_size=false, double offset=0.0, bool test=false, bool two_offset=true) const;

    // ******** Graph building ********
    void                            update_brothers_inclusion(Border *border_pave1, Border *border_pave2);
    void                            remove_inclusion(int indice);
    void                            remove_inclusion(Inclusion *inclusion);
    void                            remove_inclusion_receving(int indice);
    void                            remove_inclusion_receving(Inclusion *inclusion);

    // ******** Border Properties ********
    // Operations
    Border&                         operator&=(const Border &b);
    Border&                         operator|=(const Border &b);
    bool                            inter(const Border &b, bool with_bwd=false);
    bool                            diff(const Border &b);
    void                            complementaire();
    void                            copy_to_inner();
    void                            inter_inner(std::vector<Border*> border_list);

    // Setters

    void                            set_full();
    void                            set_full_all();
    void                            set_full_segment_in();
    void                            set_full_segment_out();
    void                            set_full_inner();
    void                            set_full_outer();
    void                            set_full_inner_in();
    void                            set_full_inner_out();
    void                            set_full_outer_in();
    void                            set_full_outer_out();
    void                            set_empty();
    void                            set_empty_outer();
    void                            set_empty_inner();
    void                            set_empty_inner_in();
    void                            set_empty_inner_out();
    void                            set_empty_outer_in();
    void                            set_empty_outer_out();
    void                            set_segment(bool in, bool out);
    void                            set_segment(ibex::Interval seg, bool inclusion);
    bool                            set_segment_in(ibex::Interval segment_in, bool inclusion);
    bool                            set_segment_out(ibex::Interval segment_out, bool inclusion);
    void                            set_pave(Pave* pave);
    void                            set_continuity_in(bool enable);
    void                            set_continuity_out(bool enable);
    void                            set_inner_mode(bool val);

    void                            set_inclusion(Border *border, int id_brother);
    void                            set_inclusion_receving(Border* border, int id_brother);
    void                            set_compute_inner(bool val);
    void                            reset_full_empty();
    void                            set_zone_propagation(bool val);
    void                            set_zone_function_in(const std::vector<bool> &zone_function);
    void                            set_zone_function_in_fwd(const std::vector<bool> &zone_function);
    void                            set_zone_function_in_bwd(const std::vector<bool> &zone_function);
    void                            set_zone_function_out(const std::vector<bool> &zone_function);
    void                            set_zone_function_out_fwd(const std::vector<bool> &zone_function);
    void                            set_zone_function_out_bwd(const std::vector<bool> &zone_function);

    void                            set_backward_function(bool val);

    void                            add_inclusions(const std::vector<Inclusion *> &inclusion_list);
    bool                            add_inclusion(Inclusion *inclusion);
    void                            add_inclusion_copy(Inclusion *inclusion);
    void                            add_inclusion_receving(Inclusion* inclusion);

    void                            push_back_zone_function_in(bool zone_active);
    void                            push_back_zone_function_out(bool zone_active);

    // Getters
    void                            get_points(std::vector<double> &x, std::vector<double> &y, bool complementary=false) const;
    const ibex::Interval            get_segment_in_union_out() const;
    const ibex::Interval            get_segment_in() const;
    const ibex::IntervalVector      get_segment_out_2D() const;
    const ibex::Interval            get_segment_out() const;
    const ibex::IntervalVector      get_segment_in_2D() const;

    const ibex::Interval            get_segment_in_inner() const;
    const ibex::Interval            get_segment_in_outer() const;
    const ibex::Interval            get_segment_out_inner() const;
    const ibex::Interval            get_segment_out_outer() const;

    const std::vector<Inclusion *>  get_inclusions() const ;
    const std::vector<Inclusion*>&  get_inclusions_receving() const;
    Inclusion *get_inclusion(int i);
    const ibex::IntervalVector &    get_position() const;
    Pave*                           get_pave() const;
    const ibex::Interval&           get_segment_full() const;
    int                             get_face() const;
    int                             size() const;
    Inclusion*                      operator[](int id);

    bool                            get_continuity_in() const;
    bool                            get_continuity_out() const;

    bool                            get_inner_mode() const;
    bool                            get_compute_inner() const;

    bool                            get_zone_propagation() const;
    bool                            get_zone_propagation_fwd() const;
    bool                            get_zone_propagation_bwd() const;
    bool                            reset_computation_zone();
    std::vector<bool>               get_zone_function_in() const;
    std::vector<bool>               get_zone_function_in_fwd() const;
    std::vector<bool>               get_zone_function_in_bwd() const;
    bool                            get_zone_function_in(int function_id) const;
    std::vector<bool>               get_zone_function_out() const;
    std::vector<bool>               get_zone_function_out_fwd() const;
    std::vector<bool>               get_zone_function_out_bwd() const;
    bool                            get_zone_function_out(int function_id) const;
    bool                            get_backward_function() const;

    // Tests
    bool                            is_empty_inner();
    bool                            is_empty_outer();
    bool                            is_full_inner();
    bool                            is_full_outer();
    bool                            is_empty();
    bool                            is_full();

/***************** Variables ******************/
private:
    ibex::Interval                  m_segment_in_outer, m_segment_out_outer;
    ibex::Interval                  m_segment_in_inner, m_segment_out_inner;
    ibex::Interval                  m_segment_full;

private:
    int                             m_face;                     // Number of the face (0=bottom, 1=right, ...)
    std::vector<Inclusion*>         m_inclusions;               // Pointer to brothers Borders
    std::vector<Inclusion*>         m_inclusions_receving;      // Pointer to inclusion that point to this border
    ibex::IntervalVector            m_position;                 // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *                          m_pave;                     // Pointer to its container

private:
    bool                            m_empty_inner;
    bool                            m_empty_outer;
    bool                            m_full_inner;
    bool                            m_full_outer;

private:
    bool                            m_enable_continuity_in, m_enable_continuity_out;
    bool                            m_active_in, m_active_out;

private:
    bool                            m_mode_inner;
    bool                            m_compute_inner;
    bool                            set_segment_in(ibex::Interval segment_in);
    bool                            set_segment_out(ibex::Interval segment_in);

private:
    bool                            m_zone_propagation;
    bool                            m_backward_function;
    std::vector<bool>               m_zone_function_in_fwd;
    std::vector<bool>               m_zone_function_in_bwd;
    std::vector<bool>               m_zone_function_out_fwd;
    std::vector<bool>               m_zone_function_out_bwd;

};

#endif // BORDER_H
