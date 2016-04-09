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
    void                            draw() const;

    // ******** Graph building ********
    void                            update_brothers_inclusion(Border *border_pave1, Border *border_pave2);
    void                            remove_inclusion(int indice);
    void                            remove_inclusion(Inclusion *inclusion);
    void                            remove_inclusion_receving(int indice);
    void                            remove_inclusion_receving(Inclusion *inclusion);

    // ******** Border Properties ********
    // Operations
    Border&                         operator&=(const Border &b);
    bool                            inter(const Border &b);
    bool                            diff(const Border &b);

    // Setters
    void                            set_full();
    void                            set_full_segment_in();
    void                            set_full_segment_out();
    void                            set_empty();
    void                            set_segment_in(ibex::Interval segment_in, bool inclusion);
    void                            set_segment_out(ibex::Interval segment_out, bool inclusion);
    void                            set_pave(Pave* pave);
    void                            set_continuity_in(bool enable);
    void                            set_continuity_out(bool enable);
    void                            set_continuity(bool enable);

    void                            set_inclusion(Border *border, int id_brother);
    void                            set_inclusion_receving(Border* border, int id_brother);
    void                            reset_full_empty();

    void                            add_inclusions(const std::vector<Inclusion *> &inclusion_list);
    bool                            add_inclusion(Inclusion *inclusion);
    bool                            add_inclusion_copy(Inclusion *inclusion);
    void                            add_inclusion_receving(Inclusion* inclusion);

    // Getters
    void                            get_points(std::vector<double> &x, std::vector<double> &y) const;
    const ibex::Interval            get_segment_in() const;
    const ibex::IntervalVector      get_segment_out_2D() const;
    const ibex::Interval            get_segment_out() const;
    const ibex::IntervalVector      get_segment_in_2D() const;
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

    // Tests
    bool                            is_empty();
    bool                            is_full();

/***************** Variables ******************/
private:
    ibex::Interval                  m_segment_in, m_segment_out;
    ibex::Interval                  m_segment_full;

private:
    int                             m_face;                     // Number of the face (0=bottom, 1=right, ...)
    std::vector<Inclusion*>         m_inclusions;               // Pointer to brothers Borders
    std::vector<Inclusion*>         m_inclusions_receving;      // Pointer to inclusion that point to this border
    ibex::IntervalVector            m_position;                 // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *                          m_pave;                     // Pointer to its container

private:
    bool                            m_empty;
    bool                            m_full;

private:
    bool                            m_enable_continuity_in, m_enable_continuity_out;
};

#endif // BORDER_H
