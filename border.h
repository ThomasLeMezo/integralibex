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
    Border(Border *border);
    ~Border(){}

    // ******** Drawing functions ********
    void draw() const;

    // ******** Graph building ********
    void add_to_brothers_inclusion(Border *border_pave1, Border *border_pave2);
    void remove_inclusion(int indice);

    // ******** Border Properties ********
    // Operations
    Border& operator&=(const Border &b);
    bool inter(const Border &b);
    bool diff(const Border &b);

    // Setters
    void set_full();
    void set_full_segment_in();
    void set_full_segment_out();
    void set_empty();
    void set_segment_in(ibex::Interval segment_in, bool inclusion);
    void set_segment_out(ibex::Interval segment_out, bool inclusion);
    void set_pave(Pave* pave);

    void set_inclusion(Border *border, int id_brother);
    void reset_full_empty();

    void add_inclusions(const std::vector<Inclusion> &inclusion_list);
    void add_inclusion(const Inclusion &inclusion);

    // Getters
    void                            get_points(std::vector<double> &x, std::vector<double> &y) const;
    const ibex::Interval            get_segment_in() const;
    const ibex::IntervalVector      get_segment_out_2D() const;
    const ibex::Interval            get_segment_out() const;
    const ibex::IntervalVector      get_segment_in_2D() const;
    const std::vector<Inclusion>&   get_inclusions() const ;
    Inclusion&                      get_inclusion(int i);
    const ibex::IntervalVector &    get_position() const;
    Pave*                           get_pave() const;
    const ibex::Interval&           get_segment_full() const;
    int                             get_face() const;
    int size() const;
    Inclusion& operator[](int id);

    // Tests
    bool is_empty();
    bool is_full();

/***************** Variables ******************/
private:
    ibex::Interval m_segment_in, m_segment_out;
    ibex::Interval m_segment_full;

private:
    int m_face;                               // Number of the face (0=bottom, 1=right, ...)
    std::vector<Inclusion> m_inclusions;          // Pointer to brothers Borders
    ibex::IntervalVector m_position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *m_pave;                             // Pointer to its container

private:
    bool m_empty;
    bool m_full;
};

#endif // BORDER_H
