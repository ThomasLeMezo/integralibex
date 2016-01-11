#ifndef INCLUSION_H
#define INCLUSION_H

#include <ibex.h>
#include <border.h>

class Border;
class Inclusion
{
public:
    Inclusion(Border* border, ibex::Function *f, int face);

    ibex::Interval          get_segment_in();
    ibex::Interval          get_segment_out();
    ibex::IntervalVector    get_position() const;
    int                     get_brother_face();

    Border* get_border();

    void                    set_border(Border* border);


private:
    Border* m_border;
    ibex::Function* m_f;
    int m_brother_face; // face of the brother

};

#endif // INCLUSION_H
