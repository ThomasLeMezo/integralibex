#ifndef INCLUSION_H
#define INCLUSION_H

#include <ibex.h>
#include <border.h>

class Border;
class Inclusion
{
public:
    Inclusion(Border* border, ibex::Function *f, int face, Border *owner);
    Inclusion(Border* border, int brother_face, Border *owner);
    Inclusion(const Inclusion &i);
    Inclusion(Inclusion *i);

    const ibex::Interval            get_segment_in() const;
    const ibex::Interval            get_segment_out() const;
    const ibex::IntervalVector      get_position() const;
    int                             get_brother_face() const;
    Border*                         get_border() const;
    ibex::Function*                 get_function() const;
    bool                            get_shortcut() const;
    Border*                         get_owner();

    void                            set_function(ibex::Function *f);
    void                            set_border(Border* border);

    bool                            is_empty() const;


private:
    Border* m_border;
    Border* m_owner;
    ibex::Function* m_f;
    int m_brother_face; // face of the brother

    bool m_shortcut;

};

#endif // INCLUSION_H
