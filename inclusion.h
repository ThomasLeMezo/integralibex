#ifndef INCLUSION_H
#define INCLUSION_H

#include <ibex.h>
#include <border.h>

class Border;
class Inclusion
{
public:
    Inclusion(Border* border, ibex::Function *f, int brother_face_axis, int brother_face_side);
    Inclusion(Border* border, int brother_face_axis, int brother_face_side);
    Inclusion(const Inclusion &i);
    Inclusion(Inclusion *i);

    const PPL::C_Polyhedron         get_volume_in() const;
    const PPL::C_Polyhedron         get_volume_out() const;
    const ibex::IntervalVector      get_position() const;
    int                             get_brother_face_axis() const;
    int                             get_brother_face_side() const;

    Border*                         get_border() const;
    ibex::Function*                 get_function() const;
    bool                            get_shortcut() const;
    Border*                         get_owner() const;

    void                            set_function(ibex::Function *f);
    void                            set_border(Border* border);
    void                            set_owner(Border* owner);

    bool                            is_empty() const;


private:
    Border* m_border;
    Border* m_owner;
    ibex::Function* m_f;
    int m_brother_face_axis; // face of the brother
    int m_brother_face_side; // face of the brother

    bool m_shortcut;

};

#endif // INCLUSION_H
