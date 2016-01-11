#include "inclusion.h"

using namespace std;
using namespace ibex;

Inclusion::Inclusion(Border *border, ibex::Function *f, int brother_face){
    m_border = border;
    m_f = f;
    m_brother_face = brother_face;
}

ibex::Interval Inclusion::get_segment_in(){
    // ToDo
    //return m_f->eval(IntervalVector(1, m_border->get_segment_in()));
}

ibex::Interval Inclusion::get_segment_out(){
    // ToDo
    //return m_f->eval(IntervalVector(1, m_border->get_segment_out()));
}

Border* Inclusion::get_border(){
    return m_border;
}

ibex::IntervalVector Inclusion::get_position() const{
    // To Do
    //return m_f->eval(get_border()->get_position());
}

int Inclusion::get_brother_face(){
    return m_brother_face;
}

void Inclusion::set_border(Border* border){
    m_border = border;
}

