#include "inclusion.h"

using namespace std;
using namespace ibex;

Inclusion::Inclusion(Border *border, ibex::Function *f, int brother_face){
    m_shortcut = false;
    m_border = border;
    m_f = f;
    m_brother_face = brother_face;
}

Inclusion::Inclusion(Border *border, int brother_face){
    m_shortcut = true;
    m_border = border;
    m_brother_face = brother_face;
    m_f = NULL;
}

Inclusion::Inclusion(const Inclusion &i){
    m_border = i.get_border();
    m_brother_face = i.get_brother_face();
    m_shortcut = i.get_shortcut();
    if(!i.get_shortcut()){
        m_f = i.get_function();
    }
}

void Inclusion::set_function(ibex::Function *f){
    if(f!=NULL){
        m_f = f;
        m_shortcut=false;
    }
}

ibex::Function* Inclusion::get_function() const{
    return m_f;
}

const Interval &Inclusion::get_segment_in() const{
    if(m_shortcut){
        return m_border->get_segment_in();
    }
    else{
        IntervalVector box = m_border->get_position();
        box[m_brother_face%2] = m_border->get_segment_in();
        IntervalVector box_out = m_f->eval_vector(box);
        return box_out[m_brother_face%2];
    }
}

const Interval &Inclusion::get_segment_out() const{
    if(m_shortcut){
        return m_border->get_segment_out();
    }
    else{
        IntervalVector box = m_border->get_position();
        box[m_brother_face%2] = m_border->get_segment_out();
        IntervalVector box_out(m_f->eval_vector(box));
        return box_out[m_brother_face%2];
    }
}

Border* Inclusion::get_border() const{
    return m_border;
}

const IntervalVector &Inclusion::get_position() const{
    if(m_shortcut){
        return m_border->get_position();
    }
    else{
        IntervalVector box = m_border->get_position();
        return m_f->eval_vector(box);
    }
}

int Inclusion::get_brother_face() const{
    return m_brother_face;
}

void Inclusion::set_border(Border* border){
    m_border = border;
}

bool Inclusion::get_shortcut() const{
    return m_shortcut;
}

