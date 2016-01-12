#include "inclusion.h"

using namespace std;
using namespace ibex;

Inclusion::Inclusion(Border *border, ibex::Function *f, int brother_face){
    m_border = border;
    m_f = f;
    m_brother_face = brother_face;
    m_shortcut = false;
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
    if(!i.get_shortcut()){
        m_f = i.get_function();
        m_shortcut = false;
    }
    else{
        m_shortcut = true;
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

ibex::Interval Inclusion::get_segment_in(){
    if(m_shortcut){
        return m_border->get_segment_in();
    }
    else{
        return m_border->get_segment_in();
    }
}

ibex::Interval Inclusion::get_segment_out(){
    if(m_shortcut){
        return m_border->get_segment_out();
    }
    else{
        return m_border->get_segment_out();
    }
}

Border* Inclusion::get_border() const{
    return m_border;
}

ibex::IntervalVector Inclusion::get_position() const{
    if(m_shortcut){
        return m_border->get_position();
    }
    else{
        return m_border->get_position();
    }
    // To Do
    //return m_f->eval(get_border()->get_position());
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

