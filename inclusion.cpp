#include "inclusion.h"

using namespace std;
using namespace ibex;

Inclusion::Inclusion(Border *border, ibex::Function *f, int brother_face_axis, int brother_face_side){
    m_shortcut = false;
    m_border = border;
    m_owner = NULL;
    m_f = f;
    m_brother_face_axis = brother_face_axis;
    m_brother_face_side = brother_face_side;
}

Inclusion::Inclusion(Border *border, int brother_face_axis, int brother_face_side){
    m_shortcut = true;
    m_border = border;
    m_owner = NULL;
    m_brother_face_axis = brother_face_axis;
    m_brother_face_side = brother_face_side;
    m_f = NULL;
}

Inclusion::Inclusion(const Inclusion &i){
    m_border = i.get_border();
    m_owner = i.get_owner();
    m_brother_face_axis = i.get_brother_face_axis();
    m_brother_face_side = i.get_brother_face_side();
    m_shortcut = i.get_shortcut();
    m_f = i.get_function();
    if(m_border == NULL)
        cout << "ERROR m_f=NULL" << endl;
}

Inclusion::Inclusion(Inclusion *i){
    m_border = i->get_border();
    m_owner = i->get_owner();
    m_brother_face_axis = i->get_brother_face_axis();
    m_brother_face_side = i->get_brother_face_side();
    m_shortcut = i->get_shortcut();
    m_f = i->get_function();
    if(m_shortcut == false && m_f == NULL)
        cout << "ERROR m_f=NULL" << endl;
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

const PPL::C_Polyhedron Inclusion::get_volume_in() const{
    if(m_shortcut){
        return m_border->get_volume_in();
    }
    else{
        IntervalVector box = ph_2_iv(m_border->get_volume_in());
        IntervalVector box_out = m_f->eval_vector(box);
        PPL::C_Polyhedron ph_out(iv_2_box(box_out));
        return ph_out;
    }
}

const C_Polyhedron Inclusion::get_volume_out() const{
    if(m_shortcut){
        return m_border->get_volume_out();
    }
    else{
        IntervalVector box = ph_2_iv(m_border->get_volume_out());
        IntervalVector box_out = m_f->eval_vector(box);
        PPL::C_Polyhedron ph_out(iv_2_box(box_out));
        return ph_out;
    }
}

Border* Inclusion::get_border() const{
    return m_border;
}

const IntervalVector Inclusion::get_position() const{
    // Do not return a reference :
    if(m_shortcut){
        return m_border->get_position();
    }
    else{
        IntervalVector box = m_border->get_position();
        return m_f->eval_vector(box);
    }
}

int Inclusion::get_brother_face_axis() const{
    return m_brother_face_axis;
}

int Inclusion::get_brother_face_side() const{
    return m_brother_face_side;
}

void Inclusion::set_border(Border* border){
    m_border = border;
}

void Inclusion::set_owner(Border* owner){
    m_owner = owner;
}

bool Inclusion::get_shortcut() const{
    return m_shortcut;
}

bool Inclusion::is_empty() const{
    if(m_border == NULL){
        return true;
    }
    else{
        return false;
    }
}

Border* Inclusion::get_owner() const{
    return m_owner;
}
