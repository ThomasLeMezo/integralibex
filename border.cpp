#include "border.h"
#include "ibex.h"
#include "vibes.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

using namespace ibex;
using namespace std;

Border::Border(const IntervalVector &position, const int face, Pave *pave): m_position(2)
{
    m_position = position;
    m_face = face;
    m_pave = pave;
    m_segment_in = Interval::EMPTY_SET;
    m_segment_out = Interval::EMPTY_SET;

    m_segment_full = position[face%2];

    m_empty = false;
    m_full = false;
}

Border::Border(Border *border): m_position(2)
{
    m_position = border->get_position();
    m_face = border->get_face();
    m_pave = border->get_pave();
    m_segment_in = border->get_segment_in();
    m_segment_out = border->get_segment_out();
    m_segment_full = border->get_segment_full();
    m_empty = false;
    m_full = false;
    m_brothers = border->get_brothers();
}

void Border::draw() const{
    // Create an IntervalVector (2D) from the segment (1D)
    IntervalVector segment_in = m_position;
    IntervalVector segment_out = m_position;

    // Find the non flat dimension and complete replaced it by the segment
    segment_in[m_face%2] = m_segment_in;
    segment_out[m_face%2] = m_segment_out;

    double pourcentage_in = min(m_segment_full.diam()*0.01, 0.01);
    double pourcentage_out = min(m_segment_full.diam()*0.001, 0.001);
    segment_in[(m_face+1)%2] += Interval(-pourcentage_in, pourcentage_in);
    segment_out[(m_face+1)%2] += Interval(-pourcentage_out, pourcentage_out);

    vibes::drawBox(segment_in, "r[r]");
    vibes::drawBox(segment_out, "b[b]");
}

void Border::get_points(std::vector<double> &x, std::vector<double> &y){
    Interval segment = m_segment_in | m_segment_out;

    if(!segment.is_empty()){
        if(m_face == 0){
            x.push_back(segment.lb());
            x.push_back(segment.ub());
            y.push_back(m_position[1].lb());
            y.push_back(m_position[1].ub());
        }
        else if(m_face == 1){
            x.push_back(m_position[0].lb());
            x.push_back(m_position[0].ub());
            y.push_back(segment.lb());
            y.push_back(segment.ub());
        }
        else if(m_face == 2){
            x.push_back(segment.ub());
            x.push_back(segment.lb());
            y.push_back(m_position[1].ub());
            y.push_back(m_position[1].lb());
        }
        else if(m_face == 3){
            x.push_back(m_position[0].ub());
            x.push_back(m_position[0].lb());
            y.push_back(segment.ub());
            y.push_back(segment.lb());
        }
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

// Add new brothers to the list
void Border::add_brothers(std::vector<Border*> brother_list){
    for(int i=0; i<brother_list.size(); i++){
        add_brothers(brother_list[i]);
    }
}

void Border::add_brothers(Border *brother){
    IntervalVector test = brother->m_position & m_position;
    if(!(test.is_empty()) && (test[0].is_degenerated() != test[1].is_degenerated())){
        m_brothers.push_back(brother);
    }
}

void Border::add_brothers_force(Border *brother){
    m_brothers.push_back(brother);
}

void Border::update_brothers(Border* border_pave1, Border* border_pave2){
    vector<Border*> list_tmp;
    list_tmp.push_back(border_pave1);
    list_tmp.push_back(border_pave2);

    for(int i=0; i<m_brothers.size(); i++){
        // Remove reference to this border inside brothers
        for(int j=0; j<m_brothers[i]->get_brothers().size(); j++){
            if(m_brothers[i]->get_brothers()[j]==this){
                m_brothers[i]->remove_brother(j);
                break;
            }
        }

        // Add reference of border_pave 1 and 2
        m_brothers[i]->add_brothers(list_tmp);
    }
}

// ********************************************************************************
// ****************** Other functions *********************************************

void Border::set_full(){
    m_segment_in = m_segment_full;
    m_segment_out = m_segment_full;
    m_empty = false;
    m_full = true;
}

void Border::set_empty(){
    m_segment_in = Interval::EMPTY_SET;
    m_segment_out = Interval::EMPTY_SET;
    m_empty = true;
    m_full = false;
}

bool Border::is_empty(){
    if(m_empty){
        return true;
    }
    else if(m_segment_in.is_empty() && m_segment_out.is_empty()){
        m_empty = true;
        return true;
    }
    else{
        return false;
    }
}

bool Border::is_full(){
    if(!m_full){
        return false;
    }
    else{
        if((m_segment_in | m_segment_out) != m_segment_full){
            m_full = false;
            return false;
        }
        else{
            return true;
        }
    }
}

void Border::set_segment_in(ibex::Interval segment_in, bool inclusion){
    if(inclusion)
        m_segment_in &= segment_in;
    else
        m_segment_in |= segment_in & m_segment_full;
}

void Border::set_segment_out(ibex::Interval segment_out, bool inclusion){
    if(inclusion)
        m_segment_out &= segment_out;
    else
        m_segment_out |= segment_out & m_segment_full;
}

void Border::set_pave(Pave* pave){
    m_pave = pave;
}

ibex::Interval Border::get_segment_in() const{
    return m_segment_in;
}

ibex::Interval Border::get_segment_out() const{
    return m_segment_out;
}

ibex::Interval Border::get_segment_full() const{
    return m_segment_full;
}

std::vector<Border*> Border::get_brothers(){
    return m_brothers;
}

Border* Border::get_brother(int i){
    if(i<m_brothers.size())
        return m_brothers[i];
    else
        return NULL;
}

ibex::IntervalVector Border::get_position() const{
    return m_position;
}

Pave* Border::get_pave() const{
    return m_pave;
}

int Border::get_face() const{
    return m_face;
}

Border& Border::operator&=(const Border &b){
    m_segment_in &= b.get_segment_in();
    m_segment_out &= b.get_segment_out();
    return *this;
}

bool Border::inter(const Border &b){
    bool change = false;
    if((m_segment_in & b.get_segment_in()) != m_segment_in)
        change = true;
    if((m_segment_out & b.get_segment_out()) != m_segment_out)
        change = true;

    m_segment_in &= b.get_segment_in();
    m_segment_out &= b.get_segment_out();
    return change;
}

bool Border::diff(const Border &b){
    bool change = false;
    Interval segment_in_r, segment_in_l;
    m_segment_in.diff(b.get_segment_in(), segment_in_r, segment_in_l);
    if(m_segment_in != (segment_in_r | segment_in_l))
        change = true;
    m_segment_in = segment_in_r | segment_in_l;

    Interval segment_out_r, segment_out_l;
    m_segment_out.diff(b.get_segment_out(), segment_out_r, segment_out_l);
    if(m_segment_out != (segment_out_r | segment_out_l))
        change = true;
    m_segment_out = segment_out_r | segment_out_l;

}

void Border::remove_brother(int indice){
    m_brothers.erase(m_brothers.begin() + indice);
}

void Border::set_brother(Border* brother, int id_brother){
    if(id_brother<m_brothers.size())
        m_brothers[id_brother] = brother;
}

void Border::reset_full_empty(){
    m_empty = false;
    m_full = false;
}
