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
    this->m_position = position;
    this->m_face = face;
    this->m_pave = pave;
    this->m_segment_in = Interval::EMPTY_SET;
    this->m_segment_out = Interval::EMPTY_SET;

    this->m_segment_full = position[face%2];

    this->m_empty = false;
    this->m_full = false;
}

Border::Border(const Border &border): m_position(2)
{
    this->m_position = border.position();
    this->m_face = border.face();
    this->m_pave = border.pave();
    this->m_segment_in = border.segment_in();
    this->m_segment_out = border.segment_out();
    this->m_segment_full = border.segment_full();
    this->m_empty = false;
    this->m_full = false;
}

void Border::draw() const{
    // Create an IntervalVector (2D) from the segment (1D)
    IntervalVector segment_in = this->m_position;
    IntervalVector segment_out = this->m_position;

    // Find the non flat dimension and complete replaced it by the segment
    segment_in[this->m_face%2] = this->m_segment_in;
    segment_out[this->m_face%2] = this->m_segment_out;

    double pourcentage = min(this->m_segment_full.diam()*0.01, 0.01);
    segment_in[(this->m_face+1)%2] += Interval(-pourcentage, pourcentage);

    vibes::drawBox(segment_in, "r[r]");
    vibes::drawBox(segment_out, "b[b]");
}

void Border::get_points(std::vector<double> &x, std::vector<double> &y){
    Interval segment = this->m_segment_in | this->m_segment_out;

    if(!segment.is_empty()){
        if(this->m_face == 0){
            x.push_back(segment.lb());
            x.push_back(segment.ub());
            y.push_back(this->m_position[1].lb());
            y.push_back(this->m_position[1].ub());
        }
        else if(this->m_face == 1){
            x.push_back(this->m_position[0].lb());
            x.push_back(this->m_position[0].ub());
            y.push_back(segment.lb());
            y.push_back(segment.ub());
        }
        else if(this->m_face == 2){
            x.push_back(segment.ub());
            x.push_back(segment.lb());
            y.push_back(this->m_position[1].ub());
            y.push_back(this->m_position[1].lb());
        }
        else if(this->m_face == 3){
            x.push_back(this->m_position[0].ub());
            x.push_back(this->m_position[0].lb());
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
        this->add_brothers(brother_list[i]);
    }
}

void Border::add_brothers(Border * brother){
    IntervalVector test = brother->m_position & this->m_position;
    if(!(test.is_empty()) && (test[0].is_degenerated() != test[1].is_degenerated())){
        this->m_brothers.push_back(brother);
    }
}

void Border::update_brothers(Border* border_pave1, Border* border_pave2){
    vector<Border*> list_tmp;
    list_tmp.push_back(border_pave1);
    list_tmp.push_back(border_pave2);

    for(int i=0; i<this->m_brothers.size(); i++){
        // Remove reference to this border inside brothers
        for(int j=0; j<this->m_brothers[i]->brothers().size(); j++){
            if(this->m_brothers[i]->brothers()[j]==this){
                this->m_brothers[i]->remove_brother(j);
                break;
            }
        }

        // Add reference of border_pave 1 and 2
        this->m_brothers[i]->add_brothers(list_tmp);
    }
}

// ********************************************************************************
// ****************** Other functions *********************************************

void Border::set_full(){
    this->m_segment_in = this->m_segment_full;
    this->m_segment_out = this->m_segment_full;
    this->m_empty = false;
    this->m_full = true;
}

void Border::set_empty(){
    this->m_segment_in = Interval::EMPTY_SET;
    this->m_segment_out = Interval::EMPTY_SET;
    this->m_empty = true;
    this->m_full = false;
}

bool Border::is_empty(){
    if(this->m_empty){
        return true;
    }
    else if(this->m_segment_in.is_empty() && this->m_segment_out.is_empty()){
        this->m_empty = true;
        return true;
    }
    else{
        return false;
    }
}

bool Border::is_full(){
    if(!this->m_full){
        return false;
    }
    else{
        if((this->m_segment_in | this->m_segment_out) != this->m_segment_full){
            this->m_full = false;
            return false;
        }
        else{
            return true;
        }
    }
}

void Border::set_segment_in(ibex::Interval segment_in, bool inclusion){
    if(inclusion)
        this->m_segment_in &= segment_in;
    else
        this->m_segment_in |= segment_in;
}

void Border::set_segment_out(ibex::Interval segment_out, bool inclusion){
    if(inclusion)
        this->m_segment_out &= segment_out;
    else
        this->m_segment_out |= segment_out;
}

void Border::set_pave(Pave* pave){
    this->m_pave = pave;
}

ibex::Interval Border::segment_in() const{
    return this->m_segment_in;
}

ibex::Interval Border::segment_out() const{
    return this->m_segment_out;
}

ibex::Interval Border::segment_full() const{
    return this->m_segment_full;
}

std::vector<Border*> Border::brothers(){
    return this->m_brothers;
}

ibex::IntervalVector Border::position() const{
    return this->m_position;
}

Pave* Border::pave() const{
    return this->m_pave;
}

int Border::face() const{
    return this->m_face;
}

Border& Border::operator&=(const Border &b){
    this->m_segment_in &= b.segment_in();
    this->m_segment_out &= b.segment_out();
    return *this;
}

void Border::diff(const Border &b){
    Interval segment_in_r, segment_in_l;
    this->m_segment_in.diff(b.segment_in(), segment_in_r, segment_in_l);
    this->m_segment_in = segment_in_r | segment_in_l;

    Interval segment_out_r, segment_out_l;
    this->m_segment_out.diff(b.segment_out(), segment_out_r, segment_out_l);
    this->m_segment_out = segment_out_r | segment_out_l;

}

void Border::remove_brother(int indice){
    this->m_brothers.erase(this->m_brothers.begin() + indice);
}

void Border::replace_brother(Border* brother, int id_brother){
    this->m_brothers[id_brother] = brother;
}
