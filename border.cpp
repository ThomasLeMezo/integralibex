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
    m_segment_in_outer = Interval::EMPTY_SET;
    m_segment_out_outer = Interval::EMPTY_SET;

    m_segment_in_inner = Interval::EMPTY_SET;
    m_segment_out_inner = Interval::EMPTY_SET;

    m_segment_full = position[face%2];

    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;

    m_enable_continuity_in = true;
    m_enable_continuity_out = true;

    m_mode_inner = false;
}

Border::Border(const Border *border): m_position(2)
{
    m_position = border->get_position();
    m_face = border->get_face();
    m_pave = border->get_pave();
    m_segment_in_inner = border->get_segment_in_inner();
    m_segment_in_outer = border->get_segment_in_outer();
    m_segment_out_inner = border->get_segment_out_inner();
    m_segment_out_outer = border->get_segment_out_outer();
    m_segment_full = border->get_segment_full();

    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;

    //    m_inclusions = border->get_inclusions();
    //    m_inclusions_receving = border->get_inclusions_receving();
    m_enable_continuity_in = border->get_continuity_in();
    m_enable_continuity_out = border->get_continuity_out();
    m_mode_inner = border->get_inner_mode();
}

Border::~Border(){
    for(Inclusion *i :m_inclusions){
        delete(i);
    }
}

void Border::draw(bool same_size, double offset, bool test, bool two_offset) const{
    // Create an IntervalVector (2D) from the segment (1D)
    IntervalVector segment_in = get_segment_in_2D();
    IntervalVector segment_out =get_segment_out_2D();

    double ratio = 0.005;
    if(same_size)
        ratio = 0.01;
    double pourcentage_in = min(m_segment_full.diam()*ratio, ratio);
    double pourcentage_out = min(m_segment_full.diam()*0.01, 0.01);

    segment_in[(m_face+1)%2] += Interval(-pourcentage_in, pourcentage_in) + offset;
    if(two_offset)
        segment_out[(m_face+1)%2] += Interval(-pourcentage_out, pourcentage_out) + 2*offset;
    else
        segment_out[(m_face+1)%2] += Interval(-pourcentage_out, pourcentage_out) + offset;

    vibes::drawBox(segment_out, "b[b]");
    vibes::drawBox(segment_in, "r[r]");

    if(test && get_segment_in().is_empty()){
        double center[2];
        center[(m_face+1)%2] = m_position[(m_face+1)%2].mid() + offset;
        center[(m_face)%2] = m_position[(m_face)%2].mid();

        vibes::drawCircle(center[0], center[1], 0.01*m_segment_full.diam(), "black[black]");
    }

    if(test && get_segment_out().is_empty()){
        double center[2];
        center[(m_face+1)%2] = m_position[(m_face+1)%2].mid() + 2*offset;
        center[(m_face)%2] = m_position[(m_face)%2].mid();

        vibes::drawCircle(center[0], center[1], 0.01*m_segment_full.diam(), "black[black]");
    }
}

void Border::get_points(std::vector<double> &x, std::vector<double> &y) const{
    Interval segment = get_segment_in() | get_segment_out();

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

// Add new brothers to the list
void Border::add_inclusions(const std::vector<Inclusion*>& inclusion_list){
    for(Inclusion *i:inclusion_list){
        i->get_border()->remove_inclusion_receving(i);
        add_inclusion_copy(i);
    }
}

bool Border::add_inclusion(Inclusion *inclusion){
    // ToDo : error with inclusion.get_position() if returning a reference !!
    //    if(inclusion->get_owner()->is_empty()) // Test if the border exist
    //        return false;
    IntervalVector test = get_position() & inclusion->get_position();
    if(!(test.is_empty()) && (test[0].is_degenerated() != test[1].is_degenerated())){
        m_inclusions.push_back(inclusion);
        inclusion->set_owner(this);
        inclusion->get_border()->add_inclusion_receving(inclusion); // Add a ref to inclusion (for removal purpose)
        return true;
    }
    else{
        return false;
    }
}

void Border::add_inclusion_copy(Inclusion *inclusion){
    Inclusion *i = new Inclusion(inclusion);
    if(!add_inclusion(i)){
        delete(i);
    }
}

void Border::add_inclusion_receving(Inclusion* inclusion){
    m_inclusions_receving.push_back(inclusion);
}

void Border::update_brothers_inclusion(Border* border_pave1, Border* border_pave2){
    for(Inclusion *i:m_inclusions_receving){
        // Add reference of border_pave 1 and 2
        Inclusion *inclusion_to_pave1 = new Inclusion(i);
        Inclusion *inclusion_to_pave2 = new Inclusion(i);

        inclusion_to_pave1->set_border(border_pave1);
        inclusion_to_pave2->set_border(border_pave2);

        // Add inclusion to pave 1 and pave 2, if no success delete object
        if(!inclusion_to_pave1->get_owner()->add_inclusion(inclusion_to_pave1)){
            delete(inclusion_to_pave1);
        }
        if(!inclusion_to_pave2->get_owner()->add_inclusion(inclusion_to_pave2)){
            delete(inclusion_to_pave2);
        }

        // Remove inclusion to pave
        i->get_owner()->remove_inclusion(i);
    }
}

// ********************************************************************************
// ****************** Other functions *********************************************

void Border::set_full(){
    set_segment_in(m_segment_full);
    set_segment_out(m_segment_full);

    if(m_mode_inner){
        m_empty_inner = false;
        m_full_inner = true;
    }
    else{
        m_empty_outer = false;
        m_full_outer = true;
    }
}

void Border::set_full_segment_in(){
    set_segment_in(m_segment_full);
    if(m_mode_inner)
        m_empty_inner = false;
    else
        m_empty_outer = false;
}

void Border::set_full_segment_out(){
    set_segment_out(m_segment_full);
    if(m_mode_inner)
        m_full_inner = false;
    else
        m_full_outer = false;
}

void Border::set_empty(){
    set_segment_in(Interval::EMPTY_SET);
    set_segment_out(Interval::EMPTY_SET);
    if(m_mode_inner){
        m_empty_inner = true;
        m_full_inner = false;
    }
    else{
        m_empty_outer = true;
        m_full_outer = false;
    }
}

bool Border::is_empty_inner(){
    if(m_empty_inner)
        return true;
    else if(get_segment_in_inner().is_empty() && get_segment_out_inner().is_empty()){
        m_empty_inner = true;
        return true;
    }
    else
        return false;
}

bool Border::is_empty_outer(){
    if(m_empty_outer)
        return true;
    else if(get_segment_in_outer().is_empty() && get_segment_out_outer().is_empty()){
        m_empty_outer = true;
        return true;
    }
    else
        return false;
}

bool Border::is_empty(){
    if(m_compute_inner)
        return is_empty_inner() && is_empty_outer();
    else
        return is_empty_outer();
}

bool Border::is_full_inner(){
    if(!m_full_inner)
        return false;
    else{
        if((get_segment_in_inner() | get_segment_out_inner()) != m_segment_full){
            m_full_inner = false;
            return false;
        }
        else
            return true;
    }
}

bool Border::is_full_outer(){
    if(!m_full_inner)
        return false;
    else{
        if((get_segment_in_outer() | get_segment_out_outer()) != m_segment_full){
            m_full_outer = false;
            return false;
        }
        else
            return true;
    }
}

bool Border::is_full(){
    if(m_compute_inner)
        return is_full_inner() && is_full_outer();
    else
        return is_full_outer();
}

bool Border::set_segment_in(ibex::Interval segment_in, bool inclusion){
    bool change = false;
    Interval seg_in(Interval::EMPTY_SET);
    if(inclusion)
        seg_in = get_segment_in() & segment_in & m_segment_full;
    else
        seg_in = (get_segment_in() | segment_in) & m_segment_full;

    if(seg_in != get_segment_in()){
        set_segment_in(seg_in);
        change = true;
    }

    if(m_pave->get_diseable_singelton() && get_segment_in().is_degenerated())
        set_segment_in(Interval::EMPTY_SET);

    return change;
}

bool Border::set_segment_out(ibex::Interval segment_out, bool inclusion){
    bool change = false;
    Interval i(Interval::EMPTY_SET);
    if(inclusion)
        i = get_segment_out() & segment_out & m_segment_full;
    else
        i = (get_segment_out() | segment_out) & m_segment_full;

    if(i != get_segment_out()){
        set_segment_out(i);
        change = true;
    }

    if(m_pave->get_diseable_singelton() && get_segment_out().is_degenerated())
        set_segment_out(Interval::EMPTY_SET);

    return change;
}

bool Border::set_segment_in(ibex::Interval segment_in){
    if(!m_mode_inner)
        m_segment_in_outer = segment_in & m_segment_full;
    else
        m_segment_in_inner = segment_in & m_segment_full;
}

bool Border::set_segment_out(ibex::Interval segment_out){
    if(!m_mode_inner)
        m_segment_out_outer = segment_out & m_segment_full;
    else
        m_segment_out_inner = segment_out & m_segment_full;
}

void Border::set_pave(Pave* pave){
    m_pave = pave;
}

const ibex::Interval Border::get_segment_in() const{
    if(!m_mode_inner)
        return get_segment_in_outer();
    else
        return get_segment_in_inner();
}

const ibex::Interval Border::get_segment_out() const{
    if(!m_mode_inner)
        return get_segment_out_outer();
    else
        return get_segment_out_inner();
}

const ibex::Interval Border::get_segment_in_inner() const{
    return m_segment_in_inner;
}

const ibex::Interval Border::get_segment_in_outer() const{
    return m_segment_in_outer;
}

const ibex::Interval Border::get_segment_out_inner() const{
    return m_segment_out_inner;
}

const ibex::Interval Border::get_segment_out_outer() const{
    return m_segment_out_outer;
}

const ibex::Interval& Border::get_segment_full() const{
    return m_segment_full;
}

const ibex::IntervalVector Border::get_segment_in_2D() const{
    IntervalVector segment_in = m_position;
    segment_in[m_face%2] = get_segment_in();
    return segment_in;
}

const ibex::IntervalVector Border::get_segment_out_2D() const{
    IntervalVector segment_out = m_position;
    segment_out[m_face%2] = get_segment_out();
    return segment_out;
}

const std::vector<Inclusion *> Border::get_inclusions() const{
    return m_inclusions;
}

const std::vector<Inclusion*>& Border::get_inclusions_receving() const{
    return m_inclusions_receving;
}

Inclusion* Border::get_inclusion(int i){
    return m_inclusions[i];
}

const IntervalVector& Border::get_position() const{
    return m_position;
}

Pave* Border::get_pave() const{
    return m_pave;
}

int Border::get_face() const{
    return m_face;
}

Border& Border::operator&=(const Border &b){
    set_segment_in(b.get_segment_in(), true);
    set_segment_out(b.get_segment_out(), true);
    return *this;
}

Border& Border::operator|=(const Border &b){
    set_segment_in(b.get_segment_in(), false);
    set_segment_out(b.get_segment_out(), false);
    return *this;
}

bool Border::inter(const Border &b){
    bool change = false;
    if((get_segment_in() & b.get_segment_in()) != get_segment_in())
        change = true;
    if((get_segment_out() & b.get_segment_out()) != get_segment_out())
        change = true;

    set_segment_in(b.get_segment_in(), true);
    set_segment_out(b.get_segment_out(), true);
    return change;
}

bool Border::diff(const Border &b){
    bool change = false;
    Interval segment_in_r, segment_in_l;
    get_segment_in().diff(b.get_segment_in(), segment_in_r, segment_in_l);

    if(get_segment_in() != (segment_in_r | segment_in_l))
        change = true;
    set_segment_in(segment_in_r | segment_in_l);

    Interval segment_out_r, segment_out_l;
    get_segment_out().diff(b.get_segment_out(), segment_out_r, segment_out_l);
    if(get_segment_out() != (segment_out_r | segment_out_l))
        change = true;
    set_segment_out(segment_out_r | segment_out_l);
}

void Border::remove_inclusion(int indice){
    delete(m_inclusions[indice]);
    m_inclusions.erase(m_inclusions.begin() + indice);
}

void Border::remove_inclusion(Inclusion *inclusion){
    for(int i=0; i<m_inclusions.size(); i++){
        if(m_inclusions[i] == inclusion){
            remove_inclusion(i);
            break;
        }
    }
}

void Border::remove_inclusion_receving(int indice){
    m_inclusions_receving.erase(m_inclusions_receving.begin() + indice);
}

void Border::remove_inclusion_receving(Inclusion *inclusion){
    for(int i=0; i<m_inclusions_receving.size(); i++){
        if(m_inclusions_receving[i] == inclusion){
            remove_inclusion_receving(i);
            break;
        }
    }
}

void Border::set_inclusion(Border* border, int id_brother){
    if(id_brother<m_inclusions.size())
        m_inclusions[id_brother]->set_border(border);
}

void Border::set_inclusion_receving(Border* border, int id_brother){
    if(id_brother<m_inclusions_receving.size())
        m_inclusions_receving[id_brother]->set_owner(border);
}

void Border::reset_full_empty(){
    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;
}

int Border::size() const{
    return m_inclusions.size();
}

Inclusion* Border::operator[](int id){
    return m_inclusions[id];
}

bool Border::get_continuity_in() const{
    return m_enable_continuity_in;
}

bool Border::get_continuity_out() const{
    return m_enable_continuity_out;
}

void Border::set_continuity_in(bool enable){
    m_enable_continuity_in = enable;
}

void Border::set_continuity_out(bool enable){
    m_enable_continuity_out = enable;
}

void Border::complementaire(){
    set_segment_in(get_segment_out(), false);
    set_segment_out(get_segment_in(), false);

    Interval seg_in1, seg_in2;
    m_segment_full.diff(get_segment_in(), seg_in1,seg_in2);
    set_segment_in(seg_in1 | seg_in2);

    Interval seg_out1, seg_out2;
    m_segment_full.diff(get_segment_out(), seg_out1,seg_out2);
    set_segment_out(seg_out1 | seg_out2);
}

void Border::set_segment(bool in, bool out){
    if(in)
        set_segment_in(m_segment_full);
    else
        set_segment_in(Interval::EMPTY_SET);

    if(out)
        set_segment_out(m_segment_full);
    else
        set_segment_out(Interval::EMPTY_SET);
}

void Border::set_inner_mode(bool val){
    m_mode_inner = val;
}

bool Border::get_inner_mode() const{
    return m_mode_inner;
}

void Border::set_compute_inner(bool val){
    m_compute_inner = val;
}


