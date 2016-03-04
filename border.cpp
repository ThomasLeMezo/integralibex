#include "border.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

#include <ppl.hh>

#include "conversion.h"
#include "vtkSmartPointer.h"

using namespace ibex;
using namespace std;
using namespace Parma_Polyhedra_Library::IO_Operators;

Border::Border(const IntervalVector &position, Pave *pave, int face_axis, int face_side):
    m_position(position.size()),
    m_volume_in(position.size(), PPL::EMPTY),
    m_volume_out(position.size(), PPL::EMPTY),
    m_volume_full(position.size(), PPL::EMPTY)
{
    m_position = position;
    m_pave = pave;

    m_dim = position.size();
    m_volume_full = C_Polyhedron(iv_2_box(position));

    m_empty = false;
    m_full = false;

    m_face_axis = face_axis;
    m_face_side = face_side;
}

Border::Border(const Border *border):
    m_position(border->get_dim()),
    m_volume_in(border->get_dim(), PPL::EMPTY),
    m_volume_out(border->get_dim(), PPL::EMPTY),
    m_volume_full(border->get_dim(), PPL::EMPTY)
{
    m_position = border->get_position();
    m_pave = border->get_pave();

    m_volume_in = border->get_volume_in();
    m_volume_out = border->get_volume_out();
    m_volume_full = border->get_volume_full();

    m_empty = false;
    m_full = false;

    m_face_axis = border->get_face_axis();
    m_face_side = border->get_face_side();
    //    m_inclusions = border->get_inclusions();
    //    m_inclusions_receving = border->get_inclusions_receving();
}

Border::~Border(){
    for(int i=0; i<m_inclusions.size(); i++){
        delete(m_inclusions[i]);
    }
}

// Add new brothers to the list
void Border::add_inclusions(const std::vector<Inclusion*>& inclusion_list){
    for(int i=0; i<inclusion_list.size(); i++){
        inclusion_list[i]->get_border()->remove_inclusion_receving(inclusion_list[i]);
        add_inclusion_copy(inclusion_list[i]);
    }
}

bool Border::add_inclusion(Inclusion *inclusion){
    // ToDo : error with inclusion.get_position() if returning a reference !!
    //    if(inclusion->get_owner()->is_empty()) // Test if the border exist
    //        return false;
    IntervalVector test = m_position & inclusion->get_position();
    if(!test.is_empty()){
        int dim_degenerated = 0;
        for(int dim=0; dim<test.size(); dim++){
            if(test[dim].is_degenerated())
                dim_degenerated++;
        }
        if(dim_degenerated!=1)
            return false;

        m_inclusions.push_back(inclusion);
        inclusion->set_owner(this);
        inclusion->get_border()->add_inclusion_receving(inclusion); // Add a ref to inclusion (for removal purpose)
        return true;
    }
    return false;
}

bool Border::add_inclusion_copy(Inclusion *inclusion){
    Inclusion *i = new Inclusion(inclusion);
    if(!add_inclusion(i)){
        delete(i);
    }
}

void Border::add_inclusion_receving(Inclusion* inclusion){
    m_inclusions_receving.push_back(inclusion);
}

void Border::update_brothers_inclusion(Border* border_pave1, Border* border_pave2){
    for(int inclusion_id = 0; inclusion_id<m_inclusions_receving.size(); inclusion_id++){
        // Add reference of border_pave 1 and 2
        Inclusion *inclusion_to_pave1 = new Inclusion(m_inclusions_receving[inclusion_id]);
        Inclusion *inclusion_to_pave2 = new Inclusion(m_inclusions_receving[inclusion_id]);

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
        m_inclusions_receving[inclusion_id]->get_owner()->remove_inclusion(m_inclusions_receving[inclusion_id]);
    }
}

// ********************************************************************************
// ****************** Other functions *********************************************

void Border::set_full(){
    m_volume_in = m_volume_full;
    m_volume_out = m_volume_full;

    m_empty = false;
    m_full = true;
}

void Border::set_full_volume_in(){
    m_volume_in = m_volume_full;
    m_empty = false;
}

void Border::set_full_volume_out(){
    m_volume_out = m_volume_full;
    m_empty = false;
}

void Border::set_empty(){
    m_volume_in = PPL::C_Polyhedron(m_dim, PPL::EMPTY);
    m_volume_out = PPL::C_Polyhedron(m_dim, PPL::EMPTY);
    m_empty = true;
    m_full = false;
}

bool Border::is_empty(){
    if(m_empty){
        return true;
    }
    else if(m_volume_in.is_empty() && m_volume_out.is_empty()){
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
        C_Polyhedron ph_test(m_volume_in);
        ph_test.upper_bound_assign(m_volume_out);
        if(ph_test != m_volume_full){
            m_full = false;
            return false;
        }
        else{
            return true;
        }
    }
}

void Border::set_volume_in(PPL::C_Polyhedron volume_in, bool inclusion){
    if(inclusion)
        m_volume_in.intersection_assign(volume_in);
    else{
        C_Polyhedron ph_tmp(volume_in);
        ph_tmp.intersection_assign(m_volume_full);
        m_volume_in.upper_bound_assign(ph_tmp);
    }
}

void Border::set_volume_out(PPL::C_Polyhedron volume_out, bool inclusion){
    if(inclusion)
        m_volume_out.intersection_assign(volume_out);
    else{
        C_Polyhedron ph_tmp(volume_out);
        ph_tmp.intersection_assign(m_volume_full);
        m_volume_out.upper_bound_assign(ph_tmp);
    }
}

void Border::set_pave(Pave* pave){
    m_pave = pave;
}

const PPL::C_Polyhedron Border::get_volume_in() const{
    return m_volume_in;
}

const PPL::C_Polyhedron Border::get_volume_out() const{
    return m_volume_out;
}

const PPL::C_Polyhedron Border::get_volume_full() const{
    return m_volume_full;
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

Border& Border::operator&=(const Border &b){
    m_volume_in.intersection_assign(b.get_volume_in());
    m_volume_out.intersection_assign(b.get_volume_out());
    return *this;
}

bool Border::inter(const Border &b){
    bool change = false;
    PPL::C_Polyhedron ph_inter_in(m_volume_in);
    ph_inter_in.intersection_assign(b.get_volume_in());

    if(ph_inter_in != m_volume_in)
        change = true;

    PPL::C_Polyhedron ph_inter_out(m_volume_out);
    ph_inter_in.intersection_assign(b.get_volume_out());

    if(ph_inter_out != m_volume_out)
        change = true;

    m_volume_in = ph_inter_in;
    m_volume_out = ph_inter_out;

    return change;
}

bool Border::diff(const Border &b){
    bool change = false;

    C_Polyhedron ph_in_diff(m_volume_in);
    ph_in_diff.difference_assign(b.get_volume_in());
    if(ph_in_diff != m_volume_in)
        change = true;
    m_volume_in = ph_in_diff;

    C_Polyhedron ph_out_diff(m_volume_out);
    ph_in_diff.difference_assign(b.get_volume_out());
    if(ph_out_diff != m_volume_out)
        change = true;
    m_volume_out = ph_out_diff;

    return change;
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
    m_empty = false;
    m_full = false;
}

int Border::size() const{
    return m_inclusions.size();
}

Inclusion* Border::operator[](int id){
    return m_inclusions[id];
}

int Border::get_dim() const{
    return m_dim;
}

int Border::get_face_axis() const{
    return m_face_axis;
}
int Border::get_face_side() const{
    return m_face_side;
}

int Border::get_face() const{
    return 2*m_face_axis + m_face_side;
}

void Border::draw_vtk_get_points(vtkSmartPointer<vtkPoints> &points){
    for(int ph_id = 0; ph_id < 2; ph_id ++){
        C_Polyhedron *ph;
        if(ph_id==0)
            ph = &m_volume_in;
        else
            ph = &m_volume_out;

//        cout << "volume_full = " << m_volume_full.generators() << endl;
//        cout << "volume_in = " << m_volume_in.generators() << endl;
//        cout << "volume_out = " << m_volume_out.generators() << endl;

        if(ph->space_dimension()==2){
//            cout << ph->generators() << endl;
            for(auto &g:ph->generators()){
                if(g.is_point()){
                    std::vector<double> coord;
                    for(int cpt=0; cpt<2; cpt++){
                        coord.clear();
                        for(int i=0; i<3; i++){
                            PPL::Variable x(i);
                            if(i<2){
                                coord.push_back(g.coefficient(x).get_d()/(g.divisor().get_d()*IBEX_PPL_PRECISION));
                            }
                            else{
                                coord.push_back(-cpt/100.0);
                            }
                        }
                        points->InsertNextPoint(coord[0], coord[1], coord[2]);
                    }
                }
            }
        }
        else{
            for(auto &g:ph->generators()){
                if(g.is_point()){
                    std::vector<double> coord;
                    for(int i=0; i<3; i++){
                        PPL::Variable x(i);
                        if(g.space_dimension()>i){
                            coord.push_back(g.coefficient(x).get_d()/(g.divisor().get_d()*IBEX_PPL_PRECISION));
                        }
                        else{
                            coord.push_back(0.0);
                        }
                    }
                    points->InsertNextPoint(coord[0], coord[1], coord[2]);
                }
            }
        }
    }
}
