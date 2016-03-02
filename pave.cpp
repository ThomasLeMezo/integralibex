#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "iomanip"

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAppendPolyData.h>

#include "conversion.h"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &position, ibex::Function *f, ibex::IntervalVector u): m_position(position.size())
{
    m_position = position;    // Box corresponding to the Pave
    m_dim = position.size();
    m_borders.reserve(4);
    m_f = f;
    m_u = u;

    m_in_queue = false;
    m_copy_node = NULL;

    m_first_process = false;

    // Border building
    vector< vector<IntervalVector>> faces = get_faces(position);
    for(int face=0; face<faces.size(); face++){
        for(int side=0; side<faces[face].size(); side++){
            m_borders.push_back(new Border(face, this, face, side));
        }
    }

    m_full = false;
    m_empty = false;

    // Build ray for theta & for u
    IntervalVector theta = f->eval_vector(position);
    Linear_Expression e = Linear_Expression(0);
    std::vector<Linear_Expression> linear_expression_list;
    recursive_linear_expression_from_iv(theta, theta.size(), linear_expression_list,e);
    for(auto &l:linear_expression_list){
        m_ray_vector_field.push_back(ray(l));
    }

    e = Linear_Expression(0);
    linear_expression_list.clear();
    recursive_linear_expression_from_iv(u, u.size(), linear_expression_list,e);
    for(auto &l:linear_expression_list){
        m_ray_command.push_back(ray(l));
    }
}

Pave::Pave(const Pave *p): m_position(p->get_dim())
{
    m_position = p->get_position();    // Box corresponding to the Pave
    m_dim = p->get_dim();
    m_f = p->get_f();
    m_full = true; // Force to recompute results
    m_empty = false;
    m_in_queue = false;
    m_first_process = false;

    m_ray_command = p->get_ray_command();
    m_ray_vector_field = p->get_ray_vector_field();

    for(int face = 0; face < m_borders.size(); face++){
        Border *b_cpy = new Border(p->get_border_const(face));
        b_cpy->set_pave(this);
        m_borders.push_back(b_cpy);
    }
    m_copy_node = NULL;
}

Pave::~Pave(){
    for(auto &b:m_borders){
        delete(b);
    }
}

Pave& Pave::operator&=(const Pave &p){    
    for(int face = 0; face < m_borders.size(); face++){
        *(m_borders[face]) &= *(p.get_border_const(face));
    }
    return *this;
}

bool Pave::inter(const Pave &p){
    bool change = false;
    for(int face = 0; face <m_borders.size(); face++){
        if(m_borders[face]->inter(*(p.get_border_const(face))))
            change = true;
    }
    return change;
}

bool Pave::diff(const Pave &p){
    bool change = false;
    for(int face = 0; face<m_borders.size(); face++){
        if(m_borders[face]->diff(*(p.get_border_const(face)))){
            change = true;
        }
    }
    m_empty=false; // forces to recompute the value
    m_full=true;
}

void Pave::set_theta(IntervalVector theta){
    Linear_Expression e = Linear_Expression(0);
    std::vector<Linear_Expression> linear_expression_list;
    recursive_linear_expression_from_iv(theta, theta.size(), linear_expression_list,e);

    m_ray_vector_field.clear();
    for(auto &l:linear_expression_list){
        m_ray_vector_field.push_back(ray(l));
    }
}

void Pave::set_full(){
    for(int face=0; face<m_borders.size(); face++){
        m_borders[face]->set_full();
    }
    m_full = true;
}

void Pave::set_empty(){
    for(int face=0; face<m_borders.size(); face++){
        m_borders[face]->set_empty();
    }
    m_empty = true;
    m_full = false;
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

vtkPolyData draw_vtk(){
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();

    for(auto &ph:ph_list){
        for(auto &g:ph.generators()){
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

    vtkSmartPointer< vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    //polydata->SetVerts(vertices);

    // Create the convex hull of the pointcloud (delaunay + outer surface)
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer< vtkDelaunay3D >::New();
    delaunay->SetInputData(polydata);
    delaunay->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
    surfaceFilter->Update();

    return surfaceFilter->GetOutput();
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result){
    // Create 2 new paves
    ibex::LargestFirst bisector(0.0, 0.5);

    std::pair<ibex::IntervalVector, IntervalVector> result_boxes = bisector.bisect(m_position);

    Pave *pave1 = new Pave(result_boxes.first, m_f, m_u); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, m_f, m_u); // Right or Down

    int indice1, indice2;

    if(pave1->m_position[0] == m_position[0]){
        // Case UP/DOWN bisection
        indice1 = 2;
        indice2 = 0;
    }
    else{
        // Case LEFT/RIGHT bisection
        indice1 = 1;
        indice2 = 3;
    }

    // The order of tasks is important !
    // 1) Update pave brothers with pave1 & pave2
    for(int face=0; face<m_borders.size(); face++){
        m_borders[face]->update_brothers_inclusion(pave1->get_border(face), pave2->get_border(face));
    }

    // 2) Copy brothers Pave (this) to pave1 and pave2
    for(int face=0; face<m_borders.size(); face++){
        if(m_borders[face]->get_inclusions().size()!=0){
            if(face!=indice1){
                pave1->get_border(face)->add_inclusions(m_borders[face]->get_inclusions());
            }
            if(face!=indice2){
                pave2->get_border(face)->add_inclusions(m_borders[face]->get_inclusions());
            }
        }
    }

    // 3) Add each other to its brother list (pave1 <-> pave2)
    Inclusion *inclusion_to_pave2 = new Inclusion(pave2->get_border(indice2), indice2);
    Inclusion *inclusion_to_pave1 = new Inclusion(pave1->get_border(indice1), indice1);

    pave1->get_border(indice1)->add_inclusion(inclusion_to_pave2);
    pave2->get_border(indice2)->add_inclusion(inclusion_to_pave1);

    if(is_full()){
        pave1->set_full();
        pave2->set_full();
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

// ********************************************************************************
// ****************** UTILS functions *********************************************

void Pave::remove_brothers(Pave* p, int face){
    for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
        if(m_borders[face]->get_inclusion(i)->get_border()->get_pave() == p){
            m_borders[face]->remove_inclusion(i);
            return;
        }
    }
}

void Pave::remove_from_brothers(){
    for(int face=0; face<4; face++){
        for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
            m_borders[face]->get_inclusion(i)->get_border()->get_pave()->remove_brothers(this, m_borders[face]->get_inclusion(i)->get_brother_face());
        }
    }
}

bool Pave::is_empty(){
    if(m_empty){
        return true;
    }
    else{
        for(int i=0; i<m_borders.size(); i++){
            if(!m_borders[i]->is_empty()){
                return false;
            }
        }

        m_empty = true;
        return true;
    }
}

bool Pave::is_full(){
    if(!m_full){
        return false;
    }
    else{
        for(int face=0; face<m_borders.size(); face++){
            if(!m_borders[face]->is_full()){
                m_full = false;
                return false;
            }
        }
        m_full = true;
        return true;
    }
}

const std::vector<Pave *> Pave::get_brothers(int face){
    vector<Pave*> brothers_list;
    for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
        brothers_list.push_back(m_borders[face]->get_inclusion(i)->get_border()->get_pave());
    }
    return brothers_list;
}

void Pave::reset_full_empty(){
    m_empty = false;
    m_full = true;
    for(auto &border: m_borders){
        border->reset_full_empty();
    }
}

const vector<PPL::Generator> &Pave::get_ray_vector_field() const{
    return m_ray_vector_field;
}

const vector<PPL::Generator> &Pave::get_ray_command() const{
    return m_ray_command;
}

const IntervalVector &Pave::get_position() const{
    return m_position;
}

const std::vector<Border*> &Pave::get_borders(){
    return m_borders;
}

Border* Pave::get_border(int face){
    assert(face >=0 && face < m_borders.size());
    return m_borders[face];
}

Border* Pave::operator[](int face){
    return m_borders[face];
}

const Border* Pave::get_border_const(int face) const{
    return m_borders[face];
}

bool Pave::is_in_queue() const{
    return m_in_queue;
}

void Pave::set_in_queue(bool flag){
    m_in_queue = flag;
}

ibex::Function* Pave::get_f() const{
    return m_f;
}

void Pave::set_copy_node(Pave *p){
    m_copy_node = p;
}

Pave* Pave::get_copy_node(){
    return m_copy_node;
}

//void Pave::print(){
//    cout << "********" << endl;
//    cout << "PAVE x=" << m_position[0] << " y= " << m_position[1] << endl;
//    cout << this << endl;
//    for(int face = 0; face < m_borders.size(); face++){
//        if(m_borders[face]->get_inclusions().size()==0){
//            cout << "border=" << face << " " << &(m_borders[face])
//                 << " segment_in=" << m_borders[face]->get_segment_in()
//                 << " segment_out=" << m_borders[face]->get_segment_out()
//                 << endl;
//        }
//        else{
//            for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
//                cout << "border=" << face << " " << &(m_borders[face])
//                     << " segment_in=" << m_borders[face]->get_segment_in()
//                     << " segment_out=" << m_borders[face]->get_segment_out()
//                     << " inclusion=" << i
//                     << " *border=" << m_borders[face]->get_inclusion(i)->get_border()
//                     << " segment_full=" << m_borders[face]->get_inclusion(i)->get_border()->get_segment_full()
//                     << endl;
//            }
//        }
//    }
}

bool Pave::get_first_process() const{
    return m_first_process;
}

void Pave::set_first_process_true(){
    m_first_process = true;
}

void Pave::set_first_process_false(){
    m_first_process = false;
}

int Pave::get_dim() const{
    return m_dim;
}

int Pave::get_size() const{
    return m_borders.size();
}
