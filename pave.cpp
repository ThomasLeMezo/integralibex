#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "iomanip"

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVertex.h>
#include <vtkVertexGlyphFilter.h>

#include <vtkCellArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkCubeSource.h>
#include <vtkFloatArray.h>

#include "conversion.h"

using namespace std;
using namespace ibex;
using namespace Parma_Polyhedra_Library::IO_Operators;

Pave::Pave(const IntervalVector &position, ibex::Function *f, ibex::IntervalVector u):
    m_position(position.size()),
    m_u(position.size()),
    m_theta(position.size())
{
    m_position = position;    // Box corresponding to the Pave
    m_dim = position.size();
    m_borders.reserve(4);
    m_f = f;
    m_u = u;

    m_in_queue = false;
    m_copy_node = NULL;

    m_first_process = false;
    m_continuity = true;

    // Border building
    vector< vector<IntervalVector>> faces = get_faces(position);
    for(int face=0; face<faces.size(); face++){
        for(int side=0; side<faces[face].size(); side++){
            m_borders.push_back(new Border(faces[face][side], this, face, side));
        }
    }

    m_full = false;
    m_empty = false;

    // Build ray for theta & for u + theta
    m_theta = f->eval_vector(position);
    Linear_Expression e = Linear_Expression(0);
    std::vector<Linear_Expression> linear_expression_list;
    recursive_linear_expression_from_iv(m_theta, m_theta.size(), linear_expression_list,e);
    for(auto &l:linear_expression_list){
        if(!l.all_homogeneous_terms_are_zero()) // Case {0, 0, ...}
            m_ray_vector_field.push_back(ray(l));
    }

    Linear_Expression e_bwd = Linear_Expression(0);
    std::vector<Linear_Expression> linear_expression_list_backward;
    recursive_linear_expression_from_iv(-m_theta, m_theta.size(), linear_expression_list_backward, e_bwd);
    for(auto &l:linear_expression_list_backward){
        if(!l.all_homogeneous_terms_are_zero()) // Case {0, 0, ...}
            m_ray_vector_field_backward.push_back(ray(l));
    }

    //    e = Linear_Expression(0);
    //    linear_expression_list.clear();
    //    recursive_linear_expression_from_iv(u + theta, theta.size(), linear_expression_list,e);
    //    for(auto &l:linear_expression_list){
    //        m_ray_command.push_back(ray(l));
    //    }
}

Pave::Pave(const Pave *p):
    m_position(p->get_dim()),
    m_u(p->get_dim()),
    m_theta(p->get_dim())
{
    m_position = p->get_position();    // Box corresponding to the Pave
    m_dim = p->get_dim();
    m_f = p->get_f();
    m_u = p->get_u();
    m_full = true; // Force to recompute results
    m_empty = false;
    m_in_queue = false;
    m_first_process = false;
    m_continuity = true;

    m_ray_command = p->get_ray_command();
    m_ray_vector_field = p->get_ray_vector_field();

    for(int face = 0; face < 2*p->get_dim(); face++){
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
    m_empty = false;
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

void Pave::draw_vtk(vtkSmartPointer<vtkAppendPolyData> &polyData){
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(auto &border:m_borders){
        border->draw_vtk_get_points(points);
    }
    if(points->GetNumberOfPoints() == 0)
        return;

    // ********** Points **************
    vtkSmartPointer< vtkPolyData> pointsCollection = vtkSmartPointer<vtkPolyData>::New();
    pointsCollection->SetPoints(points);
    //polydata->SetVerts(vertices);

    //    vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    //    Colors->SetNumberOfComponents ( 3 );
    //    Colors->SetName ( "RGB" );
    //    for(int i=0; i<points->GetNumberOfPoints(); i++)
    //        Colors->InsertNextTuple3(i,i,0);
    //    pointsCollection->GetPointData()->SetVectors(Colors);

    // ********** Surface **************
    // Create the convex hull of the pointcloud (delaunay + outer surface)
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer< vtkDelaunay3D >::New();
    delaunay->SetInputData(pointsCollection);
    delaunay->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
    surfaceFilter->Update();

    // ********** Append results **************
    polyData->AddInputData(surfaceFilter->GetOutput());

    //    polyData->GetOutput()->GetPointData()->SetScalars(Colors);
}

void Pave::draw_box(vtkSmartPointer<vtkAppendPolyData> &polyData){
    vtkSmartPointer<vtkCubeSource> cubedata = vtkSmartPointer<vtkCubeSource>::New();
    double bounds[6];
    for(int i=0; i<3; i++){
        if(m_position.size()>i){
            bounds[2*i] = m_position[i].lb();
            bounds[2*i+1] = m_position[i].ub();
        }
        else{
            bounds[2*i] = 0.0;
            bounds[2*i+1] = 0.0;
        }
    }

    cubedata->SetBounds(bounds);
    cubedata->Update();
    polyData->AddInputData(cubedata->GetOutput());
}

void Pave::draw_vector_field(vtkSmartPointer<vtkAppendPolyData> &polyData){
    ibex::IntervalVector position(3);
    ibex::IntervalVector theta(3);

    if(m_position.size()>=3){
        position[0] = m_position[0];
        position[1] = m_position[1];
        position[2] = m_position[2];
        theta[0] = m_theta[0];
        theta[1] = m_theta[1];
        theta[2] = m_theta[2];
    }
    if(m_position.size()==2){
        position[0] = m_position[0];
        position[1] = m_position[1];
        position[2] = ibex::Interval(0.0, 1.0);
        theta[0] = m_theta[0];
        theta[1] = m_theta[1];
        theta[2] = ibex::Interval(0.0);
    }

    vector< vector<double>> list_point = get_points_from_iv(theta);

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer<vtkFloatArray> field = vtkSmartPointer<vtkFloatArray>::New();
    field->SetNumberOfComponents(3);
    field->SetName("Glyph");

    for(auto &vect:list_point){
        points->InsertNextPoint(position[0].mid(), position[1].mid(), position[2].mid());
        field->InsertNextTuple3(vect[0], vect[1], vect[2]);
    }

    vtkSmartPointer< vtkPolyData> dataObject = vtkSmartPointer<vtkPolyData>::New();
    dataObject->SetPoints(points);
    dataObject->GetPointData()->SetVectors(field);

    // Transform points (array) in vertex objects
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(dataObject);
    vertexFilter->Update();

    // ********** Append results **************
    polyData->AddInputData(vertexFilter->GetOutput());
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result){
    // Create 2 new paves
    ibex::LargestFirst bisector(0.0, 0.5);

    std::pair<ibex::IntervalVector, IntervalVector> result_boxes = bisector.bisect(m_position);

    Pave *pave1 = new Pave(result_boxes.first, m_f, m_u); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, m_f, m_u); // Right or Down

    // Find the axe of bissection
    int bisect_axis = 0;
    for(int i=0; i<m_position.size(); i++){
        if(result_boxes.first[i] != m_position[i]){
            bisect_axis = i;
            break;
        }
    }

    int indice1, indice2;
    indice1 = 2*bisect_axis + 1;    // Indice of the common face with pave2
    indice2 = 2*bisect_axis;

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
    Inclusion *inclusion_to_pave2 = new Inclusion(pave2->get_border(indice2), bisect_axis, 1);
    Inclusion *inclusion_to_pave1 = new Inclusion(pave1->get_border(indice1), bisect_axis, 0);

    pave1->get_border(indice1)->add_inclusion(inclusion_to_pave2);
    pave2->get_border(indice2)->add_inclusion(inclusion_to_pave1);


    // Intersect initial pave Polyhedron with pave1 & pave2
    // Union of volume IN & OUT is MANDATORY to work !!!
    if(is_full()){
        pave1->set_full();
        pave2->set_full();
    }
    else{
        PPL::C_Polyhedron volume_in_out(get_volume_in_out());
        pave1->set_volume_in_out(volume_in_out);
        pave2->set_volume_in_out(volume_in_out);
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

void Pave::disable_continuity(){
    m_continuity = false;
}

const vector<PPL::Generator> &Pave::get_ray_vector_field() const{
    return m_ray_vector_field;
}

const vector<PPL::Generator> &Pave::get_ray_vector_backward_field() const{
    return m_ray_vector_field_backward;
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
//}

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

IntervalVector Pave::get_u() const{
    return m_u;
}

IntervalVector Pave::get_theta() const{
    return m_theta;
}

bool Pave::get_continuity() const{
    return m_continuity;
}

const C_Polyhedron Pave::get_volume_in() const{
    PPL::C_Polyhedron volume_in(m_dim, PPL::EMPTY);
    for(auto &b:m_borders){
        volume_in.poly_hull_assign(b->get_volume_in());
    }
    return volume_in;
}

const C_Polyhedron Pave::get_volume_out() const{
    PPL::C_Polyhedron volume_out(m_dim, PPL::EMPTY);
    for(auto &b:m_borders){
        volume_out.poly_hull_assign(b->get_volume_out());
    }
    return volume_out;
}

const C_Polyhedron Pave::get_volume_in_out() const{
    PPL::C_Polyhedron volume_in_out(m_dim, PPL::EMPTY);
    volume_in_out.poly_hull_assign(get_volume_in());
    volume_in_out.poly_hull_assign(get_volume_out());
    return volume_in_out;
}

void Pave::set_volume_in(const C_Polyhedron &volume_in){
    for(auto &b:m_borders){
        b->set_volume_in(volume_in, false);
    }
}

void Pave::set_volume_out(const C_Polyhedron &volume_out){
    for(auto &b:m_borders){
        b->set_volume_out(volume_out, false);
    }
}

void Pave::set_volume_in_out(const PPL::C_Polyhedron& volume_in_out){
    set_volume_in(volume_in_out);
    set_volume_out(volume_in_out);
}


