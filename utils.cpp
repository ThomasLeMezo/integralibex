#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;
using namespace Parma_Polyhedra_Library::IO_Operators;

// cout << setprecision(80) << "..." << endl;

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void CtcPropagateSegment(const PPL::C_Polyhedron &volume_in, Pave *pave, vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator>& ray_vector_field_list, const std::vector<PPL::Generator>& ray_command_list){
    C_Polyhedron ph_projection(volume_in);
    if(!volume_in.is_discrete()){
        for(auto &ray:ray_vector_field_list){
            if(ray.is_ray())
                ph_projection.add_generator(ray);
            else
                cout << "RAY IS NOT A RAY : " << ray.type() << endl;
        }
    }

    for(int i=0; i<list_volume_out.size(); i++)
        list_volume_out[i].intersection_assign(ph_projection);
}

void CtcPaveForward(Pave *p, bool inclusion, bool inner){
    int nb_face = 2*p->get_dim();

    vector<PPL::C_Polyhedron> list_volume_out;
    for(int face=0; face<nb_face; face++){
        list_volume_out.push_back(PPL::C_Polyhedron(p->get_dim(), PPL::EMPTY));
    }

    for(int face = 0; face<nb_face; face++){
        PPL::C_Polyhedron volume_in(p->get_border(face)->get_volume_in());
        vector<PPL::C_Polyhedron> list_volume_out_tmp;
        for(int i=0; i<nb_face; i++){
            list_volume_out_tmp.push_back(p->get_border(i)->get_volume_full());
        }

        CtcPropagateSegment(volume_in, p, list_volume_out_tmp, p->get_ray_vector_field(), p->get_ray_command());

        for(int i=0; i<nb_face; i++){
            if(i!=face)
                list_volume_out[i].poly_hull_assign(list_volume_out_tmp[i]);
        }
    }

    for(int face = 0; face<nb_face; face++){
        p->get_border(face)->set_volume_out(list_volume_out[face], inclusion);
    }
}

///
/// \brief CtcPropagateSegmentBackward
/// \param volume_in
/// \param pave
/// \param list_volume_out
/// \param ray_vector_field_list
/// \param ray_command_list
///
void CtcPropagateSegmentBackward(PPL::C_Polyhedron &volume_in, Pave *pave, const vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator>& ray_vector_field_backward_list, const std::vector<PPL::Generator>& ray_command_list){
    if(volume_in.is_empty())
        return;
    C_Polyhedron volume_in_tmp(pave->get_dim(), PPL::EMPTY);

    for(auto &ph:list_volume_out){
        if(!ph.is_empty()){
            C_Polyhedron ph_projection(ph);

            for(auto &ray:ray_vector_field_backward_list){
                if(ray.is_ray())
                    ph_projection.add_generator(ray);
                else
                    cout << "RAY IS NOT A RAY : " << ray.type() << endl;
            }
            ph_projection.intersection_assign(volume_in);
            volume_in_tmp.poly_hull_assign(ph_projection);
        }
    }
    volume_in = volume_in_tmp;
}

void CtcPaveBackward(Pave *p, bool inclusion, bool inner){
    vector<PPL::C_Polyhedron> list_volume_in;

    for(auto &b:p->get_borders()){ // For each volume IN
        PPL::C_Polyhedron volume_in(b->get_volume_in());
        if(!volume_in.is_discrete()){ // Avoid single point IN
            vector<PPL::C_Polyhedron> list_volume_out;

            // Make a list of volume OUT
            for(auto &b_tmp:p->get_borders()){
                if(b_tmp != b){ // Do not project on the same border
                    if(!b_tmp->get_volume_out().is_discrete()) // Do not keep single point
                        list_volume_out.push_back(b_tmp->get_volume_out());
                }
            }

            CtcPropagateSegmentBackward(volume_in, p, list_volume_out, p->get_ray_vector_backward_field(), p->get_ray_command());
        }
        if(!volume_in.is_discrete()) // Do not add a single point volume IN
            list_volume_in.push_back(volume_in); // Save temporary the result
        else
            list_volume_in.push_back(PPL::C_Polyhedron(p->get_dim(), PPL::EMPTY));
    }

    // Write the results to the pave's borders
    for(auto &b:p->get_borders()){
        b->set_volume_in(list_volume_in[b->get_face()], inclusion);
    }
}

// ********************************************************************************
// ****************** Algorithm functions      ************************************

void CtcPaveConsistency(Pave *p, bool backward, bool inner){
    if(backward){
        IntervalVector position(2);
        position[0] = ibex::Interval(-3, -2);
        position[1] = ibex::Interval(-1, 0);

        if(p->get_position() == position){
            for(int i=0; i<4; i++){
                cout << "i = " << i << endl;
                cout << "volume_in = " << p->get_border(i)->get_volume_in().generators() << endl;
                cout << "volume_out = " << p->get_border(i)->get_volume_out().generators() << endl;
                cout << endl;
            }
        }

        CtcPaveBackward(p, true, inner);
        Pave p2(p);
        CtcPaveForward(&p2, true, inner);
        *p &= p2;

        if(p->get_position() == position){
            for(int i=0; i<4; i++){
                cout << "i = " << i << endl;
                cout << "volume_in = " << p->get_border(i)->get_volume_in().generators() << endl;
                cout << "volume_out = " << p->get_border(i)->get_volume_out().generators() << endl;
                cout << endl;
            }
        }

        if(inner){
            Pave p3(p);
            CtcPaveBackward(&p3, false, inner);
            *p &= p3;
        }
    }
    else{
        CtcPaveForward(p, backward, inner);
    }
}

bool CtcContinuity(Pave *p, bool backward){
    bool change = false;

    if(p->get_continuity()){
        for(int face = 0; face < p->get_borders().size(); face++){
            PPL::C_Polyhedron volume_in(p->get_dim(), PPL::EMPTY);
            PPL::C_Polyhedron volume_out(p->get_dim(), PPL::EMPTY);

            for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
                volume_in.poly_hull_assign(p->get_border(face)->get_inclusion(b)->get_volume_in());
                volume_out.poly_hull_assign(p->get_border(face)->get_inclusion(b)->get_volume_out());
            }

            if(backward){
                PPL::C_Polyhedron volume_inter_out(volume_out), volume_inter_in(volume_in);
                volume_inter_out.intersection_assign(p->get_border(face)->get_volume_in());
                volume_inter_in.intersection_assign(p->get_border(face)->get_volume_out());

                if(p->get_border(face)->get_volume_in() != volume_inter_out || p->get_border(face)->get_volume_out() != volume_inter_in){
                    change = true;
                    p->get_border(face)->set_volume_in(volume_inter_out, backward);
                    p->get_border(face)->set_volume_out(volume_inter_in, backward);
                }
            }
            else{
                PPL::C_Polyhedron volume_test(p->get_border(face)->get_volume_in());
                volume_test.poly_hull_assign(volume_out);
                volume_test.intersection_assign(p->get_border(face)->get_volume_full());

                if(p->get_border(face)->get_volume_in() != volume_test){
                    change = true;
                    p->get_border(face)->set_volume_in(volume_out, backward);
                }
            }
        }
    }
    else{
        change = true;
    }

    return change;
}
