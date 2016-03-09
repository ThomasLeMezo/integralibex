#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void CtcPropagateSegment(PPL::C_Polyhedron &volume_in, Pave *pave, vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator>& ray_vector_field_list, const std::vector<PPL::Generator>& ray_command_list){
    C_Polyhedron ph_projection(volume_in);
    if(volume_in.is_empty())
        return;
    for(auto &ray:ray_vector_field_list){
        if(ray.is_ray())
            ph_projection.add_generator(ray);
        else
            cout << "RAY IS NOT A RAY : " << ray.type() << endl;
    }

    list_volume_out.clear();
    for(int i=0; i<2*pave->get_dim(); i++){
        C_Polyhedron ph_face(ph_projection);
        ph_face.intersection_assign(pave->get_border(i)->get_volume_full());
        list_volume_out.push_back(ph_face);
    }    
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
        CtcPropagateSegment(volume_in, p, list_volume_out_tmp, p->get_ray_vector_field(), p->get_ray_command());
        if(list_volume_out_tmp.size() != 0){
            for(int face_update = 0; face_update < nb_face; face_update++){
                list_volume_out[face_update].poly_hull_assign(list_volume_out_tmp[face_update]);
            }
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
void CtcPropagateSegmentBackward(PPL::C_Polyhedron &volume_in, Pave *pave, vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator>& ray_vector_field_backward_list, const std::vector<PPL::Generator>& ray_command_list){
    if(volume_in.is_empty())
        return;
    C_Polyhedron ph_projection(pave->get_dim(), PPL::EMPTY);

    for(auto &ph:list_volume_out){
        ph_projection.poly_hull_assign(ph);
    }

    for(auto &ray:ray_vector_field_backward_list){
        if(ray.is_ray())
            ph_projection.add_generator(ray);
        else
            cout << "RAY IS NOT A RAY : " << ray.type() << endl;
    }

    volume_in.intersection_assign(ph_projection);
}

void CtcPaveBackward(Pave *p, bool inclusion, bool inner){
    int nb_face = 2*p->get_dim();

    vector<PPL::C_Polyhedron> list_volume_in;
    for(int face=0; face<nb_face; face++){
        list_volume_in.push_back(PPL::C_Polyhedron(p->get_dim(), PPL::EMPTY));
    }

    for(int face = 0; face<nb_face; face++){
        vector<PPL::C_Polyhedron> list_volume_out_tmp;
        for(int face=0; face<nb_face; face++){
            list_volume_out_tmp.push_back(p->get_border(face)->get_volume_out());
        }
        PPL::C_Polyhedron volume_in(p->get_border(face)->get_volume_in());

        CtcPropagateSegment(volume_in, p, list_volume_out_tmp, p->get_ray_vector_field(), p->get_ray_command());
        list_volume_in.push_back(volume_in);
    }
    for(int face = 0; face<nb_face; face++){
        p->get_border(face)->set_volume_in(list_volume_in[face], inclusion);
    }
}

// ********************************************************************************
// ****************** Algorithm functions      ************************************

void CtcPaveConsistency(Pave *p, bool backward, bool inner){
    if(backward){
        CtcPaveBackward(p, backward, inner);
        Pave p2(p);
        CtcPaveForward(&p2, backward, inner);
        *p &= p2;
        if(inner){
            Pave p3(p);
            CtcPaveBackward(&p3, backward, inner);
            *p &= p3;
        }
    }
    else{
        CtcPaveForward(p, backward, inner);
    }
}

bool CtcContinuity(Pave *p, bool backward){
    bool change = false;
    int nb_face = 2*p->get_dim();

    for(int face = 0; face < nb_face; face++){
        if(p->get_border(face)->get_inclusions().size() > 0){
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
                    p->get_border(face)->set_volume_in(volume_out, backward);
                    p->get_border(face)->set_volume_out(volume_in, backward);
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

    return change;
}
