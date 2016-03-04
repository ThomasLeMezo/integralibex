#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void CtcPropagateSegment(const PPL::C_Polyhedron &volume_in, Pave *pave, vector<PPL::C_Polyhedron> list_volume_out, const std::vector<PPL::Generator>& ray_theta_list, const std::vector<PPL::Generator>& ray_command_list){
    C_Polyhedron ph_projection(volume_in);
    if(volume_in.is_empty())
        return;
    for(auto &ray:ray_theta_list){
        if(ray.is_ray())
            ph_projection.add_generator(ray);
        else
            cout << "RAY IS NOT A RAY : " << ray.type() << endl;
    }
//    for(auto &ray:ray_command_list){
//        ph_projection.add_generator(ray);
//    }

    list_volume_out.clear();
    for(int i=0; i<2*pave->get_dim(); i++){
        C_Polyhedron ph_face(ph_projection);
        ph_face.intersection_assign(pave->get_border(i)->get_volume_full());
        list_volume_out.push_back(ph_face);
    }
}

void CtcPaveBackward(Pave *p, bool inclusion, bool inner){

    //    ibex::Interval segment_in[4] = {ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET};

    //    for(int face = 0; face < 4; face++){
    //        ibex::Interval seg_in = p->get_border(face)->get_segment_in();

    //        vector<ibex::Interval> seg_out;
    //        for(int j=(face+1)%4; j!=face; j=(j+1)%4){
    //            seg_out.push_back(p->get_border(j)->get_segment_out());
    //        }

    //        if(!inner)
    //            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), false, false);
    //        else
    //            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), true, true);

    //        segment_in[face] = seg_in;
    //    }

    //    for(int face = 0; face<4; face++){
    //        p->get_border(face)->set_segment_in(segment_in[face], inclusion);
    //    }
}

void CtcPaveForward(Pave *p, bool inclusion, bool inner){
    int nb_face = 2*p->get_dim();

    vector<PPL::C_Polyhedron> list_volume_out;
    for(int face=0; face<nb_face; face++){
        list_volume_out.push_back(PPL::C_Polyhedron(p->get_dim(), PPL::EMPTY));
    }

    for(int face = 0; face<nb_face; face++){
        PPL::C_Polyhedron volume_in(p->get_border(face)->get_volume_in());
        vector<PPL::C_Polyhedron> list_volume_out_tmp(list_volume_out);
        if(!inner){
            CtcPropagateSegment(volume_in, p, list_volume_out_tmp, p->get_ray_vector_field(), p->get_ray_command());
        }
        else{
            // ToDo !!!
            CtcPropagateSegment(volume_in, p, list_volume_out_tmp, p->get_ray_vector_field(), p->get_ray_command());
        }

        for(int face_update = 0; face_update < nb_face; face_update++){
            list_volume_out[face_update].poly_hull_assign(list_volume_out_tmp[face_update]);
        }
    }

    for(int face = 0; face<nb_face; face++){
        p->get_border(face)->set_volume_out(list_volume_out[face], inclusion);
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

    for(int face = 0; face < 4; face++){
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

    return change;
}
