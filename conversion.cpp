#include "conversion.h"
#include <iostream>
#include "ibex.h"
#include "ppl.hh"

using namespace std;
namespace PPL = Parma_Polyhedra_Library;
using namespace ibex;

PPL::Rational_Box iv_2_box(const ibex::IntervalVector &iv){
    Rational_Box box(iv.size());
    for(int i=0; i<box.space_dimension(); i++){
        PPL::Variable x(i);

        if(!iv[i].is_empty()){
            box.add_constraint(x >= floor(iv[i].lb()*IBEX_PPL_PRECISION));
            box.add_constraint(x <= ceil(iv[i].ub()*IBEX_PPL_PRECISION));
        }
    }
    return box;
}

ibex::IntervalVector ph_2_iv(const PPL::C_Polyhedron &ph){
    PPL::Rational_Box rb(ph);
    IntervalVector box(ph.space_dimension());
    for(int dim=0; dim<box.size(); dim++){
        PPL::Variable x(dim);
        box[dim] = ibex::Interval(rb.get_interval(x).lower().get_d()/IBEX_PPL_PRECISION,
                                  rb.get_interval(x).upper().get_d()/IBEX_PPL_PRECISION);
    }
    return box;
}

std::vector< vector<IntervalVector>> get_faces(ibex::IntervalVector pave){
    std::vector< vector<IntervalVector>> face_list;
    for(int face=0; face<pave.size(); face++){
        std::vector<IntervalVector> side_list;
        for(int side = 0; side<2; side++){
            IntervalVector tmp_face(pave.size());
            for(int i=0; i<pave.size(); i++){
                if(i==face){
                    if(side==0){
                        tmp_face[i] = ibex::Interval(pave[i].lb());
                    }
                    else{
                        tmp_face[i] = ibex::Interval(pave[i].ub());
                    }
                }
                else{
                    tmp_face[i] = pave[i];
                }
            }
            side_list.push_back(tmp_face);
        }
        face_list.push_back(side_list);
    }
    return face_list;
}

void recursive_linear_expression_from_iv(const ibex::IntervalVector &theta,
                                         int dim,
                                         std::vector<Linear_Expression> &linear_expression_list,
                                         Linear_Expression &local_linear_expression){
    if(dim > 0){
        PPL::Variable x(dim-1);
        Linear_Expression l_b = local_linear_expression;
        Linear_Expression l_u = local_linear_expression;

        // ToDo: case theta[dim] -> lb=+oo | ub=-oo
        if(std::isinf(theta[dim-1].ub())){
            linear_expression_list.push_back(Linear_Expression(x));
            cout << "INFINITY" << endl;
        }
        else
            l_u += x*ceil(theta[dim-1].ub()*IBEX_PPL_PRECISION);

        if(std::isinf(theta[dim-1].lb())){
            linear_expression_list.push_back(Linear_Expression(-x));
            cout << "INFINITY" << endl;
        }
        else
            l_b += x*floor(theta[dim-1].lb()*IBEX_PPL_PRECISION);

        recursive_linear_expression_from_iv(theta, dim-1, linear_expression_list, l_u);
        recursive_linear_expression_from_iv(theta, dim-1, linear_expression_list, l_b);
    }
    else{
        linear_expression_list.push_back(local_linear_expression);
    }
}

void recursive_get_points(int dim,
                          const std::vector<double> &pt,
                          std::vector< std::vector<double>> &list,
                          const ibex::IntervalVector &iv){
    if(dim < iv.size()){
        std::vector<double> pt_lb(pt), pt_ub(pt);
        pt_lb.push_back(iv[dim].lb());
        pt_ub.push_back(iv[dim].ub());
        recursive_get_points(dim+1, pt_lb, list,iv);
        recursive_get_points(dim+1, pt_ub, list,iv);
    }
    else{
        list.push_back(pt);
    }
}

std::vector< std::vector<double>> get_points_from_iv(const ibex::IntervalVector &iv){

    std::vector< std::vector<double>> list;
    std::vector<double> point;
    recursive_get_points(0, point, list,iv);
    return list;
}
