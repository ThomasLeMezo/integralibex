#include <iostream>
#include "ibex.h"
#include "ppl.hh"

using namespace std;
namespace PPL = Parma_Polyhedra_Library;
using namespace ibex;

#define IBEX_PPL_PRECISION 1e16

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
        box[dim] = Interval(rb.get_interval(x).lower().get_d()/IBEX_PPL_PRECISION,
                            rb.get_interval(x).upper().get_d()/IBEX_PPL_PRECISION);
    }
    return box;
}

std::vector< vector<IntervalVector>> get_faces(ibex::IntervalVector pave){
    std::vector< vector<IntervalVector>> face_list;
    for(int face=0; face<pave.size(); face++){
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
            face_list.push_back(tmp_face);
        }
    }
    return face_list;
}

void recursive_linear_expression_from_iv(const ibex::IntervalVector &theta,
                                         int dim,
                                         std::vector<Linear_Expression> &linear_expression_list,
                                         Linear_Expression &local_linear_expression){
    if(dim > 0){
        PPL::Variable x(dim-1);

        // ToDo: case theta[dim]= +oo / -oo
        Linear_Expression l_u = x*ceil(theta[dim-1].ub()*IBEX_PPL_PRECISION) + local_linear_expression;
        Linear_Expression l_l = x*floor(theta[dim-1].lb()*IBEX_PPL_PRECISION) + local_linear_expression;

        recursive_linear_expression_from_iv(theta, dim-1, linear_expression_list, l_u);
        recursive_linear_expression_from_iv(theta, dim-1, linear_expression_list, l_l);
    }
    else{
        linear_expression_list.push_back(local_linear_expression);
    }
}
