#include <iostream>

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
