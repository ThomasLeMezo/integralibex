#ifndef CONVERSION
#define CONVERSION

#include <ppl.hh>
#include <ibex.h>

#define IBEX_PPL_PRECISION 1e14

Parma_Polyhedra_Library::Rational_Box iv_2_box(const ibex::IntervalVector &iv);
std::vector< vector<ibex::IntervalVector>> get_faces(ibex::IntervalVector pave);

void recursive_linear_expression_from_iv(const ibex::IntervalVector &theta,
                                         int dim,
                                         std::vector<Linear_Expression> &linear_expression_list,
                                         Linear_Expression &local_linear_expression);

ibex::IntervalVector ph_2_iv(const Parma_Polyhedra_Library::C_Polyhedron &ph);

void recursive_get_points(int dim,
                          const std::vector<double> &pt,
                          std::vector< std::vector<double>> &list,
                          const ibex::IntervalVector &iv);

std::vector< std::vector<double>> get_points_from_iv(const ibex::IntervalVector &iv);

#endif // CONVERSION
