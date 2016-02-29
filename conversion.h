#ifndef CONVERSION
#define CONVERSION

#include <ppl.hh>
#include <ibex.h>

#define IBEX_PPL_PRECISION 1e16

Parma_Polyhedra_Library::Rational_Box iv_2_box(const ibex::IntervalVector &iv);


#endif // CONVERSION
