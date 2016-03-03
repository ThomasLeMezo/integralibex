#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"
#include "pave.h"

void CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> &theta, const ibex::IntervalVector &box_pave, const ibex::Interval &u, bool inner=false, bool inner_backward=false);

void CtcPaveForward(Pave *p, bool inclusion, bool inner);
void CtcPaveBackward(Pave *p, bool inclusion, bool inner);
void CtcPaveConsistency(Pave *p, bool backward, bool inner);
bool CtcContinuity(Pave *p, bool backward);


#endif // UTILS_H
