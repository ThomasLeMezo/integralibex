#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"
#include "pave.h"

void CtcPropagateSegment(const C_Polyhedron &volume_in, Pave *pave, vector<PPL::C_Polyhedron> list_volume_out, const std::vector<Generator> &ray_theta_list, const std::vector<Generator> &ray_command_list);

void CtcPaveForward(Pave *p, bool inclusion, bool inner);
void CtcPaveBackward(Pave *p, bool inclusion, bool inner);
void CtcPaveConsistency(Pave *p, bool backward, bool inner);
bool CtcContinuity(Pave *p, bool backward);


#endif // UTILS_H
