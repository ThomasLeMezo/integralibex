#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"
#include "pave.h"

void CtcPropagateSegment(const PPL::C_Polyhedron &volume_in, Pave *p, std::vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator> &ray_vector_field_list, const std::vector<PPL::Generator> &ray_command_list);
void CtcPropagateSegmentBackward(PPL::C_Polyhedron &volume_in, Pave *p, const std::vector<PPL::C_Polyhedron> &list_volume_out, const std::vector<PPL::Generator>& ray_vector_field_backward_list, const std::vector<PPL::Generator>& ray_command_list);

void CtcPaveForward(Pave *p, bool inclusion, bool inner);
void CtcPaveBackward(Pave *p, bool inclusion, bool inner);
void CtcPaveConsistency(Pave *p, bool backward, bool inner);
bool CtcContinuity(Pave *p, bool backward);


#endif // UTILS_H
