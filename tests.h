#ifndef TESTS_H
#define TESTS_H

#include "pave.h"

void test_draw(Pave *p, std::string drawing_name="test", bool full=true);
void test_display_flow(Pave *p);

void testTranslate();
void testRotate();

void test_CtcPropagateLeftSide();
void test_CtcPropagateRightSide();
void test_CtcPropagateFront();

void test_CtcPropagateSegment();

void test_CtcPaveForward();
void test_CtcConsistency();
void test_CtcConsistency2();
void test_CtcConsistency3();

void test_CtcConsistency_Kernel();
void test_CtcConsistency_Kernel2();
void test_CtcConsistency_Kernel3();
void test_CtcConsistency_Kernel4();

void test_Newton();
void test_rotation();

void test_diff();

void test_copy_graph();
void test_contractor_polar();

void test_imageIntegral();

void sandbox();

void test_car_on_hill();

void test_inter_pave_perimeter();
void test_possible_path();

void test_infinity();
void test_chi_function();

void test_diff_infinity();

void test_contractor();
void test_contractor2();

//void ctcNorm(ibex::IntervalVector &x, ibex::IntervalVector &u, ibex::Interval& d);
void test_contractor_n();

#endif // TESTS_H
