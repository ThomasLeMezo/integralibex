#ifndef TESTS_H
#define TESTS_H

#include "pave.h"

void test_draw(Pave *p, std::string drawing_name="test", bool full=false);
void test_display_flow(Pave *p);

void testTranslate();
void testRotate();

void test_CtcPropagateLeftSide();
void test_CtcPropagateRightSide();
void test_CtcPropagateFront();

void test_CtcPropagateSegment();

void test_CtcPaveForward();
void test_CtcPaveConsistency();

void test_Newton();
void test_rotation();

void test_diff();

void test_copy_graph();
void test_contractor_polar();

void sandbox();

#endif // TESTS_H

