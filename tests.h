#ifndef TESTS_H
#define TESTS_H

#include "pave.h"

void test_draw(Pave *p, std::string drawing_name="test");

void testTranslate();
void testRotate();

void test_CtcPropagateLeftSide();
void test_CtcPropagateRightSide();
void test_CtcPropagateFront();

void test_CtcPaveForward();
void test_CtcPaveBackward();

void test_Newton();

#endif // TESTS_H

