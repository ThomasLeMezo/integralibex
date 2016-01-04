#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"
#include "pave.h"


class Utils
{
public:

    // ***********************************************************
    // ************************ Functions ************************
    // ***********************************************************

    Utils();
    ~Utils();

    std::vector<ibex::Interval> rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y);
    void rotate_segment_and_box(ibex::IntervalVector &Sk, const ibex::Interval &theta, ibex::IntervalVector &box, bool modifyBox);
    void translate_segment_and_box(ibex::IntervalVector &Sk, ibex::IntervalVector &box, bool toZero, bool modifyBox);
    ibex::IntervalVector segment2IntervalVector(const ibex::Interval &seg, const int &face, const ibex::IntervalVector &box);

    // Contractor
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const ibex::Interval theta[], const ibex::IntervalVector &box_pave);

    ibex::CtcPolar contract_polar;
    ibex::CtcNewton* contract_newton;
    ibex::Function* vector_field_function;

    vector<bool> CtcPaveForward(Pave *p);
    vector<bool> CtcPaveBackward(Pave *p);
    vector<bool> CtcPaveConsistency(Pave *p);
    bool CtcContinuity(Pave *p);

//    void CtcPaveFlow(Pave *p);
    bool CtcNetwonPave(Pave *p);

    // ***********************************************************
    // ************************ Variables ************************
    // ***********************************************************

    ibex::Interval tab_rotation[4] = {ibex::Interval::ZERO, -ibex::Interval::HALF_PI, ibex::Interval::PI, ibex::Interval::HALF_PI};

};

#endif // UTILS_H
