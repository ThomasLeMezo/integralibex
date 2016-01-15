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

    // Contractor Outer
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> theta, const ibex::IntervalVector &box_pave);

    void CtcPaveForward(Pave *p, bool inclusion, bool inner);
    void CtcPaveBackward(Pave *p, bool inclusion, bool inner);
    void CtcPaveConsistency(Pave *p, bool backward, bool inner);
    bool CtcContinuity(Pave *p, bool backward);

    // Contractor Inner
    void CtcPropagateLeftSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u, bool final=false);
    void CtcPropagateLeftSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u, bool final=false);

    void CtcPropagateFrontInner(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u);
    void CtcPropagateFrontInner(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u);

    void CtcPropagateRightSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u, bool final=false);
    void CtcPropagateRightSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u, bool final=false);

    void CtcPropagateSegmentInner(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> theta, const ibex::IntervalVector &box_pave, const ibex::Interval &u);

    ibex::CtcPolar contract_polar;

    // ***********************************************************
    // ************************ Variables ************************
    // ***********************************************************

    ibex::Interval tab_rotation[4] = {ibex::Interval::ZERO, -ibex::Interval::HALF_PI, ibex::Interval::PI, ibex::Interval::HALF_PI};

};

#endif // UTILS_H
