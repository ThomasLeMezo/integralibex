#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"


class Utils
{
public:
    Utils();

    std::vector<ibex::Interval> rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y);
    void rotate_segment_and_box(ibex::IntervalVector &Sk, const double &theta, ibex::IntervalVector &box);
    void translate_segment_and_box(ibex::IntervalVector &Sk, ibex::IntervalVector &box, bool toZero, bool modifyBox);


    // Contractor
    void CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dy);
    void CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const ibex::IntervalVector &box);

    void CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy);
    void CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const ibex::IntervalVector &box);

    ibex::CtcPolar contract_polar;

};

#endif // UTILS_H
