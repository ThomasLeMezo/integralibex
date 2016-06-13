#ifndef UTILS_H
#define UTILS_H

#include "ibex.h"
#include "pave.h"
#include "imageintegral.h"
#include "ibex-robotics/ibex_CtcPolar.h"

class Utils
{
public:

    // ***********************************************************
    // ************************ Functions ************************
    // ***********************************************************

    Utils();
    ~Utils();

    void rotate_segment_and_box(ibex::IntervalVector &Sk, const ibex::Interval &theta, ibex::IntervalVector &box, bool modifyBox);
    void translate_segment_and_box(ibex::IntervalVector &Sk, ibex::IntervalVector &box, bool toZero, bool modifyBox);
    ibex::IntervalVector segment2IntervalVector(const ibex::Interval &seg, const int &face, const ibex::IntervalVector &box);

    std::vector<ibex::IntervalVector> get_segment_from_box(const ibex::IntervalVector &box, const double size_border=0.0);

    // Contractor Outer
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);
    void CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);

    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);
    void CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);

    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);
    void CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box, bool inner=false, const ibex::Interval &rotation=ibex::Interval::EMPTY_SET);

    void CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box_pave, bool inner=false, Border *b=NULL);

    void CtcPaveForward(Pave *p, bool inclusion, std::vector<bool> &change_tab);
    void CtcPaveBackward(Pave *p, bool inclusion, std::vector<bool> &change_tab);
    void CtcConsistency(Pave *p, bool backward, std::vector<bool> &change_tab, bool enable_function_iteration=true);
    bool CtcContinuity(Pave *p, bool backward);

    void CtcPolarCorrection(ibex::Interval &x, ibex::Interval &y, ibex::Interval &rho, ibex::Interval &theta);
    ibex::CtcPolar contract_polar;
    std::vector<ibex::IntervalVector> diff(const ibex::IntervalVector &box_initial, const ibex::IntervalVector &box_remove);
    bool test_discontinuity(const ibex::Interval &theta1, const ibex::Interval &theta2, const ibex::Interval &rotation);

    // ***********************************************************
    // ************************ Variables ************************
    // ***********************************************************

    ibex::Interval tab_rotation[4] = {ibex::Interval::ZERO, -ibex::Interval::HALF_PI, ibex::Interval::PI, ibex::Interval::HALF_PI};

    imageIntegral* m_imageIntegral;
    bool           m_imageIntegral_activated;

};

#endif // UTILS_H
