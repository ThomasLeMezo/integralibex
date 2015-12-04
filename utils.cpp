#include "utils.h"

using namespace std;
using namespace ibex;

Utils::Utils()
{

}

std::vector<ibex::Interval> Utils::rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y){
    Interval xR = cos(theta)*x -sin(theta)*y;
    Interval yR = sin(theta)*x + cos(theta)*y;
    vector<Interval> list;
    list.push_back(xR);
    list.push_back(yR);
    return list;
}

/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Utils::CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy){
    Interval X(0.0, dx);

    Interval x = Interval(dy);
    Interval y = Interval(-dx, dx);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = theta;

    contract_polar.contract(x, y, rho, theta2);
    Sk = (Sk + x ) & X;
}

void Utils::CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateFront(Sk, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dy){
    Interval x = Sk;
    Interval y = Interval(0.0, dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = (Interval::PI - theta) & (-Interval::ZERO | Interval::HALF_PI);

    this->contract_polar.contract(x, y, rho, theta2);

    Sk = y;
}

void Utils::CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateLeftSide(Sk, theta, box[1].ub());
}

void Utils::CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> pi -theta
    */

    this->CtcPropagateLeftSide(Sk, Interval::PI-theta, dy);
}

void Utils::CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateRightSide(Sk, theta, box[0].ub(), box[1].ub());
}

/**
 * @brief Utils::rotateSegment
 * @param Sk
 * @param theta
 * @param box
 *
 * Rotate a Segment from the center of the "box"
 * The box is supposed (down,left) centered to (0,0)
 */
void Utils::rotate_segment_and_box(ibex::IntervalVector &Sk, const ibex::Interval &theta, IntervalVector &box, bool modifyBox){
    IntervalVector Sk_(Sk);
    IntervalVector box_(box);
    IntervalVector box_2(2);
    Vector center = box.mid();

    Sk_ -= center;
    box_ -= center;

    Sk[0] = cos(theta)*Sk_[0] -sin(theta)*Sk_[1];
    Sk[1] = sin(theta)*Sk_[0] + cos(theta)*Sk_[1];

    box_2[0] = cos(theta)*box_[0] -sin(theta)*box_[1];
    box_2[1] = sin(theta)*box_[0] + cos(theta)*box_[1];

    Vector lb = box_2.lb();
    Sk -= lb;
    box_2 -= lb;

    if(modifyBox)
        box=box_2;

//    Sk = Sk & box;
}

/**
 * @brief Utils::translateSegment
 * @param Sk
 * @param box
 * @param toZero
 *
 * Translate a segment such as the box is down/left centered on (0,0) when toZero == true
 * Backward translation when toZero == false
 */
void Utils::translate_segment_and_box(ibex::IntervalVector &Sk, IntervalVector &box, bool toZero, bool modifyBox){
    Vector lb = box.lb();

    if(toZero){
        Sk -= lb;
        if(modifyBox)
            box -= lb;
    }
    else{
        Sk += lb;
        if(modifyBox)
            box += lb;
    }
}

