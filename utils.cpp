#include "utils.h"

using namespace std;
using namespace ibex;

Utils::Utils()
{
}

// ********************************************************************************
// ****************** Contractors functions ***************************************
/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy){
    Interval X = x & Interval(0.0, dx);

    Interval Dx = Interval(-dx, dx);
    Interval Dy = Interval(dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = theta;

    contract_polar.contract(Dx, Dy, rho, theta2);
    x_front = (x + Dx ) & X;    // A écrire sous forme de contracteur !!
    x = (x_front - Dx) & X;
}

void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateFront(x, x_front, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy){
    x = x & Interval(0.0, dx);
    y = y & Interval(0.0, dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = (Interval::PI - theta);

    this->contract_polar.contract(x, y, rho, theta2);
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateLeftSide(x, y, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> pi - theta
     ** x -> dx - x
    */

    x = Interval(dx) - x;
    this->CtcPropagateLeftSide(x, y, Interval::PI-theta, dx, dy);
    x = Interval(dx) - x;
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateRightSide(x, y, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const ibex::Interval theta[], const ibex::IntervalVector &box_pave){
    // Translate and rotate the Segment
    IntervalVector box(box_pave);
    IntervalVector box_in(box_pave);
    IntervalVector segment_in = segment2IntervalVector(seg_in, face, box);
    IntervalVector segment_out[3] = IntervalVector(2);
    for(int i=0; i<3; i++){
        segment_out[i] = segment2IntervalVector(seg_out[i], (face+1+i)%4, box);
    }

    this->translate_segment_and_box(segment_in, box, true, true);
    IntervalVector box_translate(box);
    this->rotate_segment_and_box(segment_in, this->tab_rotation[face], box, true);

    for(int i=0; i<3; i++){
        this->translate_segment_and_box(segment_out[i], box_in, true, false);
        this->rotate_segment_and_box(segment_out[i], this->tab_rotation[face], box_translate, false);
    }

    // Compute the propagation
    Interval segment_norm_in[3][2], segment_norm_out[3][2];

    for(int i=0; i<3; i++){
        for(int j=0; j<2; j++){
            segment_norm_in[i][j] = segment_in[0];
            segment_norm_out[i][j] = (segment_out[i][0].diam() > segment_out[i][1].diam()) ? segment_out[i][0] : segment_out[i][1];
        }
    }

    for(int i=0; i<2; i++){
        this->CtcPropagateRightSide(segment_norm_in[0][i], segment_norm_out[0][i], theta[i] + tab_rotation[face], box);
        this->CtcPropagateFront(segment_norm_in[1][i], segment_norm_out[1][i], theta[i] + tab_rotation[face], box);
        this->CtcPropagateLeftSide(segment_norm_in[2][i], segment_norm_out[2][i], theta[i] + tab_rotation[face], box);
    }

    // Translate and rotate back the Segment
    IntervalVector segment_contracted_out[3] = IntervalVector(2);
    segment_contracted_out[0][1] = segment_norm_out[0][0] | segment_norm_out[0][1]; segment_contracted_out[0][0] = Interval(box[0].ub());
    segment_contracted_out[1][0] = segment_norm_out[1][0] | segment_norm_out[1][1]; segment_contracted_out[1][1] = Interval(box[1].ub());
    segment_contracted_out[2][1] = segment_norm_out[2][0] | segment_norm_out[2][1]; segment_contracted_out[2][0] = Interval(box[0].lb());

    // Rotate and translate back with the initial box
    seg_out.clear();
    for(int i=0; i<3; i++){
        this->rotate_segment_and_box(segment_contracted_out[i], -tab_rotation[face], box, false);
        this->translate_segment_and_box(segment_contracted_out[i], box_in, false, false);
        // Add segment to seg_out list
        seg_out.push_back( (segment_contracted_out[i][0].diam() > segment_contracted_out[i][1].diam()) ? segment_contracted_out[i][0] : segment_contracted_out[i][1] );
    }

    // Segment in (backward)
    Interval segment_contracted_in[2] = Interval::EMPTY_SET;
    for(int j=0; j<2; j++){ // theta
        for(int i=0; i<3; i++){ // right - front - left
            segment_contracted_in[j] = segment_contracted_in[j] | segment_norm_in[i][j];
            // Union car, ce n'est pas pcq segment_norm_in peut être vide à cause de theta et non
            // des conditions des frontières !
            // Mais pas asser efficace !!
            // A REGARDER !!!
        }
    }

    segment_in = segment2IntervalVector(segment_contracted_in[0] | segment_contracted_in[1], face, box);
    this->rotate_segment_and_box(segment_in, -tab_rotation[face], box, false);
    this->translate_segment_and_box(segment_in, box_in, false, false);
    seg_in = segment_in[face%2];
}

// ********************************************************************************
// ****************** Transformation functions ************************************

std::vector<ibex::Interval> Utils::rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y){
    Interval xR = cos(theta)*x -sin(theta)*y;
    Interval yR = sin(theta)*x + cos(theta)*y;
    vector<Interval> list;
    list.push_back(xR);
    list.push_back(yR);
    return list;
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

ibex::IntervalVector Utils::segment2IntervalVector(const ibex::Interval &seg, const int &face, const ibex::IntervalVector &box){
    IntervalVector intervalVectorSegment(2);
    intervalVectorSegment[face%2] = seg;
    intervalVectorSegment[(face+1)%2] = (((face == 1) | (face == 2) ) ? Interval(box[(face+1)%2].ub()) : Interval(box[(face+1)%2].lb()));
    return intervalVectorSegment;
}





