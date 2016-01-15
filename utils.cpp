#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

Utils::Utils()
{
}

Utils::~Utils(){
}

// ********************************************************************************
// ****************** Contractors functions ***************************************
/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy){
    if(x_front.is_empty()){
        x=Interval::EMPTY_SET;
        return;
    }

    Interval X = Interval(0.0, dx);

    Interval Dx = Interval(-dx, dx);
    Interval Dy = Interval(dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = theta;

    contract_polar.contract(Dx, Dy, rho, theta2);
    x_front &= (x + Dx ) & X;

    if(!(X & (x_front - Dx)).is_empty()){
        x &= (x_front - Dx);
    }
    else{
        x |= (x_front - Dx);
    }

    if(!(X & (x_front + Dx)).is_empty()){
        x &= (x_front + Dx);
    }
    else{
        x |= (x_front + Dx);
    }

    x &= X;
}

void Utils::CtcPropagateFrontInner(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u){
    if(x_front.is_empty()){
        x=Interval::EMPTY_SET;
        return;
    }

    Interval X = Interval(0.0, dx);

    Interval Dx = Interval(-dx, dx);
    Interval Dy = Interval(dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = theta + u;

    contract_polar.contract(Dx, Dy, rho, theta2);
    x_front &= ((Interval(x.lb()) + Dx) & (Interval(x.ub()) + Dx)) & X;

    // Bwd : contraction des bornes sup et inf : idem
    x = ((Interval(x_front.lb()) + Dx) & (Interval(x_front.ub()) + Dx)) & X;
}

void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateFront(x, x_front, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateFrontInner(ibex::Interval &x, ibex::Interval &x_front, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u){
    this->CtcPropagateFrontInner(x, x_front, theta, box[0].ub(), box[1].ub(), u);
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy){
    x = x & Interval(0.0, dx);
    y = y & Interval(0.0, dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = (Interval::PI - theta);

    this->contract_polar.contract(x, y, rho, theta2);
}

void Utils::CtcPropagateLeftSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u){
    Interval x_lb = Interval(x.lb()) & Interval(0.0, dx);
    Interval y_lb = y & Interval(0.0, dy);
    Interval rho_lb = Interval::POS_REALS;
    Interval theta2_lb = (Interval::PI - theta);
    this->contract_polar.contract(x_lb, y_lb, rho_lb, theta2_lb);

    Interval x_ub = Interval(x.ub()) & Interval(0.0, dx);
    Interval y_ub = y & Interval(0.0, dy);
    Interval rho_ub = Interval::POS_REALS;
    Interval theta2_ub = (Interval::PI - theta);
    this->contract_polar.contract(x_ub, y_ub, rho_ub, theta2_ub);

    x &= (x_lb & x_ub);

    x_lb = x;
    y_lb = Interval(y.lb()) & Interval(0.0, dy);
    rho_lb = Interval::POS_REALS;
    theta2_lb = (Interval::PI - theta);
    this->contract_polar.contract(x_lb, y_lb, rho_lb, theta2_lb);

    x_ub = x;
    y_ub = Interval(y.ub()) & Interval(0.0, dy);
    rho_ub = Interval::POS_REALS;
    theta2_ub = (Interval::PI - theta);
    this->contract_polar.contract(x_ub, y_ub, rho_ub, theta2_ub);

    y &= (y_lb & y_ub);
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateLeftSide(x, y, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateLeftSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u){
    this->CtcPropagateLeftSideInner(x, y, theta, box[0].ub(), box[1].ub(), u);
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> pi - theta
     ** x -> dx - x
    */

    x = Interval(dx) - x;
    Interval theta2(Interval::PI-theta);
    this->CtcPropagateLeftSide(x, y, theta2, dx, dy);
    x = Interval(dx) - x;
}

void Utils::CtcPropagateRightSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const double &dx, const double &dy, const ibex::Interval &u){
    x = Interval(dx) - x;
    Interval theta2(Interval::PI-theta);
    Interval u2(Interval::PI-u);
    this->CtcPropagateLeftSideInner(x, y, theta2, dx, dy, u);
    x = Interval(dx) - x;
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const IntervalVector &box){
    this->CtcPropagateRightSide(x, y, theta, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateRightSideInner(ibex::Interval &x, ibex::Interval &y, const ibex::Interval &theta, const ibex::IntervalVector &box, const ibex::Interval &u){
    this->CtcPropagateRightSideInner(x, y, theta, box[0].ub(), box[1].ub(), u);
}

void Utils::CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> theta, const ibex::IntervalVector &box_pave){
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

    // **************************** OUTPUT ****************************
    // Translate and rotate back the Segment
    IntervalVector segment_contracted_out[3] = IntervalVector(2);
    segment_contracted_out[0][1] = segment_norm_out[0][0] | segment_norm_out[0][1]; segment_contracted_out[0][0] = Interval(box[0].ub());
    segment_contracted_out[1][0] = segment_norm_out[1][0] | segment_norm_out[1][1]; segment_contracted_out[1][1] = Interval(box[1].ub());
    segment_contracted_out[2][1] = segment_norm_out[2][0] | segment_norm_out[2][1]; segment_contracted_out[2][0] = Interval(box[0].lb());

    // Rotate and translate back with the initial box
    vector<Interval> seg_out_tmp;
    for(int i=0; i<3; i++){
        this->rotate_segment_and_box(segment_contracted_out[i], -tab_rotation[face], box, false);
        this->translate_segment_and_box(segment_contracted_out[i], box_in, false, false);
        // Add segment to seg_out list
        seg_out_tmp.push_back( (seg_out[i] & ((segment_contracted_out[i][0].diam() > segment_contracted_out[i][1].diam()) ? segment_contracted_out[i][0] : segment_contracted_out[i][1])));
    }
    seg_out = seg_out_tmp;

    // **************************** INPUT ****************************
    IntervalVector segment_contracted_in = IntervalVector(2);
    segment_contracted_in[0] = Interval::EMPTY_SET;
    segment_contracted_in[1] = Interval(box[1].lb());
    for(int i=0; i<3; i++){
        for(int j=0; j<2; j++){
            segment_contracted_in[0] |= segment_norm_in[i][j];
        }
    }
    // Rotate and translate back with the initial box
    this->rotate_segment_and_box(segment_contracted_in, -tab_rotation[face], box, false);
    this->translate_segment_and_box(segment_contracted_in, box_in, false, false);
    seg_in &= ((segment_contracted_in[0].diam() > segment_contracted_in[1].diam()) ? segment_contracted_in[0] : segment_contracted_in[1]);
}

void Utils::CtcPropagateSegmentInner(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> theta, const ibex::IntervalVector &box_pave, const ibex::Interval &u){
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
        this->CtcPropagateRightSideInner(segment_norm_in[0][i], segment_norm_out[0][i], theta[i] + tab_rotation[face], box, u);
        this->CtcPropagateFrontInner(segment_norm_in[1][i], segment_norm_out[1][i], theta[i] + tab_rotation[face], box, u);
        this->CtcPropagateLeftSideInner(segment_norm_in[2][i], segment_norm_out[2][i], theta[i] + tab_rotation[face], box, u);
    }

    // **************************** OUTPUT ****************************
    // Translate and rotate back the Segment
    IntervalVector segment_contracted_out[3] = IntervalVector(2);
    segment_contracted_out[0][1] = segment_norm_out[0][0] | segment_norm_out[0][1]; segment_contracted_out[0][0] = Interval(box[0].ub());
    segment_contracted_out[1][0] = segment_norm_out[1][0] | segment_norm_out[1][1]; segment_contracted_out[1][1] = Interval(box[1].ub());
    segment_contracted_out[2][1] = segment_norm_out[2][0] | segment_norm_out[2][1]; segment_contracted_out[2][0] = Interval(box[0].lb());

    // Rotate and translate back with the initial box
    vector<Interval> seg_out_tmp;
    for(int i=0; i<3; i++){
        this->rotate_segment_and_box(segment_contracted_out[i], -tab_rotation[face], box, false);
        this->translate_segment_and_box(segment_contracted_out[i], box_in, false, false);
        // Add segment to seg_out list
        seg_out_tmp.push_back( (seg_out[i] & ((segment_contracted_out[i][0].diam() > segment_contracted_out[i][1].diam()) ? segment_contracted_out[i][0] : segment_contracted_out[i][1])));
    }
    seg_out = seg_out_tmp;

    // **************************** INPUT ****************************
    IntervalVector segment_contracted_in = IntervalVector(2);
    segment_contracted_in[0] = Interval::EMPTY_SET;
    segment_contracted_in[1] = Interval(box[1].lb());
    for(int i=0; i<3; i++){
        for(int j=0; j<2; j++){
            segment_contracted_in[0] |= segment_norm_in[i][j];
        }
    }
    // Rotate and translate back with the initial box
    this->rotate_segment_and_box(segment_contracted_in, -tab_rotation[face], box, false);
    this->translate_segment_and_box(segment_contracted_in, box_in, false, false);
    seg_in &= ((segment_contracted_in[0].diam() > segment_contracted_in[1].diam()) ? segment_contracted_in[0] : segment_contracted_in[1]);
}

// ********************************************************************************
// ****************** Transformation functions ************************************

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

    if(Sk[0].is_empty() || Sk[1].is_empty()){
        Sk_[0] = Interval::EMPTY_SET;
        Sk_[1] = Interval::EMPTY_SET;
    }

    Sk_ -= center;
    box_ -= center;

    // 2nd method
    if(theta == Interval::ZERO){
        Sk[0] = Sk_[0];
        Sk[1] = Sk_[1];
        box_2[0] = box_[0];
        box_2[1] = box_[1];
    }
    else if(theta == Interval::HALF_PI){
        Sk[0] = -Sk_[1];
        Sk[1] = Sk_[0];
        box_2[0] = -box_[1];
        box_2[1] = box_[0];
    }
    else if(theta == Interval::PI){
        Sk[0] = -Sk_[0];
        Sk[1] = -Sk_[1];
        box_2[0] = -box_[0];
        box_2[1] = -box_[1];
    }
    else if(theta == -Interval::HALF_PI){
        Sk[0] = Sk_[1];
        Sk[1] = -Sk_[0];
        box_2[0] = box_[1];
        box_2[1] = -box_[0];
    }
    else{
        Sk[0] = cos(theta)*Sk_[0] - sin(theta)*Sk_[1];
        Sk[1] = sin(theta)*Sk_[0] + cos(theta)*Sk_[1];

        box_2[0] = cos(theta)*box_[0] - sin(theta)*box_[1];
        box_2[1] = sin(theta)*box_[0] + cos(theta)*box_[1];
    }

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

void Utils::CtcPaveBackward(Pave *p, bool inclusion, bool inner){

    Interval segment_in[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        Interval seg_in = p->get_border(face)->get_segment_in();

        vector<Interval> seg_out;
        for(int j=(face+1)%4; j!=face; j=(j+1)%4){
            seg_out.push_back(p->get_border(j)->get_segment_out());
        }

        if(!inner)
            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position());
        else
            this->CtcPropagateSegmentInner(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u());

        segment_in[face] = seg_in;
    }

    for(int face = 0; face<4; face++){
        p->get_border(face)->set_segment_in(segment_in[face], inclusion);
    }
}

void Utils::CtcPaveForward(Pave *p, bool inclusion, bool inner){
    Interval segment_out[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        Interval seg_in = p->get_border(face)->get_segment_in();

        vector<Interval> seg_out;
        for(int j=0; j<3; j++){
            seg_out.push_back(Interval::ALL_REALS);
        }

        if(!inner)
            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position());
        else
            this->CtcPropagateSegmentInner(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u());

        int k=0;
        for(int i=(face+1)%4; i!=face; i=(i+1)%4){
            segment_out[i] |= seg_out[k];
            k++;
        }
    }

    for(int face = 0; face<4; face++){
        p->get_border(face)->set_segment_out(segment_out[face], inclusion);
    }
}

void Utils::CtcPaveConsistency(Pave *p, bool backward, bool inner){
    if(backward){
        this->CtcPaveBackward(p, backward, inner);
        Pave p2(p);
        this->CtcPaveForward(&p2, backward, inner);
        *p &= p2;
    }
    else{
        this->CtcPaveForward(p, backward, inner);
    }
}

bool Utils::CtcContinuity(Pave *p, bool backward){
    bool change = false;

    for(int face = 0; face < 4; face++){
        Interval segment_in = Interval::EMPTY_SET;
        Interval segment_out = Interval::EMPTY_SET;

        for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
            segment_in |= p->get_border(face)->get_inclusion(b).get_segment_in();
            segment_out |= p->get_border(face)->get_inclusion(b).get_segment_out();
        }

        if(backward){
            if(p->get_border(face)->get_segment_in() != (segment_out & p->get_border(face)->get_segment_in())
                    || p->get_border(face)->get_segment_out() != (segment_in & p->get_border(face)->get_segment_out())){
                change = true;
                p->get_border(face)->set_segment_in(segment_out, backward);
                p->get_border(face)->set_segment_out(segment_in, backward);
            }
        }
        else{
            if(p->get_border(face)->get_segment_in() != (p->get_border(face)->get_segment_in() | segment_out & p->get_border(face)->get_segment_full())){
                change = true;
                p->get_border(face)->set_segment_in(segment_out, backward);
            }
        }

    }

    return change;
}

