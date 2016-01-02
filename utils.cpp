#include "utils.h"
#include "pave.h"

using namespace std;
using namespace ibex;

Utils::Utils()
{
    Variable x, y;
    vector_field_function = new Function(x, y, Return(y, 1.0*(1-pow(x, 2))*y-x));
    contract_newton = new CtcNewton(*vector_field_function, 10.0,10.0);
}

Utils::~Utils(){
    delete vector_field_function;
    delete contract_newton;
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
    Interval theta2(Interval::PI-theta);
    this->CtcPropagateLeftSide(x, y, theta2, dx, dy);
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

    // **************************** OUTPUT ****************************
    // Translate and rotate back the Segment
    IntervalVector segment_contracted_out[3] = IntervalVector(2);
    segment_contracted_out[0][1] = segment_norm_out[0][0] | segment_norm_out[0][1]; segment_contracted_out[0][0] = Interval(box[0].ub());
    segment_contracted_out[1][0] = segment_norm_out[1][0] | segment_norm_out[1][1]; segment_contracted_out[1][1] = Interval(box[1].ub());
    segment_contracted_out[2][1] = segment_norm_out[2][0] | segment_norm_out[2][1]; segment_contracted_out[2][0] = Interval(box[0].lb());

    // Rotate and translate back with the initial box
    seg_out.clear();
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

    Sk_ -= center;
    box_ -= center;

    Sk[0] = cos(theta)*Sk_[0] -sin(theta)*Sk_[1];
    Sk[1] = sin(theta)*Sk_[0] + cos(theta)*Sk_[1];

    box_2[0] = cos(theta)*box_[0] - sin(theta)*box_[1];
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

vector<bool> Utils::CtcPaveBackward(Pave *p){

    Interval segment_in[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        if(p->m_borders[face].flow_in == true){

            Interval seg_in = p->m_borders[face].segment;

            vector<Interval> seg_out;
            for(int j=(face+1)%4; j!=face; j=(j+1)%4){
                seg_out.push_back(p->m_borders[j].segment);
            }

            this->CtcPropagateSegment(seg_in, seg_out, face, p->m_theta, p->m_box);

            segment_in[face] = seg_in;
        }
    }

    vector<bool> tab_change;
    tab_change.reserve(4);
    for(int i=0; i<4; i++){
        tab_change.push_back(false);
    }

    for(int face = 0; face<4; face++){
        if(p->m_borders[face].flow_in==true){
            if((p->m_borders[face].segment & segment_in[face]) != segment_in[face]){
                tab_change[face] = true;
            }
            p->m_borders[face].segment &= segment_in[face];
        }
    }

    Pave p2(*p);
    vector<bool> tab_change_forward = this->CtcPaveForward(&p2);
    *p &= p2;

    for(int face = 0; face<4; face++){
        tab_change[face] = tab_change[face] || tab_change_forward[face];
    }

    return tab_change;
}

vector<bool> Utils::CtcPaveForward(Pave *p){
    Interval segment_out[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        if(p->m_borders[face].flow_in == true){

            Interval seg_in = p->m_borders[face].segment;
            segment_out[face] |= seg_in; // usefull ?

            vector<Interval> seg_out;
            for(int j=0; j<3; j++){
                seg_out.push_back(Interval::ALL_REALS);
            }

            this->CtcPropagateSegment(seg_in, seg_out, face, p->m_theta, p->m_box);

            int k=0;
            for(int i=(face+1)%4; i!=face; i=(i+1)%4){
                segment_out[i] |= seg_out[k];
                k++;
            }

        }
    }

    vector<bool> tab_change;
    tab_change.reserve(4);
    for(int i=0; i<4; i++){
        tab_change.push_back(false);
    }

    for(int face = 0; face<4; face++){
        if(p->m_borders[face].segment != segment_out[face]){
            tab_change[face] = true;
            p->m_borders[face].segment = segment_out[face];
        }
    }

    return tab_change;
}

bool Utils::CtcContinuity(Pave *p){
    bool change = false;

    for(int face = 0; face < 4; face++){
        Interval segment = Interval::EMPTY_SET;
        for(int b = 0; b < p->m_borders[face].brothers.size(); b++){
            segment |= p->m_borders[face].brothers[b]->segment;
        }

        if(p->m_borders[face].segment != segment)
            change = true;
        p->m_borders[face].segment = segment;
    }

    return change;
}

bool Utils::CtcNetwonPave(Pave *p){
    return false;

    if(p->is_all_brothers_full()){
        IntervalVector box_tmp = p->m_box;
        this->contract_newton->contract(box_tmp);
        if(!box_tmp.is_empty())
            this->contract_newton->contract(box_tmp);

        if(!box_tmp.is_empty() && box_tmp[0].is_degenerated() && box_tmp[1].is_degenerated()){
            cout << "-->" << p->m_box << "<>" << box_tmp << endl;

            for(int face=0; face<4; face++){
                p->m_borders[face].segment=Interval::EMPTY_SET;
            }

            return true;
        }
    }
    return false;
}

void Utils::CtcPaveFlow(Pave *p){
    bool flow_in[4] = {false, false, false, false};
    bool flow_out[4] = {false, false, false, false};

    for(int face = 0; face < 4; face++){
//        Interval seg_in = p->borders[face].position[face%2];
        Interval seg_in = Interval::ALL_REALS;
        vector<Interval> seg_out;
        for(int i=0; i<3; i++){
            seg_out.push_back(p->m_borders[(face+i+1)%4].segment);
//            seg_out.push_back(Interval::ALL_REALS);
        }
        this->CtcPropagateSegment(seg_in, seg_out, face, p->m_theta, p->m_box);

        if(!seg_out[0].is_degenerated() || !seg_out[1].is_degenerated() || !seg_out[2].is_degenerated()){
            flow_in[face] = true;
        }
        for(int i = 0; i<3; i++){
            if(!seg_out[i].is_degenerated()){
                flow_out[(face+i+1)%4] = true;
                p->m_borders[(face+i+1)%4].flow_out[face] = true;
            }
            else{
                p->m_borders[(face+i+1)%4].flow_out[face] = false;
            }
        }
    }

    for(int face = 0; face <4; face++){
        p->m_borders[face].flow_in = flow_in[face];
    }
}

