#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

Utils::Utils()
{
    m_imageIntegral_activated = false;
}

Utils::~Utils(){
    if(m_imageIntegral_activated)
        delete(m_imageIntegral);
}

// ********************************************************************************
// ****************** Contractors functions ***************************************
/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner){
    if(x_front.is_empty() || x.is_empty() || theta_list.size()==0){
        x = Interval::EMPTY_SET;
        x_front = Interval::EMPTY_SET;
        return;
    }
    Interval X = Interval(0.0, dx);
    std::vector<ibex::Interval> x_list, x_front_list;

    for(auto &theta:theta_list){
        Interval Dx = Interval(-dx, dx);
        Interval Dy = Interval(dy);
        Interval rho = Interval::POS_REALS;
        Interval theta2(theta);

        contract_polar.contract(Dx, Dy, rho, theta2);

        ibex::Interval x_front_tmp =  (x + Dx) & X & x_front;
        x_front_list.push_back(x_front_tmp);
        if(Dx.is_empty() && Dx.is_empty())
            x_list.push_back(Interval::EMPTY_SET);
        else
            x_list.push_back(x & Interval(x_front_tmp.lb()-Dx.lb(), x_front_tmp.ub()-Dx.ub()) | Interval(x_front_tmp.lb()-Dx.ub(), x_front_tmp.ub()-Dx.lb()));
    }

    Interval x_out(Interval::EMPTY_SET), x_front_out(Interval::EMPTY_SET);
    for(int i=0; i<x_front_list.size(); i++){
        x_out |= x_list[i];
        x_front_out |= x_front_list[i];
    }
    //    x_out &= X & x;
    //    x_front_out &= X & x_front;

    if(inner){
                Interval x_out_inner(Interval::ALL_REALS);
                for(int i=0; i<x_list.size(); i++){
                    if(x_list.size()>i+1){
                        if((x_list[i] & x_list[i+1]).is_degenerated()){
                            x_list[i] += x_list[i+1];
                            x_list.erase(x_list.begin()+i+1);
                            i--;
                        }
                    }
                }
                for(int i=0; i<x_list.size(); i++){
                    if(!x_list[i].is_empty()){
                        x_out_inner &= x_list[i];
                    }
                }
                x_out_inner &= X & x;
                if(!x_out_inner.is_empty())
                    x_out = x_out_inner;
    }

    x = x_out;
    x_front = x_front_out;

}

void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box, bool inner){
    this->CtcPropagateFront(x, x_front, theta_list, box[0].ub(), box[1].ub(), inner);
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner){
    if(x.is_empty() || y.is_empty() || theta_list.size()==0){
        x = Interval::EMPTY_SET;
        y = Interval::EMPTY_SET;
        return;
    }

    y &= Interval(0.0, dy);
    x &= Interval(0.0, dx);
    vector<Interval> theta2_list, x_out, y_out;
    for(auto &theta:theta_list){
        theta2_list.push_back(Interval::PI - theta);
    }

    for(auto &theta:theta2_list){
        Interval x_tmp(x), y_tmp(y);
        Interval rho(Interval::POS_REALS);
        this->contract_polar.contract(x_tmp, y_tmp, rho, theta);
        x_out.push_back(x_tmp);
        y_out.push_back(y_tmp);
    }

    x = Interval::EMPTY_SET; y = Interval::EMPTY_SET;
    for(int i=0; i<x_out.size(); i++){
        x |= x_out[i];
        y |= y_out[i];
    }

    if(inner){
                Interval x_inner(x);
                for(int i=0; i<x_out.size(); i++){
                    if(x_out.size()>i+1){
                        if((x_out[i] & x_out[i+1]).is_degenerated()){
                            x_out[i] += x_out[i+1];
                            x_out.erase(x_out.begin()+i+1);
                            i--;
                        }
                    }
                }
                for(int i=0; i<x_out.size(); i++){
                    x_inner &= x_out[i];
                }
                x_inner &= Interval(0.0, dx);

                if(!x_inner.is_empty())
                    x = x_inner;
    }

    //    x &= Interval(0.0, dx);
    //    y &= Interval(0.0, dy);
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box, bool inner){
    this->CtcPropagateLeftSide(x, y, theta_list, box[0].ub(), box[1].ub(), inner);
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy, bool inner){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> pi - theta
     ** x -> dx - x
    */

    x = Interval(dx) - x;
    vector<Interval> theta2_list;
    for(auto &theta:theta_list){
        theta2_list.push_back(Interval::PI - theta);
    }
    this->CtcPropagateLeftSide(x, y, theta2_list, dx, dy, inner);
    x = Interval(dx) - x;
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box, bool inner){
    this->CtcPropagateRightSide(x, y, theta_list, box[0].ub(), box[1].ub(), inner);
}

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void Utils::CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box_pave, bool inner){
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
    Interval segment_norm_in[3], segment_norm_out[3];

    for(int i=0; i<3; i++){
        segment_norm_in[i] = segment_in[0];
        segment_norm_out[i] = (segment_out[i][0].diam() > segment_out[i][1].diam()) ? segment_out[i][0] : segment_out[i][1];
    }

    std::vector<ibex::Interval> theta_list_rotate;
    for(auto &theta:theta_list){
        theta_list_rotate.push_back(theta + tab_rotation[face]);
    }

    this->CtcPropagateRightSide(segment_norm_in[0], segment_norm_out[0], theta_list_rotate, box, inner);
    this->CtcPropagateFront(segment_norm_in[1], segment_norm_out[1], theta_list_rotate, box, inner);
    this->CtcPropagateLeftSide(segment_norm_in[2], segment_norm_out[2], theta_list_rotate, box, inner);

    // **************************** OUTPUT ****************************
    // Translate and rotate back the Segment
    IntervalVector segment_contracted_out[3] = {IntervalVector(2, Interval::EMPTY_SET), IntervalVector(2, Interval::EMPTY_SET), IntervalVector(2, Interval::EMPTY_SET)};

    segment_contracted_out[0][1] |= segment_norm_out[0];
    segment_contracted_out[1][0] |= segment_norm_out[1];
    segment_contracted_out[2][1] |= segment_norm_out[2];

    segment_contracted_out[0][0] = Interval(box[0].ub());
    segment_contracted_out[1][1] = Interval(box[1].ub());
    segment_contracted_out[2][0] = Interval(box[0].lb());

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
        segment_contracted_in[0] |= segment_norm_in[i];
    }

    // Rotate and translate back with the initial box
    this->rotate_segment_and_box(segment_contracted_in, -tab_rotation[face], box, false);
    this->translate_segment_and_box(segment_contracted_in, box_in, false, false);
    seg_in &= ((segment_contracted_in[0].diam() > segment_contracted_in[1].diam()) ? segment_contracted_in[0] : segment_contracted_in[1]);
}

void Utils::CtcPaveBackward(Pave *p, bool inclusion, std::vector<bool> &change_tab, bool inner){
    for(int face = 0; face < 4; face++){
        Interval seg_in = p->get_border(face)->get_segment_in();

        vector<Interval> seg_out;
        for(int j=(face+1)%4; j!=face; j=(j+1)%4){
            seg_out.push_back(p->get_border(j)->get_segment_out());
        }

        this->CtcPropagateSegment(seg_in, seg_out, face, p->get_all_theta(), p->get_position(), inner);

        change_tab[face] = p->get_border(face)->set_segment_in(seg_in, inclusion) || change_tab[face];
    }
}

void Utils::CtcPaveForward(Pave *p, bool inclusion, std::vector<bool> &change_tab, bool inner){
    Interval segment_out[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        Interval seg_in;
        seg_in = p->get_border(face)->get_segment_in();

        vector<Interval> seg_out;
        for(int j=0; j<3; j++){
            seg_out.push_back(Interval::ALL_REALS);
        }

        this->CtcPropagateSegment(seg_in, seg_out, face, p->get_all_theta(), p->get_position(), inner);

        int k=0;
        for(int i=(face+1)%4; i!=face; i=(i+1)%4){
            segment_out[i] |= seg_out[k];
            k++;
        }
    }

    for(int face = 0; face<4; face++){
        change_tab[face] = p->get_border(face)->set_segment_out(segment_out[face], inclusion) || change_tab[face];
    }
}

// ********************************************************************************
// ****************** Algorithm functions      ************************************

void Utils::CtcPaveConsistency(Pave *p, bool backward, std::vector<bool> &change_tab, bool enable_function_iteration, bool inner){
//    int nb_f = p->get_f_list().size();
//    if(!enable_function_iteration)
//        nb_f=1;

//    for(int i=0; i<nb_f; i++){ //to reach fix point (more iteration might be necessary)
        if(backward){
            this->CtcPaveBackward(p, true, change_tab, inner);
            Pave *p2 = new Pave(p);
            this->CtcPaveForward(p2, true, change_tab, inner);
            *p &= *(p2);
        }
        else{
            this->CtcPaveForward(p, false, change_tab, inner);
        }
//    }

    // Test if only one border is not empty
    int nb_not_empty = 0;
    for(auto &b:p->get_borders()){
        if(!b->is_empty()){
            nb_not_empty++;
        }
    }
    if(nb_not_empty==1)
        p->set_empty();

    // Reduce impact of change when backward (mandatory)
    if(backward){
        for(int face = 0; face<4; face++){
            if(p->get_border(face)->get_segment_full() == (p->get_border(face)->get_segment_in() | p->get_border(face)->get_segment_out()))
                change_tab[face] = false;
        }
    }
}

bool Utils::CtcContinuity(Pave *p, bool backward){
    bool change = false;

    for(int face = 0; face < 4; face++){
        // **************** OUT CONTINUITY *************
        if(p->get_border(face)->get_continuity_out()){
            Interval segment_in = Interval::EMPTY_SET;

            for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
                segment_in |= p->get_border(face)->get_inclusion(b)->get_segment_in();
            }

            if(backward && (p->get_border(face)->get_segment_out() != (segment_in & p->get_border(face)->get_segment_out()))){
                change = true;
                p->get_border(face)->set_segment_out(segment_in, backward);
            }
        }

        // **************** IN CONTINUITY *************
        if(p->get_border(face)->get_continuity_in()){
            Interval segment_out = Interval::EMPTY_SET;

            for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
                segment_out |= p->get_border(face)->get_inclusion(b)->get_segment_out();
            }

            if(backward){
                if(p->get_border(face)->get_segment_in() != (segment_out & p->get_border(face)->get_segment_in())){
                    change = true;
                    p->get_border(face)->set_segment_in(segment_out, backward);
                }
            }
            else{
                if(p->get_border(face)->get_segment_in() != (p->get_border(face)->get_segment_in() | segment_out & p->get_border(face)->get_segment_full())){
                    change = true;
                    p->get_border(face)->set_segment_in(segment_out, backward);
                }
            }
        }
    }

    return change;
}

// ********************************************************************************
// ****************** Utils functions ************************************

void Utils::CtcPolarCorrection(Interval &x, Interval &y, Interval &rho, Interval &theta){
    if(x == Interval::ZERO || y == Interval::ZERO){
        Interval x_r, y_r, theta_r;
        x_r = sqrt(2)/2*(x - y);
        y_r = sqrt(2)/2*(x + y);
        theta_r += Interval::PI/4.0;
        contract_polar.contract(x_r, y_r, rho, theta_r);
        x &= sqrt(2)/2*(x_r + y_r);
        y &= sqrt(2)/2*(-x_r + y_r);
        theta -= Interval::PI/4.0;
    }
    else{
        contract_polar.contract(x, y, rho, theta);
    }
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
    else if(theta == Interval::PI || theta == -Interval::PI){
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

std::vector<IntervalVector> Utils::diff(const IntervalVector &box_initial, const IntervalVector &box_remove){
    std::vector<IntervalVector> box_result;

    IntervalVector box_diff = box_initial & box_remove;

    if(!box_diff.is_empty()){
        IntervalVector box(2);
        box[0] = Interval(box_initial[0].lb(), box_diff[0].lb());
        box[1] = Interval(box_diff[1].lb(), box_initial[1].ub());
        if(!box.is_flat())
            box_result.push_back(box);

        box[0] = Interval(box_initial[0].lb(), box_diff[0].ub());
        box[1] = Interval(box_initial[1].lb(), box_diff[1].lb());
        if(!box.is_flat())
            box_result.push_back(box);

        box[0] = Interval(box_diff[0].ub(), box_initial[0].ub());
        box[1] = Interval(box_initial[1].lb(), box_diff[1].ub());
        if(!box.is_flat())
            box_result.push_back(box);

        box[0] = Interval(box_diff[0].lb(), box_initial[0].ub());
        box[1] = Interval(box_diff[1].ub(), box_initial[1].ub());
        if(!box.is_flat())
            box_result.push_back(box);

    }
    else{
        box_result.push_back(box_initial);
    }

    return box_result;
}

std::vector<IntervalVector> Utils::get_segment_from_box(const IntervalVector &box, const double size_border){
    vector<IntervalVector> list_border;
    IntervalVector test(2);
    test[0] = Interval(box[0].lb(), box[0].ub());
    test[1] = Interval(box[1].lb()-size_border, box[1].lb());
    list_border.push_back(test);
    test[0] = Interval(box[0].ub(), box[0].ub()+size_border);
    test[1] = Interval(box[1].lb(), box[1].ub());
    list_border.push_back(test);
    test[0] = Interval(box[0].lb(), box[0].ub());
    test[1] = Interval(box[1].ub(), box[1].ub()+size_border);
    list_border.push_back(test);
    test[0] = Interval(box[0].lb()-size_border, box[0].lb());
    test[1] = Interval(box[1].lb(), box[1].ub());
    list_border.push_back(test);

    return list_border;
}

