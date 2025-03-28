#include "utils.h"
#include "pave.h"
#include "iomanip"
#include "vibes.h"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

// ********************************************************************************
// ****************** Contractors functions ***************************************
/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy){
    if(y.is_empty() || x.is_empty() || theta_list.size()==0){
        x = Interval::EMPTY_SET;
        y = Interval::EMPTY_SET;
        return;
    }
    Interval X = Interval(0.0, dx);
    std::vector<ibex::Interval> x_list, y_list;

    for(Interval theta:theta_list){
        Interval Dx = Interval(-dx, dx);
        Interval Dy = Interval(dy);
        Interval rho = Interval::POS_REALS;
        Interval theta2(theta);

        contract_polar.contract(Dx, Dy, rho, theta2);
        //                CtcPolarCorrection(Dx, Dy, rho, theta2);

        // Compute x_front
        y_list.push_back(x + Dx);

        // Compute X
        if(Dx.is_empty() && Dx.is_empty())
            x_list.push_back(Interval::EMPTY_SET);
        else{
            ibex::Interval y_tmp = (x + Dx) & X & y;
            x_list.push_back(Interval(y_tmp.lb()-Dx.lb(), y_tmp.ub()-Dx.ub()) | Interval(y_tmp.lb()-Dx.ub(), y_tmp.ub()-Dx.lb()));
        }
    }

    Interval x_out(Interval::EMPTY_SET), y_out(Interval::EMPTY_SET);
    for(int i=0; i<y_list.size(); i++){
        x_out |= x_list[i];
        y_out |= y_list[i];
    }

    x = x_out & x & X;
    y = y_out & y & X;
}

void Utils::CtcPropagateFront(ibex::Interval &x, ibex::Interval &x_front, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box){
    this->CtcPropagateFront(x, x_front, theta_list, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy){
    if(x.is_empty() || y.is_empty() || theta_list.size()==0){
        x = Interval::EMPTY_SET;
        y = Interval::EMPTY_SET;
        return;
    }

    y &= Interval(0.0, dy);
    x &= Interval(0.0, dx);
    vector<Interval> theta2_list, x_list, y_list;
    for(Interval theta:theta_list){
        theta2_list.push_back(Interval::PI - theta);
    }

    for(Interval theta:theta2_list){
        Interval x_tmp(x), y_tmp(y);
        Interval rho(Interval::POS_REALS);
        this->contract_polar.contract(x_tmp, y_tmp, rho, theta);
        //        CtcPolarCorrection(x_tmp, y_tmp, rho, theta);
        x_list.push_back(x_tmp);
        y_list.push_back(y_tmp);
    }

    Interval x_out(Interval::EMPTY_SET);
    Interval y_out(Interval::EMPTY_SET);
    for(int i=0; i<x_list.size(); i++){
        x_out |= x_list[i];
        y_out |= y_list[i];
    }

    x &= x_out;
    y &= y_out;
}

void Utils::CtcPropagateLeftSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box){
    this->CtcPropagateLeftSide(x, y, theta_list, box[0].ub(), box[1].ub());
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const double &dx, const double &dy){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> pi - theta
     ** x -> dx - x
    */

    x = Interval(dx) - x;
    vector<Interval> theta2_list;
    for(Interval theta:theta_list){
        theta2_list.push_back(Interval::PI - theta);
    }
    this->CtcPropagateLeftSide(x, y, theta2_list, dx, dy);
    x = Interval(dx) - x;
}

void Utils::CtcPropagateRightSide(ibex::Interval &x, ibex::Interval &y, const std::vector<ibex::Interval> &theta_list, const IntervalVector &box){
    this->CtcPropagateRightSide(x, y, theta_list, box[0].ub(), box[1].ub());
}

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void Utils::CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> &theta_list, const ibex::IntervalVector &box_pave){
    // Translate and rotate the Segment
    IntervalVector box(box_pave);
    IntervalVector box_in(box_pave);
    IntervalVector segment_in = segment2IntervalVector(seg_in, face, box);
    IntervalVector segment_out[3] = {IntervalVector(2), IntervalVector(2), IntervalVector(2)};
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
        segment_norm_out[i] = segment_out[i][(i+1)%2];
    }

    std::vector<ibex::Interval> theta_list_rotate;
    for(Interval theta:theta_list){
        theta_list_rotate.push_back(theta + tab_rotation[face]);
    }

    this->CtcPropagateRightSide(segment_norm_in[0], segment_norm_out[0], theta_list_rotate, box);
    this->CtcPropagateFront(segment_norm_in[1], segment_norm_out[1], theta_list_rotate, box);
    this->CtcPropagateLeftSide(segment_norm_in[2], segment_norm_out[2], theta_list_rotate, box);

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
        seg_out_tmp.push_back(seg_out[i] & segment_contracted_out[i][(face+i+1)%2]);
    }

    // **************************** INPUT ****************************

    IntervalVector segment_contracted_in = IntervalVector(2);
    segment_contracted_in[1] = Interval(box[1].lb());

    segment_contracted_in[0] = Interval::EMPTY_SET;
    for(int i=0; i<3; i++)
        segment_contracted_in[0] |= segment_norm_in[i];

    // Rotate and translate back with the initial box
    this->rotate_segment_and_box(segment_contracted_in, -tab_rotation[face], box, false);
    this->translate_segment_and_box(segment_contracted_in, box_in, false, false);

    // Write data
    seg_in &= segment_contracted_in[face%2];
    seg_out = seg_out_tmp;
}

void Utils::CtcPaveBackward(Pave *p, bool inclusion, std::vector<bool> &change_tab){
    for(int face = 0; face < 4; face++){
        Interval seg_in = p->get_border(face)->get_segment_in();
        if(!seg_in.is_empty()){

            vector<Interval> seg_out;
            for(int j=(face+1)%4; j!=face; j=(j+1)%4){
                seg_out.push_back(p->get_border(j)->get_segment_out());
            }

            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_all_theta(), p->get_position());

            change_tab[face] = p->get_border(face)->set_segment_in(seg_in, inclusion) || change_tab[face];
        }
    }
}

void Utils::CtcPaveForward(Pave *p, bool inclusion, std::vector<bool> &change_tab, bool union_functions){
    Interval segment_out[4] = {Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET, Interval::EMPTY_SET};

    for(int face = 0; face < 4; face++){
        Interval seg_in;
        seg_in = p->get_border(face)->get_segment_in();
        if(!seg_in.is_empty()){

            vector<Interval> seg_out;
            for(int j=0; j<3; j++){
                seg_out.push_back(Interval::ALL_REALS);
            }

            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_all_theta(union_functions), p->get_position());

            int k=0;
            for(int i=(face+1)%4; i!=face; i=(i+1)%4){
                segment_out[i] |= seg_out[k];
                k++;
            }
        }
    }

    for(int face = 0; face<4; face++){
        change_tab[face] = p->get_border(face)->set_segment_out(segment_out[face], inclusion) || change_tab[face];
    }
}

// ********************************************************************************
// ****************** Algorithm functions      ************************************

void Utils::CtcConsistency(Pave *p, bool backward, std::vector<bool> &change_tab, bool union_functions){

    if(backward){
        if(p->get_compute_inner() && p->get_inner_mode() && p->get_f_list().size()>1){
//        if(/*p->get_compute_inner() && p->get_inner_mode() &&*/ p->get_f_list().size()>1){
            // Case several cones
            vector<Pave*> list_pave;
            for(int i=0; i<p->get_f_list().size(); i++){
                Pave *p_tmp = new Pave(p);
                p_tmp->set_active_function(i);
                this->CtcPaveBackward(p_tmp, true, change_tab);
                list_pave.push_back(p_tmp);
            }
            // Intersect paves
            p->inter_inner(list_pave);

            // Delete pave tmp
            for(Pave *p_tmp:list_pave){
                delete(p_tmp);
            }
        }
        else{
            this->CtcPaveBackward(p, true, change_tab);
//            this->CtcPaveBackward2(p, true, change_tab);
        }
        Pave *p2 = new Pave(p);
        this->CtcPaveForward(p2, true, change_tab, union_functions); // Test ? union_functions
        *p &= *(p2);
        delete(p2);

    }
    else{
        this->CtcPaveForward(p, false, change_tab, union_functions);
    }

    // Test if only one border is not empty
    int nb_not_empty = 0;
    for(Border *b:p->get_borders()){
        if(!b->is_empty()){
            nb_not_empty++;
        }
    }
    if(nb_not_empty==1)
        p->set_empty();

    //     Reduce impact of change when backward (mandatory)
    //            if(backward && !p->get_inner_mode()){
    //                for(int face = 0; face<4; face++){
    //                    if((p->get_border(face)->get_segment_full() == (p->get_border(face)->get_segment_in() | p->get_border(face)->get_segment_out())))
    //                        change_tab[face] = false;
    //                }
    //            }
}

bool Utils::CtcContinuity(Pave *p, bool backward){
    bool change = false;

    for(int face = 0; face < 4; face++){
        // **************** OUT CONTINUITY *************
        //        if(p->get_border(face)->get_continuity_out()){
        if(backward){
            Interval segment_in = Interval::EMPTY_SET;

            for(int b = 0; b < (int)p->get_border(face)->get_inclusions().size(); b++){
//                p->get_border(face)->get_inclusion(b)->get_border()->lock_read();
                segment_in |= p->get_border(face)->get_inclusion(b)->get_segment_in();
//                p->get_border(face)->get_inclusion(b)->get_border()->unlock_read();
            }

            if(p->get_border(face)->get_segment_out() != (segment_in & p->get_border(face)->get_segment_out())){
                change = true;
//                p->get_border(face)->lock_read();
                p->get_border(face)->set_segment_out(segment_in, true);
//                p->get_border(face)->unlock_read();
            }
        }
        //        }

        // **************** IN CONTINUITY *************
        //        if(p->get_border(face)->get_continuity_in()){
        Interval segment_out = Interval::EMPTY_SET;

        for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
            p->get_border(face)->get_inclusion(b)->get_border()->lock_read();
            segment_out |= p->get_border(face)->get_inclusion(b)->get_segment_out();
            p->get_border(face)->get_inclusion(b)->get_border()->unlock_read();
        }

        if(backward){
            if(!p->get_inner_mode()){
                if(p->get_border(face)->get_segment_in() != (segment_out & p->get_border(face)->get_segment_in())){
                    change = true;
                    p->get_border(face)->lock_read();
                    p->get_border(face)->set_segment_in(segment_out, true);
                    p->get_border(face)->unlock_read();
                }
            }
        }
        else{
            if(p->get_border(face)->get_segment_in() != (p->get_border(face)->get_segment_in() | segment_out & p->get_border(face)->get_segment_full())){
                change = true;
                p->get_border(face)->lock_read();
                p->get_border(face)->set_segment_in(segment_out, false);
                p->get_border(face)->unlock_read();
            }
        }
        //        }
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

bool Utils::test_discontinuity(const Interval &theta1, const Interval &theta2, const Interval &rotation){
    bool test_cos, test_sin, test_cos1, test_sin1;

    if(rotation == Interval::ZERO){
        test_cos = !(cos(theta1) & Interval(-1)).is_empty();
        test_sin = !(sin(theta1) & Interval::ZERO).is_empty();
        test_cos1 = !(cos(theta2) & Interval(-1)).is_empty();
        test_sin1 = !(sin(theta2) & Interval::ZERO).is_empty();
    }
    else if(rotation == Interval::HALF_PI){
        test_cos = !(sin(theta1) & Interval(-1)).is_empty();
        test_sin = !(cos(theta1) & Interval::ZERO).is_empty();
        test_cos1 = !(sin(theta2) & Interval(-1)).is_empty();
        test_sin1 = !(cos(theta2) & Interval::ZERO).is_empty();
    }
    else if(rotation == Interval::PI || rotation == -Interval::PI){
        test_cos = !(cos(theta1) & Interval(1)).is_empty();
        test_sin = !(sin(theta1) & Interval::ZERO).is_empty();
        test_cos1 = !(cos(theta2) & Interval(1)).is_empty();
        test_sin1 = !(sin(theta2) & Interval::ZERO).is_empty();
    }
    else if(rotation == -Interval::HALF_PI){
        test_cos = !(sin(theta1) & Interval(1)).is_empty();
        test_sin = !(cos(theta1) & Interval::ZERO).is_empty();
        test_cos1 = !(sin(theta2) & Interval(1)).is_empty();
        test_sin1 = !(cos(theta2) & Interval::ZERO).is_empty();
    }
    else{
        cout << "ROTATION ERROR, theta = " << rotation << endl;
    }
    return (test_cos && test_sin && test_cos1 && test_sin1);
}

void Utils::CtcPaveBackward2(Pave *p, bool inclusion, std::vector<bool> &change_tab){
    IntervalVector zero(2, Interval::ZERO);
    if(zero.is_subset(p->get_vector_field()))
        return;
    IntervalVector test(2);
    test[0] = Interval(1.5, 1.625);
    test[1] = Interval(-1, -0.75);
//    if(test == p->get_position())
//        p->draw_test(512, "test_before");
    vector<IntervalVector> seg_out_list;
    for(int i=0; i<4; i++)
        seg_out_list.push_back(IntervalVector(2, Interval::EMPTY_SET));

    for(int face = 0; face < 4; face++){        
        IntervalVector in(2, Interval::EMPTY_SET);
        for(int j=(face+1)%4; j!=face; j=(j+1)%4){
            IntervalVector seg_out(p->get_border(j)->get_segment_out_2D());
            IntervalVector seg_in(p->get_border(face)->get_segment_in_2D());
            this->CtcFlow(seg_in, seg_out, p->get_backward_function()? -p->get_vector_field():p->get_vector_field());
            in |= seg_in;
            seg_out_list[j] |= seg_out;
        }

        if(p->get_border(face)->get_segment_in_2D() != in)
            change_tab[face] = true;
        p->get_border(face)->set_segment_in(in[face%2], true);
    }
    for(int face=0; face <4; face++){
        if(p->get_border(face)->get_segment_out_2D() != seg_out_list[face])
            change_tab[face] = true;
        p->get_border(face)->set_segment_out(seg_out_list[face][face%2], true);
    }
//    if(test == p->get_position())
//        p->draw_test(512, "test_after", 200);
}

void Utils::CtcFlow(ibex::IntervalVector &in, ibex::IntervalVector &out, const ibex::IntervalVector &vect){
    // assert 0 not in v.
    IntervalVector c(out-in);
    IntervalVector v(vect);
    Interval alpha(Interval::POS_REALS);

    for(int i=0; i<v.size(); i++){
        // Issue with division by interval containing zero => need to separate cases
        alpha &= ((c[i]/(v[i] & Interval::POS_REALS)) & Interval::POS_REALS) | ((c[i]/(v[i] & Interval::NEG_REALS)) & Interval::POS_REALS);
    }

    c &= alpha*v;
    out &= c+in;
    in &= out-c;
}

void Utils::CtcVect(ibex::IntervalVector &seg, const ibex::IntervalVector &vect, int face, bool in){
    if(seg.is_empty())
        return;
    IntervalVector zero(2, Interval::ZERO);
    if(zero.is_subset(vect))
        return;

    IntervalVector bool_seg(2), bool_vect(2, Interval::EMPTY_SET);
    IntervalVector v(in?-vect:vect);

    bool_seg[face%2] = Interval(0.0, 1.0);
    if(face==0)
        bool_seg[(face+1)%2] = Interval(0);
    else if(face==1)
        bool_seg[(face+1)%2] = Interval(1);
    else if(face==2)
        bool_seg[(face+1)%2] = Interval(1);
    else if(face==3)
        bool_seg[(face+1)%2] = Interval(0);

    if(!(v[0] & Interval::POS_REALS).is_empty())
        bool_vect[0] |= Interval(1);
    if(!(v[0] & Interval::NEG_REALS).is_empty())
        bool_vect[0] |= Interval(0);
    if(!(v[1] & Interval::POS_REALS).is_empty())
        bool_vect[1] |= Interval(1);
    if(!(v[1] & Interval::NEG_REALS).is_empty())
        bool_vect[1] |= Interval(0);

    if((bool_seg & bool_vect).is_empty())
        seg.set_empty();
}



