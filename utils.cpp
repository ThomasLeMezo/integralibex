#include "utils.h"
#include "pave.h"
#include "iomanip"

using namespace std;
using namespace ibex;

// cout << setprecision(80) << "..." << endl;

// ********************************************************************************
// ****************** Contractors Global functions ********************************

void CtcPropagateSegment(ibex::Interval &seg_in, std::vector<ibex::Interval> &seg_out, const int &face, const std::vector<ibex::Interval> &theta, const ibex::IntervalVector &box_pave, const ibex::Interval &u, bool inner, bool inner_backward){
//    // Translate and rotate the Segment
//    IntervalVector box(box_pave);
//    IntervalVector box_in(box_pave);
//    IntervalVector segment_in = segment2IntervalVector(seg_in, face, box);
//    IntervalVector segment_out[3] = IntervalVector(2);
//    for(int i=0; i<3; i++){
//        segment_out[i] = segment2IntervalVector(seg_out[i], (face+1+i)%4, box);
//    }

//    this->translate_segment_and_box(segment_in, box, true, true);
//    IntervalVector box_translate(box);
//    this->rotate_segment_and_box(segment_in, this->tab_rotation[face], box, true);

//    for(int i=0; i<3; i++){
//        this->translate_segment_and_box(segment_out[i], box_in, true, false);
//        this->rotate_segment_and_box(segment_out[i], this->tab_rotation[face], box_translate, false);
//    }

//    // Compute the propagation
//    ibex::Interval segment_norm_in[3][2], segment_norm_out[3][2];

//    for(int i=0; i<3; i++){
//        for(int j=0; j<2; j++){
//            segment_norm_in[i][j] = segment_in[0];
//            segment_norm_out[i][j] = (segment_out[i][0].diam() > segment_out[i][1].diam()) ? segment_out[i][0] : segment_out[i][1];
//        }
//    }

//    if(!inner){
//        for(int i=0; i<2; i++){
//            this->CtcPropagateRightSide(segment_norm_in[0][i], segment_norm_out[0][i], theta[i] + tab_rotation[face] + u, box);
//            this->CtcPropagateFront(segment_norm_in[1][i], segment_norm_out[1][i], theta[i] + tab_rotation[face] + u, box);
//            this->CtcPropagateLeftSide(segment_norm_in[2][i], segment_norm_out[2][i], theta[i] + tab_rotation[face] + u, box);
//        }
//    }
//    else{
//        vector<ibex::Interval> theta_list = {theta[0] + tab_rotation[face], theta[1] + tab_rotation[face]};
//        this->CtcPropagateRightSideInner(segment_norm_in[0][0], segment_norm_out[0][0], theta_list, box, u ,false, inner_backward);
//        this->CtcPropagateFrontInner(segment_norm_in[1][0], segment_norm_out[1][0], theta_list, box, u, inner_backward);
//        this->CtcPropagateLeftSideInner(segment_norm_in[2][0], segment_norm_out[2][0], theta_list, box, u, false, inner_backward);
//        for(int i=0; i<3; i++){
//            segment_norm_out[i][1] = ibex::Interval::EMPTY_SET;
//            segment_norm_in[i][1] = ibex::Interval::EMPTY_SET;
//        }
//    }

//    // **************************** OUTPUT ****************************
//    // Translate and rotate back the Segment
//    IntervalVector segment_contracted_out[3] = IntervalVector(2);
//    segment_contracted_out[0][1] = segment_norm_out[0][0] | segment_norm_out[0][1]; segment_contracted_out[0][0] = ibex::Interval(box[0].ub());
//    segment_contracted_out[1][0] = segment_norm_out[1][0] | segment_norm_out[1][1]; segment_contracted_out[1][1] = ibex::Interval(box[1].ub());
//    segment_contracted_out[2][1] = segment_norm_out[2][0] | segment_norm_out[2][1]; segment_contracted_out[2][0] = ibex::Interval(box[0].lb());

//    // Rotate and translate back with the initial box
//    vector<ibex::Interval> seg_out_tmp;
//    for(int i=0; i<3; i++){
//        this->rotate_segment_and_box(segment_contracted_out[i], -tab_rotation[face], box, false);
//        this->translate_segment_and_box(segment_contracted_out[i], box_in, false, false);
//        // Add segment to seg_out list
//        seg_out_tmp.push_back( (seg_out[i] & ((segment_contracted_out[i][0].diam() > segment_contracted_out[i][1].diam()) ? segment_contracted_out[i][0] : segment_contracted_out[i][1])));
//    }
//    seg_out = seg_out_tmp;

//    // **************************** INPUT ****************************
//    IntervalVector segment_contracted_in = IntervalVector(2);
//    segment_contracted_in[0] = ibex::Interval::EMPTY_SET;
//    segment_contracted_in[1] = ibex::Interval(box[1].lb());
//    for(int i=0; i<3; i++){
//        for(int j=0; j<2; j++){
//            segment_contracted_in[0] |= segment_norm_in[i][j];
//        }
//    }
//    // Rotate and translate back with the initial box
//    this->rotate_segment_and_box(segment_contracted_in, -tab_rotation[face], box, false);
//    this->translate_segment_and_box(segment_contracted_in, box_in, false, false);
//    seg_in &= ((segment_contracted_in[0].diam() > segment_contracted_in[1].diam()) ? segment_contracted_in[0] : segment_contracted_in[1]);
}

void CtcPaveBackward(Pave *p, bool inclusion, bool inner){

//    ibex::Interval segment_in[4] = {ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET};

//    for(int face = 0; face < 4; face++){
//        ibex::Interval seg_in = p->get_border(face)->get_segment_in();

//        vector<ibex::Interval> seg_out;
//        for(int j=(face+1)%4; j!=face; j=(j+1)%4){
//            seg_out.push_back(p->get_border(j)->get_segment_out());
//        }

//        if(!inner)
//            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), false, false);
//        else
//            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), true, true);

//        segment_in[face] = seg_in;
//    }

//    for(int face = 0; face<4; face++){
//        p->get_border(face)->set_segment_in(segment_in[face], inclusion);
//    }
}

void CtcPaveForward(Pave *p, bool inclusion, bool inner){
//    ibex::Interval segment_out[4] = {ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET, ibex::Interval::EMPTY_SET};

//    for(int face = 0; face < 4; face++){
//        ibex::Interval seg_in = p->get_border(face)->get_segment_in();

//        vector<ibex::Interval> seg_out;
//        for(int j=0; j<3; j++){
//            seg_out.push_back(ibex::Interval::ALL_REALS);
//        }

//        if(!inner){
//            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), false, false);
//        }
//        else{
//            this->CtcPropagateSegment(seg_in, seg_out, face, p->get_theta(), p->get_position(), p->get_u(), true, false);
//        }

//        int k=0;
//        for(int i=(face+1)%4; i!=face; i=(i+1)%4){
//            segment_out[i] |= seg_out[k];
//            k++;
//        }
//    }

//    for(int face = 0; face<4; face++){
//        p->get_border(face)->set_segment_out(segment_out[face], inclusion);
//    }
}

// ********************************************************************************
// ****************** Algorithm functions      ************************************

void CtcPaveConsistency(Pave *p, bool backward, bool inner){
    if(backward){
        CtcPaveBackward(p, backward, inner);
        Pave p2(p);
        CtcPaveForward(&p2, backward, inner);
        *p &= p2;
        if(inner){
            Pave p3(p);
            CtcPaveBackward(&p3, backward, inner);
            *p &= p3;
        }
    }
    else{
        CtcPaveForward(p, backward, inner);
    }
}

bool CtcContinuity(Pave *p, bool backward){
//    bool change = false;

//    for(int face = 0; face < 4; face++){
//        ibex::Interval segment_in = ibex::Interval::EMPTY_SET;
//        ibex::Interval segment_out = ibex::Interval::EMPTY_SET;

//        for(int b = 0; b < p->get_border(face)->get_inclusions().size(); b++){
//            segment_in |= p->get_border(face)->get_inclusion(b)->get_segment_in();
//            segment_out |= p->get_border(face)->get_inclusion(b)->get_segment_out();
//        }

//        if(backward){
//            if(p->get_border(face)->get_segment_in() != (segment_out & p->get_border(face)->get_segment_in())
//                    || p->get_border(face)->get_segment_out() != (segment_in & p->get_border(face)->get_segment_out())){
//                change = true;
//                p->get_border(face)->set_segment_in(segment_out, backward);
//                p->get_border(face)->set_segment_out(segment_in, backward);
//            }
//        }
//        else{
//            if(p->get_border(face)->get_segment_in() != (p->get_border(face)->get_segment_in() | segment_out & p->get_border(face)->get_segment_full())){
//                change = true;
//                p->get_border(face)->set_segment_in(segment_out, backward);
//            }
//        }

//    }

//    return change;
    return false;
}
