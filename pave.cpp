#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"
#include <ctime>

#include "utils.h"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &box, Scheduler *scheduler): box(2)
{
    this->box = box;    // Box corresponding to the Pave
    this->scheduler = scheduler;
    this->borders.reserve(4);

    // Border creation
    IntervalVector coordinate(2);
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].lb()); this->borders.push_back(Border(coordinate, 0, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].ub()); this->borders.push_back(Border(coordinate, 1, this));
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].ub()); this->borders.push_back(Border(coordinate, 2, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].lb()); this->borders.push_back(Border(coordinate, 3, this));

    this->speed = Interval(1.0);
    this->theta[0] = Interval::EMPTY_SET;
    this->theta[1] = Interval::EMPTY_SET;

    double mu = 0.01;

    Interval dx = box[1];
    Interval dy = 1.0*(1-pow(box[0], 2))*box[1]-box[0];

    //    Interval rho = Interval::POS_REALS;
    //    Interval t = Interval::ZERO | 2.0*Interval::PI;
    //    this->scheduler->contract_polar.contract(dx, dy,  rho, t);
    //    this->theta = t;

    Interval theta = atan2(dy, dx);
    if(theta==(-Interval::PI|Interval::PI)){
        vector<Interval> dR = this->scheduler->utils.rotate(Interval::PI, dx, dy);
        Interval thetaR = atan2(dR[1], dR[0]);
        if(thetaR.diam()<theta.diam()){
            this->theta[0] = (thetaR+Interval::PI) & (-Interval::PI | Interval::PI);
            this->theta[1] = (thetaR-Interval::PI) & (-Interval::PI | Interval::PI);
        }
        else{
            this->theta[0] = theta;
        }
    }
    else{
        this->theta[0] = theta;
    }
}

void Pave::set_theta(ibex::Interval theta){
    this->theta[0] = Interval::EMPTY_SET;
    this->theta[1] = Interval::EMPTY_SET;

    if(theta.is_subset(-Interval::PI | Interval::PI)){
        this->theta[0] = theta;
    }
    else{
        this->theta[0] = (theta & (-Interval::PI | Interval::PI));

        if(!((theta + 2*Interval::PI) & (-Interval::PI | Interval::PI) ).is_empty())
            this->theta[1] =(theta + 2*Interval::PI);
        else if (!((theta - 2*Interval::PI) & (-Interval::PI | Interval::PI)).is_empty())
            this->theta[1] = (theta - 2*Interval::PI);
    }
}

void Pave::draw(){
    // Draw the pave
    vibes::drawBox(this->box, "b[]");

    // Draw the impacted segment (in option)
    for(int i=0; i<this->borders.size(); i++){
        this->borders[i].draw();
    }

    // Draw the inside impact -> polygone (in option) ?

    // Draw theta
    double size = 0.8*min(this->box[0].diam(), this->box[1].diam())/2.0;
    for(int i=0; i<2; i++){
         vibes::drawSector(this->box[0].mid(), this->box[1].mid(), size, size, (-this->theta[i].lb())*180.0/M_PI, (-this->theta[i].ub())*180.0/M_PI, "r[]");
    }

    // Draw Polygone
//    vector<double> x, y;
//    for(int i=0; i<this->borders.size(); i++){
//        this->borders[i].get_points(x, y);
//    }
//    vibes::drawPolygon(x, y, "g[g]");
}

void Pave::draw_borders(){
    for(int i=0; i<this->borders.size(); i++){
        this->borders[i].draw();
    }
}

void Pave::bisect(vector<Pave*> &result){
    // Create 4 new paves
    ibex::LargestFirst bisector(0.0, 0.5);
    std::pair<IntervalVector, IntervalVector> result_boxes = bisector.bisect(this->box);

    Pave *pave1 = new Pave(result_boxes.first, this->scheduler); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, this->scheduler); // Right or Down

    int indice1, indice2;

    if(pave1->box[0] == this->box[0]){
        // Case UP/DOWN bisection
        indice1 = 2;
        indice2 = 0;
    }
    else{
        // Case LEFT/RIGHT bisection
        indice1 = 1;
        indice2 = 3;
    }

    // Copy brothers Pave (this) to pave1 and pave2
    for(int i=0; i<4; i++){
        if(this->borders[i].brothers.size()!=0){
            if(i!=indice1)
                pave1->borders[i].add_brothers(this->borders[i].brothers);
            if(i!=indice2)
                pave2->borders[i].add_brothers(this->borders[i].brothers);
        }
    }

    // Add each other to its brother list (pave1 <-> pave2)
    pave1->borders[indice1].brothers.push_back(&pave2->borders[indice2]);
    pave2->borders[indice2].brothers.push_back(&pave1->borders[indice1]);

    // Remove
    for(int i=0; i<4; i++){
        this->borders[i].update_brothers(&(pave1->borders[i]), &(pave2->borders[i]));
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

void Pave::process(){
    // Process all new incoming valid segment (represents as borders)

    // Only take the first box in the list because the Pave is called for each new segment in the scheduler
    Border segment = queue.front();
    queue.erase(queue.begin());

    // Add the new segment & Test if the border interesect the segment of the pave
    vector<Interval> seg_in_list = this->borders[segment.face].add_segment(segment.segment);

    for(int i=0; i<seg_in_list.size(); i++){
        computePropagation(seg_in_list[i], segment.face);
    }

//    this->draw_borders();
}

void Pave::computePropagation(Interval seg_in, int face){
    // Translate and rotate the Segment
    IntervalVector box = this->box;
    IntervalVector segment(2);
    segment[face%2] = seg_in;
    segment[(face+1)%2] = (face == 1 | face == 2 ) ? Interval(box[(face+1)%2].ub()) : Interval(box[(face+1)%2].lb());

    Interval tab_theta[4] = {Interval::ZERO, -Interval::HALF_PI, Interval::PI, Interval::HALF_PI};

    this->scheduler->utils.translate_segment_and_box(segment, box, true, true);
    this->scheduler->utils.rotate_segment_and_box(segment, tab_theta[face], box, true);

    // Compute the propagation
    Interval seg_out[3];

    Interval segment_Rside[2] = {segment[0], segment[0]};
    Interval segment_front[2] = {segment[0], segment[0]};
    Interval segment_Lside[2] = {segment[0], segment[0]};

    for(int i=0; i<2; i++){
        this->scheduler->utils.CtcPropagateRightSide(segment_Rside[i], this->theta[i] + tab_theta[face], box);
        this->scheduler->utils.CtcPropagateFront(segment_front[i], this->theta[i] + tab_theta[face], box);
        this->scheduler->utils.CtcPropagateLeftSide(segment_Lside[i], this->theta[i] + tab_theta[face], box);
    }

    // Translate and rotate back the Segment
    IntervalVector segment_Rside_v(2);
    IntervalVector segment_front_v(2);
    IntervalVector segment_Lside_v(2);

    segment_Rside_v[1] = segment_Rside[0] | segment_Rside[1]; segment_Rside_v[0] = Interval(box[0].ub());
    segment_front_v[0] = segment_front[0] | segment_front[1]; segment_front_v[1] = Interval(box[1].ub());
    segment_Lside_v[1] = segment_Lside[0] | segment_Lside[1]; segment_Lside_v[0] = Interval(box[0].lb());

    this->scheduler->utils.rotate_segment_and_box(segment_Rside_v, -tab_theta[face], box, false);
    this->scheduler->utils.rotate_segment_and_box(segment_front_v, -tab_theta[face], box, false);
    this->scheduler->utils.rotate_segment_and_box(segment_Lside_v, -tab_theta[face], box, false);

    this->scheduler->utils.translate_segment_and_box(segment_Rside_v, this->box, false, false);
    this->scheduler->utils.translate_segment_and_box(segment_front_v, this->box, false, false); // Translate back with the initial box
    this->scheduler->utils.translate_segment_and_box(segment_Lside_v, this->box, false, false);

    seg_out[0] = (segment_Rside_v[0].diam() > segment_Rside_v[1].diam()) ? segment_Rside_v[0] : segment_Rside_v[1];
    seg_out[1] = (segment_front_v[0].diam() > segment_front_v[1].diam()) ? segment_front_v[0] : segment_front_v[1];
    seg_out[2] = (segment_Lside_v[0].diam() > segment_Lside_v[1].diam()) ? segment_Lside_v[0] : segment_Lside_v[1];

    for(int i=0; i<3; i++){
        if(!seg_out[i].is_empty()){
            vector<Interval> output = this->borders[(i+1+face)%4].add_segment(seg_out[i]);

            if(output.size()!=0){
                this->borders[(i+1+face)%4].publish_to_borthers(output[0]);
                if(output.size()==1){
                    this->borders[(i+1+face)%4].publish_to_borthers(output[1]);
                }
            }
        }
    }
}

void Pave::add_new_segment(Border &b){
    this->queue.push_back(b);
}

void Pave::warn_scheduler(){
    this->scheduler->add_to_queue(this);
}

void Pave::activate_pave(){
    Border b0(this->box[0], 0);
    Border b1(this->box[1], 1);
    Border b2(this->box[0], 2);
    Border b3(this->box[1], 3);

    this->queue.push_back(b0); this->warn_scheduler();
    this->queue.push_back(b1); this->warn_scheduler();
    this->queue.push_back(b2); this->warn_scheduler();
    this->queue.push_back(b3); this->warn_scheduler();
}


