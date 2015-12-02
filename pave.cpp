#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"
#include <ctime>

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

    double mu = 0.01;

    Interval dx = box[1];
    Interval dy = 1.0*(1-pow(box[0], 2))*box[1]-box[0];

    //    Interval rho = Interval::POS_REALS;
    //    Interval t = Interval::ZERO | 2.0*Interval::PI;
    //    this->scheduler->contract_polar.contract(dx, dy,  rho, t);
    //    this->theta = t;

    Interval theta = atan2(dy, dx);
    if(theta==(-Interval::PI|Interval::PI)){
        vector<Interval> dR = scheduler->rotate(Interval::PI, dx, dy);
        Interval thetaR = atan2(dR[1], dR[0]);
        if(thetaR.diam()<theta.diam()){
            this->theta.push_back((thetaR+Interval::PI) & (-Interval::PI | Interval::PI));
            this->theta.push_back((thetaR-Interval::PI) & (-Interval::PI | Interval::PI));
        }
        else{
            this->theta.push_back(theta);
        }
    }
    else{
        this->theta.push_back(theta);
    }
}

void Pave::draw() const{
    // Draw the pave
    //    vibes::drawBox(this->box, "b[]");

    // Draw the impacted segment (in option)
    for(int i=0; i<this->borders.size(); i++){
        this->borders[i].draw();
    }

    // Draw the inside impact -> polygone (in option) ?

    // Draw theta
    double size = 0.8*min(this->box[0].diam(), this->box[1].diam())/2.0;
    for(int i=0; i<this->theta.size(); i++){
        vibes::drawSector(this->box[0].mid(), this->box[1].mid(), size, size, (-this->theta[i].lb())*180.0/M_PI, (-this->theta[i].ub())*180.0/M_PI, "r[]");
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

    for(int i=0; i<4; i++){
        if(this->borders[i].brothers.size()!=0){
            if(i!=indice1)
                pave1->borders[i].add_brothers(this->borders[i].brothers);
            if(i!=indice2)
                pave2->borders[i].add_brothers(this->borders[i].brothers);
        }
    }
    pave1->borders[indice1].brothers.push_back(&pave2->borders[indice2]);
    pave2->borders[indice2].brothers.push_back(&pave1->borders[indice1]);

    for(int i=0; i<4; i++){
        this->borders[i].update_brothers(&(pave1->borders[i]), &(pave2->borders[i]));
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

void Pave::process(){
    // Process all new incoming valid segment (represents as borders)

    // Only take the first box in the list because the Pave is called for each new segment in the scheduler
    Border segment = queue.back();
    queue.pop_back();

    // Add the new segment & Test if the border interesect the segment of the pave
    vector<Interval> seg_in_list = this->borders[segment.face].add_segment(segment.segments[0]);
    //    cout << seg_in_list.size() << endl;

    for(int i=0; i<seg_in_list.size(); i++){
        computePropagation(seg_in_list[i], segment.face);
    }
}

void Pave::computePropagation(Interval seg_in, int face){
    vector<Interval> seg_out;

    // Passage dans le repère local
    double offset_x = box[0].lb();
    double offset_y = box[1].lb();

    Interval seg_in_local;
    switch(face){
    case 0:
        seg_in_local = seg_in - box[0].lb();
        break;
    case 1:
        seg_in_local = seg_in - box[1].lb();
        break;
    case 2:
        seg_in_local = box[0].ub() - seg_in;
        break;
    case 3:
        seg_in_local = box[1].ub() - seg_in;
        break;
    }

    Interval c0 = box[face % 2] - box[face % 2].lb();
    Interval c1 = box[(face + 1) % 2] - box[(face + 1) % 2].lb();
    vector<Interval> theta;
    for(int i=0; i<this->theta.size(); i++){
        theta.push_back(this->table_rotation[face] + this->theta[i]);
    }

    // ****** RIGHT Border ******
    vector<Interval> x_right;
    vector<Interval> y_right;
    vector<Interval> rho_right;
    vector<Interval> theta_c_right;

    for(int i=0; i<this->theta.size(); i++){
        theta_c_right.push_back(Interval::PI/2.0 + theta[i]);
        x_right.push_back(c0.ub() - seg_in_local);
        y_right.push_back(c1);
        rho_right.push_back(Interval::POS_REALS);
        this->scheduler->contract_polar.contract(x_right[i], y_right[i], rho_right[i], theta_c_right[i]);
    }

    // ****** FRONT Border ******
    vector<Interval> x_front;
    vector<Interval> y_front;
    vector<Interval> rho_front;
    vector<Interval> theta_c_front;
    vector<Interval> front;

    for(int i=0; i<this->theta.size(); i++){
        theta_c_front.push_back(theta[i]);
        x_front.push_back(Interval(c1.diam()));
        y_front.push_back(Interval::ALL_REALS);
        rho_front.push_back(Interval::POS_REALS);
        this->scheduler->contract_polar.contract(x_front[i], y_front[i], rho_front[i], theta_c_front[i]);
        front.push_back((seg_in_local - y_front[i]) & c0);
    }

    // ****** LEFT Border ******
    vector<Interval> x_left;
    vector<Interval> y_left;
    vector<Interval> rho_left;
    vector<Interval> theta_c_left;

    for(int i=0; i<this->theta.size(); i++){
        theta_c_left.push_back(Interval::PI/2.0 - theta[i]);
        x_left.push_back(seg_in_local);
        y_left.push_back(c1);
        rho_left.push_back(Interval::POS_REALS);
        this->scheduler->contract_polar.contract(x_left[i], y_left[i], rho_left[i], theta_c_left[i]);
    }

    // Passage dans le repère global (et rotation)

    for(int i=0; i<this->theta.size(); i++){
        switch(face){
        case 0:
            seg_out.push_back(y_right[i] + offset_y);
            seg_out.push_back(front[i] + c0.lb() + offset_x);
            seg_out.push_back(y_left[i] + offset_y);
            break;
        case 1:
            seg_out.push_back(c1.ub() - y_right[i] + offset_x);
            seg_out.push_back(front[i] + c1.lb() + offset_y);
            seg_out.push_back(c1.ub() - y_left[i] + offset_x);
            break;
        case 2:
            seg_out.push_back(c1.ub() - y_right[i] + offset_y);
            seg_out.push_back(c0.ub() - front[i] + offset_x);
            seg_out.push_back(c1.ub() - y_left[i] + offset_y);
            break;
        case 3:
            seg_out.push_back(y_right[i] + offset_x);
            seg_out.push_back(c1.ub() - front[i] + offset_y);
            seg_out.push_back(y_left[i] + offset_x);
            break;
        }
    }

    // ******* Publish new segments *******
    for(int i=0; i<seg_out.size(); i++){
        if(!seg_out[i].is_empty()){
            this->borders[(i%3+1+face)%4].add_segment(seg_out[i]);
            this->borders[(i%3+1+face)%4].publish_to_borthers(seg_out[i]);
        }
    }
}

void Pave::push_queue(Border &b){
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

    this->queue.push_back(b0);
    this->queue.push_back(b1);
    this->queue.push_back(b2);
    this->queue.push_back(b3);

    this->warn_scheduler();
}
