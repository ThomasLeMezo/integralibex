#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &box): box(2)
{
    this->box = box;    // Box corresponding to the Pave

    // Border creation
    IntervalVector coordinate(2);
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].lb()); this->borders.push_back(Border(coordinate, 0, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].ub()); this->borders.push_back(Border(coordinate, 1, this));
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].ub()); this->borders.push_back(Border(coordinate, 2, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].lb()); this->borders.push_back(Border(coordinate, 3, this));

    this->theta = Interval(0.0);
    this->speed = Interval(1.0);

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
    double size = min(this->box[0].diam(), this->box[1].diam())/2.0;
    vibes::drawSector(this->box[0].mid(), this->box[1].mid(), size, size, (-this->theta.lb())*180.0/M_PI, (-this->theta.ub())*180.0/M_PI, "r[]");

}

void Pave::bisect(vector<Pave> *result){
    // Create 4 new paves
    std::pair<IntervalVector(2), IntervalVector(2)> result_boxes;
    ibex::LargestFirst bisector(0.0, 0.5);

    result_boxes = bisector.bisect(this->box);

    // Link the borders between the new paves & with neighbours

}

void Pave::process(){
    // Process all new incoming valid segment (represents as borders)
    while(!this->queue.empty()){
        Border segment = queue.back();
        queue.pop_back();

        // Test if the border interesect the segment of the pave
        // ToDo : Test if this is a new part of the segment !!!
        Interval seg_in = segment.segments[0] & this->box[(segment.face)%2];

        if(seg_in.diam() !=0){
            computePropagation(seg_in, segment.face);
        }

        // Then publish the impact on the neighbour paves

    }
}

void Pave::computePropagation(Interval seg_in, int face){
    vector<Interval> seg_out;
    CtcPolar contract_p0;

    // Passage dans le repère local
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

    double offset_x = box[face % 2].lb();
    double offset_y = box[(face + 1) % 2].lb();

    Interval c0 = box[face % 2] - offset_x;
    Interval c1 = box[(face + 1) % 2] - offset_y;
    Interval theta = this->table_rotation[face] + this->theta;

    // ****** RIGHT Border ******
    Interval x = c0.ub() - seg_in_local;
    Interval y = c1;
    Interval rho = Interval::POS_REALS;
    Interval theta_c = Interval::PI/2.0 + theta;

    contract_p0.contract(x, y, rho, theta_c);

    switch(face){
    case 0:
        seg_out.push_back(y + offset_y);
        break;
    case 1:
        seg_out.push_back(c1.ub() - y + offset_x);
        break;
    case 2:
        seg_out.push_back(c1.ub() - y + offset_y);
        break;
    case 3:
        seg_out.push_back(y + offset_x);
        break;
    }

    // ****** FRONT Border ******
    x = Interval(c1.diam());
    y = Interval::ALL_REALS;
    rho = Interval::POS_REALS;
    theta_c = theta;

    contract_p0.contract(x, y, rho, theta_c);
    Interval front = (seg_in_local - y) & c0;

    switch(face){
    case 0:
        seg_out.push_back(front + c0.lb() + offset_x);
        break;
    case 1:
        seg_out.push_back(front + c1.lb() + offset_y);
        break;
    case 2:
        seg_out.push_back(c0.ub() - front + offset_x);
        break;
    case 3:
        seg_out.push_back(c1.ub() - front + offset_y);
        break;
    }

    // ****** LEFT Border ******
    x = seg_in_local;
    y = c1;
    rho = Interval::POS_REALS;
    theta_c = Interval::PI/2.0 - theta;

    contract_p0.contract(x, y, rho, theta_c);

    switch(face){
    case 0:
        seg_out.push_back(y + offset_y);
        break;
    case 1:
        seg_out.push_back(c1.ub() - y + offset_x);
        break;
    case 2:
        seg_out.push_back(c1.ub() - y + offset_y);
        break;
    case 3:
        seg_out.push_back(y + offset_x);
        break;
    }

    // ******* Publish new segments *******
    for(int i=0; i<seg_out.size(); i++){
        if(!seg_out[i].is_empty()){
            this->borders[(i+1+face)%4].add_segement(seg_out[i]);
            this->borders[(i+1+face)%4].publish_to_borthers(seg_out[i]);
        }
    }
}
