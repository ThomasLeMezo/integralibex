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
    vibes::drawSector(this->box[0].mid(), this->box[1].mid(), this->speed.mid(), this->speed.mid(), (-this->theta.lb())*180.0/M_PI, (-this->theta.ub())*180.0/M_PI, "r[]");

}

void Pave::cut(vector<Pave> *result){
    // Create 4 new paves

    // Link the borders between the new paves & with neighbours

}

void Pave::process(){
    // Process all new incoming valid segment (represents as borders)
    while(!this->queue.empty()){
        Border segment = queue.back();
        queue.pop_back();

        // Test if the border interesect the segment of the pave
        // Test if this is a new part of the segment
        int border_intersected=-1;
        Interval seg_in;

        for(int i=0; i<4; i++){
            Interval inter = segment.segments[0] & this->box[segment.face%2];
            if(inter.diam() !=0){
                border_intersected=i;
                seg_in = inter;
            }
        }

        if(border_intersected!=-1){
            computePropagation(seg_in, border_intersected);
        }

        // Then publish the impact on the neighbour paves

    }
}

void Pave::computePropagation(Interval seg_in, int face){

    vector<Interval> seg_out;
    project(seg_out, seg_in, this->table_rotation[face] + this->theta, this->box[face % 2], this->box[(face +1) % 2]);

    for(int i=0; i<3; i++){
        if(!seg_out[i].is_empty()){
            this->borders[(i+1+face)%4].add_segement(seg_out[i]);
            this->borders[(i+1+face)%4].publish_to_borthers(seg_out[i]);
        }
    }
}

void Pave::project(vector<Interval> &seg_out, Interval seg_in, Interval theta, Interval c0, Interval c1){
    seg_out.push_back((c1.ub() + (c0.lb() - seg_in) * tan(theta)) & c1);
    seg_out.push_back((seg_in + tan(theta) * c1.diam()) & c0);
    seg_out.push_back((c1.lb() + (seg_in - c0.lb())* tan(Interval::PI - theta)) & c1);
}
