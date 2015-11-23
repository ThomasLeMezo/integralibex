#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &b): box(2)
{
    this->box = b;

    IntervalVector v(2);
    v[0] = b[0]; v[1] = Interval(b[1].lb()); this->borders.push_back(Border(v, 0, this));
    v[0] = b[1]; v[1] = Interval(b[0].ub()); this->borders.push_back(Border(v, 1, this));
    v[0] = b[0]; v[1] = Interval(b[1].ub()); this->borders.push_back(Border(v, 2, this));
    v[0] = b[1]; v[1] = Interval(b[0].lb()); this->borders.push_back(Border(v, 3, this));
}

void Pave::draw(){
    // Draw the pave
    vibes::drawBox(this->box, "b[]");

    // Draw the impacted segment (in option)

    // Draw the inside impact -> polygone (in option) ?
    // Difficult because of non complete interval

}

void Pave::cut(vector<Pave> *result){
    // Create 4 new paves

    // Link the borders between the new paves & with neighbours

}

void Pave::process(){
    // Process the new incoming valid segment (represents as borders)
    if(!this->queue.empty()){
        Border b = queue.back();
        queue.pop_back();

        // Test if the border interesect the segment of the pave
        // Test if this is a new part of the segment

        if(false)
            cout << "test" << endl;

        // If true validate the new part and calculate the impact on the three other borders

        // Then publish the impact on the neighbour paves

    }
}
