#include "border.h"
#include "ibex.h"
#include "vibes.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

using namespace ibex;
using namespace std;

Border::Border(const IntervalVector &position, int face, Pave *pave): position(2)
{
    this->position = position;
    this->face = face;
    this->pave = pave;
}

Border::Border(const Interval &segment, int face): position(2)
{
    this->segments.push_back(segment);
    this->face = face;
}

void Border::draw(){
    for(int i=0; i<this->segments.size(); i++){

        // Create an IntervalVector (2D) from the segment (1D)
        IntervalVector segment = this->position;

        // Find the non flat dimension and complete replaced it by the segment
        segment[this->face%2] = this->segments[i];

        vibes::drawBox(segment, "g[]");
    }
}

void Border::add_segement(Interval seg){
    // To Do : merge segments
    if(seg.is_empty())
        return;

    for(int i=0; i<this->segments.size(); i++){
        Interval inter = seg & this->segments[i];
        Interval u = seg | this->segments[i];
        if(!inter.is_empty() || u.diam() == (seg.diam() + this->segments[i].diam())){
            this->segments[i] = inter;
            return;
        }
    }
    this->segments.push_back(seg);
}

void Border::publish_to_borthers(ibex::Interval seg){
    for(int i=0; i<this->borthers.size(); i++){
        Border new_segment(seg, (this->face+2)%2);
        this->borthers[i]->pave->queue.push_back(new_segment);
    }
}

