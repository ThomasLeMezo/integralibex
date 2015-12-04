#include "border.h"
#include "ibex.h"
#include "vibes.h"

#include "iostream"
#include "stdlib.h"
#include "stdio.h"

using namespace ibex;
using namespace std;

Border::Border(const IntervalVector &position, const int face, Pave *pave): position(2)
{
    this->position = position;
    this->face = face;
    this->pave = pave;
    this->segment = Interval::EMPTY_SET;
}

Border::Border(const Interval &segment, const int face): position(2)
{
    this->segment = segment;
    this->face = face;
}

void Border::draw() const{
    // Create an IntervalVector (2D) from the segment (1D)
    IntervalVector segment = this->position;

    // Find the non flat dimension and complete replaced it by the segment
    segment[this->face%2] = this->segment;

    vibes::drawBox(segment, "g[]");
}

// ********************************************************************************
// ****************** Segment Propagation *****************************************

// Add a new segment to the list of segment (check wheter there is overlapping)
// ToDo : Check whole overlapping
// ToDo : return the added segment (for the process function)
vector<ibex::Interval> Border::add_segment(Interval seg){
    // To Do : merge segments
    vector<Interval> list_segments;

    if(seg.is_empty()){
        return list_segments;
    }

    if((seg & this->position[this->face%2]).is_empty()){
        return list_segments;
    }

    Interval new_segment = this->segment | seg;

    // Compute the new segment(s) to propagate
    Interval left, right;
    new_segment.diff(this->segment, left, right);

    if(!left.is_empty())
      list_segments.push_back(left);
    if(!right.is_empty())
      list_segments.push_back(right);

    this->segment = new_segment;
    return list_segments;
}

// Publish a new segment to a list of brothers
void Border::publish_to_borthers(ibex::Interval seg){
    Border new_segment(seg, (this->face+2)%4);

    for(int i=0; i<this->brothers.size(); i++){
        this->brothers[i]->pave->add_new_segment(new_segment);
        this->brothers[i]->pave->warn_scheduler();
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

// Add new brothers to the list
void Border::add_brothers(std::vector<Border*> brother_list){
    for(int i=0; i<brother_list.size(); i++){
        if(!(brother_list[i]->position & this->position).is_empty()){
            this->brothers.push_back(brother_list[i]);
        }
    }
}

void Border::update_brothers(Border* border_pave1, Border* border_pave2){
    for(int i=0; i<brothers.size(); i++){
        // Remove reference to this border inside brothers
        for(int j=0; j<this->brothers[i]->brothers.size(); j++){
            if(this->brothers[i]->brothers[j]==this){
                this->brothers[i]->brothers.erase(this->brothers[i]->brothers.begin()+j);
            }
        }
        // Add reference of border_pave 1 and 2
        vector<Border*> list_tmp;
        list_tmp.push_back(border_pave1);
        list_tmp.push_back(border_pave2);
        this->brothers[i]->add_brothers(list_tmp);
    }
}
