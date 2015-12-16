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

    this->flow_in = true;
    for(int i=0; i<4; i++){
        this->flow_out[i] = false;
    }
}

Border::Border(const ibex::IntervalVector& position, const int face, const ibex::Interval &segment): position(2)
{
    this->position = position;
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

void Border::get_points(std::vector<double> &x, std::vector<double> &y){
    if(!this->segment.is_empty()){
        if(this->face == 0){
            x.push_back(this->segment.lb());
            x.push_back(this->segment.ub());
            y.push_back(this->position[1].lb());
            y.push_back(this->position[1].ub());
        }
        else if(this->face == 1){
            x.push_back(this->position[0].lb());
            x.push_back(this->position[0].ub());
            y.push_back(this->segment.lb());
            y.push_back(this->segment.ub());
        }
        else if(this->face == 2){
            x.push_back(this->segment.ub());
            x.push_back(this->segment.lb());
            y.push_back(this->position[1].ub());
            y.push_back(this->position[1].lb());
        }
        else if(this->face == 3){
            x.push_back(this->position[0].ub());
            x.push_back(this->position[0].lb());
            y.push_back(this->segment.ub());
            y.push_back(this->segment.lb());
        }
    }
}

// ********************************************************************************
// ****************** Segment Propagation *****************************************

/**
 * @brief Border::add_segment
 * @param seg
 * @return
 *
 * // Add a new segment to the list of segment (check wheter there is overlapping)
 * // ToDo : Check whole overlapping
 */
vector<ibex::Interval> Border::add_segment(const Interval &seg){
    vector<Interval> list_segments;

    if(seg.is_empty())
        return list_segments;

    Interval seg_box = seg & this->position[this->face%2];
    if(seg_box.is_empty())
        return list_segments;

    Interval new_segment = this->segment | seg_box;

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

/**
 * @brief Border::plug_segment
 * @param input
 * @return true if change / false if no change
 */
bool Border::plug_segment(ibex::Interval &input, ibex::Interval &position, bool modify){

    // Intersect the input with this border segment
    input &= this->segment;

    // Test if something is removed
    if(input == this->segment){
        return false;
    }
    else{
        // Do the union with all brothers
        Interval seg_out = input;

        if(position.is_strict_subset(this->position[this->face%2])){ // New incoming segment smaller than this
            // Union of all brothers segment of the same face
            for(int i=0; i<this->brothers.size(); i++){
                if(this->brothers[i]->position[this->brothers[i]->face%2]!=position){
                    // Case input is not inside the brothers
                    seg_out |= this->brothers[i]->segment;
                }
            }
        }

        if(this->segment == (seg_out & this->segment)){
            return false;
        }
        else{
            if(modify)
                this->segment &= seg_out;
            input = (this->segment & seg_out);
            return true;
        }
    }
}

/**
 * @brief Publish a new segment to a list of brothers
 * @param seg
 */
void Border::publish_to_borthers(ibex::Interval seg, bool forward){
    Border new_segment(this->position, (this->face+2)%4, seg);

    for(int i=0; i<this->brothers.size(); i++){
        this->brothers[i]->pave->add_new_segment(new_segment, forward);
        this->brothers[i]->pave->warn_scheduler(forward);
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

// Add new brothers to the list
void Border::add_brothers(std::vector<Border*> brother_list){
    for(int i=0; i<brother_list.size(); i++){
        IntervalVector test = brother_list[i]->position & this->position;
        if(!(test.is_empty()) && (test[0].is_degenerated() != test[1].is_degenerated())){
            this->brothers.push_back(brother_list[i]);
        }
    }
}

void Border::update_brothers(Border* border_pave1, Border* border_pave2){
    vector<Border*> list_tmp;
    list_tmp.push_back(border_pave1);
    list_tmp.push_back(border_pave2);

    for(int i=0; i<brothers.size(); i++){
        // Remove reference to this border inside brothers
        for(int j=0; j<this->brothers[i]->brothers.size(); j++){
            if(this->brothers[i]->brothers[j]==this){
                this->brothers[i]->brothers.erase(this->brothers[i]->brothers.begin()+j);
                break;
            }
        }

        // Add reference of border_pave 1 and 2
        this->brothers[i]->add_brothers(list_tmp);
    }
}

void Border::set_full(){
    this->segment = this->position[this->face%2];
}
