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
}

Border::Border(const Interval &segment, const int face): position(2)
{
  this->segments.push_back(segment);
  this->face = face;
}

void Border::draw() const{
  for(int i=0; i<this->segments.size(); i++){

      // Create an IntervalVector (2D) from the segment (1D)
      IntervalVector segment = this->position;

      // Find the non flat dimension and complete replaced it by the segment
      segment[this->face%2] = this->segments[i];

      vibes::drawBox(segment, "g[]");
    }
}

// Add a new segment to the list of segment (check wheter there is overlapping)
// ToDo : Check whole overlapping
// ToDo : return the added segment (for the process function)
ibex::Interval Border::add_segment(Interval seg){
  // To Do : merge segments
  if(seg.is_empty()){
      return Interval::EMPTY_SET;
    }

  if((seg & this->position[this->face%2]).is_empty()){
      return Interval::EMPTY_SET;
    }

  for(int i=0; i<this->segments.size(); i++){
      Interval inter = seg & this->segments[i];
      Interval u = seg | this->segments[i];
      if(!inter.is_empty() || u.diam() == (seg.diam() + this->segments[i].diam())){
          this->segments[i] = inter;
          return inter;
        }
    }
  this->segments.push_back(seg);
  return seg;
}

// Publish a new segment to a list of brothers
void Border::publish_to_borthers(ibex::Interval seg){
  Border new_segment(seg, (this->face+2)%4);

  for(int i=0; i<this->brothers.size(); i++){
      this->brothers[i]->pave->push_queue(new_segment);
      this->brothers[i]->pave->warn_scheduler();
  }
}

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
