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

    visited_node = false;
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
void Pave::activate_pave(){
    for(int face=0; face<4; face++){
        Border b(this->get_border_position(face), face, this->box[face%2]);
        this->queue_forward.push_back(b);this->warn_scheduler(true);
        this->borders[face].publish_to_borthers(this->box[face%2], true);
    }
}

void Pave::set_full_continuity(){
    for(int face=0; face < 4; face++){
        this->borders[face].set_full();

        if(this->borders[face].brothers.size()==0){
            this->borders[face].set_full();
            Border b(this->borders[face].position, face, Interval::EMPTY_SET);
            this->add_new_segment(b, false);
            this->warn_scheduler(false);
        }
    }
}

IntervalVector Pave::get_border_position(int face){
    IntervalVector position(2);
    position[0] = this->box[face%2];
    position[1] = this->box[(face+1)%2];
    return position;
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Pave::draw(bool filled){
    // Draw the pave
    vibes::drawBox(this->box, "b[]");

    this->draw_borders(filled);

    // Draw theta
    double size = 0.8*min(this->box[0].diam(), this->box[1].diam())/2.0;
    for(int i=0; i<2; i++){
        vibes::drawSector(this->box[0].mid(), this->box[1].mid(), size, size, (-this->theta[i].lb())*180.0/M_PI, (-this->theta[i].ub())*180.0/M_PI, "r[]");
    }
}
void Pave::draw_borders(bool filled){
    if(!filled){
        // Draw Segments
        for(int i=0; i<this->borders.size(); i++){
            this->borders[i].draw();
        }
    }
    else{
        // Draw Polygone
        vector<double> x, y;
        for(int i=0; i<this->borders.size(); i++){
            this->borders[i].get_points(x, y);
        }
        vibes::drawPolygon(x, y, "g[g]");
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

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
            if(i!=indice1){
                pave1->borders[i].add_brothers(this->borders[i].brothers);
                //                pave1->queue_forward.push_back(this->borders[i]);this->warn_scheduler(true); // Add segments of this to pave1
            }
            if(i!=indice2){
                pave2->borders[i].add_brothers(this->borders[i].brothers);
                //                pave2->queue_forward.push_back(this->borders[i]);this->warn_scheduler(true); // Add segments of this to pave2
            }
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

// ********************************************************************************
// ****************** Segment Propagation *****************************************

void Pave::process_forward(){
    // Process all new incoming valid segment (represents as borders)
    if(queue_forward.size()==0)
        return;
    // Only take the first box in the list because the Pave is called for each new segment in the scheduler
    Border border = queue_forward.front();
    queue_forward.erase(queue_forward.begin());

    // Add the new segment & Test if the border interesect the segment of the pave
    vector<Interval> seg_in_list = this->borders[border.face].add_segment(border.segment);

    for(int i=0; i<seg_in_list.size(); i++){
        vector<Interval> seg_out;
        for(int j=0; j<3; j++){
            seg_out.push_back(Interval::ALL_REALS);
        }
        // Compute the propagation to the 3 other face
        this->scheduler->utils.CtcPropagateSegment(seg_in_list[i], seg_out, border.face, this->theta, this->box);
        // Apply the principle of continuity by sending results to neighbours

        for(int j=0; j<3; j++){
            if(!seg_out[j].is_empty()){
                vector<Interval> output = this->borders[(j+1+border.face)%4].add_segment(seg_out[j]);
                for(int k=0; k<output.size(); k++){
                    this->borders[(j+1+border.face)%4].publish_to_borthers(output[k], true);
                }
            }
        }
    }

    //    this->draw_borders();
}

void Pave::process_backward(){
    // Process all new incoming valid segment (represents as borders)
    if(queue_backward.size()==0)
        return;
    // Only take the first box in the list because the Pave is called for each new segment in the scheduler
    Border border = queue_backward.front();
    queue_backward.erase(queue_backward.begin());

    bool test_cout = false;
    IntervalVector position(2);

    position[0] = Interval(-9.42, -9.4);
    position[1] = Interval(-1.2, -1);

    if(!(position & this->box).is_empty()){
        //        this->scheduler->draw(1024, false);
        cout << "***********" << endl << "TEST - face=" << border.face << endl;
        test_cout = true;
        //this->scheduler->draw(1024, false);
    }

    if(this->borders[border.face].plug_segment(border.segment, true)){
        Interval seg_out[4] = Interval::EMPTY_SET;

        // Propagation of faces
        for(int face = 0; face<4; face++){
            Border *border_propagation = &(this->borders[face]);
            for(int brother=0; brother < border_propagation->brothers.size(); brother++){
                //Définition du segment
                Interval seg_in = border_propagation->brothers[brother]->segment;
//                Interval seg_in = this->borders[face].segment; // Converge plus lentement !!

                vector<Interval> seg_out_brother;
                for(int j=(face+1)%4; j!=((face+3)%4); j=(j+1)%4){
                    seg_out_brother.push_back(this->borders[j].segment);
                }
                // Propagation
                this->scheduler->utils.CtcPropagateSegment(seg_in, seg_out_brother, face, this->theta, this->box);

                // Ajout au seg_out
                int k=0;
                for(int j=(face+1)%4 ; j!=(face+4)%4; j=((j+1)%4)){
                    seg_out[j] |= seg_out_brother[k];
                    k++;
                }
            }
        }

        // Publish results to neighbours + modify this pave borders
        for(int j=0; j<4; j++){
            if(j != border.face){
                if((this->borders[j].flow_out[border.face]==true) && (this->borders[j].flow_in == false)){
                    if(this->borders[j].plug_segment(seg_out[j], false)){ // modify segment
                        this->borders[j].publish_to_borthers(seg_out[j], false); // publish to neighbour
                    }
                }
            }
        }
    }
}

void Pave::compute_flow(){
    bool flow_in[4] = {false, false, false, false};
    bool flow_out[4] = {false, false, false, false};

    for(int face = 0; face < 4; face++){
        Interval seg_in = this->borders[face].position[face%2];
        vector<Interval> seg_out;
        for(int i=0; i<3; i++){
            seg_out.push_back(this->borders[(face+i+1)%4].segment);
        }
        this->scheduler->utils.CtcPropagateSegment(seg_in, seg_out, face, this->theta, this->box);

        if(!seg_out[0].is_empty() || !seg_out[1].is_empty() || !seg_out[2].is_empty()){
            flow_in[face] = true;
        }
        for(int i = 0; i<3; i++){
            if(!seg_out[i].is_empty()){
                flow_out[(face+i+1)%4] = true;
                this->borders[(face+i+1)%4].flow_out[face] = true;
            }
            else{
                this->borders[(face+i+1)%4].flow_out[face] = false;
            }
        }
    }

    for(int face = 0; face <4; face++){
        this->borders[face].flow_in = flow_in[face];
    }
}


/**
 * @brief Pave::warn_scheduler
 * Warn the scheduler that this Pave has new segment to forward/backward process
 */
void Pave::warn_scheduler(bool forward){
    this->scheduler->add_to_queue(this, forward);
}

// ********************************************************************************
// ****************** Graph functions *********************************************

void Pave::compute_successors(){
    for(int face=0; face<4; face++){
        Border test_border(this->get_border_position(face), face, this->box[face%2]);
        vector<Interval> output;
        for(int i=0; i<3; i++)
            output.push_back(Interval::ALL_REALS);

        this->scheduler->utils.CtcPropagateSegment(test_border.segment, output, test_border.face, this->theta, this->box);

        for(int i=0; i<3; i++){
            if(!output[i].is_empty()){
                for(int j=0; j< this->borders[(i+1+face)%4].brothers.size(); j++){
                    this->borders[(i+1+face)%4].brothers[j]->pave->add_precursors(this);  // Add precursors
                    this->add_successors(this->borders[(i+1+face)%4].brothers[j]->pave);  // Add successors
                }
            }
        }
    }
}
bool Pave::test_cycle(Pave* p_test, int depth, int depth_max){
    if(this->precursors.size()==0){
        return false;
    }

    if(depth!=0 && this==p_test){
        return true;
    }
    else if(depth < depth_max && !this->visited_node){
        this->visited_node = true;
        bool test = false;
        for(int i=0; i<this->successors.size(); i++){
            test = test || this->successors[i]->test_cycle(p_test, depth+1, depth_max);
            if(test)
                return true;
        }
        return false;
    }
    else{
        return false;
    }
}
void Pave::add_precursors(Pave* p){
    IntervalVector intersection = this->box & p->box;
    if(!intersection.is_empty() && (intersection[0].is_degenerated() != intersection[1].is_degenerated())){
        bool flag = true;
        for(int i=0; i<this->precursors.size(); i++){
            if(this->precursors[i] == p)
                flag=false;
        }
        if(flag)
            this->precursors.push_back(p);
    }
}
void Pave::add_successors(Pave* p){
    IntervalVector intersection = this->box & p->box;
    if(!intersection.is_empty() && (intersection[0].is_degenerated() != intersection[1].is_degenerated())){
        bool flag = true;
        for(int i=0; i<this->successors.size(); i++){
            if(this->successors[i] == p)
                flag=false;
        }
        if(flag)
            this->successors.push_back(p);
    }
}
void Pave::clear_graph(){
    this->successors.clear();
    this->precursors.clear();
}

// ********************************************************************************
// ****************** TEST functions *********************************************
/**
 * @brief Pave::add_new_segment
 * @param b
 * Add new segment to queue_forward list of Pave
 * Used by
 */
void Pave::add_new_segment(Border &b, bool forward){
    if(forward){
        this->queue_forward.push_back(b);
    }
    else{
        this->queue_backward.push_back(b);
    }
}
