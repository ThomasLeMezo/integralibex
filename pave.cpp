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

Pave::Pave(const IntervalVector &box, ibex::Function f): box(2)
{
    this->box = box;    // Box corresponding to the Pave
    this->borders.reserve(4);

    // Border building
    IntervalVector coordinate(2);
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].lb()); this->borders.push_back(Border(coordinate, 0, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].ub()); this->borders.push_back(Border(coordinate, 1, this));
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].ub()); this->borders.push_back(Border(coordinate, 2, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].lb()); this->borders.push_back(Border(coordinate, 3, this));

    for(int i=0; i<2; i++)
        this->theta[i] = Interval::EMPTY_SET;

    Interval dx = box[1];
    Interval dy = 1.0*(1-pow(box[0], 2))*box[1]-box[0];

    Interval theta = atan2(dy, dx);
    if(theta==(-Interval::PI|Interval::PI)){
        vector<Interval> dR = this->rotate(Interval::PI, dx, dy);
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

    full = false;
    empty = true;
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
        this->borders[face].set_full();
    }
}

void Pave::set_full_continuity(){
    this->full = true;
    this->empty = false;

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

void Pave::draw(bool filled, string color){
    // Draw the pave
    vibes::drawBox(this->box, color);

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

    bool full = this->is_full();
    bool empty = this->is_empty();

    pave1->set_full(full);
    pave1->set_empty(empty);

    pave2->set_full(full);
    pave2->set_empty(empty);

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
            }
            if(i!=indice2){
                pave2->borders[i].add_brothers(this->borders[i].brothers);
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

        if(!seg_out[0].is_degenerated() || !seg_out[1].is_degenerated() || !seg_out[2].is_degenerated()){
            flow_in[face] = true;
        }
        for(int i = 0; i<3; i++){
            if(!seg_out[i].is_degenerated()){
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

// ********************************************************************************
// ****************** UTILS functions *********************************************

double Pave::get_theta_diam(){
    double diam = 0.0;
    for(int i=0; i<2; i++){
        if(!this->theta[i].is_empty())
            diam += this->theta[i].diam();
    }
    return diam;
}

bool Pave::get_brother_empty(int level){
    // ToDo : improve function by adding visited_node option
    if(!this->is_empty())
        return false;

    for(int face=0; face<4; face++){
        for(int i=0; i<this->borders[face].brothers.size(); i++){
            if(level == 0){
                if(!this->borders[face].brothers[i]->pave->is_empty()){
                    return false;
                }
            }
            else{
                if(!this->borders[face].brothers[i]->pave->get_brother_empty(level-1)){
                    return false;
                }
            }
        }
    }

    return true;
}

bool Pave::all_brothers_full(int level){
    if(!this->is_full())
        return false;

    for(int face=0; face<4; face++){
        for(int i=0; i<this->borders[face].brothers.size(); i++){
            if(level == 0){
                if(!this->borders[face].brothers[i]->pave->is_full()){
                    return false;
                }
            }
            else{
                if(!this->borders[face].brothers[i]->pave->all_brothers_full(level-1)){
                    return false;
                }
            }
        }
    }

    return true;
}

void Pave::remove_brothers(Pave* p, int face){
    for(int i=0; i<this->borders[face].brothers.size(); i++){
        if(this->borders[face].brothers[i]->pave == p){
            this->borders[face].brothers.erase(this->borders[face].brothers.begin()+i);
            return;
        }
    }
}

void Pave::remove_from_brothers(){
    for(int face=0; face<4; face++){
        for(int i=0; i<this->borders[face].brothers.size(); i++){
            this->borders[face].brothers[i]->pave->remove_brothers(this, (face+2)%4);
        }
    }
}

bool Pave::is_empty(){
    if(this->empty)
        return true;

    for(int i=0; i<4; i++){
        if(!this->borders[i].is_empty()){
            this->empty = false; // Normaly useless
            return false;
        }
    }

    this->empty = true;
    return true;
}

bool Pave::is_full(){
    if(!this->full){
        return false;
    }
    else{
        for(int face=0; face<4; face++){
            if(!this->borders[face].is_full()){
                this->full = false;
                return false;
            }
        }
        return true;
    }
}

void Pave::set_empty(bool val){
    this->empty = val;
}

bool Pave::set_full(bool val){
    this->full = val;
}

bool Pave::netwon_test(){
    return false;

    if(this->all_brothers_full()){
        IntervalVector box_tmp = this->box;
        this->scheduler->utils.contract_newton->contract(box_tmp);
        if(!box_tmp.is_empty())
            this->scheduler->utils.contract_newton->contract(box_tmp);

        if(!box_tmp.is_empty() && box_tmp[0].is_degenerated() && box_tmp[1].is_degenerated()){
            cout << "-->" << this->box << "<>" << box_tmp << endl;

            for(int face=0; face<4; face++){
                this->borders[face].segment=Interval::EMPTY_SET;
            }

            return true;
        }
    }
    return false;
}

bool Pave::copy_segment(Pave *p){
    if(this->box == p->box){
        for(int face=0; face < 4; face++){
            this->borders[face].segment = p->borders[face].segment;
        }
        return true;
    }
    else{
        return false;
    }
}

bool Pave::equal_segment(Pave *p){
    for(int face = 0; face < 4; face++){
        if(this->borders[face].segment != p->borders[face].segment)
            return false;
    }
    return true;
}

vector<Pave*> Pave::get_brothers(int face){
    vector<Pave*> brothers_list;
    for(int i=0; i<this->borders[face].brothers.size(); i++){
        brothers_list.push_back(this->borders[face].brothers[i]->pave);
    }
    return brothers_list;
}

std::vector<ibex::Interval> Pave::rotate(const ibex::Interval &theta, const ibex::Interval &x, const ibex::Interval &y){
    Interval xR = x*cos(theta) - y*sin(theta);
    Interval yR = x*sin(theta) + y*cos(theta);
    vector<Interval> list;
    list.push_back(xR);
    list.push_back(yR);
    return list;
}
