#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &box, ibex::Function *f): m_box(2)
{
    this->m_box = box;    // Box corresponding to the Pave
    this->m_borders.reserve(4);
    this->m_f = f;

    this->m_in_queue = false;

    // Border building
    IntervalVector coordinate(2);
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].lb()); this->m_borders.push_back(Border(coordinate, 0, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].ub()); this->m_borders.push_back(Border(coordinate, 1, this));
    coordinate[0] = box[0]; coordinate[1] = Interval(box[1].ub()); this->m_borders.push_back(Border(coordinate, 2, this));
    coordinate[1] = box[1]; coordinate[0] = Interval(box[0].lb()); this->m_borders.push_back(Border(coordinate, 3, this));

    m_full = false;
    m_empty = false;

    for(int i=0; i<2; i++)
        this->m_theta[i] = Interval::EMPTY_SET;

    Interval dx = box[1];
    Interval dy = 1.0*(1-pow(box[0], 2))*box[1]-box[0];

    Interval theta = atan2(dy, dx);
    if(theta==(-Interval::PI|Interval::PI)){
        vector<Interval> dR = this->rotate(Interval::PI, dx, dy);
        Interval thetaR = atan2(dR[1], dR[0]);
        if(thetaR.diam()<theta.diam()){
            this->m_theta[0] = (thetaR+Interval::PI) & (-Interval::PI | Interval::PI);
            this->m_theta[1] = (thetaR-Interval::PI) & (-Interval::PI | Interval::PI);
        }
        else{
            this->m_theta[0] = theta;
        }
    }
    else{
        this->m_theta[0] = theta;
    }
}

Pave::Pave(const Pave &p): m_box(2)
{
    this->m_box = p.m_box;    // Box corresponding to the Pave
    this->m_f = p.m_f;
    this->m_full = p.m_full;
    this->m_empty = p.m_empty;
    this->m_theta[0] = p.m_theta[0];
    this->m_theta[2] = p.m_theta[1];
    this->m_in_queue = false;

    for(int face = 0; face< 4; face++){
        this->m_borders.push_back(p.m_borders[face]); // Copy the border !
        this->m_borders[face].set_pave(this);
    }
}

Pave::~Pave(){
}

Pave& Pave::operator&=(const Pave &p){
    for(int face = 0; face <4; face++){
        this->m_borders[face] &= p.m_borders[face];
    }
    return *this;
}

void Pave::diff(const Pave &p){
    for(int face = 0; face<4; face++){
        this->m_borders[face].diff(p.m_borders[face]);
    }
    this->m_empty=false; // forces to recompute the value
    this->m_full=true;
}

void Pave::set_theta(ibex::Interval theta){
    this->m_theta[0] = Interval::EMPTY_SET;
    this->m_theta[1] = Interval::EMPTY_SET;

    if(theta.is_subset(-Interval::PI | Interval::PI)){
        this->m_theta[0] = theta;
    }
    else{
        this->m_theta[0] = (theta & (-Interval::PI | Interval::PI));

        if(!((theta + 2*Interval::PI) & (-Interval::PI | Interval::PI) ).is_empty())
            this->m_theta[1] =(theta + 2*Interval::PI);
        else if (!((theta - 2*Interval::PI) & (-Interval::PI | Interval::PI)).is_empty())
            this->m_theta[1] = (theta - 2*Interval::PI);
    }
}

void Pave::set_full(){
    for(int face=0; face<4; face++){
        this->m_borders[face].set_full();
    }
    this->m_full = true;
}

void Pave::set_empty(){
    for(int face=0; face<4; face++){
        this->m_borders[face].set_empty();
    }
    this->m_empty = true;
}

IntervalVector Pave::get_border_position(int face){
    IntervalVector position(2);
    position[0] = this->m_box[face%2];
    position[1] = this->m_box[(face+1)%2];
    return position;
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Pave::draw(bool filled, string color){
    // Draw the pave
    vibes::drawBox(this->m_box, color);

    this->draw_borders(filled);

    // Draw theta
    double size = 0.8*min(this->m_box[0].diam(), this->m_box[1].diam())/2.0;
    for(int i=0; i<2; i++){
        vibes::drawSector(this->m_box[0].mid(), this->m_box[1].mid(), size, size, (-this->m_theta[i].lb())*180.0/M_PI, (-this->m_theta[i].ub())*180.0/M_PI, "r[]");
    }
}
void Pave::draw_borders(bool filled){
    if(!filled){
        // Draw Segments
        for(int i=0; i<this->m_borders.size(); i++){
            this->m_borders[i].draw();
        }
    }
    else{
        // Draw Polygone
        vector<double> x, y;
        for(int i=0; i<this->m_borders.size(); i++){
            this->m_borders[i].get_points(x, y);
        }
        vibes::drawPolygon(x, y, "g[g]");
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result){
    // Create 4 new paves
    ibex::LargestFirst bisector(0.0, 0.5);
    std::pair<IntervalVector, IntervalVector> result_boxes = bisector.bisect(this->m_box);

    Pave *pave1 = new Pave(result_boxes.first, this->m_f); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, this->m_f); // Right or Down

    int indice1, indice2;

    if(pave1->m_box[0] == this->m_box[0]){
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
        if(this->m_borders[i].brothers().size()!=0){
            if(i!=indice1){
                pave1->m_borders[i].add_brothers(this->m_borders[i].brothers());
            }
            if(i!=indice2){
                pave2->m_borders[i].add_brothers(this->m_borders[i].brothers());
            }
        }
    }

    // Add each other to its brother list (pave1 <-> pave2)
    pave1->m_borders[indice1].add_brothers(&pave2->m_borders[indice2]);
    pave2->m_borders[indice2].add_brothers(&pave1->m_borders[indice1]);

    // Remove
    for(int i=0; i<4; i++){
        this->m_borders[i].update_brothers(&(pave1->m_borders[i]), &(pave2->m_borders[i]));
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

// ********************************************************************************
// ****************** UTILS functions *********************************************

double Pave::get_theta_diam(){
    double diam = 0.0;
    for(int i=0; i<2; i++){
        if(!this->m_theta[i].is_empty())
            diam += this->m_theta[i].diam();
    }
    return diam;
}

void Pave::remove_brothers(Pave* p, int face){
    for(int i=0; i<this->m_borders[face].brothers().size(); i++){
        if(this->m_borders[face].brothers()[i]->pave() == p){
            this->m_borders[face].remove_brother(i);
            return;
        }
    }
}

void Pave::remove_from_brothers(){
    for(int face=0; face<4; face++){
        for(int i=0; i<this->m_borders[face].brothers().size(); i++){
            this->m_borders[face].brothers()[i]->pave()->remove_brothers(this, (face+2)%4);
        }
    }
}

bool Pave::is_empty(){
    if(this->m_empty){
        return true;
    }
    else{
        for(int i=0; i<4; i++){
            if(!this->m_borders[i].is_empty()){
                return false;
            }
        }

        this->m_empty = true;
        return true;
    }
}

bool Pave::is_full(){
    if(!this->m_full){
        return false;
    }
    else{
        for(int face=0; face<4; face++){
            if(!this->m_borders[face].is_full()){
                this->m_full = false;
                return false;
            }
        }
        this->m_full = true;
        return true;
    }
}

vector<Pave*> Pave::get_brothers(int face){
    vector<Pave*> brothers_list;
    for(int i=0; i<this->m_borders[face].brothers().size(); i++){
        brothers_list.push_back(this->m_borders[face].brothers()[i]->pave());
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

//void Pave::tarjan_compute_successors(){
//    this->m_tarjan_index = 0;
//    this->m_tarjan_lowlink = 0;
//    this->m_tarjan_on_stack = false;
//    this->m_tarjan_successors.clear();

//    for(int face=0; face<4; face++){
//        for(int b=0; b<this->m_borders[face].brothers().size(); b++){
//            if(!(this->m_borders[face].brothers()[b]->segment_in() & this->m_borders[face].segment_out()).is_empty()){
//                this->m_tarjan_successors.push_back(this->m_borders[face].brothers()[b]->pave());
//            }
//        }
//    }
//    if(get_theta_diam()>=2*M_PI){
//        this->m_tarjan_successors.push_back(this);
//    }
//}

//void Pave::strongconnect(int &index, std::vector<Pave*> *S,std::vector<std::vector<Pave*>> *SCC){
//    this->m_tarjan_index = index;
//    this->m_tarjan_lowlink = index;
//    index++;
//    S->push_back(this);
//    this->m_tarjan_on_stack = true;
//    Pave *v = this;

//    for(int i=0; i<this->m_tarjan_successors.size(); i++){
//        Pave *w = this->m_tarjan_successors[i];
//        if(w->m_tarjan_index == 0){
//            w->strongconnect(index, S, SCC);
//            v->m_tarjan_lowlink = min(v->m_tarjan_lowlink, w->m_tarjan_lowlink);
//        }
//        else if(w->m_tarjan_on_stack){
//            v->m_tarjan_lowlink = min(v->m_tarjan_lowlink, w->m_tarjan_index);
//        }
//    }

//    if(v->m_tarjan_index == v->m_tarjan_lowlink){
//        std::vector<Pave *> C;
//        Pave *w;
////        cout << "NEW LIST C : v = " << v->m_tarjan_index << endl;
//        do{
//            w = S->back();
////            cout << "Node (" << w->m_box << ") " << '\t' << "= " << w->m_tarjan_index << "," << w->m_tarjan_lowlink << " size = " << w->m_tarjan_successors.size() << endl;
//            S->pop_back();
//            w->m_tarjan_on_stack = false;
//            C.push_back(w);
//        } while(w != v);

//        SCC->push_back(C);
//    }
//}
