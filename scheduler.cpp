#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(){

}

void Scheduler::set_initial_pave(const IntervalVector &box){
    Pave *p = new Pave(box, this);
    this->pave_list.push_back(p);
}

void Scheduler::add_to_queue(Pave* pave){
    this->pave_queue.push_back(pave);
}

void Scheduler::SIVIA(double epsilon_theta, int iterations_max){

    int iterations = 0;
    vector<Pave *> tmp_pave_list(this->pave_list);

    this->pave_list.clear();

    while(tmp_pave_list.size()!=0 & iterations<iterations_max){
        iterations++;

        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        double diam = tmp->theta[0].diam();
        if(tmp->theta.size()==2){
            diam += tmp->theta[1].diam();
        }

        if(diam < epsilon_theta){
            this->pave_list.push_back(tmp);
        }
        else{
            tmp->bisect(tmp_pave_list);
        }
    }

    for(int i=0; i<tmp_pave_list.size(); i++){
        this->pave_list.push_back(tmp_pave_list[i]);
    }
}

void Scheduler::process(int max_iterations){
    int iterations=0;
    while(this->pave_queue.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = this->pave_queue.back();
        this->pave_queue.pop_back();
        pave->process();

        if(iterations%100 == 0){
            cout << iterations << endl;
        }
    }

    cout << "queue size = " << this->pave_queue.size() << endl;
}

void Scheduler::draw(){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->draw();
    }
}

std::vector<ibex::Interval> Scheduler::rotate(ibex::Interval theta, ibex::Interval x, ibex::Interval y){
    Interval xR = cos(theta)*x -sin(theta)*y;
    Interval yR = sin(theta)*x + cos(theta)*y;
    vector<Interval> list;
    list.push_back(xR);
    list.push_back(yR);
    return list;
}

void Scheduler::print_pave_info(double x, double y){

    Pave* p = this->get_pave(x, y);
    for(int i= 0; i<p->borders.size(); i++){
        cout << "border " << i << " " << p->borders[i].position << p->borders[i].segment << endl;
    }
}

Pave* Scheduler::get_pave(double x, double y){
    IntervalVector position(2);
    position[0] = Interval(x);
    position[1] = Interval(y);

    for(int i=0; i<this->pave_list.size(); i++){
        if(!(position & this->pave_list[i]->box).is_empty()){
            return this->pave_list[i];
        }
    }
}

void Scheduler::add_segment(int id_box){
    this->pave_list[id_box]->activate_pave();
}

void Scheduler::add_segment(double x, double y){
    this->get_pave(x, y)->activate_pave();
}

/**
 ** CtcPropagateFront supposed that the down left box corner is (0,0)
 **
*/
void Scheduler::CtcPropagateFront(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy){
    Interval Y(0.0, dy);

    Interval x = Interval(-dx, dx);
    Interval y = Interval(dx);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = Interval::HALF_PI - theta;

    contract_polar.contract(x, y, rho, theta2);

    Sk = (Sk + x ) & Y;
}

void Scheduler::CtcPropagateLeftSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dy){
    Interval x = Sk;
    Interval y = Interval(0.0, dy);
    Interval rho = Interval::POS_REALS;
    Interval theta2 = Interval::PI - theta;

    contract_polar.contract(x, y, rho, theta2);

    Sk = y;
}

void Scheduler::CtcPropagateRightSide(ibex::Interval &Sk, const ibex::Interval &theta, const double &dx, const double &dy){
    /** Apply a symetry to CtcPropagateLeftSide
     ** theta -> -theta
     ** [X] -> dx - [X]
    */

    CtcPropagateLeftSide(Sk, -theta, dy);
    Sk = dx - Sk;
}
