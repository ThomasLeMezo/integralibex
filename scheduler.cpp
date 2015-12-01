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
    cout << "Border size = " << this->pave_list[0]->borders.size() << endl;
}

void Scheduler::add_to_queue(Pave* pave){
    this->pave_queue.push_back(pave);
}

void Scheduler::SIVIA(double epsilon_theta, int iterations_max){

    int iterations = 0;

    for(int i=0; i<1000; i++){
        Pave* tmp = this->pave_list.front();
        this->pave_list.erase(this->pave_list.begin());
        tmp->bisect(this->pave_list);
    }
}

void Scheduler::add_segment(const IntervalVector &box){

}

void Scheduler::process(){
    while(this->pave_queue.size() != 0){
        Pave *pave = this->pave_queue.back();
        this->pave_queue.pop_back();
        pave->process();
    }
}

void Scheduler::draw(){
    cout << "Border size = " << this->pave_list[0]->borders.size() << endl;
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->draw();
    }
}

