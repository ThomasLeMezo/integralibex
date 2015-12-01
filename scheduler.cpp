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

        if(tmp->theta.diam() < epsilon_theta){
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

void Scheduler::add_segment(){
    this->pave_list[1000]->activate_pave();
}

void Scheduler::process(int max_iterations){
    int iterations=0;
    while(this->pave_queue.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = this->pave_queue.back();
        this->pave_queue.pop_back();
        pave->process();
    }

    cout << "queue size = " << this->pave_queue.size() << endl;
}

void Scheduler::draw(){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->draw();
    }
}

