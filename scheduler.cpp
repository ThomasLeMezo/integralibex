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

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<iterations_max){


        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        double diam = 0.0;
        if(!(tmp->theta[0].is_empty()))
            diam += tmp->theta[0].diam();
        if(!(tmp->theta[1].is_empty()))
            diam += tmp->theta[1].diam();

        if(diam < epsilon_theta){
            this->pave_list.push_back(tmp);
            iterations++;
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
        Pave *pave = this->pave_queue.front();
        this->pave_queue.erase(this->pave_queue.begin());

        pave->process();

        if(iterations%100 == 0){
            cout << iterations << endl;
        }
    }
//    cout << "queue size = " << this->pave_queue.size() << endl;
}

void Scheduler::draw(){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->draw(true);
    }
}

void Scheduler::print_pave_info(double x, double y, string color){

    Pave* p = this->get_pave(x, y);
    cout << "BOX = " << p->box << endl;
    cout << p << endl;
    for(int i= 0; i<p->borders.size(); i++){
        cout << "border " << i << '\t' << p->borders[i].position << '\t' << p->borders[i].segment << endl;
    }
    cout << "theta " << p->theta[0] << " " << p->theta[1] << endl;

    for(int i=0; i<p->borders.size(); i++){
        for(int j = 0; j<p->borders[i].brothers.size(); j++){
            cout << "border=" << i << " brother=" << j << " " << p->borders[i].brothers[j]->pave << endl;
        }
    }

    double r=0.5*min(p->box[0].diam(), p->box[1].diam())/2.0;

    vibes::drawCircle(p->box[0].mid(), p->box[1].mid(), r, color);

    cout << endl;
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


