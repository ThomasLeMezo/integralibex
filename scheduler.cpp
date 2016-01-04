#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(){
    this->draw_nb = 0;
}

Scheduler::~Scheduler(){
    for(int i=0; i<this->pave_list.size(); i++){
        delete(this->pave_list[i]);
    }
    for(int i=0; i<this->pave_list_empty.size(); i++){
        delete(this->pave_list_empty[i]);
    }
    this->pave_list.clear();
    this->pave_list_empty.clear();
}

// ********************************************************************************
// ****************** Setup Init functions ****************************************

void Scheduler::set_initial_pave(const IntervalVector &box, ibex::Function *f){
    Pave *p = new Pave(box, f);
    this->pave_list.push_back(p);
}

Pave* Scheduler::get_pave(double x, double y){
    IntervalVector position(2);
    position[0] = Interval(x);
    position[1] = Interval(y);

    for(int i=0; i<this->pave_list.size(); i++){
        if(!(position & this->pave_list[i]->m_box).is_empty()){
            return this->pave_list[i];
        }
    }
}

void Scheduler::set_full(){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->set_full();
    }
}

// ********************************************************************************
// ****************** Propagation functions ***************************************

void Scheduler::SIVIA(double epsilon_theta, int iterations_max, bool not_full_test){

    int iterations = 0;
    vector<Pave *> tmp_pave_list(this->pave_list);

    this->pave_list.clear();

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<iterations_max){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        double diam = 0.0;
        if(!(tmp->m_theta[0].is_empty()))
            diam += tmp->m_theta[0].diam();
        if(!(tmp->m_theta[1].is_empty()))
            diam += tmp->m_theta[1].diam();

        if(diam < epsilon_theta || (not_full_test && tmp->is_full())){
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
    int iterations = 0;
    while(this->pave_queue.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = this->pave_queue.front();
        this->pave_queue.erase(this->pave_queue.begin());

        bool change = this->utils.CtcContinuity(pave);
        if(change){
            vector<bool> changeBackward = this->utils.CtcPaveBackward(pave);

            // Warn scheduler to process new pave
            for(int face=0; face<4; face++){
                if(/*changeForward[face] ||*/ changeBackward[face]){
                    vector<Pave*> brothers_pave = pave->get_brothers(face);
                    for(int i=0; i<brothers_pave.size(); i++){
                        this->pave_queue.push_back(brothers_pave[i]);
                    }
                }
            }
        }

        if(iterations%100 == 0){
            cout << iterations << endl;
        }
    }
}

void Scheduler::process_SIVIA_cycle(int iterations_max, int pave_max, int process_iterations_max){
    if(this->pave_list.size()!=1)
        return;

    int iterations = 0;
    this->set_full();

    if(iterations < iterations_max){
        this->SIVIA(0.0, 4, false); // Start with 4 boxes
        this->set_full();
        for(int i=0; i<this->pave_list.size(); i++){
            this->pave_queue.push_back(this->pave_list[i]);
        }

        this->process(process_iterations_max);
        iterations++;
    }

    while(this->pave_list.size()<pave_max && this->pave_list.size()!=0 && iterations < iterations_max){
        this->SIVIA(0.0, 2*this->pave_list.size(), true);
        // Set full continuity
        this->set_full();
        for(int i=0; i<this->pave_list.size(); i++){
            this->pave_queue.push_back(this->pave_list[i]);
        }

        // Process the backward with the subpaving
        this->process(process_iterations_max);

        // Remove empty pave
        for(int i=0; i<this->pave_list.size(); i++){
            // || (this->utils.CtcNetwonPave(this->pave_list[i]))
            if(this->pave_list[i]->is_all_brother_empty()){
                this->pave_list[i]->remove_from_brothers();
                this->pave_list_empty.push_back(this->pave_list[i]);
                this->pave_list.erase(this->pave_list.begin() + i);
                if(i!=0)
                    i--;
            }
        }

        cout << "****** " << iterations <<  " ******" << endl;

        iterations++;
    }
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Scheduler::draw(int size, bool filled){

    stringstream ss;
    ss << "integralIbex" << this->draw_nb;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    for(int i=0; i<this->pave_list_empty.size(); i++){
        this->pave_list_empty[i]->draw(filled, "gray[]");
    }

    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->draw(filled);
    }

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    this->draw_nb++;
}

// ********************************************************************************
// ****************** TEST functions **********************************************

void Scheduler::print_pave_info(double x, double y, string color){

    Pave* p = this->get_pave(x, y);
    cout << "BOX = " << p->m_box << endl;
    cout << p << endl;
    for(int i= 0; i<p->m_borders.size(); i++){
        cout << "border " << i << '\t' << p->m_borders[i].position() << '\t' << p->m_borders[i].segment_in() << p->m_borders[i].segment_out()<< endl;
    }
    cout << "theta " << p->m_theta[0] << " " << p->m_theta[1] << endl;

    for(int i=0; i<p->m_borders.size(); i++){
        for(int j = 0; j<p->m_borders[i].brothers().size(); j++){
            cout << "border=" << i << " brother=" << j << " " << p->m_borders[i].brothers()[j]->pave() << endl;
        }
    }

    double r=0.5*min(p->m_box[0].diam(), p->m_box[1].diam())/2.0;

    vibes::drawCircle(p->m_box[0].mid(), p->m_box[1].mid(), r, color);

    cout << endl;
}
