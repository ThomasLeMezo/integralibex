#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(){
    this->draw_nb = 0;
}

// ********************************************************************************
// ****************** Setup Init functions ****************************************

void Scheduler::set_initial_pave(const IntervalVector &box){
    Pave *p = new Pave(box, this);
    this->pave_list.push_back(p);
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

void Scheduler::set_full_continuity(){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->set_full_continuity();
    }
}

// ********************************************************************************
// ****************** Propagation functions ***************************************

void Scheduler::add_to_queue(Pave* pave, bool forward){
    if(forward){
        this->pave_queue_forward.push_back(pave);
    }
    else{
        this->pave_queue_backward.push_back(pave);
    }
}

void Scheduler::SIVIA(double epsilon_theta, int iterations_max, bool not_full_test){

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
    int iterations=0;
    while(this->pave_queue_forward.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = this->pave_queue_forward.front();
        this->pave_queue_forward.erase(this->pave_queue_forward.begin());

        pave->process_forward();

        if(iterations%100 == 0){
            cout << iterations << endl;
        }
    }
    //    cout << "queue size = " << this->pave_queue_forward.size() << endl;
}

void Scheduler::process_backward(int max_iterations){
    for(int i=0; i<this->pave_list.size(); i++){
        this->pave_list[i]->compute_flow();
    }

    int iterations=0;
    while(this->pave_queue_backward.size() != 0 & iterations < max_iterations){
        // ToDo : Backward/forward combination !!
        iterations++;
        Pave *pave = this->pave_queue_backward.front();
        this->pave_queue_backward.erase(this->pave_queue_backward.begin());

        pave->process_backward();

        if(iterations%100 == 0){
            cout << iterations << endl;
        }
    }
}

void Scheduler::process_graph(int iterations_max, int pave_max){

    vector<Pave*> pave_list_tmp;
    int nb_graph_node_before = this->pave_list.size();
    int nb_graph_node_after = 0;
    int iterations_fix_pt = 0;

    while(this->pave_list.size()<pave_max & this->pave_list.size()!=0){
        this->SIVIA(M_PI/20.0, 2*this->pave_list.size(), false);
        nb_graph_node_before = this->pave_list.size();
        nb_graph_node_after = 0;

        while(nb_graph_node_before != nb_graph_node_after & iterations_fix_pt < iterations_max){
            nb_graph_node_before = this->pave_list.size();
            nb_graph_node_after = 0;
            for(int i=0; i<this->pave_list.size(); i++){
                this->pave_list[i]->compute_successors();
            }

            for(int i=0; i<this->pave_list.size(); i++){
                for(int j=0; j<this->pave_list.size(); j++)
                    this->pave_list[j]->visited_node = false;

                if(this->pave_list[i]->test_cycle(this->pave_list[i], 0, this->pave_list.size())){
                    pave_list_tmp.push_back(this->pave_list[i]);
                    nb_graph_node_after++;
                }
            }

            this->pave_list = pave_list_tmp;
            pave_list_tmp.clear();
            for(int i=0; i<this->pave_list.size(); i++){
                this->pave_list[i]->clear_graph();
            }
            iterations_fix_pt ++;
        }
        cout << "iterations = " << iterations_fix_pt << " Nb pave = " << this->pave_list.size() << endl;
    }
}

void Scheduler::process_SIVIA_cycle(int iterations_max, int pave_max, int backward_iterations_max){

    if(this->pave_list.size()!=1)
        return;

    int iterations = 1;

    this->SIVIA(0.0, 4, false); // Start with 4 boxes
    this->set_full_continuity();
    this->process_backward(backward_iterations_max);

    while(this->pave_list.size()<pave_max && this->pave_list.size()!=0 && iterations < iterations_max){
        //cout << "SIVIA" << endl;
        this->SIVIA(0.0, 2*this->pave_list.size(), true);
        // Set full continuity
        //cout << "CONTINUITY" << endl;
        this->set_full_continuity();

        // Process the backward with the subpaving
        //cout << "PROCESS" << endl;
        this->process_backward(backward_iterations_max);

        // Remove empty pave
        //cout << "REMOVE" << endl;
        for(int i=0; i<this->pave_list.size(); i++){
            if(this->pave_list[i]->get_brother_empty()){
                this->pave_list[i]->remove_from_brothers();
                this->empty_pave_list.push_back(this->pave_list[i]);
                this->pave_list.erase(this->pave_list.begin() + i);
                if(i!=0)
                    i--;
            }
        }

        if(iterations%100==0){
            cout << "****** " << iterations <<  " ******" << endl;
        }

        iterations++;
    }
}

// ********************************************************************************
// ****************** Draw functions **********************************************

void Scheduler::draw(int size, bool filled){

    stringstream ss;
    ss << "integralIbex" << this->draw_nb;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    for(int i=0; i<this->empty_pave_list.size(); i++){
        this->empty_pave_list[i]->draw(filled, "gray[]");
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
