#include "scheduler.h"
#include "vibes.h"
#include "ibex.h"

using namespace std;
using namespace ibex;

Scheduler::Scheduler(){
    this->m_draw_nb = 0;
}

Scheduler::~Scheduler(){    
    for(int i=0; i<this->m_global_pave_list.size(); i++){
        for(int j=0; j<this->m_global_pave_list[j].size(); j++){
            delete(this->m_global_pave_list[i][j]);
        }
    }

    for(int i=0; i<this->m_global_pave_list_empty.size(); i++){
        for(int j=0; j<this->m_global_pave_list_empty[j].size(); j++){
            delete(this->m_global_pave_list_empty[i][j]);
        }
    }
}

// ********************************************************************************
// ****************** Setup Init functions ****************************************

void Scheduler::set_initial_pave(const IntervalVector &box, ibex::Function *f){
    this->m_global_pave_list.clear();

    Pave *p = new Pave(box, f);
    vector<Pave*> pave_list, pave_queue, pave_empty;
    pave_list.push_back(p);

    this->m_global_pave_list.push_back(pave_list);
    this->m_global_pave_queue.push_back(pave_queue);
    this->m_global_pave_list_empty.push_back(pave_empty);
}

Pave* Scheduler::get_pave(std::vector<Pave*> &pave_list, double x, double y){
    IntervalVector position(2);
    position[0] = Interval(x);
    position[1] = Interval(y);

    for(auto &pave : pave_list){
        if(!(position & pave->m_box).is_empty()){
            return pave;
        }
    }
    return NULL;
}

void Scheduler::activate_pave(std::vector<Pave*> &pave_list, std::vector<Pave*> &pave_queue, double x, double y){
    Pave *pave = get_pave(pave_list, x, y);
    pave->set_full();
    for(int face=0; face<4; face++){
        vector<Pave*> brothers_pave = pave->get_brothers(face);
        for(auto &pave : brothers_pave){
            pave_queue.push_back(pave);
        }
    }
}

void Scheduler::set_full(std::vector<Pave*> &pave_list){
    for(auto &pave : pave_list){
        pave->set_full();
    }
}

// ********************************************************************************
// ****************** Propagation functions ***************************************

void Scheduler::SIVIA(std::vector<Pave*> &pave_list, std::vector<Pave*> &pave_queue, double epsilon_theta, int iterations_max, bool backward){

    int iterations = 0;
    vector<Pave *> tmp_pave_list(pave_list);

    pave_list.clear();

    while(tmp_pave_list.size()!=0 & (iterations+tmp_pave_list.size())<iterations_max){
        Pave* tmp = tmp_pave_list.front();
        tmp_pave_list.erase(tmp_pave_list.begin());

        double diam = 0.0;
        if(!(tmp->m_theta[0].is_empty()))
            diam += tmp->m_theta[0].diam();
        if(!(tmp->m_theta[1].is_empty()))
            diam += tmp->m_theta[1].diam();

        if(diam < epsilon_theta){// || (not_full_test && tmp->is_full() && diam < M_PI)){
            pave_list.push_back(tmp);
            iterations++;
        }
        else{
            tmp->bisect(tmp_pave_list);
        }
    }

    for(int i=0; i<tmp_pave_list.size(); i++){
        if(backward){
            tmp_pave_list[i]->set_full();
            pave_queue.push_back(tmp_pave_list[i]);
        }
        pave_list.push_back(tmp_pave_list[i]);
    }
}

void Scheduler::process(std::vector<Pave*> &pave_queue, int max_iterations, bool backward){
    int iterations = 0;
    while(pave_queue.size() != 0 & iterations < max_iterations){
        iterations++;
        Pave *pave = pave_queue.front();
        pave_queue.erase(pave_queue.begin());
        pave->m_in_queue = false;

        bool change = this->m_utils.CtcContinuity(pave, backward);
        if(change){
            this->m_utils.CtcPaveConsistency(pave, backward);

            // Warn scheduler to process new pave
            for(int face=0; face<4; face++){
                vector<Pave*> brothers_pave = pave->get_brothers(face);
                for(int i=0; i<brothers_pave.size(); i++){
                    if(brothers_pave[i]->m_in_queue == false){
                        pave_queue.push_back(brothers_pave[i]);
                        brothers_pave[i]->m_in_queue = true;
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
    if(this->m_global_pave_list.size()!=1 && this->m_global_pave_list[0].size() !=1)
        return;

    int iterations = 0;
    this->set_full(this->m_global_pave_list[0]);

    if(iterations < iterations_max){
        this->SIVIA(this->m_global_pave_list[0], this->m_global_pave_queue[0], 0.0, 4, true); // Start with 4 boxes
        this->process(this->m_global_pave_queue[0], process_iterations_max,  true);
        iterations++;
    }

    for(int nb_graph=0; nb_graph<this->m_global_pave_list.size(); nb_graph++){
        while(this->m_global_pave_list[nb_graph].size()<pave_max && this->m_global_pave_list[nb_graph].size()!=0 && iterations < iterations_max){
            this->SIVIA(this->m_global_pave_list[nb_graph], this->m_global_pave_queue[nb_graph], 0.0, 2*this->m_global_pave_list[nb_graph].size(), true);

            // Process the backward with the subpaving
            cout << "****** " << iterations <<  " ****** (" << this->m_global_pave_list[nb_graph].size() << ")" << endl;
            this->process(this->m_global_pave_queue[nb_graph], process_iterations_max, true);

            // Remove empty pave
            for(int i=0; i<this->m_global_pave_list[nb_graph].size(); i++){
                if(this->m_global_pave_list[nb_graph][i]->is_empty()){
                    this->m_global_pave_list[nb_graph][i]->remove_from_brothers();
                    this->m_global_pave_list_empty[nb_graph].push_back(this->m_global_pave_list[nb_graph][i]);
                    this->m_global_pave_list[nb_graph].erase(this->m_global_pave_list[nb_graph].begin() + i);
                    if(i!=0)
                        i--;
                }
            }

            if(this->m_global_pave_list[nb_graph].size()==0)
                break;

            // ***************************************************
            // Copy graph & propagate one Pave + intersect with cycle
            // Find the first non-full & non-empty Pave
#if 0
            cout << "REMOVE INSIDE" << endl;
            Pave *pave_start;
            for(int i=0; i<this->m_global_pave_list[nb_graph].size(); i++){
                if(!this->m_global_pave_list[nb_graph][i]->is_empty() && !this->m_global_pave_list[nb_graph][i]->is_full()){
                    pave_start = this->m_global_pave_list[nb_graph][i];
                    break;
                }
            }

            vector<Pave*> pave_list = copy_graph(nb_graph, true, pave_start); //
            vector<Pave*> pave_queue;

            for(int face=0; face<4; face++){    // Warn scheduler
                vector<Pave*> brothers_pave = pave_start->m_copy_node->get_brothers(face);
                for(auto &pave : brothers_pave){
                    pave_queue.push_back(pave);
                }
            }
            this->process(pave_queue, process_iterations_max, false);


            // Make the difference of the two pave_list
            vector<Pave*> pave_list2 = copy_graph(nb_graph, false, NULL);

            for(int i=0; i<pave_list.size(); i++){
                *(this->m_global_pave_list[nb_graph][i]) &= *(pave_list[i]);
                pave_list2[i]->diff(*(this->m_global_pave_list[nb_graph][i]));
            }

            this->m_global_pave_list.push_back(pave_list2);
            vector<Pave*> pave_queue2, pave_empty2;
            this->m_global_pave_list_empty.push_back(pave_empty2);
            this->m_global_pave_queue.push_back(pave_queue2);

            // delete pave_list
            for(int i=0; i<pave_list.size(); i++){
                delete(pave_list[i]);
            }
#endif

            iterations++;
        }
        iterations = 0;
    }
}

vector<Pave*> Scheduler::copy_graph(int graph, bool empty, Pave* pave_keep){
    vector<Pave*> pave_list;
    for(int i=0; i<this->m_global_pave_list[graph].size(); i++){
        Pave *p = new Pave(*(this->m_global_pave_list[graph][i]));
        if(empty & pave_keep != this->m_global_pave_list[graph][i])
            p->set_empty();
        pave_list.push_back(p);
        this->m_global_pave_list[graph][i]->m_copy_node = p;
    }

    for(int i=0; i<this->m_global_pave_list[graph].size(); i++){
        Pave* pave_root = this->m_global_pave_list[graph][i];
        Pave* pave_copy = pave_list[i];

        for(int face = 0; face<4; face++){
            for(int j=0; j<pave_root->m_borders[face].brothers().size(); j++){
                pave_copy->m_borders[face].brothers()[j]->set_pave(pave_root->m_borders[face].brothers()[j]->pave()->m_copy_node);
            }
        }
    }

    return pave_list;
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Scheduler::draw(int size, bool filled){

    for(int graph=0; graph<this->m_global_pave_list.size(); graph++){
        stringstream ss;
        ss << "integralIbex" << this->m_draw_nb << "-" << graph;
        vibes::newFigure(ss.str());
        vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

        for(int i=0; i<this->m_global_pave_list_empty[graph].size(); i++){
            this->m_global_pave_list_empty[graph][i]->draw(filled, "gray[]");
        }

        for(int i=0; i<this->m_global_pave_list[graph].size(); i++){
            this->m_global_pave_list[graph][i]->draw(filled);
        }

        vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    }
    this->m_draw_nb++;
}

// ********************************************************************************
// ****************** TEST functions **********************************************

void Scheduler::print_pave_info(int graph, double x, double y, string color){

    Pave* p = this->get_pave(this->m_global_pave_list[graph],x, y);
    if(p==NULL){
        p = this->get_pave(this->m_global_pave_list_empty[graph],x, y);
        if(p==NULL){
            cout << "PAVE NOT FOUND" << endl;
            return;
        }
    }

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
