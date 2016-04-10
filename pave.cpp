#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "iomanip"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &position, const std::vector<ibex::Function*> &f_list, const ibex::IntervalVector &u, bool diseable_singeleton, bool active):
    m_position(2),
    m_u_iv(2)
{
    m_position = position;    // Box corresponding to the Pave
    m_borders.reserve(4);
    m_f_list = f_list;
    m_active_function = 0;
    m_active = active;
    m_diseable_singeleton = diseable_singeleton;

    m_u_iv = u;

    m_in_queue = false;
    m_copy_node = NULL;

    m_first_process = false;

    m_inner = false;

    // Border building
    IntervalVector coordinate(2);
    coordinate[0] = position[0]; coordinate[1] = Interval(position[1].lb()); m_borders.push_back(new Border(coordinate, 0, this));
    coordinate[1] = position[1]; coordinate[0] = Interval(position[0].ub()); m_borders.push_back(new Border(coordinate, 1, this));
    coordinate[0] = position[0]; coordinate[1] = Interval(position[1].ub()); m_borders.push_back(new Border(coordinate, 2, this));
    coordinate[1] = position[1]; coordinate[0] = Interval(position[0].lb()); m_borders.push_back(new Border(coordinate, 3, this));

    m_full = true;
    m_empty = false;

    /////////////////////////////// CASE THETA ///////////////////////////////
    for(int i=0; i<f_list.size(); i++){
        std::vector<ibex::Interval> theta = compute_theta(f_list[i], u);
        m_theta_list.push_back(theta);
    }


    /////////////////////////////// CASE CMD U ///////////////////////////////
//    if(u[0] == Interval::ZERO && u[1] == Interval::ZERO){
        m_u.push_back(Interval::ZERO);
//    }
//    else{
//        Interval theta_u = atan2(u[1], u[0]);
//        if(theta_u==(-Interval::PI|Interval::PI)){
//            Interval thetaR_u = atan2(-u[1], -u[0]); // PI rotation ({dx, dy} -> {-dx, -dy})
//            if(thetaR_u.diam()<theta_u.diam()){
//                if(thetaR_u.is_subset(-Interval::PI | Interval::PI)){
//                    m_u.push_back((thetaR_u + Interval::PI) & (Interval::ZERO | Interval::PI)); // theta[0] in [0, pi]
//                    m_u.push_back(((thetaR_u + Interval::PI) & (Interval::PI | Interval::TWO_PI)) - Interval::TWO_PI); // theta[1] in [-pi, 0]
//                }
//                else{
//                    cout << "****************** ERROR ******************" << endl;
//                }

//            }
//            else{
//                m_u.push_back(theta_u);
//            }
//        }
//        else if(theta_u.is_empty()){
//            m_u.push_back(-Interval::PI|Interval::PI);
//        }
//        else{
//            m_u.push_back(theta_u);
//        }
//        if(m_u[0].is_empty()){
//            cout << "ERROR - Pave "<< theta_u << u[0] << u[1] << m_position << endl;
//        }
//    }
}

const std::vector<ibex::Interval> Pave::compute_theta(ibex::Function *f, ibex::IntervalVector u){
    std::vector<ibex::Interval> theta_list;

    for(int i=0; i<2; i++)
        theta_list.push_back(Interval::EMPTY_SET);
    if(f!=NULL){
        IntervalVector dposition = f->eval_vector(m_position) + u;

        Interval dx = dposition[0];
        Interval dy = dposition[1];

        Interval theta = atan2(dy, dx);

        if(theta==(-Interval::PI|Interval::PI)){
            Interval thetaR = atan2(-dy, -dx); // PI rotation ({dx, dy} -> {-dx, -dy})
            if(thetaR.diam()<theta.diam()){
                if(thetaR.is_subset(-Interval::PI | Interval::PI)){
                    theta_list[0] = (thetaR + Interval::PI) & (Interval::ZERO | Interval::PI); // theta[0] in [0, pi]
                    theta_list[1] = ((thetaR + Interval::PI) & (Interval::PI | Interval::TWO_PI)) - Interval::TWO_PI; // theta[1] in [-pi, 0]
                }
                else{
                    cout << "****************** ERROR ******************" << endl;
                }

            }
            else{
                theta_list[0] = theta;
            }
        }
        else if(theta.is_empty()){
            theta_list[0] = -Interval::PI|Interval::PI;
        }
        else{
            theta_list[0] = theta;
        }
        if(theta_list[0].is_empty()){
            cout << "ERROR - Pave "<< theta << dx << dy << m_position << endl;
        }
    }

    return theta_list;
}

Pave::Pave(const Pave *p):
    m_position(2),
    m_u_iv(2)
{
    m_position = p->get_position();    // Box corresponding to the Pave
    m_f_list = p->get_f_list();
    m_active_function = p->get_active_function();
    m_u_iv = p->get_u_iv();
    m_full = true; // Force to recompute results
    m_empty = false;
    m_in_queue = false;
    m_first_process = false;
    m_inner = p->get_inner();
    m_active = p->is_active();
    m_diseable_singeleton = p->get_diseable_singelton();

    m_theta_list = p->get_theta_list();
    m_u = p->get_u();

    for(int face = 0; face < 4; face++){
        Border *b = new Border(p->get_border_const(face));
        m_borders.push_back(b); // Copy the border !
        m_borders[face]->set_pave(this);
    }
    m_copy_node = NULL;
}

Pave::~Pave(){
    for(int face=0; face<4; face++){
        delete(m_borders[face]);
    }
}

Pave& Pave::operator&=(const Pave &p){
    for(int face = 0; face <4; face++){
        *(m_borders[face]) &= *(p.get_border_const(face));
    }
    return *this;
}

Pave& Pave::operator|=(const Pave &p){
    for(int face = 0; face <4; face++){
        *(m_borders[face]) |= *(p.get_border_const(face));
    }
    return *this;
}

bool Pave::inter(const Pave &p){
    bool change = false;
    for(int face = 0; face <4; face++){
        if(m_borders[face]->inter(*(p.get_border_const(face))))
            change = true;
    }
    return change;
}

bool Pave::diff(const Pave &p){
    bool change = false;
    for(int face = 0; face<4; face++){
        if(m_borders[face]->diff(*(p.get_border_const(face)))){
            change = true;
        }
    }
    m_empty=false; // forces to recompute the value
    m_full=true;
}

void Pave::set_theta(ibex::Interval theta){
    m_theta_list.clear();

    std::vector<ibex::Interval> thetas;
    for(int i=0; i<2; i++)
        thetas.push_back(Interval::EMPTY_SET);

    if(theta.is_subset(-Interval::PI | Interval::PI)){
        thetas[0] = theta;
    }
    else{
        thetas[0] = (theta & (-Interval::PI | Interval::PI));

        if(!((theta + 2*Interval::PI) & (-Interval::PI | Interval::PI) ).is_empty())
            thetas[1] =(theta + 2*Interval::PI);
        else if (!((theta - 2*Interval::PI) & (-Interval::PI | Interval::PI)).is_empty())
            thetas[1] = (theta - 2*Interval::PI);
    }

    m_theta_list.push_back(thetas);
}

void Pave::set_full(){
    for(int face=0; face<4; face++){
        m_borders[face]->set_full();
    }
    m_full = true;
}

void Pave::set_full_in(){
    for(int face=0; face<4; face++){
        m_borders[face]->set_full_segment_in();
    }
    m_full = true;
}

void Pave::set_full_out(){
    for(int face=0; face<4; face++){
        m_borders[face]->set_full_segment_out();
    }
    m_full = true;
}

void Pave::set_empty(){
    for(int face=0; face<4; face++){
        m_borders[face]->set_empty();
    }
    m_empty = true;
    m_full = false;
}

// ********************************************************************************
// ****************** Drawing functions *******************************************

void Pave::draw_position(){
    double size = 0.5*min(m_position[0].diam(), m_position[1].diam())/2.0;
    vibes::drawCircle(m_position[0].mid(), m_position[1].mid(), size, "b[b]");
}

void Pave::draw(bool filled, string color, bool borders_only, bool cmd_u){
    // Draw the pave
    if(borders_only){
        draw_borders(filled, "[#00FF00AA]");
    }
    else{

        vibes::drawBox(m_position, color);

        if(m_inner)
            draw_borders(filled, "g[g]");
        else
            draw_borders(filled, "y[y]");
        // Draw theta

        double size = 0.8*min(m_position[0].diam(), m_position[1].diam())/2.0;

        if(cmd_u){
            // ToDo: check m_u input m_u[0] & m_u[1] add without check
//            Interval theta_u = ((m_theta[0] | m_theta[1]).lb() + m_u[0]) & ((m_theta[0] | m_theta[1]).ub() + m_u[1]);
//            Interval theta_u_bwd = Interval::PI + ((Interval((m_theta[0] | m_theta[1]).lb()) + m_u[0]) & (Interval((m_theta[0] | m_theta[1]).ub()) + m_u[1]));

//            for(int face =0; face<4; face++){
//                if(!get_border(face)->get_segment_in().is_empty()){

//                    IntervalVector segment_in = get_border(face)->get_segment_in_2D();

//                    vibes::drawSector(segment_in[0].lb(), segment_in[1].lb(), size*0.5, size*0.5, (-theta_u.lb())*180.0/M_PI, (-theta_u.ub())*180.0/M_PI, "r[]");
//                    vibes::drawSector(segment_in[0].ub(), segment_in[1].ub(), size*0.5, size*0.5, (-theta_u.lb())*180.0/M_PI, (-theta_u.ub())*180.0/M_PI, "r[]");
//                }
//                if(!get_border(face)->get_segment_out().is_empty()){
//                    IntervalVector segment_out = get_border(face)->get_segment_out_2D();

//                    vibes::drawSector(segment_out[0].lb(), segment_out[1].lb(), size*0.3, size*0.3, (-theta_u_bwd.lb())*180.0/M_PI, (-theta_u_bwd.ub())*180.0/M_PI, "b[]");
//                    vibes::drawSector(segment_out[0].ub(), segment_out[1].ub(), size*0.3, size*0.3, (-theta_u_bwd.lb())*180.0/M_PI, (-theta_u_bwd.ub())*180.0/M_PI, "b[]");
//                }
//            }

        }

        for(int i=0; i<2; i++){
            for(int k=0; k<m_theta_list.size(); k++)
                vibes::drawSector(m_position[0].mid(), m_position[1].mid(), size, size, (-m_theta_list[k][i].lb())*180.0/M_PI, (-m_theta_list[k][i].ub())*180.0/M_PI, "r[]");
        }


    }
}
void Pave::draw_borders(bool filled, string color_polygon){
    if(!filled){
        // Draw Segments
        for(int i=0; i<m_borders.size(); i++){
            m_borders[i]->draw();
        }
    }
    else{
        // Draw Polygone
        vector<double> x, y;
        for(int i=0; i<m_borders.size(); i++){
            m_borders[i]->get_points(x, y);
        }
        vibes::drawPolygon(x, y, color_polygon);
    }
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result){
    // Create 4 new paves
    ibex::LargestFirst bisector(0.0, 0.5);

    std::pair<IntervalVector, IntervalVector> result_boxes = bisector.bisect(m_position);

    Pave *pave1 = new Pave(result_boxes.first, m_f_list, m_u_iv, m_diseable_singeleton, m_active); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, m_f_list, m_u_iv, m_diseable_singeleton, m_active); // Right or Down

    int indice1, indice2;

    if(pave1->m_position[0] == m_position[0]){
        // Case UP/DOWN bisection
        indice1 = 2;
        indice2 = 0;
    }
    else{
        // Case LEFT/RIGHT bisection
        indice1 = 1;
        indice2 = 3;
    }

    // The order of tasks is important !
    // 1) Update pave brothers with pave1 & pave2
    for(int face=0; face<4; face++){
        m_borders[face]->update_brothers_inclusion(pave1->get_border(face), pave2->get_border(face));
    }

    // 2) Copy brothers Pave (this) to pave1 and pave2
    for(int face=0; face<4; face++){
        if(m_borders[face]->get_inclusions().size()!=0){
            if(face!=indice1){
                pave1->get_border(face)->add_inclusions(m_borders[face]->get_inclusions());

            }
            if(face!=indice2){
                pave2->get_border(face)->add_inclusions(m_borders[face]->get_inclusions());
            }
        }

        if(face!=indice1){
            pave1->get_border(face)->set_continuity_in(m_borders[face]->get_continuity_in());
            pave1->get_border(face)->set_continuity_out(m_borders[face]->get_continuity_out());
        }
        if(face!=indice2){
            pave2->get_border(face)->set_continuity_in(m_borders[face]->get_continuity_in());
            pave2->get_border(face)->set_continuity_out(m_borders[face]->get_continuity_out());
        }
    }

    // 3) Add each other to its brother list (pave1 <-> pave2)
    Inclusion *inclusion_to_pave2 = new Inclusion(pave2->get_border(indice2), indice2);
    Inclusion *inclusion_to_pave1 = new Inclusion(pave1->get_border(indice1), indice1);

    pave1->get_border(indice1)->add_inclusion(inclusion_to_pave2);
    pave2->get_border(indice2)->add_inclusion(inclusion_to_pave1);

    if(is_full()){
        pave1->set_full();
        pave2->set_full();
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

// ********************************************************************************
// ****************** UTILS functions *********************************************

double Pave::get_theta_diam(){
    double diam = 0.0;
    for(int i=0; i<m_theta_list[m_active_function].size(); i++){
        if(!m_theta_list[m_active_function][i].is_empty())
            diam += m_theta_list[m_active_function][i].diam();
    }
    return diam;
}

void Pave::remove_brothers(Pave* p, int face){
    for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
        if(m_borders[face]->get_inclusion(i)->get_border()->get_pave() == p){
            m_borders[face]->remove_inclusion(i);
            return;
        }
    }
}

void Pave::remove_from_brothers(){
    for(int face=0; face<4; face++){
        for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
            m_borders[face]->get_inclusion(i)->get_border()->get_pave()->remove_brothers(this, m_borders[face]->get_inclusion(i)->get_brother_face());
        }
    }
}

bool Pave::is_empty(){
    if(m_empty){
        return true;
    }
    else{
        for(int i=0; i<4; i++){
            if(!m_borders[i]->is_empty()){
                return false;
            }
        }

        m_empty = true;
        return true;
    }
}

bool Pave::is_full(){
    if(!m_full){
        return false;
    }
    else{
        for(int face=0; face<4; face++){
            if(!m_borders[face]->is_full()){
                m_full = false;
                return false;
            }
        }
        m_full = true;
        return true;
    }
}

const std::vector<Pave *> Pave::get_brothers(int face){
    vector<Pave*> brothers_list;
    for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
        brothers_list.push_back(m_borders[face]->get_inclusion(i)->get_border()->get_pave());
    }
    return brothers_list;
}

void Pave::reset_full_empty(){
    m_empty = false;
    m_full = true;
    for(auto &border: m_borders){
        border->reset_full_empty();
    }
}

const Interval &Pave::get_theta(int i) const{
    if(i==0)
        return m_theta_list[m_active_function][0];
    else if(i==1)
        return m_theta_list[m_active_function][1];
    else
        return NULL;
}

const Interval &Pave::get_u(int i) const{
    if(i==0)
        return m_u[0];
    else if(i==1)
        return m_u[1];
    else
        return NULL;
}

const std::vector<Interval> Pave::get_u() const{
    return m_u;
}

const IntervalVector &Pave::get_u_iv() const{
    return m_u_iv;
}

const IntervalVector &Pave::get_position() const{
    return m_position;
}

const std::vector<Border*> &Pave::get_borders(){
    return m_borders;
}

Border* Pave::get_border(int face){
    assert(face >=0 && face < 4);
    return m_borders[face];
}

Border* Pave::operator[](int face){
    return m_borders[face];
}

const Border* Pave::get_border_const(int face) const{
    if(face >=0 && face < 4)
        return m_borders[face];
    else
        return NULL;
}

bool Pave::is_in_queue() const{
    return m_in_queue;
}

void Pave::set_in_queue(bool flag){
    m_in_queue = flag;
}

ibex::Function* Pave::get_f() const{
    return m_f_list[m_active_function];
}

void Pave::set_copy_node(Pave *p){
    m_copy_node = p;
}

Pave* Pave::get_copy_node(){
    return m_copy_node;
}

const std::vector<Interval> Pave::get_theta() const{
    return m_theta_list[m_active_function];
}

void Pave::print(){
    cout << "********" << endl;
    cout << "PAVE x=" << m_position[0] << " y= " << m_position[1] << endl;
    cout << this << endl;
    cout << "theta[0]=" << m_theta_list[m_active_function][0] << " theta[1]=" << m_theta_list[m_active_function][1] << " u=" << m_u_iv << endl;
    for(int face = 0; face < 4; face++){
        if(m_borders[face]->get_inclusions().size()==0){
            cout << "border=" << face << " " << &(m_borders[face])
                 << " segment_in=" << m_borders[face]->get_segment_in()
                 << " segment_out=" << m_borders[face]->get_segment_out()
                 << endl;
        }
        else{
            for(int i=0; i<m_borders[face]->get_inclusions().size(); i++){
                cout << "border=" << face << " " << &(m_borders[face])
                     << " segment_in=" << m_borders[face]->get_segment_in()
                     << " segment_out=" << m_borders[face]->get_segment_out()
                     << " inclusion=" << i
                     << " *border=" << m_borders[face]->get_inclusion(i)->get_border()
                     << " segment_full=" << m_borders[face]->get_inclusion(i)->get_border()->get_segment_full()
                     << endl;
            }
        }
    }
}

bool Pave::get_first_process() const{
    return m_first_process;
}

void Pave::set_first_process_true(){
    m_first_process = true;
}

void Pave::set_first_process_false(){
    m_first_process = false;
}

bool Pave::get_inner() const{
    return m_inner;
}

void Pave::set_inner(bool inner){
    m_inner = inner;
}

std::vector<ibex::Function *> Pave::get_f_list() const{
    return m_f_list;
}

void Pave::set_active_function(int id){
    if(id>=0 && id < m_f_list.size())
        m_active_function = id;
}

int Pave::get_active_function() const{
    return m_active_function;
}

std::vector<std::vector<ibex::Interval>> Pave::get_theta_list() const{
    return m_theta_list;
}

bool Pave::is_active() const{
    return m_active;
}

void Pave::set_continuity_out(bool enable){
    for(int face=0; face<4; face++)
        m_borders[face]->set_continuity_out(enable);
}

void Pave::set_continuity_in(bool enable){
    for(int face=0; face<4; face++)
        m_borders[face]->set_continuity_in(enable);
}

bool Pave::get_diseable_singelton() const{
    return m_diseable_singeleton;
}
