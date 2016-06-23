#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "iomanip"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &position, const std::vector<ibex::Function*> &f_list, bool diseable_singeleton, bool active):
    m_position(2)
{
    m_position = position;    // Box corresponding to the Pave
    m_borders.reserve(4);
    m_f_list = f_list;
    m_active_function = 0;
    m_active = active;
    m_diseable_singeleton = diseable_singeleton;

    // Graph markers
    m_in_queue_inner = false;
    m_in_queue_outer = false;
    m_copy_node = NULL;
    m_first_process = false;
    m_inner_mode = false;
    m_compute_inner = false;

    m_marker_attractor = false;
    m_external_border = false;
    m_removed_pave_inner = false;
    m_removed_pave_outer = false;

    // Border building
    IntervalVector coordinate(2);
    coordinate[0] = position[0]; coordinate[1] = Interval(position[1].lb()); m_borders.push_back(new Border(coordinate, 0, this));
    coordinate[1] = position[1]; coordinate[0] = Interval(position[0].ub()); m_borders.push_back(new Border(coordinate, 1, this));
    coordinate[0] = position[0]; coordinate[1] = Interval(position[1].ub()); m_borders.push_back(new Border(coordinate, 2, this));
    coordinate[1] = position[1]; coordinate[0] = Interval(position[0].lb()); m_borders.push_back(new Border(coordinate, 3, this));

    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;

    m_zone_propagation = false;
    m_backward_function = false;

    /////////////////////////////// THETA ///////////////////////////////
    for(ibex::Function* f:f_list){
        std::vector<ibex::Interval> theta = compute_theta(f, false);
        m_theta_list.push_back(theta);
        for(Interval t:theta){
            if(!t.is_empty())
                m_theta.push_back(t);
        }

        std::vector<ibex::Interval> theta_bwd = compute_theta(f, true);
        m_theta_list_bwd.push_back(theta_bwd);
        for(Interval t:theta_bwd){
            if(!t.is_empty())
                m_theta_bwd.push_back(t);
        }
    }

    if(m_theta.size()==0)
        m_theta.push_back(Interval::EMPTY_SET);
    if(m_theta_bwd.size()==0)
        m_theta_bwd.push_back(Interval::EMPTY_SET);
}

const std::vector<ibex::Interval> Pave::compute_theta(ibex::Function *f, bool backward_function){
    std::vector<ibex::Interval> theta_list;

    for(int i=0; i<2; i++)
        theta_list.push_back(Interval::EMPTY_SET);
    if(f!=NULL){
        IntervalVector dposition = f->eval_vector(m_position);

        Interval dx = dposition[0];
        Interval dy = dposition[1];
        if(backward_function){
            dx = -dposition[0];
            dy = -dposition[1];
        }

        Interval theta = atan2(dy, dx);

        if(theta==(-Interval::PI|Interval::PI)){
            Interval thetaR = atan2(-dy, -dx); // PI rotation ({dx, dy} -> {-dx, -dy})
            if(thetaR.diam()<theta.diam()){
                if(thetaR.is_subset(-Interval::PI | Interval::PI)){
                    theta_list[0] = (thetaR + Interval::PI) & (Interval::ZERO | Interval::PI); // theta[0] in [0, pi]
                    theta_list[1] = ((thetaR + Interval::PI) & (Interval::PI | Interval::TWO_PI)) - Interval::TWO_PI; // theta[1] in [-pi, 0]
                }
                else
                    cout << "****************** ERROR ******************" << endl;
            }
            else
                theta_list[0] = theta;
        }
        else if(theta.is_empty())
            theta_list[0] = -Interval::PI|Interval::PI;
        else
            theta_list[0] = theta;
        if(theta_list[0].is_empty())
            cout << "ERROR - Pave "<< theta << dx << dy << m_position << endl;
    }

    if(theta_list[1].is_empty())
        theta_list.erase(theta_list.begin()+1);

    return theta_list;
}

Pave::Pave(const Pave *p):
    m_position(2)
{
    m_position = p->get_position();    // Box corresponding to the Pave
    m_f_list = p->get_f_list();
    m_active_function = p->get_active_function();
    // Force to recompute results
    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;
    m_in_queue_inner = p->is_in_queue_inner();
    m_in_queue_outer = p->is_in_queue_outer();
    m_first_process = p->get_first_process();
    m_active = p->is_active();
    m_diseable_singeleton = p->get_diseable_singelton();
    m_external_border = p->is_external_border();
    m_removed_pave_inner = p->is_removed_pave_inner();
    m_removed_pave_outer = p->is_removed_pave_outer();
    m_inner_mode = p->get_inner_mode();
    m_compute_inner = p->get_compute_inner();

    m_theta_list = p->get_theta_list_fwd();
    m_theta_list_bwd = p->get_theta_list_bwd();
    m_theta = p->get_all_theta_fwd();
    m_theta_bwd = p->get_all_theta_bwd();
    m_backward_function = p->get_backward_function();

    for(int face = 0; face < 4; face++){
        Border *b = new Border(p->get_border_const(face));
        m_borders.push_back(b); // Copy the border !
        get_border(face)->set_pave(this);
    }
    m_copy_node = NULL;
    m_marker_attractor = p->is_marked_attractor();
}

Pave::~Pave(){
    for(int face=0; face<4; face++){
        delete(m_borders[face]);
    }
}

Pave& Pave::operator&=(const Pave &p){
    for(int face = 0; face <4; face++){
        *(get_border(face)) &= *(p.get_border_const(face));
    }
    return *this;
}

Pave& Pave::operator|=(const Pave &p){
    for(int face = 0; face <4; face++){
        *(get_border(face)) |= *(p.get_border_const(face));
    }
    return *this;
}

bool Pave::inter(const Pave &p, bool with_bwd){
    bool change = false;
    for(int face = 0; face <4; face++){
        if(get_border(face)->inter(*(p.get_border_const(face)), with_bwd))
            change = true;
    }
    return change;
}

bool Pave::inter_inner(const std::vector<Pave*> pave_list){
    for(int face = 0; face <4; face++){
        vector<Border*> border_list;
        for(Pave *p:pave_list)
            border_list.push_back(p->get_border(face));
        get_border(face)->inter_inner(border_list);
    }
}

bool Pave::diff(const Pave &p){
    bool change = false;
    for(int face = 0; face<4; face++){
        if(get_border(face)->diff(*(p.get_border_const(face)))){
            change = true;
        }
    }
    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;
    return change;
}

void Pave::set_theta(ibex::Interval theta){
    std::vector<ibex::Interval> theta_list;
    theta_list.push_back(theta);
    set_theta(theta_list);
}

void Pave::set_theta(std::vector<ibex::Interval> theta_list){
//    m_theta_list.clear();

//    for(Interval &theta:theta_list){
//        std::vector<ibex::Interval> thetas;
//        for(int i=0; i<2; i++)
//            thetas.push_back(Interval::EMPTY_SET);

//        if(theta.is_subset(-Interval::PI | Interval::PI)){
//            thetas[0] = theta;
//        }
//        else{
//            thetas[0] = (theta & (-Interval::PI | Interval::PI));

//            if(!((theta + 2*Interval::PI) & (-Interval::PI | Interval::PI) ).is_empty())
//                thetas[1] =(theta + 2*Interval::PI);
//            else if (!((theta - 2*Interval::PI) & (-Interval::PI | Interval::PI)).is_empty())
//                thetas[1] = (theta - 2*Interval::PI);
//        }

//        m_theta_list.push_back(thetas);
//    }
    cout << "DEPRECATED THETA FUNCTION" << endl;
}

void Pave::set_full(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full();
    }
    if(m_inner_mode)
        m_full_inner = true;
    else
        m_full_outer = true;
}

void Pave::set_full_all(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_all();
    }
    m_full_inner = true;
    m_full_outer = true;
}

void Pave::set_full_in(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_segment_in();
    }
    if(m_inner_mode)
        m_full_inner = true;
    else
        m_full_outer = true;
}

void Pave::set_full_out(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_segment_out();
    }
    if(m_inner_mode)
        m_full_inner = true;
    else
        m_full_outer = true;
}

void Pave::set_empty(){
    for(int face=0; face<4; face++){
        get_border(face)->set_empty();
    }
    if(m_inner_mode){
        m_empty_inner = true;
        m_full_inner = false;
    }
    else{
        m_empty_outer = true;
        m_full_outer = false;
    }
}

void Pave::set_empty_outer(){
    for(int face=0; face<4; face++){
        get_border(face)->set_empty_outer();
    }
    m_empty_outer = true;
    m_full_outer = false;
}

void Pave::set_empty_inner(){
    for(int face=0; face<4; face++){
        get_border(face)->set_empty_inner();
    }
    m_empty_inner = true;
    m_full_inner = false;
}

void Pave::set_full_outer(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_outer();
    }
    m_empty_outer = false;
    m_full_outer = true;
}

void Pave::set_full_inner(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_inner();
    }
    m_empty_inner = false;
    m_full_inner = true;
}

void Pave::set_segment(bool in, bool out){
    for(int face=0; face<4; face++){
        get_border(face)->set_segment(in, out);
    }
    if(m_inner_mode){
        m_empty_inner = false;
        m_full_inner = true;
    }
    else{
        m_empty_outer = false;
        m_full_outer = true;
    }
}
// ********************************************************************************
// ****************** Drawing functions *******************************************

void Pave::draw(bool filled, bool borders_only){
    // Magenta = #FF00FF
    // Gray light =  #D3D3D3
    // Blue = #4C4CFF

    // Draw the pave
    if(borders_only){
        draw_borders(filled, "#00FF00AA[#00FF00AA]");
    }
    else{
        vibes::drawBox(m_position, "#D3D3D3[]");
        if(m_compute_inner){
            if(!m_external_border){
                bool mode = get_inner_mode();
                /// OUTER
                set_inner_mode(false);
                draw_borders(true, "#4C4CFF[#4C4CFF]", true); // blue

                /// INNER
                set_inner_mode(true);
                draw_borders(true, "#FF00FF[#FF00FF]", true); // magenta

                Pave *p_polygon = new Pave(this);
                p_polygon->set_inner_mode(true);
                set_inner_mode(false);
                p_polygon->inter(*this);
                p_polygon->draw_borders(true, "y[y]");
                delete(p_polygon);

                set_inner_mode(mode);
            }
        }
        else{
            if(!m_external_border){
                vibes::drawBox(m_position, "#4C4CFF[#4C4CFF]");
                draw_borders(filled, "y[y]"); // yellow
            }
        }

        // Draw theta
        draw_theta();
    }
}

void Pave::draw_theta() const{
    double size = 0.8*min(m_position[0].diam(), m_position[1].diam())/2.0;

    std::vector<std::string> color_map;
    color_map.push_back("black[gray]");
    color_map.push_back("black[#A8A8A8]");
    color_map.push_back("g[]");

    for(int k=0; k<(int)get_theta_list().size(); k++){
        double size_ratio = size * (1-0.1*k);
        for(Interval i:get_theta_list(k)){
            vibes::drawSector(m_position[0].mid(), m_position[1].mid(), size_ratio, size_ratio, (-i.lb())*180.0/M_PI, (-i.ub())*180.0/M_PI, color_map[k%color_map.size()]);
        }
    }
}

void Pave::draw_borders(bool filled, string color_polygon, bool complementary) const{
    if(!filled){
        // Draw Segments
        for(int face=0; face<(int)get_borders_const().size(); face++){
            //bool same_size, double offset, bool test
            get_border_const(face)->draw(false, -0.01*m_position[face%2].diam(), false, false);
        }
    }
    else{
        // Draw Polygone
        vector<double> x, y;
        for(Border *b:m_borders){
            b->get_points(x, y, complementary);
        }
        if(x.size()>0)
            vibes::drawPolygon(x, y, color_polygon);

        //        int alone_points;
        //        int nb_point = x.size();
        //        if(nb_point>0){
        //            vector<int> starting_point;
        //            for(int i=0; i<(int)x.size(); i++){
        //                if(!((x[i]==x[(i+1)%nb_point] && y[i]==y[(i+1)%nb_point]) || !(x[i]==x[(i-1+nb_point)%nb_point] && y[i]==y[(i-1+nb_point)%nb_point])))
        //                    alone_points++;
        //                else
        //                    starting_point.push_back(i);
        //            }

        //            if(nb_point == 8 && alone_points == 4){
        //                for(int start:starting_point){
        //                    vector<double> xp, yp;
        //                    for(int i=(start-1+nb_point)%nb_point; i<=(start+1+nb_point)%nb_point; i++){
        //                        xp.push_back(x[i]);
        //                        yp.push_back(y[i]);
        //                    }
        //                    vibes::drawPolygon(xp, yp, color_polygon);
        //                }

        //            }
        //            else{
        //                vibes::drawPolygon(x, y, color_polygon);
        //            }
        //        }
    }
}

void Pave::draw_test(int size, string comment) const{
    vibes::beginDrawing();
    stringstream ss;
    ss << "integralIbex - pave=" << m_position << comment;
    vibes::newFigure(ss.str());
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",size,"height",size));

    vibes::drawBox(m_position, "black[]");
    draw_borders(true, "y[y]");
    double offset[2] = {m_position[0].diam()*0.1, m_position[1].diam()*0.1};

    for(int i=0; i<(int)m_borders.size(); i++){
        m_borders[i]->draw(true, offset[i%2]*(((i+3)%4>1)?-1:1), true);
    }
    draw_theta();

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result, bool backward){
    // Create 4 new paves
    ibex::LargestFirst bisector(0.0, 0.5);

    std::pair<IntervalVector, IntervalVector> result_boxes = bisector.bisect(m_position);

    Pave *pave1 = new Pave(result_boxes.first, m_f_list, m_diseable_singeleton, m_active); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, m_f_list, m_diseable_singeleton, m_active); // Right or Down

    pave1->set_active_function(get_active_function());
    pave2->set_active_function(get_active_function());

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

    if(backward){
        pave1->set_full_all();
        pave2->set_full_all();

#if 1
        if(m_borders[(indice1+1)%4]->is_empty() || m_borders[(indice1+3)%4]->is_empty()){
            bool theta_inside = false;
            bool theta_outside = false;
            for(vector<Interval> &theta_list:get_theta_list()){

                for(Interval &theta:theta_list){
                    if(theta.is_empty())
                        break;

                    switch(indice1){
                    case 1:
                        // Case LEFT/RIGHT bisection
                        if(!(theta & -Interval::HALF_PI).is_empty()){
                            Interval theta_centered = theta + Interval::HALF_PI;
//                            if(m_position[1].diam()*(fabs(atan(theta_centered.lb()))+fabs(atan(theta_centered.ub())))<m_position[0].diam()/2.0)
                            if(m_position[1].diam()*(fabs(atan(theta_centered.lb()))) < m_position[0].diam()/2.0
                                    && m_position[1].diam()*(fabs(atan(theta_centered.ub())))<m_position[0].diam()/2.0)
                                theta_inside = true;
                            else{
                                theta_inside = false;
                                break;
                            }
                        }
                        if(!(theta & Interval::HALF_PI).is_empty()){
                            Interval theta_centered = theta - Interval::HALF_PI;
                            if(m_position[1].diam()*(fabs(atan(theta_centered.lb()))) < m_position[0].diam()/2.0
                                    && m_position[1].diam()*(fabs(atan(theta_centered.ub())))<m_position[0].diam()/2.0)
                                theta_inside = true;
                            else
                                theta_outside = true;
                        }
                        break;
                    case 2:
                        // Case UP/DOWN bisection
                        if(!(theta & Interval::ZERO).is_empty()){
                            Interval theta_centered = theta - Interval::ZERO;
                            if(m_position[0].diam()*(fabs(atan(theta_centered.lb()))) < m_position[1].diam()/2.0
                                    && m_position[0].diam()*(fabs(atan(theta_centered.ub())))<m_position[1].diam()/2.0)
                                theta_inside = true;
                            else
                                theta_outside = true;
                        }
                        if(!(theta & Interval::PI).is_empty()){
                            Interval theta_centered = theta - Interval::PI;
                            if(m_position[0].diam()*(fabs(atan(theta_centered.lb()))) < m_position[1].diam()/2.0
                                    && m_position[0].diam()*(fabs(atan(theta_centered.ub())))<m_position[1].diam()/2.0)
                                theta_inside = true;
                            else
                                theta_outside = true;
                        }
                        break;
                    }
                }
            }
            if(theta_inside && !theta_outside){
                pave1->get_border(indice1)->set_empty();
                pave2->get_border(indice2)->set_empty();

                if(m_borders[(indice1+1)%4]->is_empty()){
                    pave1->get_border((indice1+1)%4)->set_empty();
                    pave2->get_border((indice2+1)%4)->set_empty();

                }
                if(m_borders[(indice1+3)%4]->is_empty()){
                    pave1->get_border((indice1+3)%4)->set_empty();
                    pave2->get_border((indice2+3)%4)->set_empty();
                }
            }

            if(!is_border()){
                for(int face=0; face<4; face++){
                    if(face!=indice1){
                        if((get_border(face)->get_segment_in() & pave1->get_border(face)->get_segment_in()).is_empty()
                                && (get_border(face)->get_segment_out() & pave1->get_border(face)->get_segment_out()).is_empty())
                            pave1->get_border(face)->set_empty();

                    }
                    if(face!=indice2){
                        if((get_border(face)->get_segment_in() & pave2->get_border(face)->get_segment_in()).is_empty()
                                && (get_border(face)->get_segment_out() & pave2->get_border(face)->get_segment_out()).is_empty())
                            pave2->get_border(face)->set_empty();
                    }
                }
            }
        }
#endif
    }

    result.push_back(pave1);
    result.push_back(pave2);
}

// ********************************************************************************
// ****************** UTILS functions *********************************************

double Pave::get_theta_diam(int active_function) const{
    if((active_function ==-1)){
        double diam_max = 0.0;
        for(int active_function = 0; active_function<(int)m_f_list.size(); active_function++){
            double diam = 0.0;
            for(Interval theta:get_theta_list(active_function)){
                if(!theta.is_empty())
                    diam += theta.diam();
            }
            diam_max=max(diam_max, diam);
        }
        return diam_max;
    }
    else{
        double diam = 0.0;
        for(Interval theta:get_theta_list(active_function)){
            if(!theta.is_empty())
                diam += theta.diam();
        }
        return diam;
    }
}

double Pave::get_theta_diam_max() const{
    double diam_max = 0.0;
    for(int active_function = 0; active_function<(int)m_f_list.size(); active_function++){
        double diam = 0.0;
        for(Interval theta:get_theta_list(active_function)){
            if(!theta.is_empty())
                diam += theta.diam();
        }
        diam_max=max(diam_max, diam);
    }
    return diam_max;
}

double Pave::get_theta_diam_min(){
    double diam_min = 2*M_PI;
    for(int active_function = 0; active_function<(int)m_f_list.size(); active_function++){
        double diam = 0.0;
        for(Interval theta:get_theta_list()[active_function]){
            if(!theta.is_empty())
                diam += theta.diam();
        }
        diam_min=min(diam_min, diam);
    }
    return diam_min;
}

void Pave::remove_brothers(Pave* p, int face){
    for(int i=0; i<(int)m_borders[face]->get_inclusions().size(); i++){
        if(m_borders[face]->get_inclusion(i)->get_border()->get_pave() == p){
            m_borders[face]->remove_inclusion(i);
            return;
        }
    }
}

void Pave::remove_from_brothers(){
    for(int face=0; face<4; face++){
        for(int i=0; i<(int)m_borders[face]->get_inclusions().size(); i++){
            m_borders[face]->get_inclusion(i)->get_border()->get_pave()->remove_brothers(this, m_borders[face]->get_inclusion(i)->get_brother_face());
        }
    }
}

bool Pave::is_empty(){
    if(!m_compute_inner)
        return is_empty_outer();
    else{
        if(m_inner_mode)
            return is_empty_inner();
        else
            return is_empty_outer();
    }
}

bool Pave::is_empty_inter(){
    if(!m_compute_inner)
        return is_empty_outer();
    else{
        if(is_empty_inner() && is_empty_outer())
            return true;
        else
            return false;
    }
}

bool Pave::is_empty_inner(){
    if(m_empty_inner)
        return true;
    else{
        for(Border *b:m_borders){
            if(!b->is_empty_inner())
                return false;
        }
        m_empty_inner = true;
        return true;
    }
}

bool Pave::is_empty_outer(){
    if(m_empty_outer)
        return true;
    else{
        for(Border *b:m_borders){
            if(!b->is_empty_outer())
                return false;
        }
        m_empty_outer = true;
        return true;
    }
}

bool Pave::is_full_inner(){
    if(!m_full_inner)
        return false;
    else{
        for(Border *b:m_borders){
            if(!b->is_full_inner()){
                m_full_inner = false;
                return false;
            }
        }
    }
    //m_full_inner = true;
    return true;
}

bool Pave::is_full_outer(){
    if(!m_full_outer)
        return false;
    else{
        for(Border *b:m_borders){
            if(!b->is_full_outer()){
                m_full_outer = false;
                return false;
            }
        }
    }
    m_full_outer = true;
    return true;
}

bool Pave::is_full(){
    if(!m_compute_inner)
        return is_full_outer();
    else{
        if(m_inner_mode)
            return is_full_inner();
        else
            return is_full_outer();
    }
}

bool Pave::is_full_inter(){
    if(!m_compute_inner)
        return is_full_outer();
    else{
        if(is_full_inner() && is_full_outer())
            return true;
        else
            return false;
    }
}

bool Pave::is_full_geometricaly() const{
    IntervalVector box(2, Interval::EMPTY_SET);
    int non_empty_face = 0;
    for(Border *b:m_borders){
        box |= b->get_segment_out_2D();
        box |= b->get_segment_in_2D();
        if(!b->is_empty())
            non_empty_face++;
    }
    if(box == m_position && non_empty_face>2)
        return true;
    else
        return false;
}

const std::vector<Pave *> Pave::get_brothers(int face){
    vector<Pave*> brothers_list;
    for(int i=0; i<(int)m_borders[face]->get_inclusions().size(); i++){
        if(!m_borders[face]->get_inclusion(i)->get_border()->get_pave()->is_removed_pave())
            brothers_list.push_back(m_borders[face]->get_inclusion(i)->get_border()->get_pave());
    }
    return brothers_list;
}

const std::vector<Pave *> Pave::get_all_brothers(){
    vector<Pave*> brothers_list;
    for(int face = 0; face<4; face++){
        for(Inclusion *i:m_borders[face]->get_inclusions()){
            brothers_list.push_back(i->get_border()->get_pave());
        }
    }
    return brothers_list;
}

void Pave::reset_full_empty(){
    m_empty_inner = false;
    m_empty_outer = false;
    m_full_inner = true;
    m_full_outer = true;
    for(Border *border: m_borders){
        border->reset_full_empty();
    }
}

const Interval &Pave::get_theta(int i) const{
    if(i==0){
        if(m_backward_function)
            return m_theta_list_bwd[m_active_function][0];
        else
            return m_theta_list[m_active_function][0];
    }
    else if(i==1){
        if(m_backward_function)
            return m_theta_list_bwd[m_active_function][1];
        else
            return m_theta_list[m_active_function][1];
    }
    else
        return NULL;
}

const IntervalVector &Pave::get_position() const{
    return m_position;
}

std::vector<Border*> &Pave::get_borders(){
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
    if(!m_compute_inner)
        return is_in_queue_outer();
    else{
        if(m_inner_mode)
            return is_in_queue_inner();
        else
            return is_in_queue_outer();
    }
}

void Pave::set_in_queue(bool flag){
    if(!m_compute_inner)
        set_in_queue_outer(flag);
    else{
        if(m_inner_mode)
            set_in_queue_inner(flag);
        else
            set_in_queue_outer(flag);
    }
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
    if(m_backward_function)
        return m_theta_list_bwd[m_active_function];
    else
        return m_theta_list[m_active_function];
}

const std::vector<Interval> Pave::get_all_theta(bool all) const{
    if(all || m_active_function==-1){
        if(m_backward_function)
            return m_theta_bwd;
        else
            return m_theta;
    }
    else
        return get_theta();
}

std::vector<Interval> Pave::get_all_theta_fwd() const{
    return m_theta;
}
std::vector<Interval> Pave::get_all_theta_bwd() const{
    return m_theta_bwd;
}

void Pave::print(){
    cout << "********" << endl;
    cout << "PAVE x=" << get_position()[0] << " y= " << m_position[1] << endl;
    cout << this << endl;
    cout << "theta[0]=" << get_theta(0) << " theta[1]=" << get_theta(1) << endl;
    for(int face = 0; face < 4; face++){
        if(get_border(face)->get_inclusions().size()==0){
            cout << "border=" << face << " " << (get_border(face))
                 << " segment_in=" << get_border(face)->get_segment_in()
                 << " segment_out=" << get_border(face)->get_segment_out()
                 << endl;
        }
        else{
            for(int i=0; i<get_border(face)->get_inclusions().size(); i++){
                cout << "border=" << face << " " << (get_border(face))
                     << " segment_in=" << get_border(face)->get_segment_in()
                     << " segment_out=" << get_border(face)->get_segment_out()
                     << " inclusion=" << i
                     << " *border=" << get_border(face)->get_inclusion(i)->get_border()
                     << " segment_full=" << get_border(face)->get_inclusion(i)->get_border()->get_segment_full()
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

std::vector<ibex::Function *> Pave::get_f_list() const{
    return m_f_list;
}

void Pave::set_active_function(int id){
    m_active_function = id;
}

int Pave::get_active_function() const{
    return m_active_function;
}

std::vector<std::vector<ibex::Interval>> Pave::get_theta_list() const{
    if(m_backward_function)
        return m_theta_list_bwd;
    else
        return m_theta_list;
}

std::vector<std::vector<ibex::Interval>> Pave::get_theta_list_fwd() const{
    return m_theta_list;
}

std::vector<ibex::Interval> Pave::get_theta_list(int function_id) const{
    if(m_backward_function)
        return m_theta_list_bwd[function_id];
    else
        return m_theta_list[function_id];
}

std::vector<std::vector<ibex::Interval>> Pave::get_theta_list_bwd() const{
    return m_theta_list_bwd;
}

bool Pave::is_active() const{
    return m_active;
}

void Pave::set_continuity_out(bool enable){
    for(Border *b:get_borders())
        b->set_continuity_out(enable);
}

void Pave::set_continuity_in(bool enable){
    for(Border *b:get_borders())
        b->set_continuity_in(enable);
}

bool Pave::get_diseable_singelton() const{
    return m_diseable_singeleton;
}

bool Pave::is_border() const{
    for(Border *b:get_borders_const()){
        if(b->get_inclusions().size() == 0){
            return true;
        }
    }
    for(Border *b:get_borders_const()){

        Interval segment_border = Interval::EMPTY_SET;
        for(Inclusion *i:b->get_inclusions()){
            segment_border |= i->get_segment_full();

            if(!i->get_border()->get_pave()->is_active() || i->get_border()->get_pave()->is_external_border()){
                return true;
            }
        }
        if(segment_border != b->get_segment_full()){
            return true;
        }
    }
    return false;
}

void Pave::complementaire(){
    for(Border *b:get_borders()){
        b->complementaire();
    }
}

IntervalVector Pave::get_bounding_pave() const{
    IntervalVector box(2, Interval::EMPTY_SET);

    for(Border *b:get_borders_const()){
        box |= b->get_segment_in_2D();
        box |= b->get_segment_out_2D();
    }
    return box;
}

bool Pave::is_marked_attractor() const{
    return m_marker_attractor;
}

void Pave::set_marker_attractor(bool val){
    m_marker_attractor = val;
}

void Pave::set_active(bool val){
    m_active = val;
}

bool Pave::is_external_border() const{
    return m_external_border;
}

void Pave::set_external_border(bool val){
    m_external_border = val;
}

bool Pave::is_removed_pave() const{
    if(m_compute_inner){
        if(m_inner_mode)
            return is_removed_pave_inner();
        else
            return is_removed_pave_outer();
    }
    else
        return is_removed_pave_outer();
}

bool Pave::is_removed_pave_union() const{
    if(m_compute_inner){
        return is_removed_pave_outer() || is_removed_pave_inner() ;
    }
    else
        return is_removed_pave_outer();
}

void Pave::set_removed_pave(bool val){
    if(m_compute_inner){
        if(m_inner_mode)
            m_removed_pave_inner = val;
        else
            m_removed_pave_outer = val;
    }
    else
        m_removed_pave_outer = val;
}

void Pave::set_removed_pave_inner(bool val){
    m_removed_pave_inner = val;
}

void Pave::set_removed_pave_outer(bool val){
    m_removed_pave_outer = val;
}

bool Pave::is_near_empty(){
    vector<Pave*> brothers = get_all_brothers();
    for(Pave *p:brothers){
        if(!p->is_external_border() && (p->is_empty() || p->is_removed_pave())){
            return true;
        }
    }
    return false;
}

void Pave::print_theta_list(){
    int count = 0;
    cout << "theta ";
    for(vector<Interval> &theta_list:get_theta_list()){
        cout << "(" << count << ") ";
        for(Interval &theta:theta_list){
            cout << theta << " ";
        }
        count++;
    }
    cout << endl;
}

const std::vector<Border *> Pave::get_borders_const() const{
    return m_borders;
}

void Pave::set_compute_inner(bool val){
    m_compute_inner = val;
    for(Border *b:m_borders)
        b->set_compute_inner(val);
}

void Pave::set_inner_mode(bool val){
    if(!m_compute_inner)
        set_compute_inner(true);
    m_inner_mode = val;
    for(Border *b:m_borders){
        b->set_inner_mode(val);
    }
}

bool Pave::get_inner_mode() const{
    return m_inner_mode;
}

bool Pave::get_compute_inner() const{
    return m_compute_inner;
}

bool Pave::is_removed_pave_outer() const{
    return m_removed_pave_outer;
}

bool Pave::is_removed_pave_inner() const{
    return m_removed_pave_inner;
}

void Pave::set_in_queue_inner(bool flag){
    m_in_queue_inner = flag;
}

void Pave::set_in_queue_outer(bool flag){
    m_in_queue_outer = flag;
}

bool Pave::is_in_queue_inner() const{
    return m_in_queue_inner;
}

bool Pave::is_in_queue_outer() const{
    return m_in_queue_outer;
}

void Pave::copy_to_inner(){
    for(Border *b:m_borders)
        b->copy_to_inner();
}

bool Pave::get_zone_propagation() const{
    return m_zone_propagation;
}

void Pave::reset_computation_zone(){
    m_zone_propagation = false;
    for(Border *b:m_borders)
        b->reset_computation_zone();
}

void Pave::set_zone_propagation(bool val){
    m_zone_propagation = val;
    for(Border *b:m_borders)
        b->set_zone_propagation(val);
}

void Pave::set_backward_function(bool val){
    m_backward_function = val;
}

bool Pave::get_backward_function() const{
    return m_backward_function;
}
