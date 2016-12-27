#include "pave.h"
#include "vibes.h"
#include "border.h"

#include "iostream"
#include "iomanip"

using namespace std;
using namespace ibex;

Pave::Pave(const IntervalVector &position, const std::vector<ibex::Function*> &f_list, bool diseable_singeleton, bool active):
    m_position(2), m_search_box(2), m_vector_field(2)
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
    m_first_process_inner = false;
    m_first_process_outer = false;
    m_inner_mode = false;
    m_compute_inner = false;

    m_marker_attractor = false;
    m_marker = false;
    m_external_border = false;
    m_removed_pave_inner = false;
    m_removed_pave_outer = false;

    m_bassin = false;
    m_infinity_pave = false;

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

    m_vector_field = compute_vector_field(f_list);

    if(m_theta.size()==0)
        m_theta.push_back(Interval::EMPTY_SET);
    if(m_theta_bwd.size()==0)
        m_theta_bwd.push_back(Interval::EMPTY_SET);

    m_theta_union_list = compute_theta_union(f_list, false);
    m_theta_union_list_bwd = compute_theta_union(f_list, true);

    m_cpt_continuity_inner = 0;
    m_cpt_continuity_outer = 0;
    m_cpt_consistency_inner = 0;
    m_cpt_consistency_outer= 0;

    if(get_theta_diam_min()<2*M_PI)
        m_theta_more_than_two_pi = false;
    else
        m_theta_more_than_two_pi = true;
}

const std::vector<ibex::Interval> Pave::compute_theta(ibex::Function *f, bool backward_function){
    if(f!=NULL){
        IntervalVector dposition = f->eval_vector(m_position);

        Interval dx = dposition[0];
        Interval dy = dposition[1];
        if(backward_function){
            dx = -dposition[0];
            dy = -dposition[1];
        }
        return compute_theta(dx, dy);
    }
    else{
        cout << "ERROR : f is NULL" << endl;
        exit(1);
    }
}

const ibex::IntervalVector Pave::compute_vector_field(std::vector<ibex::Function*> f_list){
    IntervalVector dposition(2, Interval::EMPTY_SET);
    for(Function *f:f_list)
        dposition = dposition | f->eval_vector(m_position);
    return dposition;
}

const std::vector<ibex::Interval> Pave::compute_theta_union(std::vector<ibex::Function *> f_list, bool backward_function){
    if(f_list.size()!=0){
        Interval dx(Interval::EMPTY_SET), dy(Interval::EMPTY_SET);
        for(Function *f:f_list){
            IntervalVector dposition = f->eval_vector(m_position);
            dx |= dposition[0];
            dy |= dposition[1];
        }

        if(backward_function){
            dx = -dx;
            dy = -dy;
        }
        return compute_theta(dx, dy);
    }
    else{
        cout << "ERROR : f is NULL" << endl;
        exit(1);
    }
}

const std::vector<ibex::Interval> Pave::compute_theta(ibex::Interval dx, ibex::Interval dy) {
    std::vector<ibex::Interval> theta_list;

    for(int i=0; i<2; i++)
        theta_list.push_back(Interval::EMPTY_SET);

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


    if(theta_list[1].is_empty())
        theta_list.erase(theta_list.begin()+1);

    return theta_list;
}

Pave::Pave(const Pave *p):
    m_position(2), m_search_box(2), m_vector_field(2)
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
    m_first_process_inner = p->get_first_process_inner();
    m_first_process_outer = p->get_first_process_outer();
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
    m_theta_union_list = p->get_theta_union_list_fwd();
    m_theta_union_list_bwd = p->get_theta_union_list_bwd();
    m_backward_function = p->get_backward_function();

    m_vector_field = p->get_vector_field();

    for(int face = 0; face < 4; face++){
        Border *b = new Border(p->get_border_const(face));
        m_borders.push_back(b); // Copy the border !
        get_border(face)->set_pave(this);
    }
    m_copy_node = NULL;
    m_marker_attractor = p->is_marked_attractor();
    m_marker = p->is_marked();

    m_zone_propagation = p->get_zone_propagation();

    m_cpt_consistency_inner = p->get_cpt_consistency_inner();
    m_cpt_consistency_outer = p->get_cpt_consistency_outer();
    m_cpt_continuity_inner = p->get_cpt_continuity_inner();
    m_cpt_continuity_outer= p->get_cpt_continuity_outer();

    m_bassin = p->is_bassin();
    m_infinity_pave = p->is_infinity_pave();
    m_theta_more_than_two_pi = p->is_theta_more_than_two_pi();

    m_search_box = p->get_search_box();
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
    if(p.is_bassin())
        set_bassin(true);
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
    for(Border *b:m_borders){
        b->set_empty_outer();
    }
    m_empty_outer = true;
    m_full_outer = false;
}

void Pave::set_empty_inner(){
    for(int face=0; face<4; face++)
        get_border(face)->set_empty_inner();
    m_empty_inner = true;
    m_full_inner = false;
}

void Pave::set_empty_inner_in(){
    for(Border *b:m_borders)
        b->set_empty_inner_in();
}

void Pave::set_empty_inner_out(){
    for(Border *b:m_borders)
        b->set_empty_inner_out();
}

void Pave::set_empty_outer_in(){
    for(Border *b:m_borders)
        b->set_empty_outer_in();
}

void Pave::set_empty_outer_out(){
    for(Border *b:m_borders)
        b->set_empty_outer_out();
}

void Pave::set_full_outer_in(){
    for(Border *b:m_borders)
        b->set_full_outer_in();
}

void Pave::set_full_outer_out(){
    for(Border *b:m_borders)
        b->set_full_outer_out();
}

void Pave::set_full_inner_in(){
    for(Border *b:m_borders)
        b->set_full_inner_in();
}

void Pave::set_full_inner_out(){
    for(Border *b:m_borders)
        b->set_full_inner_out();
}

void Pave::set_full_outer(){
    for(int face=0; face<4; face++){
        get_border(face)->set_full_outer();
    }
    m_empty_outer = false;
    m_full_outer = true;
}

void Pave::set_full_inner(){
    for(Border *b:m_borders){
        b->set_full_inner();
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

void Pave::draw(bool filled, bool inner_only){
    // Magenta = #FF00FF
    // Gray light =  #D3D3D3
    // Blue = #4C4CFF

    // Draw the pave
    //    if(borders_only){
    //        draw_borders(filled, "#00FF00AA[#00FF00AA]");
    //    }
    //    else{
    IntervalVector position(m_position);
    if(!m_position.is_unbounded()){
        vibes::drawBox(m_position, "#D3D3D3[]");
    }
    else{
        IntervalVector box_draw(m_position);
        if(box_draw[0].lb() == NEG_INFINITY)
            box_draw[0] = Interval(m_search_box[0].lb()-0.1*m_search_box[0].diam(), box_draw[0].ub());
        if(box_draw[1].lb() == NEG_INFINITY)
            box_draw[1] = Interval(m_search_box[1].lb()-0.1*m_search_box[1].diam(), box_draw[1].ub());
        if(box_draw[0].ub() == POS_INFINITY)
            box_draw[0] = Interval(box_draw[0].lb(), m_search_box[0].ub()+0.1*m_search_box[0].diam());
        if(box_draw[1].ub() == POS_INFINITY)
            box_draw[1] = Interval(box_draw[1].lb(), m_search_box[1].ub()+0.1*m_search_box[1].diam());
        vibes::drawBox(box_draw, "#000000[]", vibes::Params("LineStyle", "--")); //D3D3D3
        position = box_draw;
    }
    if(m_compute_inner){
        if(!m_external_border){
            bool mode = get_inner_mode();
            bool backward = get_backward_function();
            /// OUTER
            if(!inner_only){
                set_inner_mode(false);
                draw_borders(true, "#4C4CFF[#4C4CFF]", true); // blue
            }

            /// INNER
            set_inner_mode(true);
            //            if(!is_bassin())
            draw_borders(true, "#FF00FF[#FF00FF]", true); // magenta
            //            else
            //                draw_borders(true, "#FF0000[#FF0000]", true); // red (inside bassin)

            /// POLYGON
            if(!m_position.is_unbounded()){
                Pave *p_polygon = new Pave(this);
                p_polygon->set_inner_mode(true);
                if(!inner_only){
                    set_inner_mode(false);
                    p_polygon->inter(*this, false);
                }
                p_polygon->draw_borders(true, "y[y]");
                delete(p_polygon);
            }
            else{
                if(!is_empty())
                    vibes::drawBox(position, "y[y]");
            }

            set_inner_mode(mode);
            set_backward_function(backward);
        }
        else if(is_bassin()){
            vibes::drawBox(m_position, "#FF0000[#FF0000]");
        }
        else if(m_infinity_pave){
            if(!is_empty_outer())
                vibes::drawBox(position, "orange[orange]");
        }
    }
    else{
        if(!m_external_border){
            vibes::drawBox(m_position, "#4C4CFF[#4C4CFF]");
            draw_borders(filled, "y[y]"); // yellow
        }
        else{
            reset_full_empty();
            if(!is_empty())
                vibes::drawBox(position, "y[y]");
        }
    }

    // Draw theta
    //    set_backward_function(false);
    draw_theta(position);
    //    }

    //    if(m_segment_list.size()!=0){
    //        for(vector<IntervalVector> &segment:m_segment_list){
    //            vector<double> x, y;
    //            for(IntervalVector &pt:segment){
    //                x.push_back(pt[0].mid());
    //                y.push_back(pt[1].mid());
    //            }
    //            vibes::drawLine(x, y, "g[g]",vibesParams("LineWidth",10.0));
    //        }
    //    }
}

void Pave::draw_theta(IntervalVector position) const{
    double size = 0.8*min(position[0].diam(), position[1].diam())/2.0;

    std::vector<std::string> color_map;
    color_map.push_back("black[gray]");
    color_map.push_back("black[#A8A8A8]");
    color_map.push_back("g[]");

    for(int k=0; k<(int)get_theta_list().size(); k++){
        double size_ratio = size * (1-0.1*k);
        for(Interval i:get_theta_list(k)){
            vibes::drawSector(position[0].mid(), position[1].mid(), size_ratio, size_ratio, (-i.lb())*180.0/M_PI, (-i.ub())*180.0/M_PI, color_map[k%color_map.size()]);
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

        if(!is_removed_pave() && is_theta_more_than_two_pi()){
            vibes::drawBox(get_position(), color_polygon);
        }
        else{
            vector<double> x, y;
            for(Border *b:m_borders){
                b->get_points(x, y, complementary);
            }
            if(x.size()>0)
                vibes::drawPolygon(x, y, color_polygon);
        }
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
    draw_theta(m_position);

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}

// ********************************************************************************
// ****************** Paving building *********************************************

void Pave::bisect(vector<Pave*> &result, bool backward, bool apply_heuristic){
    // Create 4 new paves
    ibex::LargestFirst bisector(0.0, 0.5);
    std::pair<IntervalVector, IntervalVector> result_boxes(IntervalVector(2), IntervalVector(2));

    if(apply_heuristic){
//        std::pair<IntervalVector, IntervalVector> bisection_0 = m_position.bisect(0);
//        std::pair<IntervalVector, IntervalVector> bisection_1 = m_position.bisect(1);

//        IntervalVector vectField_b0_first = get_f_list()[0]->eval_vector(bisection_0.first);
//        IntervalVector vectField_b0_second = get_f_list()[0]->eval_vector(bisection_0.second);
//        IntervalVector vectField_b1_first  = get_f_list()[0]->eval_vector(bisection_1.first);
//        IntervalVector vectField_b1_second = get_f_list()[0]->eval_vector(bisection_1.second);

//        double b0_first_max = vectField_b0_first.diam().max();
//        double b0_second_max = vectField_b0_second.diam().max();
//        double b1_first_max = vectField_b1_first.diam().max();
//        double b1_second_max = vectField_b1_second.diam().max();

//        if(!isinf(b0_first_max) && !isinf(b1_first_max)){
//            if(b0_first_max > b1_first_max){
//                result_boxes = bisection_1;
//            }
//            else{
//                result_boxes = bisection_0;
//            }
//        }
//        else{
            result_boxes = bisector.bisect(m_position);
//        }
        //result_boxes = m_position.bisect(dim_bisect);
        // To Do test smaller angle...
    }
    else{
        result_boxes = bisector.bisect(m_position);
    }

    Pave *pave1 = new Pave(result_boxes.first, m_f_list, m_diseable_singeleton, m_active); // Left or Up
    Pave *pave2 = new Pave(result_boxes.second, m_f_list, m_diseable_singeleton, m_active); // Right or Down

    pave1->set_active_function(get_active_function());
    pave2->set_active_function(get_active_function());

    pave1->set_backward_function(get_backward_function());
    pave2->set_backward_function(get_backward_function());

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

#if 0
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

    if(apply_heuristic && /*pave1->get_theta_diam_max()<M_PI &&*/ pave1->get_theta_diam_max()>M_PI/10)
        result.insert(result.begin(), pave1);
    else
        result.push_back(pave1);

    if(apply_heuristic && /*pave2->get_theta_diam_max()<M_PI &&*/ pave2->get_theta_diam_max()>M_PI/10)
        result.insert(result.begin(), pave2);
    else
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
        for(Interval theta:get_theta_list(active_function)){
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
    for(Inclusion *i:get_border(face)->get_inclusions()){
        //        if(i->get_border()->get_pave()->is_active())
        if(!i->get_border()->get_pave()->is_removed_pave_outer())
            brothers_list.push_back(i->get_border()->get_pave());
    }
    return brothers_list;
}

const std::vector<Pave *> Pave::get_all_brothers(){
    vector<Pave*> brothers_list;
    for(Border *b:get_borders()){
        for(Inclusion *i:b->get_inclusions()){
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
    if(i==0 || i==1){
        if(m_backward_function)
            return m_theta_list_bwd[m_active_function][i];
        else
            return m_theta_list[m_active_function][i];
    }
    else{
        cout << "ERROR get_theta != {0,1}" << endl;
        exit(1);
    }
}

const IntervalVector &Pave::get_position() const{
    return m_position;
}

std::vector<Border*> &Pave::get_borders(){
    return m_borders;
}

Border* Pave::get_border(int face){
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
    cout << "mode " << (get_inner_mode()?"inner":"outer") << endl;
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
    if(!m_compute_inner){
        return get_first_process_outer();
    }
    else{
        if(m_inner_mode)
            return get_first_process_inner();
        else
            return get_first_process_outer();
    }
}

bool Pave::get_first_process_inner() const{
    return m_first_process_inner;
}

bool Pave::get_first_process_outer() const{
    return m_first_process_outer;
}

void Pave::set_first_process_inner(bool val){
    m_first_process_inner = val;
}

void Pave::set_first_process_outer(bool val){
    m_first_process_outer = val;
}

void Pave::set_first_process(bool val){
    if(!m_compute_inner){
        return set_first_process_outer(val);
    }
    else{
        if(m_inner_mode)
            return set_first_process_inner(val);
        else
            return set_first_process_outer(val);
    }
}

void Pave::set_first_process_all(bool val){
    set_first_process_outer(val);
    set_first_process_inner(val);
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

std::vector<ibex::Interval> Pave::get_theta_union_list() const{
    if(m_backward_function)
        return m_theta_union_list_bwd;
    else
        return m_theta_union_list;
}

std::vector<ibex::Interval> Pave::get_theta_union_list_fwd() const{
    return m_theta_union_list;
}

std::vector<ibex::Interval> Pave::get_theta_union_list_bwd() const{
    return m_theta_union_list_bwd;
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

void Pave::set_diseable_singelton(bool val){
    m_diseable_singeleton = val;
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

bool Pave::is_marked() const{
    return m_marker;
}

void Pave::set_marker(bool val){
    m_marker = val;
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

bool Pave::is_near_external_border() const{
    for(Border *b:m_borders){
        for(Inclusion *i:b->get_inclusions()){
            if(i->get_border()->get_pave()->is_external_border())
                return true;
        }
    }
    return false;
}

///
/// \brief Pave::is_trajectory_external_escape
/// \return
/// Outer approximation of true (not all case covered)
///
bool Pave::is_trajectory_external_escape() const{
    for(Border *b:m_borders){
        if(m_theta_more_than_two_pi || !b->get_segment_out().is_empty()){
            for(Inclusion *i:b->get_inclusions()){
                if(i->get_border()->get_pave()->is_external_border())
                    return true;
            }
        }
    }
    return false;
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

bool Pave::is_near_removed_inner(){
    if(is_removed_pave_inner())
        return false;
    vector<Pave*> brothers = get_all_brothers();
    for(Pave *p:brothers){
        if(p->is_removed_pave_inner() && !p->is_external_border() && !p->is_empty()){
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
    for(Border *b:m_borders)
        b->set_backward_function(val);
}

bool Pave::get_backward_function() const{
    return m_backward_function;
}

int Pave::get_cpt_consistency_inner() const{
    return m_cpt_consistency_inner;
}
int Pave::get_cpt_consistency_outer() const{
    return m_cpt_consistency_outer;
}
int Pave::get_cpt_continuity_inner() const{
    return m_cpt_continuity_inner;
}
int Pave::get_cpt_continuity_outer() const{
    return m_cpt_continuity_outer;
}

void Pave::increment_cpt_consistency(){
    if(m_inner_mode)
        m_cpt_consistency_inner++;
    else
        m_cpt_consistency_outer++;
}

void Pave::increment_cpt_continuity(){
    if(!m_compute_inner){
        m_cpt_continuity_outer++;
    }
    else{
        if(m_inner_mode)
            m_cpt_continuity_inner++;
        else
            m_cpt_continuity_outer++;
    }
}

bool Pave::is_bassin() const{
    return m_bassin;
}

void Pave::set_bassin(bool val){
    m_bassin = val;
}

double Pave::get_area_outer() const{

    vector<double> x, y;
    for(Border *b:m_borders){
        b->get_points(x, y, false);
    }

    // Remove duplicate points
    for(int i=0; i<x.size()-1; i++){
        if(x[i]==x[i+1] && y[i]==y[i+1]){
            x.erase(x.begin()+i);
            y.erase(y.begin()+i);
            i--;
        }
    }

    // Compute area
    double area = 0.0;
    int nb_pt = x.size();
    for(int i=0; i<nb_pt; i++){
        area += x[i]*y[(i+1)%nb_pt]-y[i]*x[(i+1)%nb_pt];
    }
    area = 0.5*fabs(area);
    return area;
}

double Pave::get_perimeter() const{

    Pave *p1 = new Pave(this);
    p1->set_diseable_singelton(false);
    Pave *p2 = new Pave(this);
    p2->complementaire();
    p1->inter(p2);

    double perimeter = 0.0;
    vector<double> x, y;

    for(Border *b:p1->get_borders()){
        if(!b->is_empty()){
            x.push_back(b->get_segment_out_2D()[0].mid());
            y.push_back(b->get_segment_out_2D()[1].mid());
        }
    }

    if(x.size() == 2){
        perimeter = std::sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
    }
    else{
        IntervalVector box = get_position();
        perimeter = std::sqrt((box[0].diam())*(box[0].diam())+(box[1].diam())*(box[1].diam()));
    }

    delete(p1);
    delete(p2);
    return perimeter;
}

bool Pave::is_possible_path(IntervalVector ptA, IntervalVector ptB){
    Interval dx, dy;
    dx = ptB[0] - ptA[0];
    dy = ptB[1] - ptA[1];
    vector<Interval> theta_test_list = compute_theta(dx, dy);

    bool is_possible = false;
    if(theta_test_list[0].is_subset(get_theta_union_list()[0]))
        is_possible = true;
    if(theta_test_list.size()==2 && get_theta_union_list().size()==2){
        if(theta_test_list[1].is_subset(get_theta_union_list()[1]))
            is_possible = true;
    }

    return is_possible;
}

bool Pave::is_positive_invariant(){
    //    reset_full_empty();

    if(is_empty())
        return true;
    else if(is_trajectory_external_escape())
        return false;
    else if(is_full() && !is_theta_more_than_two_pi())
        return true;
    else if(is_theta_more_than_two_pi()){
        vector< vector<IntervalVector>> segment_list;
        bool one_not_empty = false;
        for(Border *b:get_borders()){
            for(Inclusion *i:b->get_inclusions()){
                Interval i_segment_in_union_out = i->get_border()->get_segment_in_union_out();
                if(!i_segment_in_union_out.is_empty())
                    one_not_empty = true;
                Interval seg_range =  b->get_segment_full() & i->get_border()->get_segment_full();
                if((i_segment_in_union_out & seg_range) != seg_range
                        || b->get_segment_in_union_out().is_empty()){

                    Interval neighbour_diff_A, neighbour_diff_B;
                    i->get_border()->get_segment_full().diff(i_segment_in_union_out, neighbour_diff_A, neighbour_diff_B);
                    Interval segment_in_out = b->get_segment_in_union_out();
                    IntervalVector segment_in_out_2D = b->get_segment_in_union_out_2D();

                    if(!(neighbour_diff_A & segment_in_out).is_empty()){
                        IntervalVector neighbour_diff_A_2D = i->get_border()->get_position();
                        neighbour_diff_A_2D[i->get_border()->get_face()%2] = neighbour_diff_A;
                        IntervalVector result = segment_in_out_2D & neighbour_diff_A_2D;

                        vector<IntervalVector> segment;
                        IntervalVector ptA(result.lb());
                        IntervalVector ptB(result.ub());
                        segment.push_back(ptA);
                        segment.push_back(ptB);
                        segment_list.push_back(segment);
                    }
                    if(!(neighbour_diff_B & segment_in_out).is_empty()){
                        IntervalVector neighbour_diff_B_2D = i->get_border()->get_position();
                        neighbour_diff_B_2D[i->get_border()->get_face()%2] = neighbour_diff_B;
                        IntervalVector result = segment_in_out_2D & neighbour_diff_B_2D;

                        vector<IntervalVector> segment;
                        IntervalVector ptA(result.lb());
                        IntervalVector ptB(result.ub());
                        segment.push_back(ptA);
                        segment.push_back(ptB);
                        segment_list.push_back(segment);
                    }

                }
            }
        }
        m_segment_list = segment_list;
        if(one_not_empty)
            return true;
        else
            return false;
    }
    else{
        // Get the list of points
        vector<double> x;
        vector<double> y;
        for(Border *b:m_borders){
            b->get_points(x, y);
        }

        // Put the first point at the last position
        x.push_back(x[0]);
        y.push_back(y[0]);
        x.erase(x.begin());
        y.erase(y.begin());

        // Compute the list of segments
        vector< vector<IntervalVector>> segment_list;
        vector< bool> segment_side;

        for(int i=0; i<x.size()-1; i+=2){
            if(x[i]!=x[i+1] || y[i]!=y[i+1]){
                IntervalVector ptA(2);
                ptA[0] = Interval(x[i]);
                ptA[1] = Interval(y[i]);
                IntervalVector ptB(2);
                ptB[0] = Interval(x[i+1]);
                ptB[1] = Interval(y[i+1]);

                vector<IntervalVector> segment;
                segment.push_back(ptA); segment.push_back(ptB);
                segment_list.push_back(segment);

                double segment1[2] = {x[i+1]-x[i],y[i+1]-y[i]};
                double segment2[2] = {x[(i+2)%x.size()]-x[i+1], y[(i+2)%x.size()]-y[i+1]};
                double v_product = segment1[0]*segment2[1]-segment2[0]*segment1[1];
                if(v_product>=0)
                    segment_side.push_back(true);
                else
                    segment_side.push_back(false);
            }
        }

        // Check if all segments are consistents with theta ?
        for(int i=0; i< segment_list.size(); i++){
            vector<ibex::Interval> theta_half_circle = compute_half_circle(segment_list[i][0], segment_list[i][1], segment_side[i]);
            bool is_inside = false;
            for(vector<Interval> &f_theta:m_theta_list){
                if(f_theta[0].is_subset(theta_half_circle[0]) || f_theta[0].is_subset(theta_half_circle[1])){
                    if(f_theta.size()==2){
                        if(theta_half_circle.size() == 2){
                            if(f_theta[1].is_subset(theta_half_circle[1])){
                                is_inside = true;
                                break;
                            }
                        }
                    }
                    else{
                        is_inside = true;
                        break;
                    }
                }
            }
            if(!is_inside)
                return false;
        }
        m_segment_list = segment_list;
        return true;
    }
}

const std::vector<ibex::Interval> Pave::compute_half_circle(const IntervalVector pt_start, const IntervalVector pt_end, bool trigo_rotation){
    IntervalVector ptA(2), ptB(2);
    vector<ibex::Interval> theta_list;

    if(trigo_rotation){
        ptA = pt_start;
        ptB = pt_end;
    }
    else{
        ptA = pt_end;
        ptB = pt_start;
    }

    IntervalVector delta(2);
    delta[0] = ptB[0] - ptA[0];
    delta[1] = ptB[1] - ptA[1];

    Interval thetaP = atan2(delta[1], delta[0]);
    Interval thetaN = atan2(-delta[1], -delta[0]);

    // Calcul de l'angle
    Interval theta_high = thetaP | Interval::PI;
    Interval theta_low = (-Interval::PI) | thetaN;

    Interval theta_inter = theta_high & theta_low;

    if(theta_inter.is_empty()){
        theta_list.push_back(theta_high + Interval(0, 1e-12)); // [0, pi] part
        theta_list.push_back(theta_low + Interval(-1e-12, 0)); // [-pi, 0] part
    }
    else{
        theta_list.push_back(theta_inter);
    }
    return theta_list;
}

vector<vector<IntervalVector> > Pave::get_segment_list(){
    return m_segment_list;
}

void Pave::reset_segment_list(){
    vector<vector<IntervalVector> > list_empty;
    m_segment_list = list_empty;
}

bool Pave::is_theta_more_than_two_pi() const{
    return m_theta_more_than_two_pi;
}

bool Pave::is_infinity_pave() const{
    return m_infinity_pave;
}

void Pave::set_infinity_pave(bool val, IntervalVector search_box){
    m_infinity_pave = val;
    m_search_box = search_box;
}

ibex::IntervalVector Pave::get_search_box() const{
    return m_search_box;
}

ibex::IntervalVector Pave::get_vector_field() const{
    return m_vector_field;
}
