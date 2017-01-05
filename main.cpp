#include <ibex.h>
#include <vibes.h>

#include <pave.h>
#include <scheduler.h>
#include <tests.h>

#include <iostream>
#include <ctime>

#include <graphdot.h>

using namespace std;
using namespace ibex;

void test(){
    //    testTranslate();
    //    testRotate();
    //    test_CtcPropagateLeftSide();
    //    test_CtcPropagateRightSide();
    //    test_CtcPropagateFront();
    //    test_CtcPropagateSegment();

    //    test_CtcPaveForward();
    //    test_CtcConsistency();
    //    test_CtcConsistency2();
    //    test_CtcConsistency3();
    //    test_CtcConsistency_Kernel();
    //    test_CtcConsistency_Kernel2();
    //    test_CtcConsistency_Kernel3();
    //    test_CtcConsistency_Kernel4();

    //    test_inter_pave_perimeter();
    //    test_possible_path();

    //    test_contractor_polar();

    //    test_rotation();

    //    test_diff();
    //    test_copy_graph();

    //    test_imageIntegral();
    //    test_car_on_hill();

    //    sandbox();

    //    test_infinity();
    //    test_chi_function();

    //    Variable x, y;
    //    try{
    //        ibex::Function f("function.txt");
    //    }
    //    catch(const ibex::SyntaxError & e){
    //        std::cerr << e << endl;
    //    }

    test_diff_infinity();
}

void van_der_pol_cycle(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,(1.0*(1.0-pow(x, 2))*y-x)));

    //    Variable q, p;
    //    ibex::Function f(q, p, Return(4*p*(q*q+p*p)+20*p,
    //                                   -(4*q*(q*q+p*p)-20*q)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4,4);
    box[1] = Interval(-4,4);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, false, false, false, false);

    //int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside, bool near_bassin, bool stop_first_pos_invariant
    s.cameleon_cycle(10, 5, 1e9, true, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    //    cout << "AREA OUTER = " << s.get_graph_list(0)->get_area_outer() << endl;
    //    vector<double> perimeters = s.get_graph_list(0)->get_perimeters();
    //    cout << "perimeters = ";
    //    for(double p:perimeters)
    //        cout << "[" << p << "]";
    //    cout << endl;
    s.draw(1024, true);
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void simon_cos(){

}

void ball(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,-10-0.1*y*abs(y)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(0.0, 20.0);
    box[1] = Interval(-20.0, 20.0);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, false);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(10.0, 16.0);
    activated_pave[1] = Interval(-4.0, 4.0);

    ibex::Function f_sym(x, y, Return(x, -y-2.0));
    s.set_symetry(&f_sym,3, 3);

    //    s.cameleon_propagation(20, 1000000, activated_pave);
    s.cameleon_propagation_with_inner(14,1e9,activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
}

void hybrid_system(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(1.0+0.0*x,1.0+0.0*x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(0, 20);
    box[1] = Interval(0, 1);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, false);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(0.0, 0.5);
    activated_pave[1] = Interval(0.0,0.1);

    ibex::Function f_sym(x, y, Return(x, y-1));
    s.set_symetry(&f_sym,2,0);

    ibex::Function f_sym2(x, y, Return(x, y+1));
    s.set_symetry(&f_sym2,0,2);

//        s.cameleon_propagation(10, 1e9, activated_pave);
    s.cameleon_propagation_with_inner(10,1e9,activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
//    s.print_pave_info(0, 1.4,0.86);
//    s.print_pave_info(0, 1.4,0.33);
}

void station_keeping_attractor(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable phi, d;
    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1.0/d-1)*sin(phi)),
                                    -cos(phi)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = -Interval::PI | Interval::PI;
    //    box[1] = Interval(1e-5, 10.0);
    box[1] = Interval(0.0, 10.0);

    //    box[0] = Interval(-2, -1.1);
    //    box[1] = Interval(0.6,1.6);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, false);

    /////////////// Symetries ///////////////
    ibex::Function f_sym23(phi, d, Return(phi-Interval::TWO_PI, d));
    s.set_symetry(&f_sym23,1, 3);

    ibex::Function f_sym32(phi, d, Return(phi+Interval::TWO_PI, d));
    s.set_symetry(&f_sym32,3, 1);

    ibex::Function f_sym00(phi, d, Return(phi+Interval::TWO_PI, d));
    s.set_symetry(&f_sym00,0, 0);

    /////////////// Inner curve ///////////////
    ibex::Function f_inside_curve(phi, d, d*(d+2*sin(phi))+1/2.0);
    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    //    s.cameleon_cycle(14, 5, 1e9, false, false);
    s.cameleon_viability(8, 1e9, true);
    //    s.cameleon_propagation(15, 1e6, activated_pave, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true, "", false);
    s.print_pave_info(0, -1.17,1.47,"b[b]");

    /// Truth
    vector<double> x, y;
    double c=1.0/sqrt(2.0);
    for(double t=-M_PI; t<=M_PI; t+=0.01){
        y.push_back(sqrt(pow(cos(t)+c, 2)+pow(sin(t), 2)));
        double phi = atan2(cos(t), -sin(t))+M_PI-atan2(sin(t), cos(t)+c);
        x.push_back(2*atan(tan(phi/2.0)));
    }
    vibes::drawPolygon(x, y, "red[]");

    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);

}

void car_on_the_hill_attractor(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    //    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 13.0);
    box[1] = Interval(-16, 16);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, true, false);
    //    vector<IntervalVector> bassin_list;
    //    Scheduler s(box, bassin_list, f_list, u, true, false, false);

    /////////////// Compute ///////////////
    s.compute_attractor(12, 1e9);
    //    s.cameleon_cycle(12, 5, 1e9, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true, "invert f");
    s.print_pave_info(0, 10.2,-4.5,"b[b]");
}

void car_on_the_hill_kernel(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 13.0);
    box[1] = Interval(-16, 16);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, MAZE_BORDER_INNER_IN_FULL, MAZE_BORDER_INNER_OUT_FULL, MAZE_BORDER_OUTER_IN_FULL, MAZE_BORDER_OUTER_OUT_EMPTY);

    /////////////// Compute ///////////////
    bool is_attractor = s.compute_attractor(20, 1e9);
    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true, "attractor"); vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);

    if(is_attractor){
        s.attractor_to_kernel();
        //    s.draw(1024, true, "invert"); vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
        s.cameleon_viability(3, 1e9);

        cout << "************************" << endl;
        cout << "TOTAL TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

        /////////////// Drawing ///////////////
        s.draw(1024, true);
        vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
        //    s.print_pave_info(0, 3.15,-3.2,"b[b]");
        //    s.print_pave_info(0, 12.5,0.2,"b[b]");
    }
}

void car_on_the_hill_bassin(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));

    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    ibex::Function f3(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2));

    std::vector<ibex::Function*> f_list;
    //    f_list.push_back(&f1);
    //    f_list.push_back(&f2);
    f_list.push_back(&f3);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 13.0);
    //    box[1] = Interval(-16, 16);
    box[1] = Interval(-8, 10);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(2.2,8.8) + Interval(-0.1, 0.1);
    box_remove[1] = Interval(-0.5,0.5);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    // int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside, bool compute_inner
    //    s.cameleon_cycle(15, 5, 1e9, false, false, true);
    s.cameleon_viability(7, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true);

    //    s.print_pave_info(0, 5.5, 0.0,"b[b]");
    //    s.print_pave_info(0, -0.5, 0.44,"b[b]");

}

void cercle_capture_bassin(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2, x, y;
    //    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    ibex::Function f(x, y, Return(x-(x+y)*(x*x+y*y),
                                  y+(x-y)*(x*x+y*y)));


    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4.0, 4.0);
    //    box[1] = Interval(-16, 16);
    box[1] = Interval(-4.0, 4.0);

    // Points d'équilibre (stable)
    // x1 = 2.311(6-7)
    // x1 = 7.78(10/09)

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);
    //    box_remove[0] = Interval(2.2,2.4) + Interval(-0.5, 0.5);
    //    box_remove[1] = Interval(-0.1,0.1);
    //    list_boxes_removed.push_back(box_remove);
    //    box_remove[0] = Interval(0.9, 1.1);
    //    box_remove[1] = Interval(-0.1, 0.1);

    box_remove[0] = Interval(1.9, 2.1);
    box_remove[1] = Interval(-0.1, 0.1);

    list_boxes_removed.push_back(box_remove);

    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_cycle(15, 5, 1e9, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);

    //    s.print_pave_info(0, -0.8, -0.7,"b[b]");
    //    s.print_pave_info(0, -0.5, 0.44,"b[b]");
}

void pendulum_capture_bassin(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x2,
                                    -9.81/1.0*sin(x1)-1.0*x2));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4, 4);
    //    box[1] = Interval(-16, 16);
    box[1] = Interval(-7, 7);

    // Points d'équilibre (stable)
    // x1 = 2.311(6-7)
    // x1 = 7.78(10/09)

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);
    box_remove[0] = Interval(-0.1, 0.1);
    box_remove[1] = Interval(-0.1,0.1);
    list_boxes_removed.push_back(box_remove);

    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_cycle(14, 5, 1e9, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
    //    vibes::axisLimits(-1,13, -8,11);

    //    s.print_pave_info(0, -0.8, -0.7,"b[b]");
    //    s.print_pave_info(0, -0.5, 0.44,"b[b]");
}

void car_on_the_hill_limit_path(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    //    ibex::Function f1(x1, x2, Return(-x2,
    //                                    -(-9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0)));
    ibex::Function f1(x1, x2, Return(-x2,
                                     -(-9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0)));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 13.0);
    box[1] = Interval(-16, 16);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, false, false, false, false);

    IntervalVector activated_pave(2);
    //    activated_pave[0] = Interval(11.028646, 11.028647); // Point limite : x0 = 11.02864(6-7)
    activated_pave[0] = Interval(-1.0);
    activated_pave[1] = Interval(0.0);

    s.cameleon_propagation(25, 1e9, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
}

void car_on_the_hill_integrator(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));

    //    ibex::Function f2(x1, x2, Return(x2,
    //                                    -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    //    ibex::Function f3(x1, x2, Return(x2,
    //                                    -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2));

    std::vector<ibex::Function*> f_list;
    //    f_list.push_back(&f2);
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-2.0, 13.0);
    box[1] = Interval(-6, 6);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, false);

    vector<IntervalVector> initial_pave_list;
    IntervalVector activated_pave(2);
    //    activated_pave[0] = Interval(11.028646, 11.028647); // Point limite : x0 = 11.02864(6-7)
    activated_pave[0] = Interval(-1.0, 1.0);
    activated_pave[1] = Interval(-1.0, 1.0);

    //    activated_pave[0] = Interval(6.5); // 0.4
    //    activated_pave[1] = Interval(4.5);
    initial_pave_list.push_back(activated_pave);
    //    activated_pave[0] = Interval(7.7809, 7.7810);
    //    initial_pave_list.push_back(activated_pave);

    //    s.cameleon_propagation(18, 1e9, activated_pave); // 25
    s.cameleon_propagation_with_inner(15, 1e9, activated_pave); // 25

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);
    vibes::drawBox(activated_pave, "red[]");
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
}

void pendulum_cycle(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x2,
                                    -9.81/1.0*sin(x1)-1.0*x2));
    //    ibex::Function f(x1, x2, Return(x2,
    //                                    -9.81/1.0*sin(x1)));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-3,3);
    box[1] = Interval(-3,3);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, true);

    /////////////// Inner curve ///////////////
    // x2^2/2 +g(1-cos(x1)) <= 9/2
    ibex::Function f_inside_curve(x1, x2, pow(x2,2)/2.0 +9.81*(1.0-cos(x1))-4.5);
    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    s.cameleon_viability(15, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void pendulum_integration(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81/1.0*sin(x1+(-Interval::PI/6.0|Interval::PI/6.0))-1.0*x2));
    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81/1.0*sin(x1+(-Interval::PI/6.0|Interval::PI/6.0))-1.0*x2+9.81*sin(x1+(-Interval::PI/6.0|Interval::PI/6.0))));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(-3,3);
    box[1] = Interval(-3,3);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, MAZE_BORDER_INNER_IN_FULL, MAZE_BORDER_INNER_OUT_FULL, MAZE_BORDER_OUTER_IN_EMPTY, MAZE_BORDER_OUTER_OUT_EMPTY);

    /////////////// Initial condition ///////////////
    IntervalVector initial_box(2);
    initial_box[0] = Interval(0.0);
    initial_box[1] = Interval(0.0);

    /////////////// Compute ///////////////
    s.cameleon_propagation_with_inner(10,1e9,initial_box);
    //    s.cameleon_propagation(13,1e9,initial_box);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    vibes::drawCircle(initial_box[0].mid(), initial_box[0].mid(), 0.1, "red[]");
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void hann_max_pos_inv(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(-x1+2*(pow(x1,2))*x2,
                                    -x2));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4,4);
    box[1] = Interval(-4,4);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, true);

    /////////////// Inner curve ///////////////
    ibex::Function f_inside_curve(x1, x2, pow(x1,2) + pow(x2,2) - 1/2.0);

    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    s.cameleon_viability(8, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void predator_prey_max_pos_inv(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    // −3x 1 + 4x 21 − 0.5x 1 x 2 − x 31
    // −2.1x 2 + x 1 x 2
    ibex::Function f(x1, x2, Return(-3.0*x1+4.0*pow(x1,2)-0.5*x1*x2-pow(x1, 3),
                                    -2.1*x2+x1*x2));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-0.01,4);
    box[1] = Interval(-0.01,5);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, true);

    /////////////// Inner curve ///////////////
    //Log[(x2/x1)] - (1/2)/(21/10) (2 x1 + x2) >= 0 & x2 <= 3 & x2 >=0 & x1 >= 0
    //    ibex::Function f_inside_curve(x1, x2, max(-x1, max(x2-3, -(log(x2/x1)-0.5/2.1*(2*x1+x2)))));
    ibex::Function f_inside_curve(x1, x2, 0.5*((-(21.0/10.0) + x1)*(505.0/147.0*(-(21.0/10.0) + x1) + 10.0/21.0*(-(99.0/50.0) + x2)) + (10.0/21.0*(-(21.0/10.0) + x1) + (2665.0*(-(99.0/50.0) + x2))/1386.0)*(-(99.0/50.0) + x2)) - 1.0/8.0);

    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    s.cameleon_viability(8, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    s.print_pave_info(0, 0.0,0.0,"b[b]");
}

void circle3_max_pos_inv(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x1-x2-pow(x1,3),
                                    x1-pow(x1,2)*x2));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4,4);
    box[1] = Interval(-4,4);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, true);

    /////////////// Inner curve ///////////////
    ibex::Function f_inside_curve(x1, x2, pow(x1,2) + pow(x2,2) - 3.0);
    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    s.cameleon_viability(11, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void dipole_max_pos_inv(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return((x1-1.0)/pow(x2*x2+pow(x1-1.0,2),1.5)-(x1+1.0)/pow(x2*x2+pow(x1+1.0,2),1.5),
                                    x2/pow(x2*x2+pow(x1-1.0,2),1.5)-x2/pow(x2*x2+pow(x1+1.0,2),1.5)));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4,4);
    box[1] = Interval(-4,4);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, true);

    /////////////// Inner curve ///////////////
    ibex::Function f_inside_curve(x1, x2, pow(x1+1,2) + pow(x2,2) - 1.0);
    s.push_back_inside_curve(&f_inside_curve);

    /////////////// Compute ///////////////
    s.cameleon_viability(15, 1e9, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true, "", false);
    //    s.print_pave_info(0, -3.8, -3.97,"b[b]");
}

void integrator(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2, q, p;
    ibex::Function f0(x1, x2, Return(1.0+0.0*x1,
                                     -sin(x2)+ Interval(-0.1, 0.1)));
    //    ibex::Function f0(q, p, Return(4*p*(q*q+p*p)+20*p,
    //                                   -(4*q*(q*q+p*p)-20*q)));
//        ibex::Function f0(x1, x2, Return(1.0+0.0*x1, cos(x2)*cos(x2)));
//        ibex::Function f1(x1, x2, Return(1.0+0.0*x1,
//                                        -sin(x2)-0.1));
//        ibex::Function f2(x1, x2, Return(1.0+0.0*x1,
//                                        -sin(x2)+0.1));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f0);
//      f_list.push_back(&f1);
//      f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(0,5);
    box[1] = Interval(-2,2);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, false, false);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(0.3,0.5);
    activated_pave[1] = Interval(-0.5,0.5);

//    s.cameleon_propagation(13, 1e9, activated_pave);
    s.cameleon_propagation_with_inner(13, 1e9, activated_pave);

    cout << endl << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);
    vibes::drawBox(activated_pave, "red[]");
//    s.print_pave_info(0,0.0,0.0);
//    s.print_pave_info(0, 5.5, 0.0);
//    s.print_pave_info(0, 2.5, 2.3);
}

void van_der_pol_integration(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f1(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-6.0, 6.0);
    box[1] = Interval(-6.0, 6.0);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, false, false);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(-4.0, -3.0);
    activated_pave[1] = Interval(3.0, 4.0);

    //    s.cameleon_propagation(20, 1e9, activated_pave); // 25
    s.cameleon_propagation_with_inner(15, 1e9, activated_pave); // 25

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
    vibes::drawBox(activated_pave, "red[]");

    //    s.print_pave_info(0, -3.99, 3.055);
}

void van_der_pol_kernel(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f1(x, y, Return(y,-(1.0*(1.0-pow(x, 2))*y-x + 1.0)));
    ibex::Function f2(x, y, Return(y,-(1.0*(1.0-pow(x, 2))*y-x - 1.0)));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(-8.0, 8.0);
    box[1] = Interval(-8.0, 8.0);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, false, false);

    /////////////// Compute ///////////////
    bool is_attractor = s.compute_attractor(14, 1e9);
    s.draw(1024, true, "attractor");
    if(is_attractor){
        s.attractor_to_kernel();
        s.draw(1024, true, "invert");
        s.cameleon_viability(3, 1e9);
    }

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    if(is_attractor)
        s.draw(1024, true);
}

void bassin_van_der_pol(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f1(x, y, Return(-y,-(1.0*(1.0-pow(x, 2))*y-x)));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-4.0, 4.0);
    box[1] = Interval(-4.0, 4.0);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    //    box_remove[0] = Interval(1.0) + Interval(-0.4, 0.4);
    //    box_remove[1] = Interval(1.0) + Interval(-0.4, 0.4);
    box_remove[0] = Interval(-0.4, 0.4);
    box_remove[1] = Interval(-0.4, 0.4);

    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    // int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside, bool compute_inner
    //    s.cameleon_cycle(15, 5, 1e9, false, false, true);
    s.cameleon_viability(7, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true);

    //    s.print_pave_info(0, 5.5, 0.0,"b[b]");
    //    s.print_pave_info(0, -0.5, 0.44,"b[b]");

}

void car_on_the_hill_trajectory(){
//    const clock_t begin_time = clock();
    struct timeval start, end;
    gettimeofday(&start, NULL);


    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
//    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(-2.0, 13.0);
    box[1] = Interval(-10, 10);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, MAZE_BORDER_INNER_IN_FULL, MAZE_BORDER_INNER_OUT_FULL,
                                                        MAZE_BORDER_OUTER_IN_EMPTY, MAZE_BORDER_OUTER_OUT_EMPTY);

    IntervalVector paveA(2);
//    paveA[0] = Interval(0.0, 1.0);
//    paveA[1] = Interval(0.0, 1.0);
    paveA[0] = Interval(-1, 2.0);
    paveA[1] = Interval(-6, 6);

    IntervalVector paveB(2);
    //    paveB[0] = Interval(2.0, 3.0);
    //    paveB[1] = Interval(2.0, 3.0);
    paveB[0] = Interval(11.5,12.0);
    paveB[1] = Interval(-6.0, 6.0);

    s.find_path(15, 1e9, paveA, paveB); // 25

    gettimeofday(&end, NULL);

    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;

//    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    cout << "TIME = " << delta << endl;
    s.draw(1024, true);
    vibes::drawBox(paveA, "g[]");
    vibes::drawBox(paveB, "g[]");
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
}

void cos_trajectory(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(1.0 + 0.0*x1,
                                     sin(x1)+Interval(-0.05, 0.05)));
//    ibex::Function f2(x1, x2, Return(x2,
//                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
//    f_list.push_back(&f2);

    IntervalVector box(2);
    box[0] = Interval(0.0, 20.0);
    box[1] = Interval(-1.0, 4.0);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_OFF, true, true, false, false);

    IntervalVector paveA(2);
    paveA[0] = Interval(0.0,0.08);
    paveA[1] = Interval(-0.078,0.078);

    IntervalVector paveB(2);
    paveB[0] = Interval(1,18.5);
    paveB[1] = Interval(2.1,3);

//    s.find_path(11, 1e9, paveA, paveB); // 25
    s.cameleon_propagation(1, 1e9, paveA);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
    vibes::drawBox(paveA, "g[]");
    vibes::drawBox(paveB, "g[]");
    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
}

void bassin_ratschan6(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2+(1-x1*x1-x2*x2)*x1,
                                     -x1+(1-x1*x1-x2*x2)*x2));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-8.0, 8.0);
    box[1] = Interval(-8.0, 8.0);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(-2, 2);
    box_remove[1] = Interval(-2, 2);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_viability(10, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
}

void bassin_parrilo(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(-x1+(1+x1)*x2,
                                     -(1+x1)*x1));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-8.0, 8.0);
    box[1] = Interval(-8.0, 8.0);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(-0.7,0.9);
    box_remove[1] = Interval(-0.7,0.9);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_viability(10, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
}

void bassin_ratschan3(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(-4*x1*x1*x1+6*x1*x1-2*x1,
                                     -2*x2));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-8,8);
    box[1] = Interval(-8,8);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(-2,2);
    box_remove[1] = Interval(-2,2);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_viability(10, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
}

void integrator_ratschan3(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(-(-4*x1*x1*x1+6*x1*x1-2*x1),
                                     -(-2*x2)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 2.0);
    box[1] = Interval(-1.0, 1.0);
    //    box[0] = Interval(-8,8);
    //    box[1] = Interval(-8,8);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, false, false);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(-0.1,0.1);
    activated_pave[1] = Interval(-0.1,0.1);
    //    activated_pave[0] = Interval(-2,2);
    //    activated_pave[1] = Interval(-2,2);

    //    s.cameleon_propagation(19, 1e9, activated_pave); // 25
    s.cameleon_propagation_with_inner(13, 1e9, activated_pave); // 25

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
    s.print_pave_info(0, 0.4,0.1);

}

void bassin_genesio(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(-x1+x2,
                                     0.1*x1-2*x2-x1*x1-0.1*x1*x1*x1));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-100,100);
    box[1] = Interval(-1000,1000);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(-0.8,0.8);
    box_remove[1] = Interval(-0.8,0.8);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_viability(10, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
    //    vibes::axisLimits(box[0].lb()-1.0,box[0].ub()+1.0, box[1].lb()-1.0,box[1].ub()+1.0);
}

void integrator_genesio(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(-(-x1+x2),
                                     -(0.1*x1-2*x2-x1*x1-0.1*x1*x1*x1)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-70,70);
    box[1] = Interval(-200,200);

    Scheduler s(box, f_list, MAZE_DISEABLE_SINGLETON_ON, true, true, false, false);

    IntervalVector activated_pave(2);
    //    activated_pave[0] = Interval(-10,10);
    //    activated_pave[1] = Interval(-10,10);
    activated_pave[0] = Interval(-3,3);
    activated_pave[1] = Interval(-3,3);

    //    s.cameleon_propagation(19, 1e9, activated_pave); // 25
    s.cameleon_propagation_with_inner(18, 1e9, activated_pave); // 25

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);

}

void bassin_bacha(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                     -0.5*x2-sin(x1+0.412)+sin(0.412)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);

    IntervalVector box(2);
    box[0] = Interval(-6,6);
    box[1] = Interval(-6,6);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);

    box_remove[0] = Interval(-0.42,0.42);
    box_remove[1] = Interval(-0.42,0.42);
    list_boxes_removed.push_back(box_remove);

    // const IntervalVector &box, const vector<IntervalVector> &bassin_boxes, const std::vector<ibex::Function *> &f_list,
    // const IntervalVector &u, bool diseable_singleton, bool border_in, bool border_out
    Scheduler s(box, list_boxes_removed, f_list, MAZE_DISEABLE_SINGLETON_ON, false, true); // diseable singleton = true

    /////////////// Compute ///////////////
    s.cameleon_viability(8, 1e9, true); // 10 = 256 s

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    /////////////// Drawing ///////////////
    s.draw(1024, true);
}


int main()
{
    /// **** BALL ***** //
//        ball();
//        hybrid_system();

    /// **** STATION KEEPING ***** //
    //    station_keeping_attractor();

    /// **** CAR ON THE HILL ***** //
    //    car_on_the_hill_attractor();
    //    car_on_the_hill_bassin();

//        car_on_the_hill_kernel();

//        car_on_the_hill_trajectory();
//        car_on_the_hill_integrator();
    //    car_on_the_hill_limit_path();

    /// **** CAPTURE BASSIN ***** //
    //    pendulum_capture_bassin();
    //    cercle_capture_bassin();

    /// **** VAN DER POL ***** //
        van_der_pol_cycle();
//        van_der_pol_integration();
    //    van_der_pol_kernel();

    /// **** INTEGRATOR ***** //
//        integrator();

    /// **** BASSIN ***** //
    //    bassin_ratschan6();
    //    bassin_ratschan3();
    //    bassin_parrilo();
    //    bassin_genesio();
    //    bassin_bacha();
    //    bassin_van_der_pol();

    //    integrator_genesio();
    //    integrator_ratschan3();

    /// **** Max POS INV ***** //
    //    hann_max_pos_inv();
    //    predator_prey_max_pos_inv();
    //    circle3_max_pos_inv();
    //    dipole_max_pos_inv();

    /// **** PENDULUM ***** //
    //  pendulum_cycle();
//      pendulum_integration();


    /// **** DUAL TUBE **** //
//    cos_trajectory();

    /// **** TEST ***** //
    //    test();
    return 0;
}
