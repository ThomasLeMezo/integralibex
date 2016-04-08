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
//    test_CtcPaveConsistency();
    test_CtcPaveConsistency2();

//    test_contractor_polar();

//    test_rotation();

//    test_diff();
//    test_copy_graph();

//    test_imageIntegral();
//    test_car_on_hill();

//    sandbox();
}

void van_der_pol_cycle(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-4.0, 4.0);
    box[1] = Interval(-4.0, 4.0);

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval::ZERO;
    Scheduler s(box, f_list, u);

    s.cameleon_cycle(13, 5, 1e9, true, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);
//    s.print_pave_info(0, -2.6,-2.5,"b[b]");
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

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval::ZERO;
    Scheduler s(box, f_list, u);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(12.0);
    activated_pave[1] = Interval(-2.0);

    ibex::Function f_sym(x, y, Return(x, -y-2.0));
    s.set_symetry(&f_sym,3, 3);

    s.cameleon_propagation(20, 1000000, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
}

void capture_attractor(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable phi, d;
    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
                                    -cos(phi)));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = -Interval::PI | Interval::PI;
    box[1] = Interval(0.01, 10.0);

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval::ZERO;
    Scheduler s(box, f_list, u);

    /////////////// Symetries ///////////////
    ibex::Function f_sym23(phi, d, Return(phi-Interval::TWO_PI, d));
    s.set_symetry(&f_sym23,1, 3);

    ibex::Function f_sym32(phi, d, Return(phi+Interval::TWO_PI, d));
    s.set_symetry(&f_sym32,3, 1);

    /////////////// Inner ///////////////
    Variable t;
    ibex::Function f_inner(t, Return(2*atan(tan((atan2(cos(t), -sin(t))+Interval::PI-atan2(sin(t), cos(t)+1.0/sqrt(2.0)))/2.0)),
                               sqrt(pow(cos(t)+1/sqrt(2.0), 2)+pow(sin(t), 2))));
    s.set_imageIntegral(box, &f_inner, Interval::ZERO | Interval::TWO_PI, 15,5000);

    /////////////// Compute ///////////////
    s.cameleon_cycle(15, 5, 1e9, false, false, false);
//    s.cameleon_propagation(15, 1e6, activated_pave, false);  

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true);
//    s.print_pave_info(0, -1.64,0.11,"b[b]");

    /// Truth
    vector<double> x, y;
    double c=1.0/sqrt(2.0);
    for(double t=-M_PI; t<=M_PI; t+=0.01){
        y.push_back(sqrt(pow(cos(t)+c, 2)+pow(sin(t), 2)));
        double phi = atan2(cos(t), -sin(t))+M_PI-atan2(sin(t), cos(t)+c);
        x.push_back(2*atan(tan(phi/2.0)));
    }
    vibes::drawPolygon(x, y, "blue[]");

}

void car_on_the_hill(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(0.0, 12.0);
    box[1] = Interval(-6.0, 6.0);

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval(-0.5, 0.5);

    Scheduler s(box, f_list, u);

    /////////////// Compute ///////////////
    s.cameleon_cycle(14, 5, 1e9, false, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true);
//    s.print_pave_info(0, -1.64,0.11,"b[b]");

}

void car_on_the_hill_v2(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));

    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f2);
    f_list.push_back(&f);

    IntervalVector box(2);
    box[0] = Interval(-1.0, 13.0);
    box[1] = Interval(-16.0, 16.0);

    std::vector<IntervalVector> list_boxes_removed;
    IntervalVector box_remove(2);
    box_remove[0] = Interval(-0.4,0.4);
    box_remove[1] = Interval(-0.3,0.3);
    list_boxes_removed.push_back(box_remove);
    box_remove[0] = Interval(2.4,3.4);
    box_remove[1] = Interval(-0.3,0.3);
    list_boxes_removed.push_back(box_remove);
    box_remove[0] = Interval(5.4,6.4);
    box_remove[1] = Interval(-0.3,0.3);
    list_boxes_removed.push_back(box_remove);
    box_remove[0] = Interval(8.4,9.4);
    box_remove[1] = Interval(-0.3,0.3);
    list_boxes_removed.push_back(box_remove);
    box_remove[0] = Interval(11.9,12.4);
    box_remove[1] = Interval(-0.3,0.3);
    list_boxes_removed.push_back(box_remove);

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval(-0.5, 0.5);

    Scheduler s(box, list_boxes_removed, f_list, u);

    /////////////// Compute ///////////////
    s.cameleon_cycle(3, 5, 1e9, false, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    /////////////// Drawing ///////////////
    s.draw(1024, true);


//    s.print_pave_info(0, 3.03, 0.14,"b[b]");
//    s.print_pave_info(0, 0.2,5.96,"b[b]");

}

int main()
{
//    ball();
//    capture_attractor();
//    car_on_the_hill();
    car_on_the_hill_v2();
//    van_der_pol_cycle();
//    test();

#if 0
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    Interval u = -Interval::HALF_PI/2 | Interval::HALF_PI/2;
    Scheduler s(box, &f, u);

    vector<IntervalVector> activated_paves;
    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(-2.2);
    activated_pave[1] = Interval(0.0);
    activated_paves.push_back(activated_pave);

    s.cameleon_propagation(15, 1000000, activated_paves, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);

#endif

#if 0
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

//    Interval u = -2*Interval::PI/3 | 2*Interval::PI/3;
    Interval u = Interval::ZERO;
    Scheduler s(box, &f, u);

    s.cameleon_cycle(10, 5, 1e9, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);

//    s.print_pave_info(-5, 6, "r[]");

#endif
    return 0;
}
