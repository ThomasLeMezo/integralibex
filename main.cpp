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

    test_CtcPaveForward();
//    test_CtcPaveConsistency();

//    test_contractor_polar();

//    test_rotation();

//    test_diff();
//    test_copy_graph();

//    sandbox();
}

void van_der_pol_cycle(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(2);
    box[0] = Interval(-4.0, 4.0);
    box[1] = Interval(-4.0, 4.0);

    Interval u = Interval::ZERO;
    Scheduler s(box, &f, u);

    s.cameleon_cycle(12, 5, 1e9, false, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);
}

void ball(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,-10-1*y));

    IntervalVector box(2);
    box[0] = Interval(0.0, 20.0);
    box[1] = Interval(-20.0, 20.0);

    Interval u = Interval::ZERO;
    Scheduler s(box, &f, u);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(12.0);
    activated_pave[1] = Interval(-2.0);

    ibex::Function f_sym(x, y, Return(x, -y));
    s.set_symetry(&f_sym,3, 3);

    s.cameleon_propagation(15, 1000000, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
}

void capture_attractor(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable phi, d;
    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
                                    -cos(phi)));

    IntervalVector box(2);
    box[0] = -Interval::PI | Interval::PI;
    box[1] = Interval(0, 10.0);

    Interval u = Interval::ZERO;
    Scheduler s(box, &f, u);

    ibex::Function f_sym23(phi, d, Return(phi-Interval::TWO_PI, d));
    s.set_symetry(&f_sym23,1, 3);

    ibex::Function f_sym32(phi, d, Return(phi+Interval::TWO_PI, d));
    s.set_symetry(&f_sym32,3, 1);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(2);
    activated_pave[1] = Interval(3.0);

    s.cameleon_cycle(2, 5, 1e9, false, false);
//    s.cameleon_propagation(15, 1e6, activated_pave, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, false);
}

int main()
{
//    ball();
    capture_attractor();
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
