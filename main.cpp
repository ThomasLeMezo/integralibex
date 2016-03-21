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
namespace PPL = Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;

void test(){

//    test_CtcPropagateSegment();
    test_CtcPropagateSegmentBackward();

//    test_copy_graph();
//    sandbox();
}

void van_der_pol_cycle(){
    const clock_t begin_time = clock();
//    vibes::beginDrawing();
    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
//    ibex::Function f(x, y, Return(-0.9*sin(y), -0.9*cos(x)));

    IntervalVector box(2);
    box[0] = ibex::Interval(-4.0, 4.0);
    box[1] = ibex::Interval(-4.0, 4.0);

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;
    Scheduler s(box, &f, u);

    s.cameleon_cycle(12, 5, 1e6, false, false);
//    s.cameleon_cycle(12, 5, 1e9, true, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

//    s.draw();
}

void lorenz_attractor(){
    const clock_t begin_time = clock();
    ibex::Variable x, y, z;
    ibex::Function f(x, y, z, Return(10.0*(y-x), 28.0*x-y-x*z, x*y-8.0/3.0*z));

    IntervalVector box(3);
    box[0] = ibex::Interval(-30.0, 30.0);
    box[1] = ibex::Interval(-30.0, 30.0);
    box[2] = ibex::Interval(-30.0, 30.0);

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;
    Scheduler s(box, &f, u);

    s.cameleon_cycle(10, 5, 1e6, false, false);
//    s.cameleon_cycle(12, 5, 1e9, true, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw();
}

void ball(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,-10-1*y));

    IntervalVector box(2);
    box[0] = ibex::Interval(0.0, 20.0);
    box[1] = ibex::Interval(-20.0, 20.0);

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;
    Scheduler s(box, &f, u);

    IntervalVector activated_pave(2);
    activated_pave[0] = ibex::Interval(12.0);
    activated_pave[1] = ibex::Interval(-2.0);

    ibex::Function f_sym(x, y, Return(x, -y));
    s.set_symetry(&f_sym,1, 0, 1, 0);

    s.cameleon_propagation(15, 1000000, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw();
}

void integration(){
    int dim = 3;
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    ibex::Variable x, y, z;
//    ibex::Function f(x, y, z, Return(10.0*(y-x), 28.0*x-y-x*z, x*y-8.0/3.0*z));
//    ibex::Function f(x, y, z, Return(1.0+0.0*x, -sin(x)+cos(y), 1.0+0.0*z));
    ibex::Function f(x, y, z, Return(1.0+0.0*x, 0.0*y, -sin(x)));
//    ibex::Function f(x, y, v, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(dim);
    box[0] = ibex::Interval(-30.0, 30.0);
    box[1] = ibex::Interval(-30.0, 30.0);
    box[2] = ibex::Interval(-30.0, 30.0);

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;
    Scheduler s(box, &f, u);

    IntervalVector activated_pave(dim);
    activated_pave[0] = ibex::Interval(0.0);
    activated_pave[1] = ibex::Interval(0.0);
    activated_pave[2] = ibex::Interval(0.0);

//    ibex::Function f_sym(x, y, Return(x, -y));
//    s.set_symetry(&f_sym,3, 3);

    s.cameleon_propagation(8, 100000, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw();
}

void capture_attractor(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    ibex::Variable phi, d;
    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
                                    -cos(phi)));

    IntervalVector box(2);
    box[0] = -ibex::Interval::PI | ibex::Interval::PI;
    box[1] = ibex::Interval(0.01, 10.0);

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;
    Scheduler s(box, &f, u);

    ibex::Function f_sym23(phi, d, Return(phi-ibex::Interval::TWO_PI, d));
    s.set_symetry(&f_sym23,1, 1, 1, 0);

    ibex::Function f_sym32(phi, d, Return(phi+ibex::Interval::TWO_PI, d));
    s.set_symetry(&f_sym32,1, 0, 1, 1);

    IntervalVector activated_pave(2);
    activated_pave[0] = ibex::Interval(2);
    activated_pave[1] = ibex::Interval(3.0);

    s.cameleon_cycle(20, 5, 1e9, true, false);
//    s.cameleon_propagation(15, 1e6, activated_pave, false);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw();
}

int main()
{
//    integration();
//    ball();
//    capture_attractor();
    van_der_pol_cycle();
//    lorenz_attractor();
//    test();

#if 0
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(2);
    box[0] = ibex::Interval(-10.0, 10.0);
    box[1] = ibex::Interval(-10.0, 10.0);

    ibex::Interval u = -ibex::Interval::HALF_PI/2 | ibex::Interval::HALF_PI/2;
    Scheduler s(box, &f, u);

    vector<IntervalVector> activated_paves;
    IntervalVector activated_pave(2);
    activated_pave[0] = ibex::Interval(-2.2);
    activated_pave[1] = ibex::Interval(0.0);
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
