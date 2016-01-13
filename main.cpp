#include <ibex.h>
#include <vibes.h>

#include <pave.h>
#include <scheduler.h>

#include <tests.h>

#include <iostream>
#include <ctime>

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
//    test_CtcPaveBackward();
//    test_CtcPaveConsistency();

//    test_Newton();

//    test_rotation();

//    test_diff();
    test_copy_graph();
}

void ball(){
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,-10-1*y));

    IntervalVector box(2);
    box[0] = Interval(0.0, 20.0);
    box[1] = Interval(-20.0, 20.0);

    Scheduler s(box, &f);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(12.0);
    activated_pave[1] = Interval(-2.0);

    ibex::Function f_sym(x, y, Return(x, -y));
    s.set_symetry(&f_sym,3);

    s.cameleon_propagation(15, 1000000, activated_pave);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);
}

int main()
{

#if 1
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    Scheduler s(box, &f);

    vector<IntervalVector> activated_paves;
    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(10.0);
    activated_pave[1] = Interval(4.0);
    activated_paves.push_back(activated_pave);

    activated_pave[0] = Interval(0.4);
    activated_pave[1] = Interval(0.4);
    activated_paves.push_back(activated_pave);

    activated_pave[0] = Interval(-0.4);
    activated_pave[1] = Interval(0.4);
    activated_paves.push_back(activated_pave);

    activated_pave[0] = Interval(0.4);
    activated_pave[1] = Interval(-0.4);
    activated_paves.push_back(activated_pave);

    activated_pave[0] = Interval(-0.4);
    activated_pave[1] = Interval(-0.4);
    activated_paves.push_back(activated_pave);

    s.cameleon_propagation(15, 1000000, activated_paves);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
    s.draw(1024, true);

#endif

#if 0
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;


    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
//    ibex::Function f(x, y, Return(y,-sin(x)-y));
//    ibex::Function f(x, y, Return(-y-10*x*x+5*x*y+y*y,x+x*x-25*x*y));
//    ibex::Function f(x, y, Return(cos(2*M_PI*y/10),cos(2*M_PI*x/10)));

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);
    Scheduler s(box, &f);

    s.cameleon_cycle(9, 5, 1000000, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);
//    s.print_pave_info(-5, 6, "r[]");

#endif

#if 0
   test();
#endif
    return 0;
}
