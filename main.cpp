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

    test_diff();
}

int main()
{

#if 1
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    ibex::Function f(x, y, Return(y,-10-1*y));

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(0.0, 20.0);
    box[1] = Interval(-20.0, 20.0);
    s.set_initial_pave(box, &f);

    IntervalVector activated_pave(2);
    activated_pave[0] = Interval(15,15);
    activated_pave[1] = Interval(0.0,0.0);

    s.cameleon_propagation(20, 1000000, activated_pave, 50);

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

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box, &f);

    s.cameleon_cycle(10, 5, 1000000, true);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);

//    s.print_pave_info(-5, 6, "r[]");

#endif

#if 0
   test();
#endif
    return 0;
}
