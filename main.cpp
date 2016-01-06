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

    test_CtcPaveForward();
//    test_CtcPaveBackward();
//    test_CtcPaveConsistency();

//    test_Newton();

//    test_rotation();
}

int main()
{

#if 0
    vibes::beginDrawing();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    ibex::Function f;
    s.set_initial_pave(box, &f);
    s.SIVIA(s.m_global_pave_list[0], s.m_global_pave_queue[0], M_PI/10.0, 3000, false);

    s.activate_pave(s.m_global_pave_list[0], s.m_global_pave_queue[0], -4.0, 6.0);

    s.process(s.m_global_pave_queue[0], 100000, false);
    s.draw(500, true);

//    s.print_pave_info(2.0,0.1, "r[]");
#endif

#if 1
    const clock_t begin_time = clock();
    vibes::beginDrawing();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    ibex::Function f;
    s.set_initial_pave(box, &f);

    s.process_SIVIA_cycle(1, 4000, 100000);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);

//    s.print_pave_info(-5, 6, "r[]");

#endif

#if 0
   test();
#endif
    return 0;
}
