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


#if 0
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box, &f);
    s.SIVIA(s.m_global_pave_list[0], s.m_global_pave_queue[0], 0.0*M_PI/10.0, 5000, false);

    // ****************************************************

    vector<Pave*> pave_list, pave_queue;
    s.copy_graph(pave_list, s.m_global_pave_list[0], true);

    s.activate_pave(pave_list, pave_queue, 3, 4);
    s.process(pave_queue, 100000, false);
    s.draw(pave_list, 500, true);

//    s.activate_pave(s.m_global_pave_list[0], s.m_global_pave_queue[0], 3, 4);
//    s.process(s.m_global_pave_queue[0], 100000, false);
//    s.draw(500, true);

#endif

#if 1
    const clock_t begin_time = clock();
    vibes::beginDrawing();
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box, &f);

    s.process_SIVIA_cycle(12, 4000, 200000);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);

//    s.print_pave_info(-5, 6, "r[]");

#endif

#if 0
   test();
#endif
    return 0;
}
