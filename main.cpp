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
    test_CtcPaveConsistency();

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

    s.set_initial_pave(box);
    s.SIVIA(M_PI/10.0, 3000, false);

    s.add_segment(-1.28, 4.0);
//    s.add_segment(1.78, -6.42);
//    s.add_segment(0.0, 0.5);

    s.process(200000);
    s.draw(500, true);

//    s.print_pave_info(2.0,0.1, "r[]");
//    s.print_pave_info(1.9,-0.1, "g[]");
//    s.print_pave_info(2.1,-0.1, "y[]");
//    s.print_pave_info(-0.56,-1.96, "y[]");
//    s.print_pave_info(0.2,3.0, "b[]");
#endif

#if 0
    vibes::beginDrawing();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.process_graph(1000, 1000);
    s.draw(500, true);

    cout << "Nb of paves = " << s.pave_list.size() << endl;
#endif

#if 0
    vibes::beginDrawing();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box);
    s.SIVIA(M_PI/10.0, 5000, false);
    s.set_full_continuity();
    s.process_backward(100000);
    s.draw(1024, true);

    s.print_pave_info(-5,-5, "b[]");

    cout << "Nb of paves = " << s.pave_list.size() << endl;
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

    s.process_SIVIA_cycle(15, 4000, 100000);

    cout << "TIME = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

    s.draw(1024, true);

//    s.print_pave_info(-5, 6);

#endif

#if 0
   test();
#endif
    return 0;
}
