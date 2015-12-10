#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

#include <tests.h>

using namespace std;
using namespace ibex;

void test(){
    //    testTranslate();
    //    testRotate();
    //    test_CtcPropagateLeftSide();
    //    test_CtcPropagateRightSide();
    //    test_CtcPropagateFront();
        test_Propagate();
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
    s.SIVIA(M_PI/10.0, 3000);

    s.add_segment(-1.28, 4.0);
//    s.add_segment(1.78, -6.42);
//    s.add_segment(0.0, 0.5);

    s.process(200000);
    s.draw(500);

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
    s.draw(500);

    cout << "Nb of paves = " << s.pave_list.size() << endl;
#endif
    return 0;
}
