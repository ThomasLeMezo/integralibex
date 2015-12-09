#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

#include <tests.h>

using namespace std;
using namespace ibex;

int main()
{
#if 0
    vibes::beginDrawing();
    vibes::newFigure("integralIBEX");
    vibes::drawSector(0.0, 0.0, 1.0, 1.0, 0.0, 359.0, "r[]");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500, "viewbox", "equal"));
    vibes::axisAuto();
#endif

    // *************************
#if 0
    Interval dx(-1.0);
    Interval dy(-0.1,0.1);
    Interval theta = atan2(dy, dx);
    cout << theta << endl;

    CtcPolar p;
    Interval rho = Interval::POS_REALS;
    Interval theta_p = Interval::ALL_REALS;
    p.contract(dx, dy, rho, theta_p);

    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "rho = " << rho << endl;
    cout << "theta_p = " << theta_p << endl;
#endif

#if 1
    vibes::beginDrawing();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box);
    //s.SIVIA(M_PI/10.0, 6);

//    s.add_segment(-1.28, 4.0);
//    s.add_segment(1.78, -6.42);
//    s.add_segment(0.0, 0.5);

    //s.process(200000);
    //s.draw(500);
    s.process_graph(1000, 1000);
    s.draw(500);

    cout << "Nb of paves = " << s.pave_list.size() << endl;

//    s.print_pave_info(2.0,0.1, "r[]");
//    s.print_pave_info(1.9,-0.1, "g[]");
//    s.print_pave_info(2.1,-0.1, "y[]");
//    s.print_pave_info(-0.56,-1.96, "y[]");
//    s.print_pave_info(0.2,3.0, "b[]");
#else

//    testTranslate();
//    testRotate();
//    test_CtcPropagateLeftSide();
//    test_CtcPropagateRightSide();
//    test_CtcPropagateFront();
    test_Propagate();
#endif
    return 0;
}
