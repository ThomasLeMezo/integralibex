#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

using namespace std;
using namespace ibex;

int main()
{
//    vibes::beginDrawing();
//    vibes::newFigure("integralIBEX");
//    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
//    vibes::axisAuto();

//    vibes::drawSector(0.0, 0.0, 1.0, 1.0, 359.0, 0.0, "r[]");

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
    vibes::newFigure("integralIBEX2");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    vibes::axisAuto();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box);

    s.SIVIA(M_PI/10.0, 20000);

    cout << s.pave_list.size() << endl;
    s.add_segment(12833);
    s.process(100000);
    s.draw();

    cout << s.pave_list[0]->box << endl;
#endif

    return 0;
}

