#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

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
    vibes::newFigure("integralIBEX2");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));


    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box);

    s.SIVIA(M_PI/10.0, 1000);

    s.add_segment(-1.28, 4.0);
    s.process(100);
    s.draw();

    vibes::axisAuto();
    vibes::setFigureProperties(vibesParams("viewbox", "equal"));

    cout << "Nb of paves = " << s.pave_list.size() << endl;

    s.print_pave_info(-0.7,3.91);
#endif

    return 0;
}

