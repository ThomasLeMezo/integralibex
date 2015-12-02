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

    vibes::beginDrawing();
    vibes::newFigure("integralIBEX2");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    vibes::axisAuto();

    Scheduler s;

    IntervalVector box(2);
    box[0] = Interval(-10.0, 10.0);
    box[1] = Interval(-10.0, 10.0);

    s.set_initial_pave(box);

    s.SIVIA(M_PI/10.0, 10000);

    s.add_segment();
    s.process(100000);
    s.draw();

    return 0;
}

