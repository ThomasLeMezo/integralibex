#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

using namespace std;
using namespace ibex;

int main()
{
    vibes::beginDrawing();
    vibes::newFigure("integralIBEX");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));

    IntervalVector box(2);
    box[0]=Interval(0,1);
    box[1]=Interval(0,1);
    Pave p1(box);
//    p1.theta = (-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0);
//    p1.theta = (-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0) + Interval::HALF_PI;
//    p1.theta = (-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0) + Interval::PI;
    p1.theta = (-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0) - Interval::HALF_PI;


    int face = 2;
    Border b = Border(Interval(0.5, 0.75), face);
    p1.borders[face].segments.push_back(Interval(0.5, 0.75));
    p1.queue.push_back(b);

    p1.process();

    p1.draw();
    vibes::axisAuto();

    return 0;
}

