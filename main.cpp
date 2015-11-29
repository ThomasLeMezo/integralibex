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
    box[1]=Interval(0,2);
    Pave p1(box);
    p1.theta = Interval(M_PI/3.0, M_PI/2.0);

//    IntervalVector box2(2);
//    box2[0] = Interval(0.5, 0.75);
//    box2[1] = Interval(0.0);

//    Border b(box2, 0, NULL);
    p1.borders[3].segments.push_back(Interval(0.5, 0.75));

//    p1.queue.push_back(segment);
//    p1.process();

    p1.draw();

    return 0;
}

