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

    Border b = Border(Interval(0.0, 0.5), 2);
    p1.queue.push_back(b);

    p1.process();

    p1.draw();
    vibes::axisAuto();

    // ************ TEST BISECT *************

    vibes::newFigure("integralIBEX2");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    vector<Pave*> bisect_pave;

    box[0]=Interval(0,2);
    box[1]=Interval(0,2);
    Pave p2(box);
    p2.bisect(bisect_pave);

    vector<Pave*> bisect_pave1;
    for(int i=0; i<bisect_pave.size(); i++){
        bisect_pave[i]->bisect(bisect_pave1);
    }

    Border b_test = Border(Interval(0.0, 0.5), 2);
    bisect_pave1[1]->queue.push_back(b_test);

    cout << "bisect_pave1[0]" << bisect_pave1[0]->box << endl;

    for(int i=0; i<bisect_pave1.size(); i++){
        bisect_pave1[i]->theta = (-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0) - Interval::HALF_PI;
    }

    for(int i=0; i<4; i++){
        for(int i=0; i<bisect_pave1.size(); i++){
            bisect_pave1[i]->process();
        }
    }

    for(int i=0; i<bisect_pave1.size(); i++){
        (bisect_pave1[i])->draw();
    }

    vibes::axisAuto();

    return 0;
}

