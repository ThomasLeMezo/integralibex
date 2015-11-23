#include <iostream>
#include <ibex.h>
#include <pave.h>
#include <vibes.h>

using namespace std;
using namespace ibex;

int main()
{
    ibex::Interval i(10.0);
    cout << "Hello World! : " << i << endl;


    vibes::beginDrawing();
    vibes::newFigure("integralIBEX");
    vibes::setFigureProperties(vibesParams("x",0,"y",40,"width",150,"height",150));

    IntervalVector box(2);
    box[0]=Interval(0,1);
    box[1]=Interval(0,2);
    Pave p1(box);
    p1.draw();

    return 0;
}

