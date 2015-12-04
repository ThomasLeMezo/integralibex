#include <utils.h>

#include <pave.h>
#include <vibes.h>

using namespace ibex;
using namespace std;

void testTranslate(){
    cout << "TEST TRANSLATE" << endl;
    Utils u;

    IntervalVector Sk(2);
    IntervalVector box(2);

    box[0] = Interval(1.0, 10.0);
    box[1] = Interval(2.0, 10.0);

    Sk[0] = Interval(5.0, 10.0);
    Sk[1] = Interval(10.0);

    u.translate_segment_and_box(Sk, box, true, true);

    cout << Sk << endl;
    cout << box << endl;
}

void testRotate(){
    cout << "TEST ROTATE" << endl;
    Utils u;

    IntervalVector Sk(2);
    IntervalVector box(2);

    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 2.0);

    Sk[0] = Interval(0.0, 0.75);
    Sk[1] = Interval(0.0);

    u.rotate_segment_and_box(Sk, M_PI, box);

    cout << Sk << endl;
//    cout << box << endl;
}

void test_CtcPropagateLeftSide(){
    cout << "TEST CtcPropagateLeftSide" << endl;
    Utils u;

    Interval Sk = Interval(0.0, 1.0);
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

//    Interval theta = Interval::PI/4.0 | Interval::HALF_PI;
    Interval theta = Interval::PI | 4*Interval::PI/5.0;

    u.CtcPropagateLeftSide(Sk, theta, box);

    cout << Sk << endl;
}

void test_CtcPropagateRightSide(){
    cout << "TEST test_CtcPropagateRightSide" << endl;
    Utils u;

    Interval Sk = Interval(0.0, 1.0);
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Interval theta = Interval::ZERO | Interval::PI/5.0;

    u.CtcPropagateRightSide(Sk, theta, box);

    cout << Sk << endl;
}

void test_CtcPropagateFront(){
    cout << "TEST test_CtcPropagateFront" << endl;
    Utils u;

    Interval Sk = Interval(0.0, 1.0);
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Interval theta = Interval::ZERO | Interval::PI/4.0;

    cout << Sk << endl;
    u.CtcPropagateFront(Sk, theta, box);

    cout << Sk << endl;
}

void test_Propagate(){
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Pave p(box, NULL);

    Border b(Interval(0.25, 0.75), 1);
    p.add_new_segment(b);
    p.set_theta( (Interval::ZERO | Interval::PI));
    p.process();

    vibes::beginDrawing();
    vibes::newFigure("test");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));

    p.draw();

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}
