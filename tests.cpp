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

    u.rotate_segment_and_box(Sk, M_PI, box, true);

    cout << Sk << endl;
//    cout << box << endl;
}

void test_CtcPropagateLeftSide(){
    cout << "TEST CtcPropagateLeftSide" << endl;
    Utils u;

    Interval x = Interval(0.0, 1.0);
    Interval y = Interval::ALL_REALS;
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

//    Interval theta = Interval::PI/4.0 | Interval::HALF_PI;
    Interval theta = Interval::PI | 4*Interval::PI/5.0;

    u.CtcPropagateLeftSide(x, y, theta, box);

    cout << x << endl;
    cout << y << endl;
}

void test_CtcPropagateRightSide(){
    cout << "TEST test_CtcPropagateRightSide" << endl;
    Utils u;

    Interval x = Interval(0.0, 1.0);
    Interval y = Interval::ALL_REALS;

    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Interval theta = Interval::ZERO | Interval::PI/5.0;

    u.CtcPropagateRightSide(x, y, theta, box);

    cout << x << endl;
    cout << y << endl;
}

void test_CtcPropagateFront(){
    cout << "TEST test_CtcPropagateFront" << endl;
    Utils u;

    Interval x = Interval(0.0, 1.0);
    Interval x_front = Interval::ALL_REALS;
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Interval theta = Interval::ZERO | Interval::PI/4.0;

    u.CtcPropagateFront(x, x_front, theta, box);

    cout << x << endl;
    cout << x_front << endl;
}

void test_Propagate(){
    IntervalVector box(2);
    box[0] = Interval(-10.0, 0.0);
    box[1] = Interval(0.0, 10.0);
    Scheduler s;
    Pave p(box, &s);

//    Border b(Interval(0.5, 1.0), 0);
    Border b(p.get_border_position(1), 1, Interval(0.0, 10.0));
    Border b2(p.get_border_position(0), 0, Interval(-10.0, 0.0));
//    Border b(Interval(1.0, 3.0), 3);

    p.add_new_segment(b, true);
    p.borders[0].segment = Interval(-10.0, 0.0);
//    p.set_theta( (-3.88*Interval::PI/8.0 | -Interval::HALF_PI));
//    p.set_theta( (6*Interval::PI/8.0 | 7*Interval::PI/8.0));
//    p.set_theta( (6*Interval::PI/8.0 | 7*Interval::PI/8.0) + Interval::HALF_PI);
//    p.set_theta( (Interval::PI/8.0 | 2*Interval::PI/8.0) - 2* Interval::PI/3.0);
//    p.set_theta(-2*Interval::PI/3.0 | -3*Interval::PI/4.0);
      p.set_theta( -Interval::HALF_PI | Interval::HALF_PI);

    p.process_forward();
//    p.process_forward();

    vibes::beginDrawing();
    vibes::newFigure("test");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));

    p.draw(false);

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}

void test_Backward(){
    IntervalVector box(2);
    box[0] = Interval(0.5, 1.5);
    box[1] = Interval(1.0, 3.0);
    Scheduler s;
    Pave p(box, &s);

    for(int i=0;i<4; i++){
        p.borders[i].set_full();
    }

    Border b(p.get_border_position(0), 0, Interval(0.5, 1.0));
//    Border b(Interval(1.0, 2.0), 1);
//    Border b2(Interval(1.0, 1.25), 2);
//    Border b(Interval(1.0, 3.0), 3);

//    p.add_new_segment(b2);
    p.add_new_segment(b, false);
//    p.set_theta( (-3.88*Interval::PI/8.0 | -Interval::HALF_PI));
//    p.set_theta( (6*Interval::PI/8.0 | 7*Interval::PI/8.0));
//    p.set_theta( (6*Interval::PI/8.0 | 7*Interval::PI/8.0) + Interval::HALF_PI);
    p.set_theta( (-Interval::HALF_PI/8.0 | Interval::HALF_PI/8.0) - 0*Interval::PI/2);
//    p.set_theta(-2*Interval::PI/3.0 | -3*Interval::PI/4.0);

    p.compute_flow();
    p.process_backward();

    vibes::beginDrawing();
    vibes::newFigure("test");
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));

    p.draw(false);

    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}
