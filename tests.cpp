#include <utils.h>

#include <scheduler.h>
#include <vibes.h>

using namespace ibex;
using namespace std;

void test_draw(Pave *p, string drawing_name, bool full=false){
    vibes::beginDrawing();
    vibes::newFigure(drawing_name);
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    p->draw(full);
    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}

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
    Interval x_front = Interval::EMPTY_SET;
    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    Interval theta = Interval::ZERO | Interval::PI/4.0;

    u.CtcPropagateFront(x, x_front, theta, box);

    cout << "x=" << x << endl;
    cout << "x_front=" << x_front << endl;
}

void test_CtcPropagateSegment(){
    cout << "TEST test_CtcPropagateSegment" << endl;
    Utils u;

    IntervalVector box(2);
    box[0] = Interval(0.0, 1.0);
    box[1] = Interval(0.0, 1.0);

    int face = 3;
    vector<Interval> theta = {Interval::HALF_PI | 5.0*Interval::HALF_PI/4.0, Interval::EMPTY_SET};
    Interval seg_in = Interval(0,1);
    vector<Interval> seg_out;
    for(int j=0; j<3; j++){
        seg_out.push_back(Interval::ALL_REALS);
    }

    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;

    u.CtcPropagateSegment(seg_in, seg_out, face, theta, box);

    cout << "----------" << endl;
    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;
}

void test_CtcPaveForward(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(0, 1);
    box[1] = Interval(0, 1);

    Function f;
    Pave p(box, &f);

//    p.set_theta(-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0);
    p.set_theta(Interval::HALF_PI | 5.0*Interval::HALF_PI/4.0);

//    p.m_borders[0].set_full();
//    p.m_borders[2].segment = Interval(0,1);
//    p.m_borders[1].segment = Interval(0,1);

    test_draw(&p, "test_before");

    u.CtcPaveForward(&p, false);

    test_draw(&p, "test_after");
}

void test_CtcPaveConsistency(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(0, 1);
    box[1] = Interval(0, 1);

    Function f;
    Pave p(box, &f);

    p.set_theta(3.0*Interval::HALF_PI/4.0 | 5.0*Interval::HALF_PI/4.0);

    p.get_border(1)->set_full();
    p.get_border(2)->set_full();
    p.get_border(3)->set_full();

    test_draw(&p, "test_before");

    u.CtcPaveConsistency(&p, true);

    test_draw(&p, "test_after");
}

void test_rotation(){
    IntervalVector box(2);
    box[0] = Interval(0,1);
    box[1] = Interval(0,1);

    IntervalVector Sk(2);
    Sk[0] = Interval(0);
    Sk[1] = Interval(0);

    Interval theta = -Interval::PI/2.0;

    cout << "Sk=" << Sk << endl;
    cout << "box=" << box << endl;

    Utils u;
    u.rotate_segment_and_box(Sk, theta, box, true);

    cout << "-----" << endl;
    cout << "Sk=" << Sk << endl;
    cout << "box=" << box << endl;
}

void test_diff(){
    IntervalVector box(2);
    box[0] = Interval(0,1);
    box[1] = Interval(0,1);
    Function f;

    Pave p1(box, &f);
    Pave p2(box, &f);

    p1.get_border(0)->set_full();
    p1.get_border(1)->set_full();

    p2.set_full();

    test_draw(&p1, "p1", true);
    test_draw(&p2, "p2_before", true);

    p2.diff(p1);
    test_draw(&p2, "p2_after", true);

}

void test_copy_graph(){
    IntervalVector box(2);
    box[0] = Interval(0,1);
    box[1] = Interval(0,1);
    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    Utils u;

    Graph g(box, &f, &u, 1);
    g.sivia(0.0, 4, false, false);

    Graph g2(&g, 2);
    Graph g3(&g, g.get_pave(1.0, 1.0), 3);

    g.print();
    g2.print();
    g3.print();
}
