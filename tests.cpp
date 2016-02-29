#include <utils.h>

#include <scheduler.h>
#include <vibes.h>
#include "iomanip"
#include "graphdot.h"

using namespace ibex;
using namespace std;

void test_draw(Pave *p, string drawing_name, bool full=false){
    vibes::beginDrawing();
    vibes::newFigure(drawing_name);
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    p->draw(full, "black[]", false, true);
    vibes::setFigureProperties(vibesParams("viewbox", "equal"));
    vibes::axisAuto();
}

void testTranslate(){
    cout << "TEST TRANSLATE" << endl;
    Utils u;

    ibex::IntervalVector Sk(2);
    ibex::IntervalVector box(2);

    box[0] = ibex::Interval(1.0, 10.0);
    box[1] = ibex::Interval(2.0, 10.0);

    Sk[0] = ibex::Interval(5.0, 10.0);
    Sk[1] = ibex::Interval(10.0);

    u.translate_segment_and_box(Sk, box, true, true);

    cout << Sk << endl;
    cout << box << endl;
}

void testRotate(){
    cout << "TEST ROTATE" << endl;
    Utils u;

    ibex::IntervalVector Sk(2);
    ibex::IntervalVector box(2);

    box[0] = ibex::Interval(0.0, 1.0);
    box[1] = ibex::Interval(0.0, 2.0);

    Sk[0] = ibex::Interval(0.0, 0.75);
    Sk[1] = ibex::Interval(0.0);

    u.rotate_segment_and_box(Sk, M_PI, box, true);

    cout << Sk << endl;
//    cout << box << endl;
}

void test_CtcPropagateLeftSide(){
    cout << "TEST CtcPropagateLeftSide" << endl;
    Utils u;

    ibex::Interval x = ibex::Interval(0.0, 1.0);
    ibex::Interval y = ibex::Interval::ALL_REALS;
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0.0, 1.0);
    box[1] = ibex::Interval(0.0, 1.0);

//    ibex::Interval theta = ibex::Interval::PI/4.0 | ibex::Interval::HALF_PI;
    ibex::Interval theta = ibex::Interval::PI | 4*ibex::Interval::PI/5.0;

    u.CtcPropagateLeftSide(x, y, theta, box);

    cout << x << endl;
    cout << y << endl;
}

void test_CtcPropagateRightSide(){
    cout << "TEST test_CtcPropagateRightSide" << endl;
    Utils u;

    ibex::Interval x = ibex::Interval(0.0, 1.0);
    ibex::Interval y = ibex::Interval::ALL_REALS;

    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0.0, 1.0);
    box[1] = ibex::Interval(0.0, 1.0);

    ibex::Interval theta = ibex::Interval::ZERO | ibex::Interval::PI/5.0;

    u.CtcPropagateRightSide(x, y, theta, box);

    cout << x << endl;
    cout << y << endl;
}

void test_CtcPropagateFront(){
    cout << "TEST test_CtcPropagateFront" << endl;
    Utils u;

    ibex::Interval x = ibex::Interval(0.0, 1.0);
    ibex::Interval x_front = ibex::Interval::EMPTY_SET;
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0.0, 1.0);
    box[1] = ibex::Interval(0.0, 1.0);

    ibex::Interval theta = ibex::Interval::ZERO | ibex::Interval::PI/4.0;

    u.CtcPropagateFront(x, x_front, theta, box);

    cout << "x=" << x << endl;
    cout << "x_front=" << x_front << endl;
}

void test_CtcPropagateSegment(){
    cout << "TEST test_CtcPropagateSegment" << endl;
    Utils u;

    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0.0, 1.0);
    box[1] = ibex::Interval(0.0, 1.0);

    int face = 3;
    vector<ibex::Interval> theta = {ibex::Interval::HALF_PI | 5.0*ibex::Interval::HALF_PI/4.0, ibex::Interval::EMPTY_SET};
    ibex::Interval seg_in = ibex::Interval(0,1);
    vector<ibex::Interval> seg_out;
    for(int j=0; j<3; j++){
        seg_out.push_back(ibex::Interval::ALL_REALS);
    }

    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;

    u.CtcPropagateSegment(seg_in, seg_out, face, theta, box, ibex::Interval::EMPTY_SET);

    cout << "----------" << endl;
    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;
}

void test_CtcPaveForward(){
    Utils u;
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(-2.03, -1.955);
    box[1] = ibex::Interval(0.47, 0.545);

//    ibex::Interval command = -ibex::Interval::PI/8 | ibex::Interval::PI/8;
    ibex::Interval command = -ibex::Interval::PI/4 | ibex::Interval::PI/4;
//    ibex::Interval command = -ibex::Interval::PI | ibex::Interval::PI;

    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    Pave p(box, &f, command);

//    p.set_theta(-ibex::Interval::HALF_PI/4.0 | ibex::Interval::HALF_PI/4.0);
//    p.set_theta((ibex::Interval::HALF_PI | 5.0*ibex::Interval::HALF_PI/4.0) + ibex::Interval::PI/3);
//    p.set_theta(-ibex::Interval::HALF_PI | ibex::Interval::HALF_PI);
    p.get_border(0)->set_full_segment_in();
    p.get_border(3)->set_full_segment_in();

//    p.get_border(0)->set_segment_in(ibex::Interval(0.5, 0.9), false);

    test_draw(&p, "test_before");

    u.CtcPaveForward(&p, false, true);

    test_draw(&p, "test_after");
}

void test_CtcPaveConsistency(){
    Utils u;
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(-1.66897, -1.5708);
    box[1] = ibex::Interval(0.0880469, 0.166094);

    ibex::Interval command = ibex::Interval::ZERO;
//    ibex::Interval command = -ibex::Interval::HALF_PI| ibex::Interval::PI;
//    ibex::Interval command = -5*ibex::Interval::HALF_PI/6.0| 5*ibex::Interval::HALF_PI/6.0;
//    ibex::Interval command = -ibex::Interval::PI/4 | ibex::Interval::PI/4;
//    ibex::Interval command = ibex::Interval(-1.0472, 1.0472);

//    Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    ibex::Variable phi, d;
    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
                                    -cos(phi)));
    Pave p(box, &f, command);
    p.set_theta(p.get_theta()[0] + (-ibex::Interval::PI/40.0 | ibex::Interval::PI/40.0) + ibex::Interval::HALF_PI);

//    p.set_theta(ibex::Interval::HALF_PI + ibex::Interval::PI/4);
//    p.set_theta((-ibex::Interval::HALF_PI/16.0 | ibex::Interval::HALF_PI/16.0)+2*ibex::Interval::PI/3);
//    p.set_theta(ibex::Interval(1.5708,2.67795));

//    p.get_border(0)->set_full_segment_in();
//    p.get_border(1)->set_full_segment_in();
    p.get_border(2)->set_full_segment_in();
    p.get_border(3)->set_full_segment_in();
    p.get_border(1)->set_segment_in(ibex::Interval(0.157751, 0.166094), false);
//    p.get_border(1)->set_segment_in(ibex::Interval(-10,0.0), false);

//    p.get_border(0)->set_full_segment_out();
    p.get_border(1)->set_full_segment_out();
    p.get_border(2)->set_full_segment_out();
    p.get_border(3)->set_full_segment_out();
//    p.get_border(2)->set_segment_out(ibex::Interval(0.1, 1.0), false);

//    p.get_border(2)->set_segment_out(ibex::Interval(-10,-5), false);

    test_draw(&p, "test_before");
    u.CtcPaveConsistency(&p, true, false);

    vibes::beginDrawing();
    vibes::newFigure("drawing_name2");
    test_draw(&p, "test_after");
    cout << setprecision(80) << endl;
    p.print();
}

void test_contractor_polar(){
    Utils u;
    ibex::Interval x_ub = ibex::Interval::ALL_REALS;
    ibex::Interval y_ub = ibex::Interval::ZERO;
    ibex::Interval rho_ub = ibex::Interval::ALL_REALS;
    ibex::Interval theta2_ub = ibex::Interval::ALL_REALS;

    cout << x_ub << y_ub << rho_ub << theta2_ub << endl;

    ibex::Interval x_r, y_r;
    x_r = sqrt(2)/2*(x_ub - y_ub);
    y_r = sqrt(2)/2*(x_ub + y_ub);
    theta2_ub += ibex::Interval::PI/4.0;
    ibex::CtcAngle ctcAngle;
    ctcAngle.contract(x_r, y_r, theta2_ub);

    x_ub &= sqrt(2)/2*(x_r + y_r);
    y_ub &= sqrt(2)/2*(-x_r + y_r);
    theta2_ub -= ibex::Interval::PI/4.0;

    cout << x_ub << y_ub << rho_ub << theta2_ub << endl;

//    ibex::Interval x = ibex::Interval::ZERO;
//    ibex::Interval y = ibex::Interval::ALL_REALS;
//    ibex::Interval theta = ibex::Interval::ALL_REALS;

//    const double d2PI   = (2*ibex::Interval::PI).ub();
//    ibex::Interval theta_tmp = atan2(y, x);
//    cout << theta_tmp << endl;
//    bwd_imod(theta, theta_tmp, d2PI);
//    cout << theta << theta_tmp << endl;
//    theta = ibex::Interval::HALF_PI | 3*ibex::Interval::HALF_PI;
//    bwd_angle(theta, y, x);

//    cout << x << y << theta << theta_tmp << endl;
}

void test_rotation(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);

    ibex::IntervalVector Sk(2);
    Sk[0] = ibex::Interval(0);
    Sk[1] = ibex::Interval(0);

    ibex::Interval theta = -ibex::Interval::PI/2.0;

    cout << "Sk=" << Sk << endl;
    cout << "box=" << box << endl;

    Utils u;
    u.rotate_segment_and_box(Sk, theta, box, true);

    cout << "-----" << endl;
    cout << "Sk=" << Sk << endl;
    cout << "box=" << box << endl;
}

void test_diff(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);
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
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);
    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    Utils u;

    Graph g(box, &f, &u, 1);
    g.sivia(0.0, 4, false, false);

    GraphDot graphDot(&g);
    graphDot.write("g.dot");

    Graph g2(&g, 2);

    GraphDot graphDot2(&g2);
    graphDot2.write("g2.dot");

    Graph g3(&g, g.get_pave(1.0, 1.0), 3);

    GraphDot graphDot3(&g3);
    graphDot3.write("g3.dot");

//    g.print();
//    g2.print();
//    g3.print();
}

void sandbox(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(5,10);
    box[1] = ibex::Interval(0,10);

    ibex::Interval test = ibex::Interval(0, 10);
    cout << test.lb() << endl;

//    Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

//    ibex::IntervalVector dposition = f.eval_vector(box);

//    ibex::Interval dx = dposition[0];
//    ibex::Interval dy = dposition[1];

//    ibex::Interval theta = atan2(dy, dx);
//    cout << setprecision(80) << theta << endl;
//    cout << ibex::Interval::HALF_PI << endl;
}
