#include <utils.h>

#include <scheduler.h>
#include <vibes.h>
#include "iomanip"
#include "graphdot.h"

using namespace ibex;
using namespace std;
//using namespace cv;

void test_draw(Pave *p, string drawing_name="test", bool full=true){
    vibes::beginDrawing();
    vibes::newFigure(drawing_name);
    vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
    p->draw(full, false);
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
    vector<Interval> theta_list; theta_list.push_back(theta);

    u.CtcPropagateLeftSide(x, y, theta_list, box);

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
    vector<Interval> theta_list; theta_list.push_back(theta);

    u.CtcPropagateRightSide(x, y, theta_list, box);

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
    vector<Interval> theta_list; theta_list.push_back(theta);

    u.CtcPropagateFront(x, x_front, theta_list, box);

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
    box[0] = Interval(-2.03, -1.955);
    box[1] = Interval(0.47, 0.545);

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p(box, f_list);

    //    p.set_theta(-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0);
    //    p.set_theta((Interval::HALF_PI | 5.0*Interval::HALF_PI/4.0) + Interval::PI/3);
    //    p.set_theta(-Interval::HALF_PI | Interval::HALF_PI);
    p.get_border(0)->set_full_segment_in();
    p.get_border(3)->set_full_segment_in();

    //    p.get_border(0)->set_segment_in(Interval(0.5, 0.9), false);

    test_draw(&p, "test_before");

std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcPaveForward(&p, false, change_tab);

    test_draw(&p, "test_after");
}

void test_CtcConsistency(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(3, 3.05);
    box[1] = Interval(0.1, 0.190625);


    //    Interval command = Interval::ZERO;
    ////    Interval command = -Interval::HALF_PI| Interval::PI;
    ////    Interval command = -5*Interval::HALF_PI/6.0| 5*Interval::HALF_PI/6.0;
    ////    Interval command = -Interval::PI/4 | Interval::PI/4;
    ////    Interval command = Interval(-1.0472, 1.0472);

    //    Variable x, y;
    //    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    //    Variable phi, d;
    //    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
    //                                    -cos(phi)));

    Variable x1, x2;
    ibex::Function f(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p(box, f_list);
    //    p.set_theta(p.get_theta()[0] + (-Interval::PI/40.0 | Interval::PI/40.0) + Interval::HALF_PI);

    //    p.set_theta(Interval::HALF_PI + Interval::PI/4);
    //    p.set_theta((-Interval::HALF_PI/16.0 | Interval::HALF_PI/16.0)+2*Interval::PI/3);
    //    p.set_theta(Interval(1.5708,2.67795));

    //    BOX = ([3, 3.05] ; [0.1, 0.190625])
    //    0xbf92b0
    //    Border ID	Position ([x], [y])	segment_in	segment_out
    //    border 0	position=([3, 3.05] ; [0.1, 0.1])    	in=[ empty ]	out=[ empty ]continuity = 1
    //    border 1	position=([3.05, 3.05] ; [0.1, 0.190625])    	in=[ empty ]	out=[0.1, 0.190625]continuity = 1
    //    border 2	position=([3, 3.05] ; [0.190625, 0.190625])    	in=[3.00253, 3.05]	out=[ empty ]continuity = 1
    //    border 3	position=([3, 3] ; [0.1, 0.190625])    	in=[0.1, 0.190625]	out=[ empty ]continuity = 1


    //    p.get_border(0)->set_full_segment_in();

    p.get_border(2)->set_segment_in(Interval(3.00253, 3.05), false);
    p.get_border(3)->set_segment_in(Interval(0.1, 0.190625), false);

    //    p.get_border(0)->set_full_segment_out();
    p.get_border(1)->set_segment_out(Interval(0.1, 0.190625), false);

    test_draw(&p, "test_before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    test_draw(&p, "test_after");
    cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency2(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(4.2, 4.85);
    box[1] = Interval(-1.09375, -0.1);

    Variable x1, x2, x, y;
    //    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    ibex::Function f1(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));

    ibex::Function f2(x1, x2, Return(x2,
                                     -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));


    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list);

    //    p.get_border(0)->set_segment_in(Interval(-3, -2), false);
    //    p.get_border(0)->set_segment_out(Interval(-0.277939, -0.2), false);
    p.get_border(0)->set_full();
    //    p.get_border(0)->set_full_segment_out();

    //    p.get_border(1)->set_segment_in(Interval(-4, -3.75), false);
    //    p.get_border(1)->set_segment_out(Interval(0.459375, 0.5125), false);
    p.get_border(1)->set_full();
    //    p.get_border(1)->set_full_segment_out();

    p.get_border(2)->set_segment_in(Interval(4.2, 4.25), false);
    p.get_border(2)->set_segment_out(Interval(4.2, 4.25), false);
    //    p.get_border(2)->set_full_segment_out();
//    p.get_border(2)->set_full();
    //    p.get_border(2)->set_full_segment_out();

    //    p.get_border(3)->set_segment_in(Interval(0.480101, 0.5125), false);
    //    p.get_border(3)->set_segment_out(Interval(-3, -2), false);
    //    p.get_border(3)->set_full_segment_out();
    p.get_border(3)->set_full();

//    test_draw(&p, "test_before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    test_draw(&p, "test_after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency3(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(0, 1);
    box[1] = Interval(0,1);

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p(box, f_list, false);

    vector<Interval> theta_list;
    theta_list.push_back(Interval::ZERO | Interval::HALF_PI/2.0);
    theta_list.push_back(Interval::ZERO | -Interval::HALF_PI/2.0);

    p.set_theta(theta_list);
    //    p.set_theta(Interval::HALF_PI | Interval::PI);
    //    p.set_theta(Interval::ZERO | Interval::HALF_PI);
    //    p.set_theta(Interval::ZERO | -Interval::HALF_PI);

    p.get_border(0)->set_full_segment_out();
    p.get_border(0)->set_full_segment_in();

    p.get_border(1)->set_full_segment_out();
    p.get_border(1)->set_full_segment_in();

    p.get_border(2)->set_full_segment_out();
    p.get_border(2)->set_full_segment_in();

    p.get_border(3)->set_full_segment_out();
    p.get_border(3)->set_full_segment_in();

    test_draw(&p, "test_before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    test_draw(&p, "test_after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency_Kernel(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(2.0625, 2.5);
    box[1] = Interval(2.5, 3);

//    BOX = ([5.83594, 5.89062] ; [-0.5625, -0.53125])

    Variable x1, x2;
    //    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    ibex::Function f1(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list);

    //    nb	in_inner	out_inner	in_outer	out_outer
//    0	[2.0625, 2.35504]	[2.0625, 2.35504]	[2.0625, 2.5]	[2.0625, 2.5]
//    1	[2.71503, 3]	[2.71503, 3]	[2.5, 3]	[2.5, 3]
//    2	[2.0625, 2.5]	[2.0625, 2.5]	[2.0625, 2.5]	[2.0625, 2.5]
//    3	[2.5, 3]	[2.5, 3]	[2.5, 3]	[2.5, 3]

//    p.set_inner_mode(false);
//    p.get_border(0)->set_segment_in(Interval(5.83594, 5.89062), false);
//    p.get_border(1)->set_segment_in(Interval(-0.5625, -0.53125), false);
//    p.get_border(2)->set_segment_in(Interval(5.83594, 5.89062), false);

//    p.get_border(0)->set_segment_out(Interval(5.83594, 5.89062), false);
//    p.get_border(2)->set_segment_out(Interval(5.83594, 5.89062), false);
//    p.get_border(3)->set_segment_out(Interval(-0.5625, -0.53125), false);

    p.set_inner_mode(true);
    p.get_border(0)->set_segment_in(Interval(2.0625, 2.35504), false);
    p.get_border(1)->set_segment_in(Interval(2.71503, 3), false);
    p.get_border(2)->set_segment_in(Interval(2.0625, 2.5), false);
    p.get_border(3)->set_segment_in(Interval(2.5, 3), false);

    p.get_border(0)->set_segment_out(Interval(2.0625, 2.35504), false);
    p.get_border(1)->set_segment_out(Interval(2.71503, 3), false);
    p.get_border(2)->set_segment_out(Interval(2.0625, 2.5), false);
    p.get_border(3)->set_segment_out(Interval(2.5, 3), false);

    Graph g(&u, 0);
    g.compute_propagation_zone(&p);

    p.draw_test(512, "before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    p.draw_test(512, "after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency_Kernel2(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(9.5, 11.25);
    box[1] = Interval(-5, -4);

    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list);

    //nb	in_inner	out_inner	in_outer	out_outer
//    0	[9.5, 11.25]	[9.5, 11.25]	[9.5, 11.25]	[9.5, 11.25]
//    1	[-5, -4]	[-5, -4]	[-5, -4]	[ empty ]
//    2	[9.5, 11.25]	[9.5, 11.25]	[9.5, 11.25]	[9.5, 11.25]
//    3	[-5, -4]	[-5, -4]	[ empty ]	[-5, -4]
    p.set_inner_mode(true);
    p.get_border(0)->set_segment_in(Interval(9.5, 11.25), false);
    p.get_border(1)->set_segment_in(Interval(-5, -4), false);
    p.get_border(2)->set_segment_in(Interval(9.5, 11.25), false);
    p.get_border(3)->set_segment_in(Interval(-5, -4), false);

    p.get_border(0)->set_segment_out(Interval(9.5, 11.25), false);
    p.get_border(1)->set_segment_out(Interval(-5, -4), false);
    p.get_border(2)->set_segment_out(Interval(9.5, 11.25), false);
    p.get_border(3)->set_segment_out(Interval(-5, -4), false);

    Graph g(&u, 0);
    g.compute_propagation_zone(&p);

    p.draw_test(512, "before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    p.draw_test(512, "after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency_Kernel3(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(12.125, 13);
    box[1] = Interval(0.5, 1);

    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list);

    //nb	in_inner	out_inner	in_outer	out_outer
    //0	[12.125, 13]	[12.125, 13]	[ empty ]	[12.125, 13]
    //1	[0.5, 1]	[0.5, 1]	[ empty ]	[ empty ]
    //2	[12.125, 13]	[12.125, 13]	[12.125, 12.5339]	[ empty ]
    //3	[0.5, 1]	[0.5, 1]	[0.5, 1]	[ empty ]
    p.set_inner_mode(true);
    p.get_border(0)->set_segment_in(Interval(12.125, 13), false);
    p.get_border(1)->set_segment_in(Interval(0.5, 1), false);
    p.get_border(2)->set_segment_in(Interval(12.125, 13), false);
    p.get_border(3)->set_segment_in(Interval(0.5, 1), false);

//    p.get_border(0)->set_segment_out(Interval(12.125, 13), false);
    p.get_border(1)->set_segment_out(Interval(0.5, 1), false);
    p.get_border(2)->set_segment_out(Interval(12.125, 13), false);
//    p.get_border(3)->set_segment_out(Interval(0.5, 1), false);

    Graph g(&u, 0);
    g.compute_propagation_zone(&p);

    p.draw_test(512, "before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    p.draw_test(512, "after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcConsistency_Kernel4(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(2.9375, 3.375);
    box[1] = Interval(-3.5, -3);

    Variable x1, x2;
    ibex::Function f1(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));
    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (-1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list);

    //nb	in_inner	out_inner	in_outer	out_outer
//    0	[ empty ]	[2.9375, 3.375]	[2.9375, 3.375]	[2.9375, 3.375]
//    1	[-3.5, -3.00493]	[ empty ]	[-3.5, -3]	[ empty ]
//    2	[ empty ]	[ empty ]	[2.9375, 3.375]	[2.9375, 3.375]
//    3	[ empty ]	[ empty ]	[ empty ]	[-3.5, -3]
    p.set_inner_mode(true);
//    p.get_border(0)->set_segment_in(Interval(12.125, 13), false);
    p.get_border(1)->set_segment_in(Interval(-3.5, -3.00493), false);
//    p.get_border(2)->set_segment_in(Interval(12.125, 13), false);
//    p.get_border(3)->set_segment_in(Interval(0.5, 1), false);

    p.get_border(0)->set_segment_out(Interval(2.9375, 3.375), false);
//    p.get_border(1)->set_segment_out(Interval(0.5, 1), false);
//    p.get_border(2)->set_segment_out(Interval(12.125, 13), false);
//    p.get_border(3)->set_segment_out(Interval(0.5, 1), false);

    Graph g(&u, 0);
    g.compute_propagation_zone(&p);

    p.draw_test(512, "before");
    std:vector<bool> change_tab;
    for(int i=0; i<4; i++)
        change_tab.push_back(false);
    u.CtcConsistency(&p, true, change_tab);

    p.draw_test(512, "after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_contractor_polar(){
    Utils u;
    Interval x_ub = Interval::ALL_REALS;
    Interval y_ub = Interval::ZERO;
    Interval rho_ub = Interval::ALL_REALS;
    Interval theta2_ub = Interval::ALL_REALS;

    cout << x_ub << y_ub << rho_ub << theta2_ub << endl;

    Interval x_r, y_r;
    x_r = sqrt(2)/2*(x_ub - y_ub);
    y_r = sqrt(2)/2*(x_ub + y_ub);
    theta2_ub += Interval::PI/4.0;
    ibex::CtcAngle ctcAngle;
    ctcAngle.contract(x_r, y_r, theta2_ub);

    x_ub &= sqrt(2)/2*(x_r + y_r);
    y_ub &= sqrt(2)/2*(-x_r + y_r);
    theta2_ub -= Interval::PI/4.0;

    cout << x_ub << y_ub << rho_ub << theta2_ub << endl;

    //    Interval x = Interval::ZERO;
    //    Interval y = Interval::ALL_REALS;
    //    Interval theta = Interval::ALL_REALS;

    //    const double d2PI   = (2*Interval::PI).ub();
    //    Interval theta_tmp = atan2(y, x);
    //    cout << theta_tmp << endl;
    //    bwd_imod(theta, theta_tmp, d2PI);
    //    cout << theta << theta_tmp << endl;
    //    theta = Interval::HALF_PI | 3*Interval::HALF_PI;
    //    bwd_angle(theta, y, x);

    //    cout << x << y << theta << theta_tmp << endl;
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

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p1(box, f_list);
    Pave p2(box, f_list);

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
    Utils utils;

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Graph g(box, f_list, &utils, 1);
    g.sivia(4, GRAPH_FORWARD, false);

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

void test_imageIntegral(){
//    int sizeX = 1000;
//    int sizeY = 1000;
//    IntervalVector range(2);
//    range[0] = -Interval::PI | Interval::PI;
//    range[1] = Interval(0.01,10.0);

//    Variable t;
//    ibex::Function f(t, Return(2*atan(tan((atan2(cos(t), -sin(t))+Interval::PI-atan2(sin(t), cos(t)+1.0/sqrt(2.0)))/2.0)),
//                               sqrt(pow(cos(t)+1/sqrt(2.0), 2)+pow(sin(t), 2))));

//    std::vector<Interval> list_t;
//    list_t.push_back(Interval::ZERO | Interval::TWO_PI);

//    for(int i=0; i<10; i++){
//        std::vector<Interval> tmp_list_t(list_t);
//        list_t.clear();
//        for(Interval &i:tmp_list_t){
//            list_t.push_back(Interval(i.lb(), i.mid()));
//            list_t.push_back(Interval(i.mid(), i.ub()));
//        }
//    }

//    double factorX = sizeX/range[0].diam();
//    double factorY = sizeY/range[1].diam();

//    Mat img(sizeX,sizeY, CV_8U, Scalar(0));
//    for(Interval &i:list_t){
//        IntervalVector tmp(2);
//        tmp[0] = i;
//        IntervalVector p(f.eval_vector(tmp));
//        if(p.diam()[0] < 1)
//            rectangle(img, Point(ceil((p[0].lb()-range[0].lb())*factorX), ceil((p[1].lb()-range[1].lb())*factorY)), Point(floor((p[0].ub()-range[0].lb())*factorX), floor((p[1].ub()-range[1].lb())*factorY)), Scalar(255), 5.0);
//    }

//    vector<vector<Point> > contours;
//    vector<Vec4i> hierarchy;

//    /// Find contours
//    findContours( img, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

//    //    /// Draw contours
//    Mat img_filled = Mat::zeros( img.size(), CV_8U);
//    for( int i = 0; i< contours.size(); i++ ){
//        drawContours(img_filled, contours, i, Scalar(1), CV_FILLED, 4, hierarchy, 0, Point() );
//    }

//    Mat img_integral = Mat::zeros(img.size(), CV_32SC1); // CV_32SC1
//    cv::integral(img_filled, img_integral, CV_32SC1);


//    IntervalVector test_box(2);
//    test_box[0] = Interval(-1.7,-1.6);
//    test_box[1] = Interval(0.9,1.2);

//    if(test_box.is_interior_subset(range)){

//        int x_min, x_max, y_min, y_max;
//        x_min = floor((test_box[0].lb()-range[0].lb())*factorX);
//        x_max = ceil((test_box[0].ub()-range[0].lb())*factorX);
//        y_min = floor((test_box[1].lb()-range[1].lb())*factorY);
//        y_max = ceil((test_box[1].ub()-range[1].lb())*factorY);


//        int val_00, val_01, val_10, val_11; // val_y_x
//        val_00 = img_integral.at<int>(y_min, x_min);
//        val_01 = img_integral.at<int>(y_min, x_max);
//        val_10 = img_integral.at<int>(y_max, x_min);
//        val_11 = img_integral.at<int>(y_max, x_max);

//        cout << val_00 << " " << val_01 << " " << val_10 << " " << val_11 << endl;

//        int areaBox = (x_max-x_min)*(y_max-y_min);
//        int areaIntegral = val_11 - val_01 - val_10 + val_00;
//        cout << "aeraBox = " << areaBox << endl;
//        cout << "aeraIntegral = " << areaIntegral << endl;

//        if(areaBox == areaIntegral){
//            cout << "TRUE" << endl;
//        }
//        else{
//            cout << "FALSE" << endl;
//        }

//        rectangle(img_integral, Point(x_min, y_min), Point(x_max, y_max), Scalar(1e6), 5.0);
//    }

//    const char* window_title = "Hello, OpenCV!";
//    namedWindow (window_title, CV_WINDOW_NORMAL);
//    imshow(window_title, img);

//    window_title = "Filled";
//    namedWindow (window_title, CV_WINDOW_NORMAL);
//    imshow(window_title, img_filled);

//    window_title = "Integral";
//    namedWindow (window_title, CV_WINDOW_NORMAL);
//    imshow(window_title, img_integral);


//    waitKey(0);
}

void test_car_on_hill(){
    IntervalVector box(2);
    box[0] = Interval(0,10);
    box[1] = Interval(0,10);

    IntervalVector box_diff(2);
    box_diff[0] = Interval(4,5);
    box_diff[1] = Interval(4,5);

    IntervalVector* box_result;
    int size = box.diff(box_diff, box_result);

    cout << "box " << box << endl;
    cout << "box_diff " << box_diff << endl;
    cout << "box_result = " << box_result << endl;
    cout << "size = " << size << endl;
    for(int i=0; i<size; i++)
        cout << "box_result = " << box_result[i] << endl;
}

void sandbox(){
        vibes::beginDrawing();
        vibes::newFigure("test");
        vibes::setFigureProperties(vibesParams("x",0,"y",0,"width",500,"height",500));
        vibes::axisLimits(-5, 5, -5, 5);
        vibes::drawBox(-4, 4, -4, 4, "b[]");

        Variable x1, x2;
        ibex::Function f(x1, x2, Return(-x2,-(1.0*(1.0-pow(x1, 2))*x2-x1)));

        double dt = 0.001;
        double t_max = 10;
        double x1_min = -4;
        double x1_max = 4;
        double x2_min = -4;
        double x2_max = 4;
        double nb_point_line = 100;

        double r = min((x1_max-x1_min)/(3*nb_point_line), (x2_max-x2_min)/(3*nb_point_line));
//        double r = 0.002;

        for(double x1 = x1_min; x1<x1_max; x1+=(x1_max-x1_min)/nb_point_line){
            for(double x2 = x2_min; x2<x2_max; x2+=(x2_max-x2_min)/nb_point_line){
                IntervalVector X(2);
                X[0] = Interval(x1);
                X[1] = Interval(x2);
                bool out = false;
                for(double t=0; t<t_max; t+=dt){
                    IntervalVector k1 = f.eval_vector(X);
                    X[0] += dt*k1[0];
                    X[1] += dt*k1[1];
                    if(X[0].mid()<x1_min || X[0].mid()>x1_max || X[1].mid()<x2_min || X[1].mid()>x2_max){
                        out = true;
                        break;
                    }
                }

                if(out){
                    vibes::drawCircle(x1, x2, r, "r[r]");
                }
                else{
                    vibes::drawCircle(x1, x2, r, "g[g]");
                }
            }
        }

//        IntervalVector dposition = f.eval_vector(box);

    //    Interval dx = dposition[0];
    //    Interval dy = dposition[1];

    //    Interval theta = atan2(dy, dx);
    //    cout << setprecision(80) << theta << endl;
    //    cout << Interval::HALF_PI << endl;

    /// Angle tests
//    IntervalVector ptA(2);
//    ptA[0] = Interval(0.0);
//    ptA[1] = Interval(0.0);

//    IntervalVector ptB(2);
//    ptB[0] = Interval(1.0);
//    ptB[1] = Interval(1.0);

//    bool normal = true;
//    if(normal){
//        IntervalVector tmp(ptA);
//        ptA = ptB;
//        ptB = tmp;
//    }

//    IntervalVector deltaP(2);
//    deltaP[0] = ptB[0] - ptA[0];
//    deltaP[1] = ptB[1] - ptA[1];

//    IntervalVector deltaN(2);
//    deltaN[0] = ptA[0] - ptB[0];
//    deltaN[1] = ptA[1] - ptB[1];

//    Interval thetaP = atan2(deltaP[1], deltaP[0]);
//    Interval thetaN = atan2(deltaN[1], deltaN[0]);

//    cout << thetaP << endl;
//    cout << thetaN << endl;


//    // Calcul de l'angle
//    Interval theta_high = thetaP | Interval::PI;
//    Interval theta_low = -Interval::PI | thetaN;

//    Interval theta_inter = theta_high & theta_low;

//    if(theta_inter.is_empty()){
//        cout << "theta = " << theta_high << " " << theta_low << endl;
//    }
//    else{
//        cout << "theta = " << theta_inter << endl;
//    }

}

void test_inter_pave_perimeter(){
    //Pave(const ibex::IntervalVector &position, const std::vector<ibex::Function *> &f_list, bool diseable_singeleton=false, bool active=true);
    IntervalVector position(2);
    position[0] = Interval(0,1);
    position[1] = Interval(0,1);

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    Pave *p = new Pave(position, f_list, false, true);
    p->get_border(0)->set_full();
    p->get_border(1)->set_full();
//    p->get_border(2)->set_segment(Interval(0.8, 1), false);
//    p->get_border(3)->set_segment(Interval(0, 0.2), false);

    cout << "perimeter = " << p->get_perimeter() << endl;


}

void test_possible_path(){
    IntervalVector position(2);
    position[0] = Interval(-0.1,1);
    position[1] = Interval(-0.1,1);

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);

    Pave *p = new Pave(position, f_list, false, true);
    p->get_border(0)->set_full();
    p->get_border(1)->set_full();
    test_draw(p, "p");

    Pave *p2 = new Pave(p);
    p2->complementaire();
    test_draw(p2, "p2");

    p->print();
    p2->print();

    p->inter(p2);
    p->print();
    test_draw(p, "inter");
}

void test_infinity(){
    Variable x, y;
    ibex::Function f(x, y, Return(y,(1.0*(1.0-pow(x, 2))*y-x)));

    IntervalVector box(2);
    box[0] = Interval::POS_REALS + Interval(4);
    box[1] = Interval::NEG_REALS - Interval(4);
    cout << box << endl;
    cout << f.eval(box) << endl;
}

void test_chi_function(){
    Variable phi, d;

    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1.0/d-1)*sin(phi)),
                                    -cos(phi)));
    IntervalVector box(2);
    box[0] = Interval(-0.785398, 0);
    box[1] = Interval(0, 1.25);
    cout << box << endl;
    IntervalVector result = f.eval_vector(box);
    cout << "result = " << result << endl;

//    result[0] = Interval(-1e300, 0.14);
    cout << "result = " << result << endl;
    result[0] = Interval::NEG_REALS;
    result[1] = Interval(-1, -0.5);

    cout << isfinite(result[0].lb()) << " -> " << result[0].lb()<< endl;
    cout << isfinite(result[0].ub()) << " -> " << result[0].ub() << endl;

    Interval theta = atan2(result[1], result[0]);
    cout << "theta = " << theta << endl;

    //
    cout << "Test equal <=0" << endl;
    Interval t1(-3, -2);
    if(t1.is_subset(Interval::NEG_REALS))
        cout << "ok" << endl;

    ibex::Variable x1, x2;
    ibex::Function f_inside_curve(x1, x2, 0.5*((-(21.0/10.0) + x1)*(505.0/147.0*(-(21.0/10.0) + x1) + 10.0/21.0*(-(99.0/50.0) + x2)) + (10.0/21.0*(-(21.0/10.0) + x1) + (2665.0*(-(99.0/50.0) + x2))/1386.0)*(-(99.0/50.0) + x2)) - 1.0/8.0);
    IntervalVector test(2);
    test[0] = Interval(-0.01,0.15);
    test[1] = Interval(-0.01,0.15);
    cout << "test inside = " << f_inside_curve.eval_vector(test) << endl;
}

void test_diff_infinity(){
    Utils u;
    IntervalVector box_R(2);
    IntervalVector box_search(2);
    box_search[0] = Interval(-1,1);
    box_search[1] = Interval(-1,1);

    vector<IntervalVector> box_diff_list = u.diff(box_R, box_search);
    for(IntervalVector &box:box_diff_list)
        cout << box << endl;
}
