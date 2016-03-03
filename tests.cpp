#include <utils.h>

#include <scheduler.h>
#include <vibes.h>
#include "iomanip"
#include "graphdot.h"

using namespace ibex;
using namespace std;

void test_CtcPropagateSegment(){
    cout << "TEST test_CtcPropagateSegment" << endl;

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

//    CtcPropagateSegment(seg_in, seg_out, face, theta, box, ibex::Interval::EMPTY_SET);

    cout << "----------" << endl;
    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;
}

void test_CtcPaveForward(){
//    Utils u;
//    ibex::IntervalVector box(2);
//    box[0] = ibex::Interval(-2.03, -1.955);
//    box[1] = ibex::Interval(0.47, 0.545);

////    ibex::Interval command = -ibex::Interval::PI/8 | ibex::Interval::PI/8;
//    ibex::Interval command = -ibex::Interval::PI/4 | ibex::Interval::PI/4;
////    ibex::Interval command = -ibex::Interval::PI | ibex::Interval::PI;

//    ibex::Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
//    Pave p(box, &f, command);

////    p.set_theta(-ibex::Interval::HALF_PI/4.0 | ibex::Interval::HALF_PI/4.0);
////    p.set_theta((ibex::Interval::HALF_PI | 5.0*ibex::Interval::HALF_PI/4.0) + ibex::Interval::PI/3);
////    p.set_theta(-ibex::Interval::HALF_PI | ibex::Interval::HALF_PI);
//    p.get_border(0)->set_full_segment_in();
//    p.get_border(3)->set_full_segment_in();

////    p.get_border(0)->set_segment_in(ibex::Interval(0.5, 0.9), false);

//    test_draw(&p, "test_before");

//    u.CtcPaveForward(&p, false, true);

//    test_draw(&p, "test_after");
}

void test_CtcPaveConsistency(){
//    Utils u;
//    ibex::IntervalVector box(2);
//    box[0] = ibex::Interval(-1.66897, -1.5708);
//    box[1] = ibex::Interval(0.0880469, 0.166094);

//    ibex::Interval command = ibex::Interval::ZERO;
////    ibex::Interval command = -ibex::Interval::HALF_PI| ibex::Interval::PI;
////    ibex::Interval command = -5*ibex::Interval::HALF_PI/6.0| 5*ibex::Interval::HALF_PI/6.0;
////    ibex::Interval command = -ibex::Interval::PI/4 | ibex::Interval::PI/4;
////    ibex::Interval command = ibex::Interval(-1.0472, 1.0472);

////    Variable x, y;
////    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

//    ibex::Variable phi, d;
//    ibex::Function f(phi, d, Return(chi(cos(phi)-sqrt(2)/2, sin(phi)/d+1, (1/d-1)*sin(phi)),
//                                    -cos(phi)));
//    Pave p(box, &f, command);
//    p.set_theta(p.get_theta()[0] + (-ibex::Interval::PI/40.0 | ibex::Interval::PI/40.0) + ibex::Interval::HALF_PI);

////    p.set_theta(ibex::Interval::HALF_PI + ibex::Interval::PI/4);
////    p.set_theta((-ibex::Interval::HALF_PI/16.0 | ibex::Interval::HALF_PI/16.0)+2*ibex::Interval::PI/3);
////    p.set_theta(ibex::Interval(1.5708,2.67795));

////    p.get_border(0)->set_full_segment_in();
////    p.get_border(1)->set_full_segment_in();
//    p.get_border(2)->set_full_segment_in();
//    p.get_border(3)->set_full_segment_in();
//    p.get_border(1)->set_segment_in(ibex::Interval(0.157751, 0.166094), false);
////    p.get_border(1)->set_segment_in(ibex::Interval(-10,0.0), false);

////    p.get_border(0)->set_full_segment_out();
//    p.get_border(1)->set_full_segment_out();
//    p.get_border(2)->set_full_segment_out();
//    p.get_border(3)->set_full_segment_out();
////    p.get_border(2)->set_segment_out(ibex::Interval(0.1, 1.0), false);

////    p.get_border(2)->set_segment_out(ibex::Interval(-10,-5), false);

//    test_draw(&p, "test_before");
//    u.CtcPaveConsistency(&p, true, false);

//    vibes::beginDrawing();
//    vibes::newFigure("drawing_name2");
//    test_draw(&p, "test_after");
//    cout << setprecision(80) << endl;
//    p.print();
}

void test_copy_graph(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);
    ibex::Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    ibex::IntervalVector u(1);
    u[0] = ibex::Interval::ZERO;

    Graph g(box, &f, u, 1);
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
