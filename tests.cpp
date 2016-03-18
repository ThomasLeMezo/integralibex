#include <utils.h>

#include <scheduler.h>
#include <vibes.h>
#include "iomanip"
#include "graphdot.h"
#include "conversion.h"

using namespace ibex;
using namespace std;

using namespace Parma_Polyhedra_Library::IO_Operators;

void test_CtcPropagateSegment(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);

    ibex::Variable x, y;
    ibex::Function f(x, y, Return(1.0+0.0*x, x*x+0.1));

    ibex::IntervalVector u(2);

    Graph g(box, &f, u);

//    g.get_node_list()[0]->get_border(0)->set_full();

    PPL::C_Polyhedron p(2, PPL::EMPTY);
    PPL::Variable x_p(0), y_p(1);
    p.add_generator(PPL::point(0.1*IBEX_PPL_PRECISION*y_p));
    p.add_generator(PPL::point(0.9*IBEX_PPL_PRECISION*y_p));
    g.get_node_list()[0]->get_border(0)->set_volume_in(p, false);

    g.add_all_to_queue();
    g.set_all_first_process();

    g.process(5, false, false);

    g.draw_vtk("test");
}

void test_CtcPropagateSegmentBackward(){
    ibex::IntervalVector box(2);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);

    ibex::Variable x, y;
    ibex::Function f(x, y, Return(1.0+0.0*x, x));

    ibex::IntervalVector u(2);

    Graph g(box, &f, u);

//    g.get_node_list()[0]->get_border(0)->set_full_volume_in();
    g.get_node_list()[0]->get_border(1)->set_full_volume_in();
    g.get_node_list()[0]->get_border(2)->set_full_volume_in();
    g.get_node_list()[0]->get_border(3)->set_full_volume_in();

//    g.get_node_list()[0]->get_border(0)->set_full_volume_out();
//    g.get_node_list()[0]->get_border(1)->set_full_volume_out();
//    g.get_node_list()[0]->get_border(2)->set_full_volume_out();
//    g.get_node_list()[0]->get_border(3)->set_full_volume_out();

    PPL::Variable x_p(0), y_p(1);
    PPL::C_Polyhedron p(2, PPL::EMPTY);
    p.add_generator(PPL::point(1.0  *IBEX_PPL_PRECISION*x_p   + 0.4 *IBEX_PPL_PRECISION*y_p));
    p.add_generator(PPL::point(1.0  *IBEX_PPL_PRECISION*x_p   + 0.5 *IBEX_PPL_PRECISION*y_p));
    g.get_node_list()[0]->get_border(1)->set_volume_out(p, false);


    g.get_node_list()[0]->disable_continuity();
    g.add_all_to_queue();
    g.set_all_first_process();

    g.process(5, true, false);

    g.draw_vtk("test");
}

void test_CtcPaveForward(){

}

void test_copy_graph(){
    ibex::IntervalVector box(3);
    box[0] = ibex::Interval(0,1);
    box[1] = ibex::Interval(0,1);
    box[2] = ibex::Interval(0,1);

    ibex::Variable x, y, z;
    ibex::Function f(x, y, z, Return(x, y, z));

    ibex::IntervalVector u(3);
    u[0] = ibex::Interval::ZERO;
    u[1] = ibex::Interval::ZERO;
    u[2] = ibex::Interval::ZERO;

    Graph g(box, &f, u, 1);
    g.sivia(0.0, 8, false, false);

    GraphDot graphDot(&g);
    graphDot.write("g.dot");

//    Graph g2(&g, 2);

//    GraphDot graphDot2(&g2);
//    graphDot2.write("g2.dot");

//    Graph g3(&g, g.get_pave(1.0, 1.0), 3);

//    GraphDot graphDot3(&g3);
//    graphDot3.write("g3.dot");

//    g.print();
//    g2.print();
//    g3.print();
}

void sandbox(){
//    ibex::IntervalVector box(2);
//    box[0] = ibex::Interval(5,10);
//    box[1] = ibex::Interval(0,10);

//    ibex::Interval test = ibex::Interval(0, 10);
//    cout << test.lb() << endl;

//    Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

//    ibex::IntervalVector dposition = f.eval_vector(box);

//    ibex::Interval dx = dposition[0];
//    ibex::Interval dy = dposition[1];

//    ibex::Interval theta = atan2(dy, dx);
//    cout << setprecision(80) << theta << endl;
//    cout << ibex::Interval::HALF_PI << endl;

//    ibex::IntervalVector box(3);
//    box[0] = ibex::Interval(0, 1);
//    box[1] = ibex::Interval(0, 1);
//    box[2] = ibex::Interval(0, 2);
//    ibex::LargestFirst bisector(0.0, 0.5);

//    std::pair<ibex::IntervalVector, IntervalVector> result_boxes = bisector.bisect(box);
//    cout << result_boxes.first << endl;
//    cout << result_boxes.second << endl;

//    PPL::C_Polyhedron test(3, PPL::EMPTY);
//    cout << test.is_empty() << endl;

//    PPL::Variable x(0), y(1), z(2);
//    test.add_generator(PPL::point());
//    test.add_generator(PPL::point(x));
//    test.add_generator(PPL::point(y));

//    cout << test.is_discrete() << endl;
//    cout << test.is_topologically_closed() << endl;

}
