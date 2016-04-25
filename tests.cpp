#include <utils.h>

#include <scheduler.h>
#include <vibes.h>
#include "iomanip"
#include "graphdot.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace ibex;
using namespace std;
using namespace cv;

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

    std::vector<Interval> command;
    command.push_back(Interval::ZERO);
    command.push_back(Interval::ZERO);
    u.CtcPropagateSegment(seg_in, seg_out, face, theta, box, command);

    cout << "----------" << endl;
    cout << "seg_in = " << seg_in << endl;
    cout << "seg_out = " << seg_out[0] << seg_out[1] << seg_out[2] << endl;
}

void test_CtcPaveForward(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(-2.03, -1.955);
    box[1] = Interval(0.47, 0.545);

    IntervalVector command(2);
    command[0] = Interval(1);
    command[1] = Interval(-1, 1);

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p(box, f_list, command);

//    p.set_theta(-Interval::HALF_PI/4.0 | Interval::HALF_PI/4.0);
//    p.set_theta((Interval::HALF_PI | 5.0*Interval::HALF_PI/4.0) + Interval::PI/3);
//    p.set_theta(-Interval::HALF_PI | Interval::HALF_PI);
    p.get_border(0)->set_full_segment_in();
    p.get_border(3)->set_full_segment_in();

//    p.get_border(0)->set_segment_in(Interval(0.5, 0.9), false);

    test_draw(&p, "test_before");

    u.CtcPaveForward(&p, false, true);

    test_draw(&p, "test_after");
}

void test_CtcPaveConsistency(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(3, 3.05);
    box[1] = Interval(0.1, 0.190625);


//    Interval command = Interval::ZERO;
////    Interval command = -Interval::HALF_PI| Interval::PI;
////    Interval command = -5*Interval::HALF_PI/6.0| 5*Interval::HALF_PI/6.0;
////    Interval command = -Interval::PI/4 | Interval::PI/4;
////    Interval command = Interval(-1.0472, 1.0472);
    IntervalVector command(2);
    command[0] = Interval::ZERO;
    command[1] = Interval::ZERO;

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
    Pave p(box, f_list, command);
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
    u.CtcPaveConsistency(&p, true, false);

    test_draw(&p, "test_after");
    cout << setprecision(80) << endl;
    p.print();
}

void test_CtcPaveConsistency2(){
//    PAVE x=[-3, -2] y= [-3, -2]
//    0x7d32e0
//    theta[0]=[1.64474, 1.92957] theta[1]=[ empty ] u=([0, 0] ; [0, 0])
//    border=0 0x7d9600 segment_in=[-3, -2] segment_out=[-3, -2] inclusion=0 *border=0x7d3110 segment_full=[-3, -2]
//    border=1 0x7d9608 segment_in=[-3, -2] segment_out=[-3, -2] inclusion=0 *border=0x7d54f0 segment_full=[-3, -2]
//    border=2 0x7d9610 segment_in=[-3, -2] segment_out=[-3, -2] inclusion=0 *border=0x7d4350 segment_full=[-3, -2]
//    border=3 0x7d9618 segment_in=[-3, -2] segment_out=[-3, -2]

//    BOX = ([-0.3, -0.2] ; [0.459375, 0.5125])
//    0x2271af0

    Utils u;
    IntervalVector box(2);
    box[0] = Interval(-0.1, 0.1);
    box[1] = Interval(-0.3485, -0.1);
    IntervalVector command(2);
    command[0] = Interval::ZERO;
    command[1] = Interval::ZERO;

    Variable x1, x2, x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));
    ibex::Function f1(x1, x2, Return(x2,
                                    -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 +2.0));

    ibex::Function f2(x1, x2, Return(x2,
                                    -9.81*sin( (1.1/1.2*sin(x1)-1.2*sin(1.1*x1))/2.0 ) -0.7*x2 -2.0));


    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f1);
    f_list.push_back(&f2);
    Pave p(box, f_list, command);

//    border 0	position=([-3.875, -3.75] ; [-4, -4])     	in=[ empty ]	out=[ empty ]	continuity_in = 1 continuity_out = 1 blocked_in = [ empty ] blocked_out = [ empty ]
//    border 1	position=([-3.75, -3.75] ; [-4, -3.75])     in=[-4, -3.75]	out=[ empty ]	continuity_in = 1 continuity_out = 1 blocked_in = [ empty ] blocked_out = [ empty ]
//    border 2	position=([-3.875, -3.75] ; [-3.75, -3.75]) in=[ empty ]	out=[-3.8277, -3.75]	continuity_in = 1 continuity_out = 1 blocked_in = [ empty ] blocked_out = [ empty ]
//    border 3	position=([-3.875, -3.875] ; [-4, -3.75])   in=[ empty ]	out=[ empty ]	continuity_in = 1 continuity_out = 1 blocked_in = [ empty ] blocked_out = [ empty ]

//    p.get_border(0)->set_segment_in(Interval(-3, -2), false);
//    p.get_border(0)->set_segment_out(Interval(-0.277939, -0.2), false);
    p.get_border(0)->set_full();
//    p.get_border(0)->set_full_segment_out();

//    p.get_border(1)->set_segment_in(Interval(-4, -3.75), false);
//    p.get_border(1)->set_segment_out(Interval(0.459375, 0.5125), false);
    p.get_border(1)->set_full();
//    p.get_border(1)->set_full_segment_out();

//    p.get_border(2)->set_segment_in(Interval(-0.293823, -0.2), false);
//    p.get_border(2)->set_segment_out(Interval(1,1.4), false);
//    p.get_border(2)->set_full_segment_out();
    p.get_border(2)->set_full();
//    p.get_border(2)->set_full_segment_out();

//    p.get_border(3)->set_segment_in(Interval(0.480101, 0.5125), false);
//    p.get_border(3)->set_segment_out(Interval(-3, -2), false);
//    p.get_border(3)->set_full_segment_out();
    p.get_border(3)->set_full();

    test_draw(&p, "test_before");
    u.CtcPaveConsistency(&p, true, false);

    test_draw(&p, "test_after");
    //cout << setprecision(80) << endl;
    p.print();
}

void test_CtcPaveConsistency3(){
    Utils u;
    IntervalVector box(2);
    box[0] = Interval(0, 1);
    box[1] = Interval(0,1);
    IntervalVector command(2);
    command[0] = Interval::ZERO;
    command[1] = Interval::ZERO;

    Variable x, y;
    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p(box, f_list, command, false);
    p.set_theta(-Interval::HALF_PI | -Interval::PI);
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
    u.CtcPaveConsistency(&p, true, false);

    test_draw(&p, "test_after");
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

    IntervalVector u(2);
    u[0] = Interval::ZERO;
    u[1] = Interval::ZERO;
    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Pave p1(box, f_list, u);
    Pave p2(box, f_list, u);

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

    IntervalVector command(2);
    command[0] = Interval::ZERO;
    command[1] = Interval::ZERO;

    std::vector<ibex::Function*> f_list;
    f_list.push_back(&f);
    Graph g(box, f_list, &u, command, 1);
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

void test_imageIntegral(){
    int sizeX = 1000;
    int sizeY = 1000;
    IntervalVector range(2);
    range[0] = -Interval::PI | Interval::PI;
    range[1] = Interval(0.01,10.0);

    Variable t;
    ibex::Function f(t, Return(2*atan(tan((atan2(cos(t), -sin(t))+Interval::PI-atan2(sin(t), cos(t)+1.0/sqrt(2.0)))/2.0)),
                                sqrt(pow(cos(t)+1/sqrt(2.0), 2)+pow(sin(t), 2))));

    std::vector<Interval> list_t;
    list_t.push_back(Interval::ZERO | Interval::TWO_PI);

    for(int i=0; i<10; i++){
        std::vector<Interval> tmp_list_t(list_t);
        list_t.clear();
        for(auto &i:tmp_list_t){
            list_t.push_back(Interval(i.lb(), i.mid()));
            list_t.push_back(Interval(i.mid(), i.ub()));
        }
    }

    double factorX = sizeX/range[0].diam();
    double factorY = sizeY/range[1].diam();

    Mat img(sizeX,sizeY, CV_8U, Scalar(0));
    for(auto &i:list_t){
        IntervalVector tmp(2);
        tmp[0] = i;
        IntervalVector p(f.eval_vector(tmp));
        if(p.diam()[0] < 1)
            rectangle(img, Point(ceil((p[0].lb()-range[0].lb())*factorX), ceil((p[1].lb()-range[1].lb())*factorY)), Point(floor((p[0].ub()-range[0].lb())*factorX), floor((p[1].ub()-range[1].lb())*factorY)), Scalar(255), 5.0);
    }

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    /// Find contours
    findContours( img, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

//    /// Draw contours
    Mat img_filled = Mat::zeros( img.size(), CV_8U);
    for( int i = 0; i< contours.size(); i++ ){
        drawContours(img_filled, contours, i, Scalar(1), CV_FILLED, 4, hierarchy, 0, Point() );
    }

    Mat img_integral = Mat::zeros(img.size(), CV_32SC1); // CV_32SC1
    cv::integral(img_filled, img_integral, CV_32SC1);


    IntervalVector test_box(2);
    test_box[0] = Interval(-1.7,-1.6);
    test_box[1] = Interval(0.9,1.2);

    if(test_box.is_interior_subset(range)){

        int x_min, x_max, y_min, y_max;
        x_min = floor((test_box[0].lb()-range[0].lb())*factorX);
        x_max = ceil((test_box[0].ub()-range[0].lb())*factorX);
        y_min = floor((test_box[1].lb()-range[1].lb())*factorY);
        y_max = ceil((test_box[1].ub()-range[1].lb())*factorY);


        int val_00, val_01, val_10, val_11; // val_y_x
        val_00 = img_integral.at<int>(y_min, x_min);
        val_01 = img_integral.at<int>(y_min, x_max);
        val_10 = img_integral.at<int>(y_max, x_min);
        val_11 = img_integral.at<int>(y_max, x_max);

        cout << val_00 << " " << val_01 << " " << val_10 << " " << val_11 << endl;

        int areaBox = (x_max-x_min)*(y_max-y_min);
        int areaIntegral = val_11 - val_01 - val_10 + val_00;
        cout << "aeraBox = " << areaBox << endl;
        cout << "aeraIntegral = " << areaIntegral << endl;

        if(areaBox == areaIntegral){
            cout << "TRUE" << endl;
        }
        else{
            cout << "FALSE" << endl;
        }

        rectangle(img_integral, Point(x_min, y_min), Point(x_max, y_max), Scalar(1e6), 5.0);
    }

    const char* window_title = "Hello, OpenCV!";
    namedWindow (window_title, CV_WINDOW_NORMAL);
    imshow(window_title, img);

    window_title = "Filled";
    namedWindow (window_title, CV_WINDOW_NORMAL);
    imshow(window_title, img_filled);

    window_title = "Integral";
    namedWindow (window_title, CV_WINDOW_NORMAL);
    imshow(window_title, img_integral);


    waitKey(0);
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
//    IntervalVector box(2);
//    box[0] = Interval(5,10);
//    box[1] = Interval(0,10);

//    Interval test = Interval(0, 10);
//    cout << test.lb() << endl;

//    Variable x, y;
//    ibex::Function f(x, y, Return(y,1.0*(1.0-pow(x, 2))*y-x));

//    IntervalVector dposition = f.eval_vector(box);

//    Interval dx = dposition[0];
//    Interval dy = dposition[1];

//    Interval theta = atan2(dy, dx);
//    cout << setprecision(80) << theta << endl;
//    cout << Interval::HALF_PI << endl;
}
