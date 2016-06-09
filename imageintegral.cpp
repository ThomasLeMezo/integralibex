#include "imageintegral.h"

#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace ibex;
using namespace cv;

imageIntegral::imageIntegral(const IntervalVector &range, ibex::Function *f, const Interval &t_range, int nbBisectionT, int nbPixelXY) :
    m_range(2)
{
    m_range = range;
    m_f = f;

    m_factor[0] = 0.0;
    m_factor[1] = 0.0;

    m_img_size[0] = nbPixelXY;
    m_img_size[1] = nbPixelXY;

    m_tRange = t_range;
    m_nbBisectionT = nbBisectionT;

    compute_image();
}

void imageIntegral::set_box(IntervalVector box){
    m_range = box;
    compute_image();
}

void imageIntegral::compute_image(){
    std::vector<Interval> list_t;
    list_t.push_back(m_tRange);

    for(int i=0; i<m_nbBisectionT; i++){
        std::vector<Interval> tmp_list_t(list_t);
        list_t.clear();
        for(Interval &i:tmp_list_t){
            list_t.push_back(Interval(i.lb(), i.mid()));
            list_t.push_back(Interval(i.mid(), i.ub()));
        }
    }

    m_factor[0] = m_img_size[0]/m_range[0].diam();
    m_factor[1] = m_img_size[1]/m_range[1].diam();

    Mat img(m_img_size[0],m_img_size[1], CV_8U, Scalar(0));
    for(Interval &i:list_t){
        IntervalVector tmp(2);
        tmp[0] = i;
        IntervalVector p(m_f->eval_vector(tmp));
        if(p.diam()[0] < 1) // TO CHANGE : BUG WITH eval vector
            rectangle(img, Point(ceil((p[0].lb()-m_range[0].lb())*m_factor[0]), ceil((p[1].lb()-m_range[1].lb())*m_factor[1])),
                           Point(floor((p[0].ub()-m_range[0].lb())*m_factor[0]), floor((p[1].ub()-m_range[1].lb())*m_factor[1])),
                           Scalar(255), 5.0);
    }

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    /// Find contours
    findContours( img, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

    /// Draw contours
    Mat img_filled_inside = Mat::zeros( img.size(), CV_8U);
    Mat img_filled_contour = Mat::zeros( img.size(), CV_8U);
    Mat img_filled_mask = Mat::zeros(img.size(), CV_8U);
    for( int i = 0; i< contours.size(); i++ ){
        drawContours(img_filled_inside, contours, i, Scalar(1), CV_FILLED, 8, hierarchy, 0, Point() );
        drawContours(img_filled_contour, contours, i, Scalar(1), 6, 8, hierarchy, 0, Point() );
    }

    img_filled_mask = cv::Scalar::all(1) - img_filled_contour;
    Mat img_filled_output = img_filled_inside.mul(img_filled_mask);

    m_img_integral = Mat::zeros(img.size(), CV_32SC1); // CV_32SC1
    cv::integral(img_filled_output, m_img_integral, CV_32SC1);
}

bool imageIntegral::testBox(const ibex::IntervalVector &test_box){
    if(test_box.is_interior_subset(m_range)){
        int x_min, x_max, y_min, y_max;
        x_min = floor((test_box[0].lb()-m_range[0].lb())*m_factor[0]);
        x_max = ceil((test_box[0].ub()-m_range[0].lb())*m_factor[0]);
        y_min = floor((test_box[1].lb()-m_range[1].lb())*m_factor[1]);
        y_max = ceil((test_box[1].ub()-m_range[1].lb())*m_factor[1]);

        int val_00, val_01, val_10, val_11; // val_y_x
        val_00 = m_img_integral.at<int>(y_min, x_min);
        val_01 = m_img_integral.at<int>(y_min, x_max);
        val_10 = m_img_integral.at<int>(y_max, x_min);
        val_11 = m_img_integral.at<int>(y_max, x_max);

        int areaBox = (x_max-x_min)*(y_max-y_min);
        int areaIntegral = val_11 - val_01 - val_10 + val_00;

        if(areaBox == areaIntegral)
            return true;
        else
            return false;
    }
    else
        return false;
}
