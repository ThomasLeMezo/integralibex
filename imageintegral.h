#ifndef IMAGEINTEGRAL_H
#define IMAGEINTEGRAL_H

#include "ibex.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"


class imageIntegral
{
public:
    imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY);

    void compute_image();
    void set_box(IntervalVector box);
    bool testBox(const ibex::IntervalVector &test_box);

private:
    ibex::IntervalVector    m_range;
    ibex::Function*         m_f;

    double                  m_factor[2];
    int                     m_img_size[2];
    cv::Mat                 m_img_integral;

    ibex::Interval          m_tRange;
    int                     m_nbBisectionT;

};

#endif // IMAGEINTEGRAL_H
