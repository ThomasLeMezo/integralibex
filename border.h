#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <pave.h>

class Pave;
class Border
{
/***************** Functions ******************/
public:
    Border(const ibex::IntervalVector& position, const int face, Pave *pave);
    ~Border(){}

    // ******** Drawing functions ********
    void draw() const;

    // ******** Graph building ********
    void add_brothers(std::vector<Border *> brother_list);
    void add_brothers(Border* brother);
    void update_brothers(Border* border_pave1, Border* border_pave2);

    // ******** Border Properties ********
    // Operations
    Border& operator&=(const Border &p);
    void remove_brother(int indice);

    // Setters
    void set_full();
    void set_empty();

    void set_segment_in(ibex::Interval segment_in);
    void set_segment_out(ibex::Interval segment_out);

    // Getters
    void get_points(std::vector<double> &x, std::vector<double> &y);
    ibex::Interval segment_in() const;
    ibex::Interval segment_out() const;
    std::vector<Border *> brothers();
    ibex::IntervalVector position();
    Pave* pave();
    ibex::Interval segment_full();

    // Tests
    bool is_empty();
    bool is_full();

/***************** Variables ******************/
private:
    ibex::Interval m_segment_in, m_segment_out;
    ibex::Interval m_segment_full;

private:
    int m_face;                               // Number of the face (0=bottom, 1=right, ...)
    std::vector<Border*> m_brothers;          // Pointer to brothers Borders
    ibex::IntervalVector m_position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave *m_pave;                             // Pointer to its container

private:
    bool m_empty;
    bool m_full;
};

#endif // BORDER_H
