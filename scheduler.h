#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <QObject>

#include "graph.h"
#include "ibex.h"
#include "utils.h"

enum MAZE_DISEABLE_SINGLETON {MAZE_DISEABLE_SINGLETON_ON=true, MAZE_DISEABLE_SINGLETON_OFF=false};
enum MAZE_BORDER_INNER_IN {MAZE_BORDER_INNER_IN_FULL=true, MAZE_BORDER_INNER_IN_EMPTY=false};
enum MAZE_BORDER_INNER_OUT {MAZE_BORDER_INNER_OUT_FULL=true, MAZE_BORDER_INNER_OUT_EMPTY=false};
enum MAZE_BORDER_OUTER_IN {MAZE_BORDER_OUTER_IN_FULL=true, MAZE_BORDER_OUTER_IN_EMPTY=false};
enum MAZE_BORDER_OUTER_OUT {MAZE_BORDER_OUTER_OUT_FULL=true, MAZE_BORDER_OUTER_OUT_EMPTY=false};

class Pave;
class Border;
class Scheduler : public QObject
{
    Q_OBJECT
/***************** Functions ******************/
private:
    Scheduler(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton);
public:
    Scheduler(const ibex::IntervalVector &box, const vector<ibex::IntervalVector> &remove_boxes, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_in=true, bool border_out=true);
    Scheduler(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, bool diseable_singleton, bool border_inner_in, bool border_inner_out, bool border_outer_in, bool border_outer_out);
    ~Scheduler();

    void cameleon_cycle(int iterations_max, int graph_max, int process_iterations_max, bool remove_inside, bool do_not_bisect_inside=false, bool stop_first_pos_invariant=false);
    void cameleon_propagation(int iterations_max, int process_iterations_max, const vector<ibex::IntervalVector> &initial_boxes);
    void cameleon_propagation(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe);
    void cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, const vector<ibex::IntervalVector> &initial_boxes);
    void cameleon_propagation_with_inner(int iterations_max, int process_iterations_max, ibex::IntervalVector &initial_boxe);
    bool compute_attractor(int iterations_max, int process_iterations_max);
    void cameleon_viability(int iterations_max, int process_iterations_max, bool border_condition=false);
    void find_path(int iterations_max, int process_iterations_max, const ibex::IntervalVector &boxA, const ibex::IntervalVector &boxB);

    void attractor_to_kernel();

    void set_symetry(ibex::Function *f, int face_in, int face_out);
    void push_back_inside_curve(ibex::Function *curve);
//    void set_imageIntegral(const ibex::IntervalVector &range, ibex::Function *f, const ibex::Interval &t_range, int nbBisectionT, int nbBisectionXY);
    void set_inner_mode(bool val);
    void set_external_boundary(bool in, bool out);

    // ******** Drawing functions ********
    void draw(int size, bool filled, string comment="", bool positive_invariant=true);
    void print_pave_info(int graph, double x, double y, string color="black[g]");

    // Getter
    Graph* get_graph_list(int i);
    int get_graph_id();


/***************** Variables ******************/
private:
    std::vector<Graph*> m_graph_list;
    int m_graph_id;

    Utils m_utils;

signals:
      void iteration_status(int iteration, int iteration_max);
      void publishLog(const QString text);
};

#endif // SCHEDULER_H
