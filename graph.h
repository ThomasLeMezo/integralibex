#ifndef GRAPH_H
#define GRAPH_H

#include <pave.h>
#include <utils.h>

enum GRAPH_BW_FW_DIRECTION {GRAPH_FORWARD=false, GRAPH_BACKWARD=true};
//enum GRAPH_COMPUTE_INNER {GRAPH_COMP_INNER_FALSE=false, GRAPH_COMP_INNER_TRUE=true};
//enum GRAPH_MODE {GRAPH_INNER=true, GRAPH_OUTER=false};

class Graph
{
public:
    Graph(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, Utils *utils, int graph_id=0, bool diseable_singleton=false);
    Graph(Utils *utils, int graph_id);
    Graph(Graph* g, int graph_id=-1);
    Graph(Graph* g, Pave* activated_node, int graph_id=-1);
    ~Graph();

    int                         process(int max_iterations, GRAPH_BW_FW_DIRECTION direction, bool union_functions=false);
    void                        sivia(int nb_node, GRAPH_BW_FW_DIRECTION direction, bool do_not_bisect_empty=false, bool do_not_bisect_full=false, bool apply_heuristic=false);
    void                        remove_empty_node();
    void                        mark_empty_node();

    void                        inter(const Graph &g, bool with_bwd=false);
    void                        diff(const Graph &g);
//    const Graph &               operator|(const Graph &g);

    void                        push_back(Pave* p);
    void                        push_back_queue(Pave *p);
    void                        push_back_external_border(Pave *p);
    void                        push_back_inside_curve(ibex::Function* curve);

    void                        copy_to_inner();

    void                        forward(int process_iterations_max);
    void                        backward(int process_iterations_max);

    // Test
    bool                        is_empty();
    bool                        is_empty_node_queue();
    bool                        identify_attractor();
    bool                        is_no_active_function();
    bool                        is_positive_invariant();
    bool                        is_sufficiently_discretized();

    // Setter
    void                        set_full();
    void                        set_active_pave(const ibex::IntervalVector &box);
    void                        initialize_queues_with_initial_condition(const std::vector<ibex::IntervalVector> &box_list);
    void                        initialize_queues_with_initial_condition(const ibex::IntervalVector &box);
    void                        initialize_queues_with_initial_condition(const ibex::Function *curve);
    void                        set_symetry(ibex::Function *f, int face_in, int face_out);
    void                        set_empty();
    void                        set_empty_outer_full_inner();
    void                        set_all_first_process();
    void                        set_active_f(int id);
    void                        set_external_boundary(bool in, bool out);
    void                        set_all_active();

    void                        set_inner_mode(bool val);
    void                        set_compute_inner(bool val);
    void                        set_positive_invariant(bool val);

    void                        set_backward_function(bool val);

    void                        update_queue(bool border_condition=true, bool empty_condition=false);
    void                        clear_node_queue();
    void                        clear_node_queue_all();
    void                        clear_node_queue_outer();
    void                        clear_node_queue_inner();
    void                        add_all_to_queue();
    void                        add_to_all_queue(Pave *p);
    void                        add_to_queue_inner(Pave *p);
    void                        add_to_queue_outer(Pave *p);
    void                        add_to_queue(Pave *p);
    void                        reset_marker_attractor();
    void                        reset_marker(vector<Pave*> list);
    void                        set_marker(vector<Pave*> list, bool val);
    void                        complementaire();
    void                        push_back_pos_attractor();
    void                        reset_pave_segment_list();

    void                        reset_queues();
    void                        reset_full_empty();

    // Getter
    Pave*                       get_pave(double x, double y) const;
    const std::vector<Pave *>   get_pave(const ibex::IntervalVector &box) const;
    std::vector<Pave *>&        get_node_list();
    std::vector<Pave *>&        get_border_list();
    const std::list<Pave *> &   get_node_queue() const;
    std::list<Pave *> &         get_node_queue_access();
    Pave*                       get_node_const(int i) const;
    Pave*                       get_semi_full_node();
    Utils*                      get_utils();
    int                         size() const;
    ibex::Function*             get_f_inclusion_std();
    int                         get_graph_id();
    ibex::IntervalVector        get_bounding_box() const;
    int                         get_f_size() const;
    ibex::IntervalVector        get_search_box() const;
    void                        get_recursive_attractor(Pave* p, vector<Pave*> &list);
    int                         get_alive_node() const;
    bool                        get_inner_mode();
    bool                        get_compute_inner();
    double                      get_area_outer();
    std::vector<double>         get_perimeters();

    void                        get_recursive_zone(Pave* p, vector<Pave*> &list);
    void                        get_recursive_contour(Pave* p, vector<Pave*> &list);
    std::vector<std::vector<Pave *> >     get_contour_nodes();
    bool                        get_positive_invariant() const;
    std::vector<std::vector<std::vector< std::vector<ibex::IntervalVector>>>> get_pos_attractor_list() const;
    std::list<ibex::Function *> get_inside_curve_list() const;


    // Other functions
    Pave&                       operator[](int id);
    void                        build_graph();
    void                        pop_front_queue();
    void                        pop_back_queue();

    void                        print_pave_info(double x, double y, string color) const;
    void                        print() const;
    void                        draw(int size, bool filled, string comment="", bool inner_only=false, int position=0, bool pos_invariant=true);
    void                        drawInner(bool filled);

    void                        compute_propagation_zone(Pave *p, bool compute_anyway=false);
    void                        compute_all_propagation_zone(bool compute_anyway=false);
    void                        reset_computation_zone();



private:
    std::vector<Pave*> m_node_list;
    std::vector<Pave*> m_node_border_list;
    std::list<Pave*> m_node_queue_inner;
    std::list<Pave*> m_node_queue_outer;

    ibex::IntervalVector m_search_box;
    std::list<ibex::Function*> m_inside_curve_list;

    Utils *m_utils;

    int m_graph_id;
    int m_count_alive;

    bool m_compute_inner;
    bool m_inner_mode;

    bool m_positive_invariant;

    std::vector<std::vector<std::vector< std::vector<ibex::IntervalVector>>>> m_pos_attractor_list;

public:
    bool debug_marker1, debug_marker2;
};

#endif // GRAPH_H
