#ifndef GRAPH_H
#define GRAPH_H

#include <pave.h>
#include <utils.h>

class Graph
{
public:
    Graph(const ibex::IntervalVector &box, const std::vector<ibex::Function *> &f_list, Utils *utils, const ibex::IntervalVector &u, int graph_id=0, bool diseable_singleton=false);
    Graph(Utils *utils, int graph_id);
    Graph(Graph* g, int graph_id=-1);
    Graph(Graph* g, Pave* activated_node, int graph_id=-1);
    ~Graph();

    int                         process(int max_iterations, bool backward, bool enable_function_iteration=true, bool inner=false);
    void                        sivia(int nb_node, bool backward, bool do_not_bisect_empty=false, bool do_not_bisect_full=false, double theta_limit=0.0);
    void                        remove_empty_node();
    void                        mark_empty_node();

    bool                        inter(const Graph &g);
    bool                        diff(const Graph &g);

    void                        push_back(Pave* p);
    void                        push_back_external_border(Pave *p);

    // Test
    bool                        is_empty();
    bool                        identify_attractor();
    bool                        is_no_active_function();

    // Setter
    void                        set_full();
    void                        set_active_pave(const ibex::IntervalVector &box);
    void                        set_symetry(ibex::Function *f, int face_in, int face_out);
    void                        set_empty();
    void                        set_all_first_process();
    void                        set_active_f(int id);
    void                        set_external_boundary(bool in, bool out);
    void                        set_all_active();
    bool                        set_all_inner(bool inner);

    void                        update_queue(bool border_condition=true, bool empty_condition=false);
    void                        clear_node_queue();
    void                        add_all_to_queue();
    void                        reset_marker_attractor();
    void                        complementaire();

    // Getter
    Pave*                       get_pave(double x, double y) const;
    const std::vector<Pave *>   get_pave(const ibex::IntervalVector &box) const;
    std::vector<Pave *>&        get_node_list();
    const std::vector<Pave *>   get_node_queue() const;
    std::vector<Pave *>&        get_node_queue_access();
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


    // Other functions
    Pave&                       operator[](int id);
    void                        build_graph();

    void                        print_pave_info(double x, double y, string color) const;
    void                        print() const;
    void                        draw(int size, bool filled, string comment="");
    void                        drawInner(bool filled);

    void                        mark_full_pave_as_inner();
    int                         get_alive_node();

    void                        propagate_inner();
    void                        recursive_mark_inner(Pave* p);

private:
    std::vector<Pave*> m_node_list;
    std::vector<Pave*> m_node_border_list;
    std::vector<Pave*> m_node_queue;

    ibex::IntervalVector m_search_box;

    Utils *m_utils;

    int m_graph_id;
    int m_drawing_cpt;
    int m_count_alive;

    bool m_compute_inner;

public:
    bool debug_marker1, debug_marker2;
};

#endif // GRAPH_H
