#ifndef GRAPH_H
#define GRAPH_H

#include <pave.h>
#include <utils.h>

class Graph
{
public:
    Graph(const ibex::IntervalVector &box, ibex::Function *f, ibex::IntervalVector u, int graph_id);
    Graph(Graph* g, int graph_id=-1);
    Graph(Graph* g, Pave* activated_node, int graph_id=-1);
    ~Graph();

    int process(int max_iterations, bool backward, bool inner);
    void sivia(double epsilon_theta, int nb_node, bool backward, bool do_not_bisect_empty);
    void remove_empty_node();

    bool inter(const Graph &g);
    bool diff(const Graph &g);

    // Test
    bool is_empty();

    // Setter
    void set_full();
    void set_active_pave(const ibex::IntervalVector &box);
    void set_symetry(ibex::Function *f, int axis_in, int side_in, int axis_out, int side_out);
    void set_empty();
    void clear_node_queue();
    void add_all_to_queue();
    void set_all_first_process();

    // Getter
    Pave*               get_pave(double x, double y) const;
    const std::vector<Pave *> get_pave(const ibex::IntervalVector &box) const;
    const std::vector<Pave *>& get_node_list() const;
    const std::vector<Pave *> get_node_queue() const;
    Pave*               get_node_const(int i) const;
    Pave*               get_semi_full_node();
    int                 size() const;
    ibex::Function*     get_f_inclusion_std();
    int                 get_graph_id();

    // Other functions
    Pave& operator[](int id);
    void draw_vtk(string filename);

private:
    std::vector<Pave*> m_node_list;
    std::vector<Pave*> m_node_empty_list;
    std::vector<Pave*> m_node_queue;

    int m_graph_id;
    int m_drawing_cpt;
};

#endif // GRAPH_H
