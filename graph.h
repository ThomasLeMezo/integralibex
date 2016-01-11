#ifndef GRAPH_H
#define GRAPH_H

#include <pave.h>
#include <utils.h>

class Graph
{
public:
    Graph(const ibex::IntervalVector &box, ibex::Function *f, Utils *utils, int graph_id);
    Graph(Graph* g, int graph_id=-1);
    Graph(Graph* g, Pave* activated_node, int graph_id=-1);
    ~Graph();

    int process(int max_iterations, bool backward);
    void sivia(double epsilon_theta, int iterations_max, bool backward, bool bisect_empty);
    void remove_empty_node();

    bool inter(const Graph &g);
    bool diff(const Graph &g);

    // Test
    bool is_empty();

    // Setter
    void set_full();
    void set_active_pave(const ibex::IntervalVector &box);
    void set_y_symetry();
    void set_empty();

    // Getter
    Pave*               get_pave(double x, double y);
    std::vector<Pave*>  get_pave(const ibex::IntervalVector &box);
    std::vector<Pave*>  get_node_list();
    Pave*               get_node_const(int i) const;
    Pave*               get_semi_full_node();
    Utils*              get_utils();
    int                 size() const;
    ibex::Function*     get_f_inclusion_std();

    // Other functions
    void print_pave_info(double x, double y, string color);
    void draw(int size, bool filled);

private:
    std::vector<Pave*> m_node_list;
    std::vector<Pave*> m_node_empty_list;
    std::vector<Pave*> m_node_queue;

    Utils *m_utils;

    int m_graph_id;
    int m_drawing_cpt;

    ibex::Function *m_f_inclusion_std;
};

#endif // GRAPH_H
