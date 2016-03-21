#ifndef PAVE_H
#define PAVE_H

#include <ibex.h>
#include <border.h>
#include <vtkAppendPolyData.h>
#include <vtkSmartPointer.h>

namespace PPL = Parma_Polyhedra_Library;

class Border;
class Pave
{

    /***************** Functions ******************/
public:
    Pave(const ibex::IntervalVector &position, ibex::Function *f, ibex::IntervalVector u);
    Pave(const Pave *p);
    ~Pave();

    Pave&                       operator&=(const Pave &p);
    bool                        inter(const Pave &p);
    bool                        diff(const Pave &p);

    // ******** Drawing functions ********
    void                        draw_vtk(vtkSmartPointer<vtkAppendPolyData> &polyData);
    void                        draw_box(vtkSmartPointer<vtkAppendPolyData> &polyData);
    void                        draw_vector_field(vtkSmartPointer<vtkAppendPolyData> &polyData);
//    void                        print();

    // ******** Graph building ********
    void                        bisect(vector<Pave *> &result);
    void                        remove_from_brothers();
    void                        remove_brothers(Pave* p, int face);

    // ******** Pave Properties ********
    // Tests
    bool                        is_empty();
    bool                        is_full();
    bool                        is_in_queue() const;

    // Setter
    void                        set_full();
    void                        set_empty();
    void                        set_theta(ibex::IntervalVector theta);
    void                        set_in_queue(bool flag);
    void                        set_copy_node(Pave *p);
    void                        set_first_process_true();
    void                        set_first_process_false();

    void                        set_volume_in(const PPL::C_Polyhedron& volume_in);
    void                        set_volume_out(const PPL::C_Polyhedron& volume_out);
    void                        set_volume_in_out(const PPL::C_Polyhedron& volume_in_out);

    void                        reset_full_empty();
    void                        disable_continuity();

    void                        update_theta();
    void                        update_ray_vector(const ibex::IntervalVector &theta, std::vector<Generator> &ray_vector_list);

    // Getters
    const std::vector<Pave*>            get_brothers(int face);
    const vector<PPL::Generator>&       get_ray_vector_field() const;
    const vector<PPL::Generator> &      get_ray_vector_backward_field() const;
    const vector<PPL::Generator>&       get_ray_command() const;
    const ibex::IntervalVector&         get_position() const;
    const C_Polyhedron                  get_volume_in() const;
    const C_Polyhedron                  get_volume_out() const;
    const C_Polyhedron                  get_volume_in_out() const;

    const std::vector<Border *> &       get_borders();
    Border*                             get_border(int face);
    const Border*                       get_border_const(int face) const;

    Pave*                               get_copy_node();
    ibex::Function*                     get_f() const;
    ibex::IntervalVector                get_u() const;
    ibex::IntervalVector                get_theta() const;

    bool                                get_first_process() const;
    int                                 get_dim() const;
    int                                 get_size() const;
    bool                                get_continuity() const;
    ibex::IntervalVector                get_compute_zone() const;

    const PPL::C_Polyhedron&            get_volume_universe() const;

    ibex::IntervalVector                get_bounding_box() const;
    bool                                is_on_border() const;

    Border* operator[](int face);

    /***************** Variables ******************/
private:
    std::vector<PPL::Generator> m_ray_vector_field;
    std::vector<PPL::Generator> m_ray_vector_field_backward;
    std::vector<PPL::Generator> m_ray_command;

    int                         m_dim;

    ibex::IntervalVector        m_position;
    ibex::IntervalVector        m_compute_zone;
    std::vector<Border*>        m_borders;

    ibex::Function              *m_f;
    ibex::IntervalVector        m_u;
    ibex::IntervalVector        m_theta;
    Pave*                       m_copy_node;

private:
    bool                        m_empty;
    bool                        m_full;

    bool                        m_in_queue;

    bool                        m_first_process;
    bool                        m_continuity;

    PPL::C_Polyhedron           m_volume_universe;
};

#endif // PAVE_H
