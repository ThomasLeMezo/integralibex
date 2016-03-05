#ifndef BORDER_H
#define BORDER_H

#include <ibex.h>
#include <ppl.hh>
#include <pave.h>
#include <inclusion.h>

#include "vtkSmartPointer.h"
#include "vtkPoints.h"

namespace PPL = Parma_Polyhedra_Library;

class Pave;
class Inclusion;
class Border
{
/***************** Functions ******************/
public:
    Border(const ibex::IntervalVector& position, Pave *pave, int face_axis, int face_side);
    Border(const Border *border);
    ~Border();

    // ******** Graph building ********
    void                            update_brothers_inclusion(Border *border_pave1, Border *border_pave2);
    void                            remove_inclusion(int indice);
    void                            remove_inclusion(Inclusion *inclusion);
    void                            remove_inclusion_receving(int indice);
    void                            remove_inclusion_receving(Inclusion *inclusion);

    void                            draw_vtk_get_points(vtkSmartPointer<vtkPoints> &points);

    // ******** Border Properties ********
    // Operations
    Border&                         operator&=(const Border &b);
    bool                            inter(const Border &b);
    bool                            diff(const Border &b);

    // Setters
    void                            set_full();
    void                            set_full_volume_in();
    void                            set_full_volume_out();
    void                            set_empty();
    void                            set_volume_in(PPL::C_Polyhedron volume_in, bool inclusion);
    void                            set_volume_out(PPL::C_Polyhedron volume_out, bool inclusion);
    void                            set_pave(Pave* pave);

    void                            set_inclusion(Border *border, int id_brother);
    void                            set_inclusion_receving(Border* border, int id_brother);
    void                            reset_full_empty();

    void                            add_inclusions(const std::vector<Inclusion *> &inclusion_list);
    bool                            add_inclusion(Inclusion *inclusion);
    bool                            add_inclusion_copy(Inclusion *inclusion);
    void                            add_inclusion_receving(Inclusion* inclusion);

    // Getters
    int                             get_dim() const;
    int                             get_face_axis() const;
    int                             get_face_side() const;
    int                             get_face() const;

    const PPL::C_Polyhedron         get_volume_in() const;
    const PPL::C_Polyhedron         get_volume_out() const;
    const PPL::C_Polyhedron         get_volume_full() const;

    const std::vector<Inclusion *>  get_inclusions() const ;
    const std::vector<Inclusion*>&  get_inclusions_receving() const;
    Inclusion*                      get_inclusion(int i);

    const ibex::IntervalVector &    get_position() const;
    Pave*                           get_pave() const;
    const ibex::Interval&           get_segment_full() const;
    int                             size() const;
    Inclusion*                      operator[](int id);

    // Tests
    bool                            is_empty();
    bool                            is_full();

/***************** Variables ******************/
private:
    PPL::C_Polyhedron       m_volume_in, m_volume_out;
    PPL::C_Polyhedron       m_volume_full;

private:
    int                     m_dim;
    int                     m_face_axis;    // x=1, y=2 etc.
    int                     m_face_side;    // lb = 0, ub = 1

    std::vector<Inclusion*> m_inclusions;          // Pointer to brothers Borders
    std::vector<Inclusion*> m_inclusions_receving;    // Pointer to inclusion that point to this border
    ibex::IntervalVector    m_position;          // Position of the border ([x], [y]) where one of the dimension is singleton

    Pave                    *m_pave;                 // Pointer to its container

private:
    bool                    m_empty;
    bool                    m_full;
};

#endif // BORDER_H
