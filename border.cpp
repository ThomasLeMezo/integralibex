#include "border.h"
#include "ibex.h"

using namespace ibex;

Border::Border(const ibex::IntervalVector &p, int face, Pave *pave): position(2)
{
    this->position = p;
    this->face = face;
    this->pave = pave;
}

