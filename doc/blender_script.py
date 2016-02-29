import bpy
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector

class Box:
    x=[-1.0, 1.0]
    y=[-1.0, 1.0]
    z=[-1.0, 1.0]

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Border:
    inf = NaN
    sup = NaN

    def __init__(self):
        self.inf = NaN
        self.sup = NaN

    def __init__(self, inf, sup):
        self.inf = inf
        self.sup = sup

    def is_nan(self):
        if inf == NaN or sup == NaN :
            return True
        else:
            return False

class Pave:
    border = [Border()]*4

    def __init__(self, b0, b1, b2, b3):
        self.border[0] = b0
        self.border[1] = b1
        self.border[2] = b2
        self.border[3] = b3

class Cube:
    pave = [Pave()]*6

    def __init__(self, p0, p1, p2, p3, p4, p5):
        self.pave[0] = p0
        self.pave[1] = p1
        self.pave[2] = p2
        self.pave[3] = p3
        self.pave[4] = p4
        self.pave[5] = p5

if __name__ == "__main__":
    ## Box mesh
    box = Box([-1, 1], [-1,1], [-1, 1])

    box_verts = []
    for pt_x in box.x:
        for pt_y in box.y:
            for pt_z in box.z:
                box_verts.append((pt_x, pt_y, pt_z))
    # 
    box_edges = []
    for x in range(2):
        for y in range(0, 2):
            for z in range(0, 2):
                box_edges.append((x*4+y*2+z, ((x+1)%2)*4 + y*2 + z))
                box_edges.append((x*4+y*2+z, x*4 + ((y+1)%2)*2 + z))
                box_edges.append((x*4+y*2+z, x*4 + y*2 + (z+1)%2))
    box_faces = []

    mesh_box = bpy.data.meshes.new(name="Node")
    mesh_box.from_pydata(box_verts, box_edges, box_faces)
    # useful for development when the mesh may be invalid.
    #mesh.validate(verbose=True)

    ## Polygone mesh
    cube = Cube(Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0)),
        Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0)),
        Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0)),
        Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0)),
        Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0)),
        Pave(Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0), Border(0.0, 1.0))
    )

    cube_verts = []
    cube_edges = []
    cube_faces = []
    for p in Cube.pave:
        for i in range(4):
            if not p.border[i].is_empty and not p.border[(i+1)%4].is_empty:
                


    bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)
    object_data_add(bpy.context, mesh_box)
