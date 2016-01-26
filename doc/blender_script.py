bl_info = {
    "name": "Node",
    "author": "Thomas Le Mezo",
    "version": (1, 0),
    "blender": (2, 75, 0),
    "location": "View3D > Add > Mesh > New Object",
    "description": "Adds a new Node Object",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
    }

import bpy
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector

def add_object(self, context):
    box_verts = []
    for pt_x in self.box.x:
        for pt_y in self.box.y:
            for pt_z in self.box.z:
                box_points.append((pt_x, pt_y, pt_z))
    # 
    box_edges = []
    for x in range(2):
        for y in range(0, 2):
            for z in range(0, 2):
                box_edges.append((x*4+y*2+z, ((x+1)%2)*4 + y*2 + z))
                box_edges.append((x*4+y*2+z, x*4 + ((y+1)%2)*2 + z))
                box_edges.append((x*4+y*2+z, x*4 + y*2 + (z+1)%2))
    box_faces = []

    mesh = bpy.data.meshes.new(name="New Object Mesh")
    mesh.from_pydata(box_verts, box_edges, box_faces)
    # useful for development when the mesh may be invalid.
    # mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)

class IntervalVector3Property(bpy.types.PropertyGroup):
    x=[0.0, 0.0]
    x[0] = bpy.props.FloatProperty() 
    x[1] = bpy.props.FloatProperty()

    y=[0.0, 0.0]
    y[0] = bpy.props.FloatProperty() 
    y[1] = bpy.props.FloatProperty()

    z=[0.0, 0.0]
    z[0] = bpy.props.FloatProperty() 
    z[1] = bpy.props.FloatProperty()

bpy.utils.register_class(IntervalVector3Property)
bpy.types.Property.IntervalVector3Property = bpy.props.PointerProperty(type=IntervalVector3Property)

class OBJECT_OT_add_object(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_node"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}

    box = IntervalVector3Property()

    def execute(self, context):
        add_object(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_object.bl_idname,
        text="Add Node",
        icon='PLUGIN')


# This allows you to right click on a button and link to the manual
def add_object_manual_map():
    url_manual_prefix = "http://wiki.blender.org/index.php/Doc:2.6/Manual/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_node", "Modeling/Objects"),
        )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_object)
    bpy.utils.register_manual_map(add_object_manual_map)
    bpy.types.INFO_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object)
    bpy.utils.unregister_manual_map(add_object_manual_map)
    bpy.types.INFO_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()

    bpy.ops.mesh.add_node(box=[(-1, 1), (-1, 1), (-1,1)])
