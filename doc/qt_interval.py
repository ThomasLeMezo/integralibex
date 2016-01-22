#!/usr/bin/python3.4
from dumper import *
import gdb

def qdump__ibex__Interval(d, value):

	# d.putItemCount(2)
	d.putNumChild(2)
	with Children(d):
		d.putSubItem("INF", value["itv"]["INF"])
		d.putSubItem("SUP", value["itv"]["SUP"])

def qdump__ibex__IntervalVector(d, value):
	size = int(value["n"])
	array = value["vec"]

	d.putItemCount(size)
	d.putNumChild(size)
	with Children(d):
		for i in range(size):
			d.putSubItem(str(i), value["vec"][i])