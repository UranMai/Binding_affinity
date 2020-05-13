import os
import sys
import math
import numpy
import pymol
import argparse
import itertools
from pymol.cgo import *

def Box(pX, pY, pZ, x, y, z):
	'''
	Sets up the search box within the protein, which is
	used in the docking protocol
	'''
	pymol.cmd.pseudoatom('Position', pos=[pX, pY, pZ])
	([X, Y, Z],[a, b, c]) = pymol.cmd.get_extent('Position')
	pymol.cmd.show('spheres', 'Position')
	minX = X+float(x)
	minY = Y+float(y)
	minZ = Z+float(z)
	maxX = X-float(x)
	maxY = Y-float(y)
	maxZ = Z-float(z)
	boundingBox = [BEGIN, LINES,
		VERTEX, minX, minY, minZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, maxX, maxY, maxZ,
		VERTEX, minX, minY, minZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, maxY, maxZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, minX, minY, minZ,
		VERTEX, minX, maxY, minZ,
		VERTEX, maxX, minY, minZ,
		VERTEX, maxX, maxY, minZ,
		VERTEX, minX, minY, maxZ,
		VERTEX, minX, maxY, maxZ,
		VERTEX, maxX, minY, maxZ,
		VERTEX, maxX, maxY, maxZ,
		END]
	boxName = 'Box'
	pymol.cmd.load_cgo(boundingBox, boxName)
	return(boxName)

parser = argparse.ArgumentParser(description='Prep ligands for AutoDock Vina')

parser.add_argument('-b',
					'--box',
					nargs='+',
					help='Draw search box')

args = parser.parse_args()

def main():
	if args.box:
		pymol.cmd.load(str(sys.argv[2]))
		pymol.cmd.extend('Box', Box)


if __name__ == '__main__': 
	main()