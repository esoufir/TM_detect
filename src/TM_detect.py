import argparse
from Bio.PDB import *
import numpy as np
import pymol
import copy

import Geometry
import sys
from Protein import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='TM_detect',
        description='This programe locates the membrane in a transmbrane protein and detects transmebrane segments')
    parser.add_argument('filename', help="A PDB file")
    
    # Optionnal arguments :
    parser.add_argument('-n', help = "Number of points to place on the sphere. (default is 15)", type=int, dest = "n", default = 15) 
    parser.add_argument('-w', help = "Initial width of the membrane. (default is 14 A)", type=float, dest = "width", default = 14) 
    parser.add_argument('-g', help = "Gap of sliding membrane along an axis. (default is 1 A)", type=float, dest = "gap", default = 1) 
    parser.add_argument('-m', help = "Gap of optimising membrane's width. (default is 1 A)", type=float, dest = "gap_membrane",default = 1) 
    
    # Parsing : 
    args = parser.parse_args()
    filename= args.filename
    # Number of points
    n = args.n
    # Chain to focus on
    chain = "A"
    # Width of the initial membrane
    width = args.width
    # Size of the sliding window when moving membrane along an axis
    gap = args.gap
    # Size of the sliding window when optimizing membrane width
    gap_membrane = args.gap_membrane
      
    # Sum-up command : 
    print(f"Command : python TM_detect.py {filename} -n {n} -w {width} -g {gap} -m {gap_membrane}")

    # Some checks on the input file : 
    check_input_file(filename)

    # Parsing the PDB file :
    protein = parse_pdb(filename, chain)
    
    # Generating the points on a demi-sphere
    directions = Geometry.find_points(2 *n, protein.mass_center)
    
    print("Calculating the planes... ")
    # For each direction...
    for d in directions:
        point = copy.deepcopy(d)
        # Find the normal vector
        normal = Geometry.find_director_vector(point=point, center_coordinate=protein.mass_center)
        # Construction of the two planes representing the membrane : 
        plane1 = Geometry.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(width) 

        axis = Geometry.Axis(p1=plane1, p2=plane2)
        best_axis_tmp =  copy.deepcopy(axis)
        
        # Looking above : 
        while axis.explore_axe(protein.amino_acid_sequence,ref = best_axis_tmp) is True : 
            # Sliding the planes if necessary
            axis.plane1.slide_plane(gap)
            axis.plane2.slide_plane(gap)         
        
        
        # Resetting start positions
        plane1 = Geometry.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(width)
        axis = Geometry.Axis(p1=plane1, p2=plane2) # pas sure de si faut le remettre

        # Looking below :
        while axis.explore_axe(protein.amino_acid_sequence, ref = best_axis_tmp):
            # Sliding the two planes
            axis.plane1.slide_plane(-gap)
            axis.plane2.slide_plane(-gap)
        
        # Saving the best position for the axis
        protein.best_positions.append(best_axis_tmp)
    
    best_axis = protein.find_best_axis()
    print("Best axis found is", best_axis)    
    
    print("Optimising membrane width...")

    # Adjusting bottom plane above   
    best_axis_tmp = optimizing_width_membrane(axis_init=best_axis, gap_membrane=gap_membrane, plane_to_consider=2, 
                              amino_acid_sequence=protein.amino_acid_sequence)   
    # Adjusting bottom plane below        
    best_axis_tmp2 = optimizing_width_membrane(axis_init=best_axis_tmp, gap_membrane=-gap_membrane, plane_to_consider=2, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    # Adjusting upper plane above        
    best_axis_tmp3 = optimizing_width_membrane(axis_init=best_axis_tmp2, gap_membrane=gap_membrane, plane_to_consider=1, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    # Adjusting upper plane below        
    best_axis_tmp4 = optimizing_width_membrane(axis_init=best_axis_tmp3, gap_membrane=-gap_membrane, plane_to_consider=1, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    print("The membrane width is ", abs(best_axis_tmp4.plane1.d - best_axis_tmp4.plane2.d))
    print("Best axis found overall is", best_axis_tmp4)
    print("Writing the TM segments in ../results/")
    best_axis_tmp4.find_tm_segment(protein)
    print("Visualising in PyMol")
    show_in_pymol(best_axis_tmp4.plane1,best_axis_tmp4.plane2, filename, mass_center=protein.mass_center)
    
    