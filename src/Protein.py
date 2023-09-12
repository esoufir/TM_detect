import argparse
from Bio.PDB import *
import numpy as np
import pymol
import copy

from AminoAcid import *  # je crois que c'est une mauvaise pratique mdr
import Vector
import sys


class Protein:
    """
    A class to represent a Protein with one chain.

    ...

    Attributes
    ----------
    name : str
        Name of the protein
    mass_center : Point
        Coordinates of the mass center
    amino_acid_sequence : list <AminoAcid>
        List of the amino acids accessible to solvant
    full_sequence : list <AminoAcid>
        Full sequence
    best_positions : list <Axis>
        List of every axis with its own best position 
    
    Methods
    -------
    find_best_axis():
        Prints the person's name and age.
    """
    def __init__(self, name):
        self.name = name
        self.mass_center = Vector.Point(0, 0, 0)
        self.amino_acid_sequence = []
        self.full_sequence = []
        self.best_positions = [] 

    def find_best_axis(self):
        """
        Finds the best axis among the list self.best_positions. 
        The best axis is the one with the best value of hydrophobicity. 

        Returns
        -------
        best_axis : Axis
            Axis with the highest value of hydrophobicity
        """
        best_axis_val = 0
        best_axis = None
        for axis in self.best_positions:
            if axis.best_hydrophobicity > best_axis_val :
                best_axis_val = axis.best_hydrophobicity
                best_axis = copy.deepcopy(axis)
        return best_axis

    def compute_mass_center(self,chain):
        """
        Compute the coordinates of the mass center.
        
        Parameters
        ----------
        chain : Chain (Bio.PDB)
            More info to be displayed (default is None)
        
        Returns
        -------
        best_axis : Point
            Coordinates of the mass center
        """
        result = chain.center_of_mass()
        x, y, z = result[0], result[1], result[2]
        return Vector.Point(x, y, z)


def check_input_file(input_file):
    if not input_file.lower().endswith(".pdb"): 
        print("Wrong extension. A PDB file is required.")
        sys.exit(0)
    try:
        # Attempt to parse the PDB file
        parser = PDBParser()
        structure = parser.get_structure("temp", input_file)
        # Check if the structure contains at least one chain labeled as "A"
        for model in structure:
            for chain in model:
                if chain.id == 'A':
                    return True    
    except Exception as e:
        print("Issue well attempting to read the PDB file. Please check the content of your input file.")
        sys.exit(0)


def caculate_solvant_accessibility(structure, input_file):
    # Using DSSP
    print("Computing solvant accessibility...")
    model = structure[0]
    dssp = DSSP(model, input_file, dssp='mkdssp')
    # Pour chaque acide aminé, on a ASA = Accessible Surface Area
    return dssp


def parse_pdb(input_file):
    print("Parsing the PDB file ...")
    p = PDBParser()
    # TODO: Adapt le X
    structure = p.get_structure("X", input_file)
    
    id_amino_acid = 1 # Counting the accessible residues
    id_full_amino_acid = 1 # Counting all residues
    protein = Protein(name=input_file[8:-4])
    
    # Compute the the solvant accessibility and set it for each amino acid of the protein 
    dssp_res = caculate_solvant_accessibility(structure, input_file=input_file)
    print("Trimming to get only exposed residues...")
    exposed_residues = []
    # Getting only the residue id and its solvant accessibility (0.3)
    i=0
    for res_id, aa, bb, asa, _, _, _, _, _,_,_,_,_,_ in dssp_res:
        if asa > 0.3: #TODO : Adapt
            exposed_residues.append(res_id)
        i+=1
    
    print(f"Found {len(exposed_residues)} exposed residues.")
    # TODO : Adapt with the condition of the exercise : 
    model = structure[0]
    print(f"Found {len(model)} chains in structure...")# TODO : ADPAPT CHAIN
    chain = model["A"] #TODO: default
    protein.compute_mass_center(chain)
    for residue in chain:
        if residue.has_id("CA"):
            atom = residue['CA']
            x, y, z = atom.get_coord()
            # Creating and setting the AminoAcid object
            new_amino_acid = AminoAcid(code=residue.get_resname(), id=residue.get_id()[1], x=x, y=y, z=z)
            # Linking the amino acid to its protein object
            protein.full_sequence.append((new_amino_acid))
            id_full_amino_acid+=1
            # If the residue is exposed, adding it to a dedicated list :
            if residue.get_id()[1] in exposed_residues:
                if residue.has_id("CA"):
                    atom = residue['CA']
                    x, y, z = atom.get_coord()
                    # Creating and setting the AminoAcid object
                    new_amino_acid = AminoAcid(code=residue.get_resname(), id=residue.get_id()[1], x=x, y=y, z=z)
                    # Linking the amino acid to its protein object
                    protein.amino_acid_sequence.append(new_amino_acid)
                    id_amino_acid += 1    
    return protein


def show_in_pymol(plane1, plane2, pdb_file, mass_center, point_x=None):
    # Run in quiet mode to avoid problem with argument parsing
    pymol.finish_launching(['pymol', '-q'])
    
    pymol.cmd.load(pdb_file, "protein")# TODO
    pymol.cmd.split_chains(pdb_file)
    pymol.cmd.remove("solvent")
   
    # Determine the dimensions of the molecule in the X-axis
    x_min, x_max = float("inf"), float("-inf")
    # Define the Y and Z range for the points
    y_min, y_max = float("inf"), float("-inf") # Adjust the Y range as needed

    for atom in pymol.cmd.get_model("protein").atom:
        x = atom.coord[0]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x
    
    for atom in pymol.cmd.get_model("protein").atom:
        y = atom.coord[1]
        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y
    # Calculate the center of plane1 based on point_x
    if point_x is not None:
        center_x = point_x.get_x()
    else:
        # If point_x is not provided, use the midpoint between x_min and x_max
        center_x = (x_min + x_max) / 2.0

    # Define the step size for sampling points
    step = 3

    # Initialize an empty list to store the points on plane 1
    points_on_plane1 = []

    # Generate points on plane 1
    for x in np.arange(center_x - (x_max - x_min) / 2, center_x + (x_max - x_min) / 2 + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 1
            z1 = (-plane1.a * x - plane1.b * y - plane1.d) / plane1.c
            points_on_plane1.append((x, y, z1))

    # Initialize an empty list to store the points on plane 2
    points_on_plane2 = []

    # Generate points on plane 2

    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 2
            z2 = (-plane2.a * x - plane2.b * y - plane2.d) /plane2.c
            points_on_plane2.append((x, y, z2))

    # Create pseudoatoms for points on plane 1 and color them
    for idx, point in enumerate(points_on_plane1):
        x, y, z = point
        atom_name = f"plane1_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="yellow")
        pymol.cmd.show("spheres", f"plane1_{idx}")

    # Create pseudoatoms for points on plane 2 and color them
    for idx, point in enumerate(points_on_plane2):
        x, y, z = point
        atom_name = f"plane2_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="blue")
        pymol.cmd.show("spheres", f"plane2_{idx}")

    # Create pseudoatoms for the provided point_x, if any
    if point_x is not None:
        atom_name = "point_x"
        pymol.cmd.pseudoatom(atom_name, pos=[point_x.get_x(), point_x.get_y(), point_x.get_z()], color="magenta")
        pymol.cmd.show("spheres", atom_name)
     
    # Mass center
    if mass_center is not None:
        atom_name = "mass_center"
        pymol.cmd.pseudoatom(atom_name, pos=[mass_center.get_x(), mass_center.get_y(), mass_center.get_z()], color="green")
        pymol.cmd.show("spheres", atom_name)

    # Show the protein structure
    pymol.cmd.show("cartoon", "protein")


def optimizing_width_membrane(gap_membrane, axis_init,amino_acid_sequence, plane_to_consider): #plane to consider is plane1 or plane 2
    # Keeping in mynid what was the best axis with the best planes yet :
    best_axis = copy.deepcopy(axis_init)
    # First slide to be able to compare two diffent axis : 
    if plane_to_consider == 1:
        # Plane 1 is slided
        axis_init.plane1.slide_plane(gap_membrane)
    else: 
        # Plane 2 is slided
        axis_init.plane2.slide_plane(gap_membrane)
    
    while axis_init.explore_axe_bis(amino_acid_sequence, ref = best_axis): 
            if plane_to_consider == 1:
                # Plane 1 is slided
                axis_init.plane1.slide_plane(gap_membrane)
            else: 
                # Plane 2 is slided
                axis_init.plane2.slide_plane(gap_membrane)
    return best_axis
           

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='TM_detect.py',
        description='This programe locates the membrane in a transmbrane protein and detects transmebrane segments')
    parser.add_argument('filename', help="A PDB file")
    
    # Optionnal arguments :
    parser.add_argument('-c', help = "Chain to analyze in the PDB file. (default is chain A)", type=str, dest = "chain", default = "A")
    parser.add_argument('-n', help = "Number of points to place on the sphere. (default is 15)", type=int, dest = "n", default = 15) # TODO : Mouais
    parser.add_argument('-w', help = "Initial width of the membrane. (default is 14 A)", type=float, dest = "width", default = 14) 
    parser.add_argument('-g', help = "Gap of sliding membrane along an axis. (default is 1 A)", type=float, dest = "gap", default = 1) 
    parser.add_argument('-m', help = "Gap of optimising membrane's width. (default is 1 A)", type=float, dest = "gap_membrane",default = 1) 
    
    # Parsing : 
    args = parser.parse_args()
    filename= args.filename
    # Number of points
    n = args.n
    # Chain to focus on
    chain = args.chain
    # Width of the initial membrane
    width = args.width
    # Size of the sliding window when moving membrane along an axis
    gap = args.gap
    # Size of the sliding window when optimizing membrane width
    gap_membrane = args.gap_membrane
      
    # Sum-up command : 
    print(f"Command : TM_detect.py {filename} -c {chain} -n {n} -w {width} -g {gap} -m {gap_membrane}")

    # Some checks on the input file : 
    check_input_file(filename)

    # Parsing the PDB file :
    protein = parse_pdb(filename)
    
    # Generating the points on a demi-sphere
    directions = Vector.find_points(n, protein.mass_center)
    
    print("Calculating the planes... ")
    # For each direction...
    for d in directions:
        point = copy.deepcopy(d)
        # Find the normal vector
        normal = Vector.find_director_vector(point=point, center_coordinate=protein.mass_center)
        # Construction of the two planes representing the membrane : 
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(width) 

        axis = Vector.Axis(p1=plane1, p2=plane2)
        best_axis_tmp =  copy.deepcopy(axis)
        
        # Looking above : 
        while axis.explore_axe(protein.amino_acid_sequence,ref = best_axis_tmp) == True : 
            # Sliding the planes if necessary
            axis.plane1.slide_plane(gap)
            axis.plane2.slide_plane(gap)         
            
        #print("BEST AXIS AFTER ABOVE EXPLORATION", best_axis_tmp)

        # Resetting start positions
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(width)
        axis = Vector.Axis(p1=plane1, p2=plane2) # pas sure de si faut le remettre

        # Looking below :
        while axis.explore_axe(protein.amino_acid_sequence, ref = best_axis_tmp):
            # Sliding the two planes
            axis.plane1.slide_plane(-gap)
            axis.plane2.slide_plane(-gap)
            
        #print("BEST AXIS AFTER below EXPLORATION", best_axis_tmp)
        
        # Saving the best position for the axis
        protein.best_positions.append(best_axis_tmp)
    
    best_axis = protein.find_best_axis()
    print("BEST AXIS BEFORE ADJUSTING IS ", best_axis)    
    
    print("OPTIMISING MEMBRANE WIDTH...")
    #print("WIDTH WAS", abs(best_axis.plane1.d - best_axis.plane2.d))
    # Adjusting the bottom plane of the best axis :

    # Keeping in mynid what was the best axis with the best planes yet :
    """best_axis_tmp = copy.deepcopy(best_axis)
    best_axis.plane2.slide_plane(gap_mebrane) 
    
    while best_axis.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp ): 
            best_axis.plane2.slide_plane(gap_mebrane)"""
    # Adjusting bottom plane above   
    best_axis_tmp = optimizing_width_membrane(axis_init=best_axis, gap_membrane=gap_membrane, plane_to_consider=2, 
                              amino_acid_sequence=protein.amino_acid_sequence)
           
    print("WIDTH 1 IS", abs(best_axis_tmp.plane1.d - best_axis_tmp.plane2.d))
   

    # Adjusting bottom plane below        
    best_axis_tmp2 = optimizing_width_membrane(axis_init=best_axis_tmp, gap_membrane=-gap_membrane, plane_to_consider=2, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    print("WIDTH 2 IS", abs(best_axis_tmp2.plane1.d - best_axis_tmp2.plane2.d))
        
    # Adjusting upper plane above        
    best_axis_tmp3 = optimizing_width_membrane(axis_init=best_axis_tmp2, gap_membrane=gap_membrane, plane_to_consider=1, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    print("WIDTH 3 IS", abs(best_axis_tmp3.plane1.d - best_axis_tmp3.plane2.d))
    
    """best_axis_tmp4 = copy.deepcopy(best_axis_tmp3)
    best_axis_tmp3.plane1.slide_plane(gap_mebrane)
    while best_axis_tmp3.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp4 ): 
            # Sliding the planes if necessary
            best_axis_tmp3.plane1.slide_plane(gap_mebrane)"""
    
    # Adjusting upper plane below        
    best_axis_tmp4 = optimizing_width_membrane(axis_init=best_axis_tmp3, gap_membrane=-gap_membrane, plane_to_consider=1, 
                              amino_acid_sequence=protein.amino_acid_sequence)
    print("WIDTH IS", abs(best_axis_tmp4.plane1.d - best_axis_tmp4.plane2.d))
    print("BEST AXIS OVERALL", best_axis_tmp4)
    
    best_axis_tmp4.find_tm_segment(protein)

    show_in_pymol(best_axis_tmp4.plane1,best_axis_tmp4.plane2, filename, mass_center=protein.mass_center)
    
    
   
            

  

    
    
    
