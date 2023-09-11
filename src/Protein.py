import argparse
from Bio.PDB import *
import numpy as np
import pymol
import copy

from AminoAcid import *  # je crois que c'est une mauvaise pratique mdr
import Vector
import sys


class Protein:
    def __init__(self, name):
        self.name = name
        self.mass_center = Vector.Point(0, 0, 0)
        self.amino_acid_sequence = []
        self.best_positions = []
        
        self.x_max = 0
        self.y_max = 0
        self.z_max = 0

        self.x_min = float('inf') # infini  
        self.y_min = float('inf') 
        self.z_min = float('inf') 

    def get_amino_acid_sequence(self):
        return self.amino_acid_sequence

    def add_to_amino_acid_sequence(self, new_amino_acid):
        self.amino_acid_sequence.append(new_amino_acid)

    def set_mass_center(self, new_mass_center):
        self.mass_center = new_mass_center

    def print_protein(self):
        for residue in self.amino_acid_sequence:
            print(residue.get_code(), end=" ")
        print()

    def find_best_axis(self):
        best_axis_val = 0
        best_number_hits = 0
        best_axis = None
        for axis in self.best_positions:
            if axis.best_hydrophobicity > best_axis_val and axis.best_number_hits > best_number_hits:
                best_axis_val = axis.best_hydrophobicity
                best_axis = copy.deepcopy(axis)
        return best_axis


# TODO: ? trouver avant les 4 points les plus éloignés et si on regarde dans des tranches en dehors de ces points =>
#  break

# TODO: Méthode qui calcule le centre de masse:,  a voir si on le calcule à la main ou pas et si on la met dans la
#  classe protéine oupas
def compute_mass_center(chain):
    result = chain.center_of_mass()
    x, y, z = result[0], result[1], result[2]
    return x, y, z


# A mettre ailleurs : 
def check_input_file(input_file):
    # TODO: Check extension, ...
    with open(input_file, "r") as handle:
        header_dict = parse_pdb_header(handle)


def caculate_solvant_accessibility(structure, input_file):
    # Using DSSP
    print("Computing solvant accessibility...")
    model = structure[0]
    dssp = DSSP(model, input_file, dssp='mkdssp')
    # Pour chaque acide aminé, on a ASA = Accessible Surface Area
    return dssp


def find_limits(self): # TODO : ca sert ? 
    for aa in self.amino_< acid_sequence:
        if aa.x > self.x_max :
            self.x_max = aa.x
        if aa.y > self.y_max :
            self.y_max = aa.y
        if aa.z > self.z_max :
            self.z_max = aa.z
        
        if aa.x < self.x_min :
            self.x_min = aa.x
        if aa.y < self.y_min :
            self.y_min = aa.y
        if aa.z < self.z_min :
            self.z_min = aa.z

def parse_pdb(input_file):
    # TODO: Parse the PDB File to get only one chain (pour le moment, cf extension) and get residues and atoms and
    #  put them in objects
    print("Parsing the PDB file ...")
    p = PDBParser()
    # TODO: Adapt le X
    structure = p.get_structure("X", input_file)
    # TODO :mieux subdiviser en fonctions
    id_amino_acid = 1
    protein = Protein(name=input_file[8:-4])
    
    
    # Compute the the solvant accessibility and set it for each amino acid of the protein 
    dssp_res = caculate_solvant_accessibility(structure, input_file=input_file)
    selected_residues = []
    for res_id, aa, ss, asa, phi, psi, sasa, area, _,_,_,_,_,_ in dssp_res:
        if asa > 0.3: # TODO : Adapt
            selected_residues.append(res_id)
    # Adapt with the condition of the exercise : 
    for model in structure:
        print(f"Found {len(model)} chains in structure...")
        for chain in model:
            xm, ym, zm = compute_mass_center(chain)
            protein.mass_center = Vector.Point(xm, ym, zm)
            for residue in chain:
                # Getting the C-alpha of the exposed residues
                if residue.get_id()[1] in selected_residues:
                    if residue.has_id("CA"):
                        atom = residue['CA']
                        x, y, z = atom.get_coord()
                        # Creating and setting the AminoAcid object
                        new_amino_acid = AminoAcid(code=residue.get_resname(), id=id_amino_acid, x=x, y=y, z=z)
                        # Linking the amino acid to its protein object
                        protein.add_to_amino_acid_sequence(new_amino_acid=new_amino_acid)
                        id_amino_acid += 1    
    return protein


def show_in_pymol(plane1, plane2, pdb_file, mass_center, point_x=None):
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, "molecule_name")

    # Determine the dimensions of the molecule in the X-axis
    x_min, x_max = float("inf"), float("-inf")
    # Define the Y and Z range for the points
    y_min, y_max = -15, 50  # Adjust the Y range as needed

    for atom in pymol.cmd.get_model("molecule_name").atom:
        x = atom.coord[0]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x
    # Calculate the center of plane1 based on point_x
    if point_x is not None:
        center_x = point_x.get_x()
        center_y = point_x.get_y()
        center_z = point_x.get_z()
    else:
        # If point_x is not provided, use the midpoint between x_min and x_max
        center_x = (x_min + x_max) / 2.0
        center_y = (y_min + y_max) / 2.0  # Adjust the Y coordinate as needed

    

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
    pymol.cmd.show("cartoon", "molecule_name")

    # Zoom to fit the view
    pymol.cmd.zoom()

    # Save an image if needed
    pymol.cmd.png("planes.png")
    # TODO XRAY pour le PNG


# TODO :All the logging file https://docs.python.org/3/howto/logging.html

if __name__ == '__main__':
    # TODO: Compléter cela
    parser = argparse.ArgumentParser(
        prog='ProgramName',
        description='What the program does',
        epilog='Text at the bottom of help')
    parser.add_argument('filename', help="A PDB file...")
    # parser.add_argument('-c', '--count')      # option that takes a value
    args = parser.parse_args()
    print(f"Command : Protein.py {args.filename}")  # to adapt

    protein = parse_pdb(args.filename)
    # Number of points
    n = 15

    directions = Vector.find_points(n, protein.mass_center)
    print("Calculating the planes... ")
    # For each direction...
    #directions.reverse() #TODO: a voir vu que le dernier point placé est le sommet (spirale)
    for d in directions:
        point = copy.deepcopy(d)
        # Find the normal vector
        normal = Vector.find_director_vector(point=point, center_coordinate=protein.mass_center)
        # Construction of the two planes representing the membrane : 
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(14)  # 14 Angstrom, à voir pour la modularité TODO


        axis = Vector.Axis(p1=plane1, p2=plane2)
        best_axis_tmp =  copy.deepcopy(axis)
        # Looking above : 
        while axis.explore_axe(protein.amino_acid_sequence,ref = best_axis_tmp) == True:
            # Keeping in mind the best axis found yet : 
            #print("TEST", axis.best_number_hits , best_axis_tmp.best_number_hits)
            """if axis.best_hydrophobicity > best_axis_tmp.best_hydrophobicity : 
                best_axis_tmp = copy.deepcopy(axis)"""
            
            # TODO: Function to "reinitialize the axis by sliding the plane + resetting to 0 the matches"
            # Sliding the planes if necessary
            axis.plane1.slide_plane(1)
            axis.plane2.slide_plane(1)  # TODO : adapte la gap je pense
        
        print("BEST AXIS AFTER ABOVE EXPLORATION", best_axis_tmp)
        
        # Resetting start positions
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(14)  # 14 Angstrom, à voir pour la modularité TODO
        axis = Vector.Axis(p1=plane1, p2=plane2) # pas sure de si faut le remettre

        # Looking below :
        while axis.explore_axe(protein.amino_acid_sequence, ref = best_axis_tmp):
            
            """if axis.best_hydrophobicity > best_axis_tmp.best_hydrophobicity  : 
                best_axis_tmp = copy.deepcopy(axis)
                print("AXIS IMPROVED BELOW")"""
            # Sliding the planes if necessary
            axis.plane1.slide_plane(-1)
            axis.plane2.slide_plane(-1)
            
            
        
        print("BEST AXIS AFTER below EXPLORATION", best_axis_tmp)
        
        # Saving the best position for the axis
        protein.best_positions.append(best_axis_tmp)

    

    best_axis = protein.find_best_axis()
    print("BEST AXIS BEFORE ADJUSTING IS ", best_axis, "WITH PLANES = ", best_axis.plane1, best_axis.plane2)    
    #show_in_pymol(best_axis.plane1,best_axis.plane2, args.filename, mass_center=protein.mass_center, point_x=best_axis.plane1.point)    
    #sys.exit(0)
    # TODO: Faire l'optimisation de la largeur de la membrane

    print("OPTIMISING MEMBRANE WIDTH...")
    print("WIDTH WAS", abs(best_axis.plane1.d - best_axis.plane2.d))
    # Adjusting the bottom plane of the best axis :
    gap = 1
    # Keeping in mynid what was the best axis with the best planes yet :
    best_axis_tmp = copy.deepcopy(best_axis)
    best_axis.plane2.slide_plane(gap) 
    
    while best_axis.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp ): 
            print("IN",best_axis.best_ratio, best_axis_tmp.best_ratio)
            # Sliding the planes if necessary
            best_axis.plane2.slide_plane(gap)
           
    print("WIDTH 1 IS", abs(best_axis_tmp.plane1.d - best_axis_tmp.plane2.d))
   
    # best_axis_tmp is the best yet
    best_axis_tmp2 = copy.deepcopy(best_axis_tmp)
    best_axis_tmp.plane2.slide_plane(-gap)
    while best_axis_tmp.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp2 ): 
            print("IN2")
            """if best_axis_tmp.best_ratio > best_axis_tmp.best_ratio : 
                best_axis_tmp2 = copy.deepcopy(best_axis_tmp)"""
                #print("BEST AXIS IMPROVED ABOVE")
            # Sliding the planes if necessary
            best_axis_tmp.plane2.slide_plane(-gap)
            
    print("WIDTH 2 IS", abs(best_axis_tmp2.plane1.d - best_axis_tmp2.plane2.d))
        
    # best_axis_tmp is the best yet
    best_axis_tmp3 = copy.deepcopy(best_axis_tmp2)
    best_axis_tmp2.plane1.slide_plane(-gap)
    while best_axis_tmp2.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp3 ): 
            print("IN3")
           
                #print("BEST AXIS IMPROVED below")
            # Sliding the planes if necessary
            best_axis_tmp2.plane1.slide_plane(-gap)
            
    print("WIDTH 3 IS", abs(best_axis_tmp3.plane1.d - best_axis_tmp3.plane2.d))
    
      # best_axis_tmp is the best yet
    best_axis_tmp4 = copy.deepcopy(best_axis_tmp3)
    best_axis_tmp3.plane1.slide_plane(gap)
    while best_axis_tmp3.explore_axe_bis(protein.amino_acid_sequence, ref =best_axis_tmp4 ): 
            print("IN4")
            # Sliding the planes if necessary
            best_axis_tmp3.plane1.slide_plane(gap)
           
    # At the end, the best axis with the good planes positions is in best_axis_tmp2
    #print("BEST AXIS OVERALL IS ", best_axis_tmp2, "WITH PLANES = ", best_axis.plane1, best_axis.plane2)
    print("WIDTH IS", abs(best_axis_tmp4.plane1.d - best_axis_tmp4.plane2.d))
    show_in_pymol(best_axis_tmp4.plane1,best_axis_tmp4.plane2, args.filename, mass_center=protein.mass_center)
    # A mettre à la toute fin : """
    #show_in_pymol(best_axis_tmp4.plane1,best_axis_tmp4.plane2, args.filename, mass_center=protein.mass_center)    

  

    
    
    
