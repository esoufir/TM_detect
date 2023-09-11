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
        best_axis = None
        for axis in self.best_positions:
            if axis.best_hydrophobicity > best_axis_val:
                best_axis_val = axis.best_hydrophobicity
                best_axis = copy.deepcopy(axis)
                print(axis)
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
                # Getting the C-alpha
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


def show_in_pymol(plane1,plane2, pdb_file, point_x = None):
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, "molecule_name")# adapty

    # TODO : Adapt
    x_min, x_max = -10, 10
    y_min, y_max = -10, 10

    # Define the step size for sampling points
    step = 1

    # Initialize an empty list to store the points
    points_on_plane = []

    # Generate points on the plane
    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation
            z = (-plane1.a * x - plane1.b * y - plane1.d) / plane1.c
            points_on_plane.append((x, y, z))

    for idx, point in enumerate(points_on_plane):
        x, y, z = point
        atom_name = f"dummy_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="yellow")
        pymol.cmd.show("spheres", f"dummy_{idx}")

    # Initialize an empty list to store the points
    points_on_plane = []

    # Generate points on the plane
    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation
            z = (-plane2.a * x - plane2.b * y - plane2.d) /plane2.c
            points_on_plane.append((x, y, z))

    for idx, point in enumerate(points_on_plane):
        x, y, z = point
        atom_name = f"dummy_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="yellow")
        pymol.cmd.show("spheres", f"dummy_{idx}")
    
     # p : 
    if point_x is not None:
        atom_name = "p"
        pymol.cmd.pseudoatom(atom_name, pos=[point_x.get_x(),point_x.get_y(),point_x.get_z()], color="blue")# TODO : Adapt
        pymol.cmd.show("spheres",atom_name)

    # Center of mass : 
    atom_name = "center_mass"
    pymol.cmd.pseudoatom(atom_name, pos=[3.19,37.1,36.2], color="magenta")# TODO : Adapt
    pymol.cmd.show("spheres",atom_name)
    
    # Protein :        
    pymol.cmd.show("cartoon", "molecule_name")

    # Zoom to fit the view
    #pymol.cmd.zoom()

    # Save an image if needed
    pymol.cmd.png("planes.png")

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
    n = 20

    directions = Vector.find_points(n, protein.mass_center)
    print("Calculating the planes... ")
    # For each direction...
    # directions.reverse() TODO: a voir vu que le dernier point placé est le sommet (spirale)
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
        while axis.explore_axe(protein.amino_acid_sequence) == True:
            # Keeping in mind the best axis found yet : 
            #print("TEST", axis.best_number_hits , best_axis_tmp.best_number_hits)
            
            if axis.best_hydrophobicity > best_axis_tmp.best_hydrophobicity : 
                best_axis_tmp = copy.deepcopy(axis)
                #print("AXIS IMPROVED ABOVE")
            
            
            #if best_axis_tmp.best_number_aa!=0 :
                #print("BEST AXIS IS", best_axis_tmp)
        
        
            # TODO: Function to "reinitialize the axis by sliding the plane + resetting to 0 the matches"
            # Sliding the planes if necessary
            axis.plane1.slide_plane(1)
            axis.plane2.slide_plane(1)  # adapte la gap je pense
            axis.best_number_hits = 0 # TODO mouais
            axis.best_number_aa = 0 # TODO : mouais
        
        #print("BEST AXIS AFTER ABOVE EXPLORATION", best_axis_tmp)
        #show_in_pymol(best_axis_tmp.plane1,best_axis_tmp.plane2, args.filename, point_x=point)    
        

        # Resetting start positions
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(14)  # 14 Angstrom, à voir pour la modularité TODO
        axis = Vector.Axis(p1=plane1, p2=plane2) # pas sure de si faut le remettre

        # Looking below :
        while axis.explore_axe(protein.amino_acid_sequence):
            
            if axis.best_hydrophobicity > best_axis_tmp.best_hydrophobicity  : 
                best_axis_tmp = copy.deepcopy(axis)
                #print("AXIS IMPROVED BELOW")
            # Sliding the planes if necessary
            axis.plane1.slide_plane(-1)
            axis.plane2.slide_plane(-1)
            axis.best_number_hits = 0 # TODO mouais
            axis.best_number_aa = 0 # TODO mouais
            
        
        #print("BEST AXIS AFTER below EXPLORATION", best_axis_tmp)
        
        # Saving the best position for the axis
        protein.best_positions.append(best_axis_tmp)

    

    best_axis = protein.find_best_axis()
    print("BEST AXIS BEFORE ADJUSTING IS ", best_axis, "WITH PLANES = ", best_axis.plane1, best_axis.plane2)    
    
    # TODO: Faire l'optimisation de la largeur de la membrane

    print("OPTIMISING MEMBRANE WIDTH...")
    print("WIDTH WAS", abs(best_axis.plane1.d - best_axis.plane2.d))
    # Adjusting the bottom plane of the best axis :
    gap = 0.1
    best_axis_tmp = copy.deepcopy(best_axis)
    best_axis.plane2.slide_plane(-gap)
    while best_axis.explore_axe_bis(protein.amino_acid_sequence): 
            print("IN",best_axis.best_ratio, best_axis_tmp.best_ratio)
            if  best_axis.best_ratio > best_axis_tmp.best_ratio : 
                best_axis_tmp = copy.deepcopy(best_axis)
                print("BEST AXIS IMPROVED BELOW")
            # Sliding the planes if necessary
            best_axis.plane2.slide_plane(-gap)
           

    # best_axis_tmp is the best yet
    best_axis_tmp2 = copy.deepcopy(best_axis_tmp)
    while best_axis_tmp.explore_axe_bis(protein.amino_acid_sequence): 
            print("IN2")
            if best_axis_tmp.best_ratio > best_axis_tmp.best_ratio : 
                best_axis_tmp2 = copy.deepcopy(best_axis_tmp)
                #print("BEST AXIS IMPROVED ABOVE")
            # Sliding the planes if necessary
            best_axis_tmp.plane2.slide_plane(gap)
            

    # best_axis_tmp is the best yet
    best_axis_tmp3 = copy.deepcopy(best_axis_tmp2)
    while best_axis_tmp2.explore_axe_bis(protein.amino_acid_sequence): 
            print("IN3")
            if best_axis_tmp2.best_ratio > best_axis_tmp2.best_ratio : 
                best_axis_tmp3 = copy.deepcopy(best_axis_tmp2)
                #print("BEST AXIS IMPROVED below")
            # Sliding the planes if necessary
            best_axis_tmp2.plane1.slide_plane(-gap)
            
    
     # best_axis_tmp is the best yet
    best_axis_tmp4 = copy.deepcopy(best_axis_tmp3)
    while best_axis_tmp2.explore_axe_bis(protein.amino_acid_sequence): 
            print("IN4")
            if best_axis_tmp3.best_ratio > best_axis_tmp3.best_ratio : 
                best_axis_tmp4 = copy.deepcopy(best_axis_tmp3)
                #print("BEST AXIS IMPROVED below")
            # Sliding the planes if necessary
            best_axis_tmp3.plane1.slide_plane(gap)
           
    # At the end, the best axis with the good planes positions is in best_axis_tmp2
    #print("BEST AXIS OVERALL IS ", best_axis_tmp2, "WITH PLANES = ", best_axis.plane1, best_axis.plane2)
    print("WIDTH IS", abs(best_axis_tmp4.plane1.d - best_axis_tmp4.plane2.d))

    # A mettre à la toute fin : 
    # show_in_pymol(best_axis_tmp4.plane1,best_axis_tmp4.plane2, args.filename)    

  

    
    
    
