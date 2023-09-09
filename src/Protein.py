import argparse
from Bio.PDB import *
import numpy as np
import pymol

import math

from AminoAcid import * # je crois que c'est une mauvaise pratique mdr
from Bio.PDB.DSSP import DSSP
import Vector

class Protein:
    def __init__(self, name):
        self.name = name
        self.mass_center = Vector.Point(0,0,0)
        self.amino_acid_sequence = []
        self.best_positions = []
    
    def get_amino_acid_sequence(self):
        return self.amino_acid_sequence
    
    def add_to_amino_acid_sequence(self, new_amino_acid):
        self.amino_acid_sequence.append(new_amino_acid)
    
    def set_mass_center(self, new_mass_center):
        self.mass_center=new_mass_center
    
    def print_protein(self):
        for residue in self.amino_acid_sequence:
            print(residue.get_code(), end=" ")
        print()
    
    def find_best_axis(self):
        best_axis_val = 0
        best_axis = None
        for axis in self.best_positions:
            if axis.best_number_hits > best_axis_val:
                best_axis_val = axis.best_number_hits
                best_axis = axis  
        return best_axis


# TODO: ? trouver avant les 4 points les plus éloignés et si on regarde dans des tranches en dehors de ces points => break   

# TODO: Méthode qui calcule le centre de masse:,  a voir si on le calcule à la main ou pas et si on la met dans la classe protéine oupas
def compute_mass_center(chain):
    result = chain.center_of_mass()
    x,y,z = result[0], result[1], result[2]
    return x,y,z


# A mettre ailleurs : 
def check_input_file(input_file):
    # TODO: Check extension, ...
    with open(input_file, "r") as handle:
        header_dict = parse_pdb_header(handle)
        print(header_dict)


def naccess(input_file):
    return Bio.PDB.NACCESS.run_naccess(input_file)


def caculate_solvant_accessibility(structure, input_file,protein): # a voir pour identifier le bon aa
    # Using DSSP
    model = structure[0]
    dssp = DSSP(model, input_file)
    # Pour chaque acide aminé, on a ASA = Accessible Surface Area
    for i in range(len(dssp.keys())):
        a_key = list(dssp.keys())[i]
        protein.get_amino_acid_sequence()[i].set_asa(dssp[a_key][3])#pas sure pour i



def parse_pdb(input_file):
    #TODO: Parse the PDB File to get only one chain (pour le moment, cf extension) and get residues and atoms and put them in objects
    print("Parsing the PDB file ...")
    p = PDBParser()
    # TODO: Adapt le X
    structure = p.get_structure("X", input_file)
    # TODO :mieux subdiviser en fonctions
    id_amino_acid = 1
    protein = Protein(name=input_file[:-4])
    # Adapt with the condition of the exercise : 
    for model in structure:
        print(f"Found {len(model)} chains in structure...")
        for chain in model:
            xm, ym, zm = compute_mass_center(chain)
            protein.mass_center = Vector.Point(xm,ym,zm)  
            for residue in chain:
                # Getting the C-alpha
                if residue.has_id("CA"):
                    atom = residue['CA']
                    x,y,z = atom.get_coord()
                    # Creating and setting the AminoAcid object
                    new_amino_acid = AminoAcid(code=residue.get_resname(), id = id_amino_acid,x=x,y=y,z=z)
                    # Linking the amino acid to its protein object
                    protein.add_to_amino_acid_sequence(new_amino_acid=new_amino_acid)
                    id_amino_acid+=1
    # Compute the the solvant accessibility and set it for each amino acid of the protein 
    
    # TODO: Ne garder que les résidus accessibles au solvant, donc ne créer que ces objets
    caculate_solvant_accessibility(structure, input_file=input_file, protein=protein)
    return protein



# TODO :All the logging file https://docs.python.org/3/howto/logging.html

if __name__ == '__main__':
    # TODO: Compléter cela
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('filename',help="A PDB file...")
    #parser.add_argument('-c', '--count')      # option that takes a value
    args = parser.parse_args()
    print(f"Command : Protein.py {args.filename}")# to adapt
    protein = parse_pdb(args.filename)
    # Number of points
    n = 40
    directions = Vector.find_points(n*2, protein.mass_center)
    print("Calculating the planes... ")
    # For each direction...
    for d in directions:
        point  = d
        # Find the normal vector
        normal = Vector.find_director_vector(point=d, center_coordinate=protein.mass_center)
        #Construction of the two planes representing the membrane : 
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(14)  # 14 Angstrom, à voir pour la modularité TODO

        #TODO : search residues in the space -> function
        # end = False
        axis = Vector.Axis(p1=plane1,p2=plane2)
        # Looking above : 
        while axis.explore_axe(protein.amino_acid_sequence) == True:
            #Sliding the planes if necessary 
            plane1.slide_plane(1) 
            plane2.slide_plane(1)# adapte la gap je pense
        # Resetting start positions
        plane1 = Vector.Plane(point=point, normal=normal)
        plane2 = plane1.complementary(14)  # 14 Angstrom, à voir pour la modularité TODO
        #Looking below : 
        while axis.explore_axe(protein.amino_acid_sequence) == True:
            #Sliding the planes if necessary 
            plane1.slide_plane(-1) 
            plane2.slide_plane(-1)
        protein.best_positions.append(axis)
        
        #print(f"With this axis, {axis.best_number_hits} hits for {axis.best_number_aa} amino acids in-between, so a ratio of : {axis.best_number_hits/axis.best_number_aa}")
    print("BEST AXIS IS ", protein.find_best_axis())

