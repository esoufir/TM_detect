import argparse
from Bio.PDB import *
import numpy as np

import math

from AminoAcid import * # je crois que c'est une mauvaise pratique mdr
from Bio.PDB.DSSP import DSSP


class Protein:
    def __init__(self, name):
        self.name = name
        self.mass_center = (0,0,0)
        self.amino_acid_sequence = []
    
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
    

    # TODO: Fonction qui détermine tous les vecteurs directeurs du plan selon le centre de masse
    def find_directions():
        pass


# TODO: Méthode qui calcule le centre de masse:,  a voir si on le calcule à la main ou pas et si on la met dans la classe protéine oupas
def compute_mass_center(chain):
    return chain.center_of_mass()


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
    id_amino_acid = 0
    protein = Protein(name=input_file[:-4])
    # Adapt with the condition of the exercise : 
    for model in structure:
        print(f"Found{len(model)} chains in structure...")
        for chain in model:
            protein.set_mass_center(compute_mass_center(chain))
            for residue in chain:
                if residue.has_id("CA"):
                    # Creating and setting the AminoAcid object
                    new_amino_acid = AminoAcid(code=residue.get_resname(), id = id_amino_acid)
                    # Getting the C-alpha
                    atom = residue['CA']
                    x,y,z = atom.get_coord()
                    # Creatingt the atom object
                    atom_object = AminoAcid.Atom(symbol=atom.get_name(),x=x,y=y,z=z, 
                                        amino_acid = new_amino_acid) 
                    # Linking the CA atom to its amino acid object
                    new_amino_acid.add_atom_to_amino_acid(atom_object)
                    # Linking the amino acid to its protein object
                    protein.add_to_amino_acid_sequence(new_amino_acid=new_amino_acid)
                    id_amino_acid+=1
            #protein.print_protein()
    
    # Compute the the solvant accessibility and set it for each amino acid of the protein 
    caculate_solvant_accessibility(structure, input_file=input_file, protein=protein)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('filename',help="A PDB file...")
    #parser.add_argument('-c', '--count')      # option that takes a value
    args = parser.parse_args()
    print(f"Command : Protein.py {args.filename}")# to adapt
    parse_pdb(args.filename)