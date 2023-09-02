import argparse
from Bio.PDB import *
import numpy as np

from AminoAcid import * # je crois que c'est une mauvaise pratique mdr

class Protein:
    def __init__(self, name):
        self.name = name
        self.mass_center = (0,0,0) # jsp
        self.amino_acid_sequence = []
    
    def get_amino_acid_sequence(self):
        return self.amino_acid_sequence
    
    def add_to_amino_acid_sequence(self, new_amino_acid):
        self.amino_acid_sequence.append(new_amino_acid)
    
    def print_protein(self):
        for residue in self.amino_acid_sequence:
            print(residue.get_code(), end=" ")
        print()

    # TODO: Méthode qui calcule le centre de masse

# A mettre ailleurs : 

def check_input_file(input_file):
    #TODO: Check extension, ...
    with open(input_file, "r") as handle:
        header_dict = parse_pdb_header(handle)
        print(header_dict)

def parse_pdb(input_file):
    #TODO: Parse the PDB File to get only one chain (pour le moment, cf extension) and get residues and atoms and put them in objects
    p = PDBParser()
    # TODO: Adapt le X
    structure = p.get_structure("X", input_file)
    id_amino_acid = 0
    protein = Protein(name=input_file[:-4])
    # Adapt with the condition of the exercise : 
    for model in structure:
        for chain in model:
            for residue in chain:
                # Creating and setting the AminoAcid object
                new_amino_acid = AminoAcid(code=residue.get_resname(), id = id_amino_acid)
                for atom in residue:
                    x,y,z = atom.get_coord()
                    atom_object = AminoAcid.Atom(symbol=atom.get_name(),x=x,y=y,z=z, 
                                       amino_acid = new_amino_acid) 
                    new_amino_acid.add_atom_to_amino_acid(atom_object)
                protein.add_to_amino_acid_sequence(new_amino_acid=new_amino_acid)
                id_amino_acid+=1
            protein.print_protein()


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