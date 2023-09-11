import Vector # changer le nom

class AminoAcid:
        hydrophobics_amino_acids = ['PHE','GLY','ILE','LEU','MET', 'VAL', 'TRP', 'TYR']
        def __init__(self, code, id,x,y,z):
            self.id = id
            self.code = code
            self.hydrophobicity  = 0
            # molecule .... ? Avec quelle structure de données la représenter
            if code in AminoAcid.hydrophobics_amino_acids:
                self.hydrophobicity = 1
            else :
                self.hydrophobicity = 0 # ou -1 ? 
            self.atoms_list = []
            self.asa = 0
            self.point = Vector.Point(x,y,z)
        
        def print_amino_acid(self):
            print("Amino Acid id", self.id)
            print("Amino Acid code", self.code)
            print("Amino Acid asa", self.asa)

        def get_code(self):
            return(self.code)
        
        def add_atom_to_amino_acid(self, new_atom):
            self.atoms_list.append(new_atom)
        
        def __str__(self):
            return f"AA n°{self.id} is {self.code} at ({self.point.get_x():.3f},{self.point.get_y():.3f},{self.point.get_z():.3f}) \n"