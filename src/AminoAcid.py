import Vector # changer le nom

class AminoAcid:
        
        hydrophobics_amino_acids = ['PHE','GLY','ILE','LEU','MET', 'VAL', 'TRP', 'TYR']
        #TODO: Que faire des OH ? 
        def __init__(self, code, id,x,y,z):
            self.id = id
            self.code = code
            self.is_hydrophobic = False
            # molecule .... ? Avec quelle structure de données la représenter
            if code in AminoAcid.hydrophobics_amino_acids:
                self.is_hydrophobic = True
            self.atoms_list = []
            self.asa = 0
            self.point = Vector.Point(x,y,z)
        
        def print_amino_acid(self):
            print("Amino Acid id", self.id)
            print("Amino Acid code", self.code)
            print("Amino Acid asa", self.asa)
            print("Amino Acid is_hydrophobic", self.is_hydrophobic)

        def get_code(self):
            return(self.code)
        def add_atom_to_amino_acid(self, new_atom):
            self.atoms_list.append(new_atom)
        
        def set_asa(self,asa):
            #Setter - Sets solvant accessibility
            self.asa=asa
        
        def __str__(self):
            return f"AA n°{self.id} is {self.code}. Is hydrophobic ? {self.is_hydrophobic}. Solvant access : {self.asa}\n"