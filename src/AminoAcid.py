class AminoAcid:
        
        hydrophobics_amino_acids = ['PHE','GLY','ILE','LEU','MET', 'VAL', 'TRP', 'TYR']
        #TODO: Que faire des OH ? 
        def __init__(self, code, id):
            self.id = id
            self.code = code
            self.is_hydrophobic = False
            # molecule .... ? Avec quelle structure de données la représenter
            if code in AminoAcid.hydrophobics_amino_acids:
                is_hydrophobic = True
            self.atoms_list = []
            self.asa = 0
        
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
        
        class Atom:
            def __init__(self, symbol,x,y,z, amino_acid):
                self.amino_acid = amino_acid # pas sure
                self.symbol=symbol
                self.x = x
                self.y = y
                self.z = z