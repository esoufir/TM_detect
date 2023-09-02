class AminoAcid:
        
        hydrophobics_amino_acids = ['PHE','GLY','ILE','LEU','MET', 'VAL', 'TRP', 'TYR']
        
        def __init__(self, code, id):
            self.id = id
            self.code = code
            is_hydrophobic = False
            # molecule .... ? Avec quelle structure de données la représenter
            if code in AminoAcid.hydrophobics_amino_acids:
                is_hydrophobic = True
            self.atoms_list = []

        def get_code(self):
            return(self.code)
        def add_atom_to_amino_acid(self, new_atom):
            self.atoms_list.append(new_atom)
        
        class Atom:
            def __init__(self, symbol,x,y,z, amino_acid):
                self.amino_acid = amino_acid # pas sure
                self.symbol=symbol
                self.x = x
                self.y = y
                self.z = z