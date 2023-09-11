import Vector # changer le nom

class AminoAcid:
        hydrophobics_amino_acids = ['PHE','GLY','ILE','LEU','MET', 'VAL', 'TRP', 'TYR']
        def __init__(self, code, id,x,y,z):
            self.id = id
            self.code = code
            self.hydrophobicity  = 0
            self.is_hydrophobic = False
            if self.code in AminoAcid.hydrophobics_amino_acids:
                self.is_hydrophobic = True
            # molecule .... ? Avec quelle structure de données la représenter
            if code in AminoAcid.hydrophobics_amino_acids:
                self.hydrophobicity = 1
            else :
                self.hydrophobicity = 0 # ou -1 ? 
            self.point = Vector.Point(x,y,z)
        
        def __str__(self):
            return f"AA n°{self.id} is {self.code} at ({self.point.get_x():.3f},{self.point.get_y():.3f},{self.point.get_z():.3f}) \n"