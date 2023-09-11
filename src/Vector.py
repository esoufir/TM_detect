import numpy as np
import math
import copy
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt

class Vector:
    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def get_coordinates(self):
        return self.coordinates

    def get_x(self):
        return self.coordinates[0]

    def get_y(self):
        return self.coordinates[1]

    def get_z(self):
        return self.coordinates[2]

    def norma(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)


class Point:
    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def __str__(self):
        return f"Point({self.coordinates[0]:.3},{self.coordinates[1]:.3},{self.coordinates[2]:.3})"

    def dot_product(self, normal):
        return -self.coordinates.dot(normal.get_coordinates())

    def get_x(self):
        return self.coordinates[0]

    def get_y(self):
        return self.coordinates[1]

    def get_z(self):
        return self.coordinates[2]

    def get_coordinates(self):
        return self.coordinates


class Plane:
    def __init__(self, point: Point, normal: Vector):
        self.point = point
        self.normal = normal
        self.a = normal.get_coordinates()[0]
        self.b = normal.get_coordinates()[1]
        self.c = normal.get_coordinates()[2]
        # self.d = -self.point.dot_product(self.normal)
        self.d = -(self.a * self.point.get_x() + self.b * self.point.get_y() + self.c * self.point.get_z())  # https://mathworld.wolfram.com/Plane.html

    def __str__(self):
        return f"{self.a:.3f}x + {self.b:.3f}y +{self.c:.3f}z + {self.d:.3f} = 0\n"

    def complementary(self, gap):
        # gap in Angstrom
        # Copying the object
        complementary_plane = copy.deepcopy(self)
        complementary_plane.d = complementary_plane.d + gap
        return complementary_plane

    def slide_plane(self, sliding_window):
        # slide the plane 
        self.d += sliding_window

    # http://mathonline.wikidot.com/point-normal-form-of-a-plane
    # TODO: Ces 2 fonctions
    def is_above(self, point):
        # Return true if the point is located over the plane (self)
        return True if (
                                   self.a * point.get_x() + self.b * point.get_y() + self.c * point.get_z() + self.d) > 0 else False

    def is_below(self, point):
        # Return true if the point is located under the plane (self)
        return True if (
                                   self.a * point.get_x() + self.b * point.get_y() + self.c * point.get_z() + self.d) < 0 else False


    def is_in(self,point):
        # Return Tru if point is in the plane self
        print("IS IN VAL", self.a * point.get_x() +self.b * point.get_y() + self.c * point.get_z() + self.d)
        return True if self.a * point.get_x() +self.b * point.get_y() + self.c * point.get_z() + self.d == 0  else False

# TODO : Error handling
class Axis:
    def __init__(self, p1, p2):
        self.plane1 = p1 # deep copy ?
        self.plane2 = p2
        self.best_number_hits = 0
        self.best_number_aa = 0
        self.best_ratio = 0
        self.best_hydrophobicity = -1000 #TODO: moauis

    def __str__(self):
        return f"AXIS with best hydro : {self.best_hydrophobicity}, {self.plane1}, {self.plane2}"

    #TODO : Potenetiellemnt à mettre dans Protéine:

    
    def explore_axe(self, amino_acid_sequence, ref): #axis de ref avec les meilleurs métriques obtenues
        in_between_planes = []
        number_atoms_in_between = 0 # en vrai sert à rien vu que c'est len de in_between_planes. 
        number_atoms_hits = 0
        for aa in (amino_acid_sequence):
            if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)) or (self.plane2.is_below(
                    aa.point) and self.plane1.is_above(aa.point)):
                number_atoms_in_between += 1
                in_between_planes.append(aa)

        number_atoms_in_between = len(in_between_planes)

        # When no more atoms between the two planes, exiting the function, we stop the exploring on this side of the axis
        if number_atoms_in_between == 0:
            print("No more atoms in between")
            return (False)
        
        # Computing the relative hydrophobicity of the selected amino_acids : to maximise
        hydrophobicity, n_hits, n_total = compute_relative_hydrophobicity(in_between_planes)

        #print(f"Hit ratio {number_atoms_hits} \t {number_atoms_in_between} = {number_atoms_hits/number_atoms_in_between}")
        # TODO : améliorer hydrophobicité seulement si le nombre de résidus hydrophobes et le total sont plus grand qu'avant ???

        if  hydrophobicity > ref.best_hydrophobicity and n_hits > ref.best_number_hits and n_total > ref.best_number_aa:
            # Updating the "best" match
            """print("NHITS was", ref.best_number_hits)
            print("N_tot was", ref.best_number_aa)
            print("HYDROPHOBICITY was",ref.best_hydrophobicity)"""
            # COPIE GENERALE ?
            ref.best_number_hits = n_hits
            ref.best_number_aa = n_total
            ref.best_hydrophobicity = hydrophobicity
            ref.plane1 = copy.deepcopy(self.plane1)
            ref.plane2 = copy.deepcopy(self.plane2)

            """print("NHITS is", ref.best_number_hits)
            print("N_tot is", ref.best_number_aa)
            print("HYDROPHOBICITY is",ref.best_hydrophobicity)"""
        return (True) if len(in_between_planes) > 0 else (False)


    def explore_axe_bis(self, amino_acid_sequence,ref):
        in_between_planes = []
        number_atoms_in_between = 0 # en vrai sert à rien vu que c'est len de in_between_planes. 
        for aa in (amino_acid_sequence):
            if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)) or (self.plane2.is_below(
                    aa.point) and self.plane1.is_above(aa.point)):
                number_atoms_in_between += 1
                in_between_planes.append(aa)

        number_atoms_in_between = len(in_between_planes)
        print("NB ATOMS CONSIDERED", number_atoms_in_between)
        # When no more atoms between the two planes, exiting the function, we stop the exploring on this side of the axis
        if number_atoms_in_between == 0:
            print("No more atoms in between")
            return (False)
        
        # Computing the relative hydrophobicity of the selected amino_acids : to maximise
        hydrophobicity, n_hits, n_total = compute_relative_hydrophobicity(in_between_planes)
        # TODO : améliorer hydrophobicité seulement si le nombre de résidus hydrophobes et le total sont plus grand qu'avant ???
        #print("NHITS was", ref.best_number_hits)
        #print("N_tot was", ref.best_number_aa)
        #print("HYDROPHOBICITY was",hydrophobicity)     
        #print("NHITS REF",n_hits)
        #print("N_tot REF", n_total)
        #print("HYDROPHOBICITY REF",ref.best_hydrophobicity)      
        if  hydrophobicity > ref.best_hydrophobicity and number_atoms_in_between > ref.best_number_aa:
            # Updating the "best" match
            ref.best_number_hits = n_hits
            ref.best_number_aa = n_total
            ref.best_hydrophobicity = hydrophobicity
            ref.plane1 = copy.deepcopy(self.plane1)
            ref.plane2 = copy.deepcopy(self.plane2)
            return True
        else:
            # If its not better
            return False

def compute_relative_hydrophobicity(amino_acid_sequence):
    relative_hydrophobicity = 0
    for aa in amino_acid_sequence:
        relative_hydrophobicity += aa.hydrophobicity
    # Number of hits and number of atoms must be significant important #TODO > 10 ? pour éviter les 1/1 ratios
    return (relative_hydrophobicity/len(amino_acid_sequence), relative_hydrophobicity, len(amino_acid_sequence))


"""def plot_plane(plane1, plane2=None, point=None):
    X, Z = np.meshgrid(range(100), range(100))
    # Calculate the corresponding y values for the plane
    # mouais, plutot faire en sorte quil ne soit jamais 0
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #mass center:

    ax.scatter3D(3.19,37.1,36.2, color="blue")

    ax.scatter3D(point.get_x(), point.get_y(), point.get_z(), color="red")
    if point is not None and plane1.d != 0:
        Y1 = (-plane1.a * X - plane1.c * Z - plane1.d) / plane1.b
        ax.plot_surface(X, Y1, Z, alpha=0.2, color="green")

    if plane2 is not None and plane2.d != 0:
        Y2 = (-plane2.a * X - plane2.c * Z - plane2.d) / plane2.b
        ax.plot_surface(X, Y2, Z, color='red', alpha=0.2)
    plt.show()"""


# TODO: faire que le demi cercle
def find_points(n_points, center_coordinates):
    points = []
    for k in range(1, n_points + 1):  # n-points pour éviter division par 0
        h = -1 + (2 * (k - 1) / (n_points - 1))
        theta = math.acos(h)
        if k == 1 or k == n_points:
            phi = 0
        else:
            phi = (phi + (3.6 / math.sqrt(n_points) * (1 / math.sqrt(1 - h * h)))) % (2 * math.pi)
        # Centering on the mass center
        x = center_coordinates.get_x() + math.sin(phi) * math.sin(theta)
        y = center_coordinates.get_y() + math.cos(theta)
        z = center_coordinates.get_z() + math.cos(phi) * math.sin(theta)
        points.append(Point(x, y, z))
    above_x_axis_points = [point for point in points if point.get_z() > center_coordinates.get_z()]
    return above_x_axis_points


def find_director_vector(point: Point, center_coordinate: Point):
    """Find the vector that goes through the two points = the normal vector of the plan"""
    # The normal vector of the plan
    return Vector(center_coordinate.get_x() - point.get_x(), center_coordinate.get_y() - point.get_y(),
                  center_coordinate.get_z() - point.get_z())


def find_normal_plan(director_vector):
    """Returns 2 normal vectors making a plan orthogonal to the director vector"""
    x, y, z = director_vector
    return (np.array([-y, x, 0]), np.array([-z, 0, x]))


def caculate_solvant_accessibility(structure, input_file):  # a voir pour identifier le bon aa
    # Using DSSP
    print("Computing solvant accessibility...")
    model = structure[0]
    dssp = DSSP(model, input_file)
    # Store ASA values for each residue
    asa_values = {}

    # Iterate over the residues in the structure
    for residue in model.get_residues():
        chain_id = residue.parent.id
        residue_id = residue.id
        key = (chain_id, residue_id)
        if key in dssp:
            asa_value = dssp[key][3]  # ASA value for the residue
            asa_values[key] = asa_value
    # Now, asa_values is a dictionary containing ASA values for each residue in the structure
    return asa_values

# TODO : Prendre le plus d'aa hits possible dans la tranche avec le meilleur ratio"""