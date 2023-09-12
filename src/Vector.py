import numpy as np
import math
import copy
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt

class Vector:
    """
    A class to represent a vector

    ...

    Attributes
    ----------
    coordinates : np.array
        The 3D coordinates of the vector. 

    Methods
    -------
    get_coordinates():
        Getter
    """

    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def get_coordinates(self):
        return self.coordinates



class Point:
    """
    A class to represent a Point

    ...

    Attributes
    ----------
    coordinates : np.array
        The 3D coordinates of the point. 

    Methods
    -------
    get_coordinates():
        Getter
    """
    
    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def __str__(self):
        """Redifining print() comportement."""
        return f"Point({self.coordinates[0]},{self.coordinates[1]},{self.coordinates[2]})"

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
    """
    A class to represent a plane in a 3D space.

    ...

    Attributes
    ----------
    point : Point
        Point on the plane. Used to define the plan with the normal vector.
    normal : Vector
        Vector normal of the plane. Used to define the plan with point. 
    a : float
        Constant in the cartesian equation. 
    b : float
        Constant in the cartesian equation. 
    c : float
        Constant in the cartesian equation. 
    d : float
        Constant in the cartesian equation. 
    Methods
    -------
    complementary(gap):
        Generates a plan parallel slided from a certain gap. 
    is_above():
        Getter
    is_below():
        Getter
    """
    def __init__(self, point: Point, normal: Vector):
        self.point = point
        self.normal = normal
        self.a = normal.get_coordinates()[0]
        self.b = normal.get_coordinates()[1]
        self.c = normal.get_coordinates()[2]
        # self.d = -self.point.dot_product(self.normal)
        self.d = -(self.a * self.point.get_x() + self.b * self.point.get_y() + self.c * self.point.get_z())  # https://mathworld.wolfram.com/Plane.html

    def __str__(self):
        """Redifining print() comportement."""
        return f"{self.a:.3f}x + {self.b:.3f}y +{self.c:.3f}z + {self.d:.3f} = 0\n"

    def complementary(self, gap):
        """ Generates a plan parallel under slided from a certain gap."""        
        # Copying the object
        complementary_plane = copy.deepcopy(self)
        complementary_plane.d = complementary_plane.d + gap
        return complementary_plane

    def slide_plane(self, sliding_window):
        """Slide the plane according to a sliding_window value"""        
        self.d += sliding_window

    # http://mathonline.wikidot.com/point-normal-form-of-a-plane

    def is_above(self, point):
        """
        Returns a bool indicating whether a point is above a plane. 

        Parameters
        ----------
        point : Point
            The point to check the location. 

        Returns
        -------
        bool
            If the point is above the plane (self)
        """
        return True if (
                                   self.a * point.get_x() + self.b * point.get_y() + self.c * point.get_z() + self.d) > 0 else False

    def is_below(self, point):
        """
        Returns a bool indicating whether a point is below a plane. 

        Parameters
        ----------
        point : Point
            The point to check the location. 

        Returns
        -------
        bool
            If the point is below the plane (self)
        """
        return True if (
                                   self.a * point.get_x() + self.b * point.get_y() + self.c * point.get_z() + self.d) < 0 else False

class Axis:
    """
    A class to represent an Axis.

    ...

    Attributes
    ----------
    plane1 : Plane
        The upper Plane. 
    plane2 : Plane
        The lower Plane.
    best_hydrophobicity : float
        Relative hydrophobicity between the two planes. 
    Methods
    -------
    get_coordinates():
        Getter
    """
    def __init__(self, p1, p2):
        self.plane1 = p1
        self.plane2 = p2
        self.best_hydrophobicity = -1000 # Low hydrophobicity, can only be improved

    def __str__(self):
        """Redifining print() comportement."""
        return f"AXIS with best hydro : {self.best_hydrophobicity}, {self.plane1}, {self.plane2}"
    
    def explore_axe(self, amino_acid_sequence, ref): #axis de ref avec les meilleurs métriques obtenues
        """
        Explore the axis and decides if it's the better one encountered yet.
        If there are no more atoms between the planes, stop the exploration of this axis.  

        Parameters
        ----------
        amino_acid_sequence : List <AminoAcid>
            The full sequence of the protein. 
        ref : Axis
            Axis of reference. The best one encountered yet.

        Returns
        -------
        bool
            Indicates whether to continue or not the exploration of the exis. 
        """
        in_between_planes = []
        n_total_hydrophobic = 0 
        n_total_hydrophile = 0
        nb_hydrophile_out_of_plan = 0
        n_hydrophobe_in_plan = 0

        # Getting only the atoms between the two planes :
        for aa in (amino_acid_sequence):
            if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)):
                in_between_planes.append(aa)
            if aa.is_hydrophobic :
                n_total_hydrophobic+=1
            else:
                n_total_hydrophile+=1
        # Hydrophile atoms not in plane : 
        for aa in amino_acid_sequence: 
            # if the atom is not in between the planes : 
            if aa not in in_between_planes and aa.is_hydrophobic == False:
                nb_hydrophile_out_of_plan+=1

        # Hydrophobic atoms in plan
        for atom_in_plan in in_between_planes:
            if atom_in_plan.is_hydrophobic : 
                n_hydrophobe_in_plan+=1
        
        # When no more atoms between the two planes, exiting the function, we stop the exploring on this side of the axis
        if len(in_between_planes) == 0:
            return (False)
        
        # Computing the relative hydrophobicity of the selected amino_acids : to maximise
        hydrophobicity = (nb_hydrophile_out_of_plan/n_total_hydrophile) + (n_hydrophobe_in_plan/n_total_hydrophobic)
        if  hydrophobicity > ref.best_hydrophobicity :
            # Updating the "best" match
            ref.best_hydrophobicity = hydrophobicity
            ref.plane1 = copy.deepcopy(self.plane1)
            ref.plane2 = copy.deepcopy(self.plane2)
        return (True) if len(in_between_planes) > 0 else (False)


    def explore_axe_bis(self, amino_acid_sequence,ref):
        """
        Similar to the upper fonction but returns False if the sliding does not improve hydriphobicity
        and so stops exploration.  
        """
        in_between_planes = []
        n_total_hydrophobic = 0 
        n_total_hydrophile = 0
        nb_hydrophile_out_of_plan = 0
        n_hydrophobe_in_plan = 0
        for aa in (amino_acid_sequence):
            if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)) or (self.plane2.is_below(
                    aa.point) and self.plane1.is_above(aa.point)):
                in_between_planes.append(aa)
            if aa.is_hydrophobic :
                n_total_hydrophobic+=1
            else:
                n_total_hydrophile+=1

        # Hydrophile atoms not in plane : 
        for aa in amino_acid_sequence: 
            # if the atom is not in between the planes : 
            if aa not in in_between_planes and aa.is_hydrophobic == False:
                nb_hydrophile_out_of_plan+=1

        # Hydrophobic atoms in plan
        for atom_in_plan in in_between_planes:
            if atom_in_plan.is_hydrophobic : 
                n_hydrophobe_in_plan+=1
        # When no more atoms between the two planes, exiting the function, we stop the exploring on this side of the axis
        if len(in_between_planes) == 0:
            return (False)
        
        # Computing the relative hydrophobicity of the selected amino_acids : to maximise
        hydrophobicity = (nb_hydrophile_out_of_plan/n_total_hydrophile) + (n_hydrophobe_in_plan/n_total_hydrophobic)
        if  hydrophobicity > ref.best_hydrophobicity :
            # Updating the "best" match
            ref.best_hydrophobicity = hydrophobicity
            ref.plane1 = copy.deepcopy(self.plane1)
            ref.plane2 = copy.deepcopy(self.plane2)
            return True
        else:
            # If its not better
            return False
    

    def find_tm_segment(self, protein):
        """
        Find and write the Transmembrane segments in an output file. 
        
        Parameters
        ----------
        protein : Protein
            Protein studied.
        
        Returns
        -------
        None
        """
        with open("../results/results_tm_segments_"+protein.name + ".txt", "w") as f_out:
            in_between_planes = [] # List of residues contained between the two planes
            for aa in (protein. amino_acid_sequence):
                if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)) or (self.plane2.is_below(
                        aa.point) and self.plane1.is_above(aa.point)):
                    in_between_planes.append(aa)
            
            tm = [] # List of transmembranes segments
            for aa in in_between_planes:
                if len(tm) ==0:
                    tm.append(aa.id)
                else:
                    # Getting the last one added
                    last_added = tm[(len(tm)-1)]
                    if last_added + 1 == aa.id : 
                        tm.append(aa.id)
                    else : 
                        # Writing in the output file
                        f_out.write(f"Transmembrane segment from residue {min(tm)} to {max(tm)}\n")
                        tm = []
                        tm.append(aa.id)

def find_points(n_points, center_coordinates):
    """
    Generates points equally distributed on a demi sphere using staff and kuijlaars method. 
        
    Parameters
    ----------
    n_points : int
        Number of points to generate. 
    center_coordinates : Point
        Coordinates of the mass center. 
        
    Returns
    -------
    points : List <Point>
        A list of points equally distributed on a demi sphere. 
    """
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
    # Get only the points from half circle :
    above_x_axis_points = [point for point in points if point.get_z() > center_coordinates.get_z()]
    return above_x_axis_points


def find_director_vector(point: Point, center_coordinate: Point):
    """Find the vector that goes through the two points = the normal vector of the plan"""
    # The normal vector of the plan
    return Vector(center_coordinate.get_x() - point.get_x(), center_coordinate.get_y() - point.get_y(),
                  center_coordinate.get_z() - point.get_z())




