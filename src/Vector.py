import matplotlib.pyplot  as plt
import numpy as np
import math
import copy 

#problème avec matplolib dans le venv
class Vector:
    def __init__(self, x,y,z):
        self.coordinates = np.array([x,y,z])
    
    def get_coordinates(self):
        return self.coordinates
    
    def get_x(self):
        return self.coordinates[0]

    def get_y(self):
        return self.coordinates[1]

    def get_z(self):
        return self.coordinates[2]
    
    def norma(self):
        return (math.sqrt(self.x**2 + self.y **2 + self.z **2))



class Point:
    def __init__(self, x,y,z):
        self.coordinates = np.array([x,y,z])
    
    def __str__(self):
        return f"Point({self.coordinates[0]:.3},{self.coordinates[1]:.3},{self.coordinates[2]:.3})"
    
    def dot_product(self, normal):
        return (-self.coordinates.dot(normal.get_coordinates()))
    def get_x(self):
        return self.coordinates[0]

    def get_y(self):
        return self.coordinates[1]

    def get_z(self):
        return self.coordinates[2]
    def get_coordinates(self):
        return self.coordinates
    


class Plane:
    def __init__(self, point:Point, normal:Vector):
        self.point=point
        self.normal=normal
        self.a = normal.get_coordinates()[0]
        self.b = normal.get_coordinates()[1]
        self.c = normal.get_coordinates()[2]
        #self.d = -self.point.dot_product(self.normal)
        self.d = -self.a * self.point.get_x() - self.b * self.point.get_y() - self.c * self.point.get_z() # https://mathworld.wolfram.com/Plane.html

    def __str__(self):
        return f"Plane equation is : {self.a}x + {self.b}y +{self.c}z + {self.d} = 0\n"
    
    
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
    #TODO: Ces 2 fonctions
    def is_above(self,point):
        # Return true if the point is located over the plane (self)
        return True if (self.a * point.get_x() +self.b * point.get_y() + self.c * point.get_z() + self.d ) > 0  else False
    
    def is_below(self, point):
        # Return true if the point is located under the plane (self)
        return True if (self.a * point.get_x() +self.b * point.get_y() + self.c * point.get_z() + self.d )   < 0 else False

    # TODO: Faire ca de manière plus efficace
    def generate_points(self):
        coordinates = []
        for x in range(10):
            for y in range(10):
                for z in range(10):
                    # SI eqation de plan == 0 : alors le point appartient au plan
                    if round(self.a * x) + round(self.b * y) + round(self.c * z) + round(self.d) == 0:
                        coordinates.append(Point(x,y,z))
        return coordinates

class Axis:
    def __init__(self,p1,p2):
        self.plane1 = p1
        self.plane2 = p2
        self.best_number_hits = 0
        self.best_number_aa = 0
    
    def __str__(self):
        return (f"AXIS with best nb hits has {self.best_number_hits}")
    
    def explore_axe(self,amino_acid_sequence): 
        # TODO: Quand la tranche qu'on regarde ne contient aucun Calpha => BREAk
        #Je crois que c'est ce qui est fait vu qu'on retourne faux des qu'on a que des trucs qui correspondent pas : a voir ce qui se passe si on tombe sur des truca mauvais de base
        in_between_planes = []
        number_atoms_in_between = 0
        number_atoms_hits = 0
        for aa in (amino_acid_sequence):
            if (self.plane1.is_below(aa.point) and self.plane2.is_above(aa.point)) or self.plane2.is_below(aa.point) and self.plane1.is_above(aa.point):
                number_atoms_in_between+=1
                if aa.asa > 0.30 :#and aa.is_hydrophobic == True
                    in_between_planes.append(aa.id)
                    number_atoms_hits+=1
                    #print("ca is under plane a and is over p2", aa)
        if number_atoms_in_between == 0:
            print("No more atoms in between")
            return False
        """print("----------------------------")
        for a in in_between_planes:
            print(a,end="+")"""
        #print(f"Hit ratio {number_atoms_hits} \t {number_atoms_in_between} = {number_atoms_hits/number_atoms_in_between}")
        if number_atoms_hits > self.best_number_hits:
            # Updating the "best" match
            self.best_number_hits=number_atoms_hits
            self.best_number_aa = number_atoms_in_between
        return True if len(in_between_planes)>0 else False


def plot_plane(plane1, plane2=None, point=None):
        X, Z = np.meshgrid(range(100), range(100))
        # Calculate the corresponding y values for the plane
        # mouais, plutot faire en sorte quil ne soit jamais 0
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter3D(point.get_x(), point.get_y(), point.get_z(), color="red")
        if point is not None and plane1.d !=0 :
            Y1 = (-plane1.a * X - plane1.c * Z - plane1.d) / plane1.b
            ax.plot_surface(X, Y1, Z, alpha=0.2, color = "green")

        if plane2 is not None and plane2.d !=0:
            Y2 = (-plane2.a * X - plane2.c * Z - plane2.d) / plane2.b
            ax.plot_surface(X, Y2, Z, color = 'red', alpha=0.2)
        plt.show()

# TODO: faire que le demi cercle
def find_points(n_points, center_coordinates):
    points = []
    for k in range(1,n_points+1):# n-points pour éviter division par 0 
        h = -1 + (2 * (k-1)/(n_points - 1))
        theta = math.acos(h)
        if k==1 or k==n_points:
            phi = 0
        else :
            phi = (phi + (3.6/math.sqrt(n_points) * (1/math.sqrt(1-h*h)))) % (2*math.pi)
        # Centering on the mass center
        x = center_coordinates.get_x() + math.sin(phi) * math.sin(theta)
        y = center_coordinates.get_y() + math.cos(theta)
        z = center_coordinates.get_z() + math.cos(phi) * math.sin(theta)
        points.append(Point(x,y,z))
    above_x_axis_points = [point for point in points if point.get_z() > center_coordinates.get_z()]
    return above_x_axis_points

def find_director_vector(point:Point, center_coordinate:Point):
    """Find the vector that goes through the two points = the normal vector of the plan"""
    # The normal vector of the plan
    return Vector(center_coordinate.get_x()-point.get_x(), center_coordinate.get_y() - point.get_y(), center_coordinate.get_z() - point.get_z())


def find_normal_plan(director_vector):
    """Returns 2 normal vectors making a plan orthogonal to the director vector"""
    x,y,z=director_vector
    return ( np.array([-y,x,0]) , np.array([-z,0,x]))

  

if __name__ == '__main__':
    # Plot points check
    mass_center = Point(3.19,37.1,36.2)
    directions = find_points(10, mass_center)
    """fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(mass_center.get_x(), mass_center.get_y(), mass_center.get_z(), color="red")
    print("Plotting the points on 3D")
    for d in directions:
        ax.scatter3D(d.get_x(), d.get_y(), d.get_z())
    plt.show()"""
    print("Calculating the planes... ")
    # Planes are defined by a point and a normal vector
    for d in directions:
        point  = d
        normal = find_director_vector(point=d,center_coordinate=mass_center)
        if normal.get_z() != 0:
            #print("Plotting the planes on 3D")
            plane = Plane(point=point, normal=normal)
            plane2 = plane.complementary(14)        
            plane.slide_plane(10) 
            #plane2.slide_plane(-60)
            # Plane = green
            # Plane2 = red
            plot_plane(plane1=plane,plane2 = plane2, point = Point(6.462 , 37.060 , 37.424))
        else:
            print("Vecteur directeur nul (z), pass --- ")
        






