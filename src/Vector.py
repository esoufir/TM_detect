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
    
    def set_d(self, d):
        self.d=d
    
    def plot_plane(self, plane2=None, point=None):
        xx, yy = np.meshgrid(range(10), range(10))
        # calculate corresponding z
        if self.normal.get_z()!=0 : # mouais, plutot faire en sorte quil ne soit jamais 0
            z = (-self.normal.get_x() * xx - self.normal.get_y() * yy - self.d) * 1. /self.normal.get_z()
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.plot_surface(xx, yy, z, alpha=0.2, color="green")
            if point is not None:
                ax.scatter3D(point.get_x(), point.get_y(), point.get_z())
            if plane2 != None:
                z2 = (-plane2.normal.get_x() * xx - plane2.normal.get_y() * yy - plane2.d) * 1. /plane2.normal.get_z()
                ax.plot_surface(xx, yy, z2, color = 'red', alpha=0.2)
            plt.show()
        else :
            print("Division par 0 imposssible : ")
        
    def complementary(self, gap):
        # gap in Angstrom
        # Copying the object
        complementary_plane = copy.deepcopy(self)
        complementary_plane.set_d(complementary_plane.d - gap)
        return complementary_plane

    def slide_plane(self, sliding_window):
        # slide the plane 
        self.d += sliding_window
    
    #TODO: Ces 2 fonctions
    def is_above(self,point):
        # Return true if the point is located over the plane (self)
        print(self,point)
        print("Result is over",self.a * point.get_x() + self.d )
        return True if self.a * point.get_x() + self.d > point.get_z() else False
    
    def is_below(self, point):
        # Return true if the point is located under the plane (self)
        print(self,point)
        print("Result is under",self.a * point.get_x() + self.b * point.get_y() + self.c * point.get_z() + self.d )
        return True if self.a * point.get_x() + self.d < 0 else False

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
    return points

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

    directions = find_points(5, mass_center)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(mass_center.get_x(), mass_center.get_y(), mass_center.get_z(), color="red")
    print("Plotting the points on 3D")
    for d in directions:
        ax.scatter3D(d.get_x(), d.get_y(), d.get_z())
    plt.show()
    print("Calculating the planes... ")
    # Planes are defined by a point and a normal vector
    for d in directions:
        point  = d
        normal = find_director_vector(point=d,center_coordinate=mass_center)
        print("Plotting the planes on 3D")
        plane = Plane(point=point, normal=normal)
        plane2 = plane.complementary(200)        
        plane.slide_plane(10) 
        plane2.slide_plane(10)
        plane.plot_plane(plane2,Point(6.462 , 37.060 , 37.424))
        print("is under ? ",plane.is_below(Point(6.462 , 37.060 , 37.424)))
        print("is over ? ",plane2.is_above(Point(6.462 , 37.060 , 37.424)))
    
        






