import matplotlib.pyplot  as plt
import numpy as np
import math

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
        self.d = -self.point.dot_product(self.normal)
    
    def get_a(self):
        return self.a
    def get_b(self):
        return self.b
    def get_c(self):
        return self.c
    def get_d(self):
        return self.d
    
    def set_d(self, d):
        self.d=d
    
    def plot_plane(self, plane2=None):
        # d of the plan equation
        d = -self.point.dot_product(self.normal)
        xx, yy = np.meshgrid(range(10), range(10))

        # calculate corresponding z
        if self.normal.get_z()!=0 : # mouais, plutot faire en sorte quil ne soit jamais 0
            z = (-self.normal.get_x() * xx - self.normal.get_y() * yy - d) * 1. /self.normal.get_z()
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.plot_surface(xx, yy, z, alpha=0.2)
            if plane2 != None:
                print("INNNNNNNNNNNNNNNnnn")
                d2 = -plane2.point.dot_product(plane2.normal)
                z2 = (-plane2.normal.get_x() * xx - plane2.normal.get_y() * yy - d2) * 1. /plane2.normal.get_z()
                ax.plot_surface(xx, yy, z2, color = 'red', alpha=0.2)
                print(plane2.d)
                print(self.d)
            plt.show()
        else :
            print("Division par 0 imposssible")
        

    
    
    def complementary(self, gap):
        # gap in Angstrom
        complementary_plane = self
        print("COMPLEMENTARY PLANE", complementary_plane.d)
        complementary_plane.set_d(complementary_plane.get_d() - gap)
        return complementary_plane




# TODO: faire que le demi cercle
# TODO : Régler le probleme du centre de la sphère = centre de masse
def find_points(n_points, center_coordinates):
    points = []
    for k in range(1,n_points+1):
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
    return Point(center_coordinate.get_x()-point.get_x(), center_coordinate.get_y() - point.get_y(), center_coordinate.get_z() - point.get_z())

"""def find_plane_equation(normal_vector, point, x,y,z):
    # When knowing a point and the normal vector, the following formula gives the plan equation : 
    xA,yA,zA = point
    # Coordinates of the normal vector : encapsule moi ca pour pouvoir self.x...
    a,b,c = normal_vector
    return (a (x-xA) + b *(y-yA) + c *(z-zA))"""

def find_complemantary_plane(plane, scale):
    return


def find_normal_plan(director_vector):
    """Returns 2 normal vectors making a plan orthogonal to the director vector"""
    x,y,z=director_vector
    return ( np.array([-y,x,0]) , np.array([-z,0,x]))

  

if __name__ == '__main__':
    print("Test vector")
    # Plot points check
    
    mass_center = Point(1,1,1)

    directions = find_points(2, mass_center)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(1, 1, 1, color="red")
    print("Plotting the points on 3D")
    for d in directions:
        print(d.get_x(), d.get_y(), d.get_z())
        ax.scatter3D(d.get_x(), d.get_y(), d.get_z())
    plt.show()
    print("Calculating the planes... ")

    # Planes are defined by a point and a normal vector

    for d in directions:
        point  = d
        normal = find_director_vector(point=d,center_coordinate=mass_center)
        print("Plotting the planes on 3D")
        plane = Plane(point=point, normal=normal)
        plane2 = plane.complementary(14)
        print("PLANES", plane.d, plane2.d)
        
        plane.plot_plane(plane2)
        






