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


class Plane:
    def __init__(self, a,b,c,d):
        self.a=a
        self.b=b
        self.c=c
        self.d=d

    def __init__(self, point:Point, normal:Vector):
        self.point=point
        self.normal=normal
    
    def plot_plane(self):
        # d of the plan equation
        d = -self.point.dot_product(self.normal)
        xx, yy = np.meshgrid(range(10), range(10))

        # calculate corresponding z
        z = (-self.normal.get_x() * xx - self.normal.get_y() * yy - d) * 1. /self.normal.get_z()
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(xx, yy, z, alpha=0.2)
        plt.show()



# TODO: faire que le demi cercle
def find_points(n_points, center_coordinates):
    points = []
    for k in range(1,n_points+1):
        h = -1 + (2 * (k-1)/(n_points - 1))
        theta = math.acos(h)
        if k==1 or k==n_points:
            phi = 0
        else :
            phi = (phi + (3.6/math.sqrt(n_points) * (1/math.sqrt(1-h*h)))) % (2*math.pi)
        x = math.sin(phi) * math.sin(theta)
        y = math.cos(theta)
        z = math.cos(phi) * math.sin(theta)
        points.append((x,y,z))
    return points

def find_director_vector(point:Point, center_coordinate:Point):
    """Find the vector that goes through the two points = the normal vector of the plan"""
    # The normal vector of the plan
    return (center_coordinate.get_x()-point.get_x(), center_coordinate.get_y() - point.get_y(), center_coordinate.get_z() - point.get_z())

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
    
    directions = find_points(20, (0,0,0))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for d in directions:
        ax.scatter3D(d[0], d[1], d[2])
    plt.show()

    # Planes are defined by a point and a normal vector
    """point  = Point(1, 2, 3)
    normal = Vector(1,1,2)

    plane = Plane(point=point, normal=normal)
    plane.plot_plane()"""






