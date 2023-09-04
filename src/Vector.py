import matplotlib.pyplot  as plt
import numpy as np
import math

#problème avec matplolib dans le venv
class Vector:
    def __init__(self, x,y,z):
        self.x=x
        self.y=y
        self.z=z

# TODO: Mieux encapsuler avec Point et tout
# TODO: faire que le demi cercle
def find_points(n_points, center_coordinates):
    x,y,z = center_coordinates
    points = []
    for k in range(1,n_points+1):
        h = -1 + (2 * (k-1)/(n_points - 1))
        theta = math.acos(h)
        if k==1 or k==n_points:
            phi = 0
        else :
            phi = (phi + (1.8/math.sqrt(n_points) * (1/math.sqrt(1-h*h)))) % (2*math.pi)
        x = math.sin(phi) * math.sin(theta)
        y = math.cos(theta)
        z = math.cos(phi) * math.sin(theta)
        points.append((x,y,z))
    return points

def find_director_vector(point, center_coordinate):
    """Find the vector that goes through the two points"""
    x1,y1,z1 = point
    x2,y2,z2 = center_coordinate
    return (x2-x1, y2 - y1, z2 - z1)

def find_plane_equation(point, center_coordinate):
    

def find_normal_plan(director_vector):
    """Returns 2 normal vectors making a plan orthogonal to the director vector"""
    x,y,z=director_vector
    return ( np.array([-y,x,0]) , np.array([-z,0,x]))
  

if __name__ == '__main__':
    print("Test vector")
    directions = find_points(20, (0,0,0))
  
    # Plot points check
    """fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for d in directions:
        ax.scatter3D(d[0], d[1], d[2])
    plt.show()"""





