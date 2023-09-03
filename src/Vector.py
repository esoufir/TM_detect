import matplotlib.pyplot  as plt
import numpy as np

class Vector:
    def __init__(self, x,y,z):
        self.x=x
        self.y=y
        self.z=z

# TODO: Mieux encapsuler avec Point et tout
def find_directions(n_directions_wanted):
    list_directions = None
    for i in range(n_directions_wanted):
        for a in range(1,10):#mouais
            for b in range(1,10):
                for c in range(1,10):
                    np.append(list_directions,np.array([a,b,c]))
    return list_directions
    
if __name__ == '__main__':
    print("Test vector")
    directions = find_directions(5)
    fig = plt.figure()
    """ax = fig.add_subplot(111, projection='3d')
    
    for d in directions:
        ax.quiver(0, 0, 0, d[0], d[1], d[2], arrow_length_ratio=0.1)  
    plt.show()"""
