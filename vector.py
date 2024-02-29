import math
import numpy as np

class Vector:
    '''Represents a 3D vector with x, y, and z components.

    This class defines a data type 'Vector' that represents a 3D vector with
    x, y, and z components, assuming the vector butt at the origin (0,0,0)

    Attributes:
        x (float): The x-component of the vector.
        y (float): The y-component of the vector.
        z (float): The z-component of the vector.
    '''
    def __init__(self, x: float, y: float, z: float) -> None:
        '''Initialize a Vector instance.

        Args:
            x (float): The x-component of the vector.
            y (float): The y-component of the vector.
            z (float): The z-component of the vector.
        '''
        self.x = x
        self.y = y
        self.z = z
        self.components = np.array([x, y, z])

    def __add__(self, other : 'Vector') -> 'Vector':
        '''Add two vectors element-wise.

        Args:
            other (Vector): Another Vector to add to this vector.

        Returns:
            Vector: A new Vector representing the element-wise sum of the two vectors.
        '''
        return Vector(*(self.components + other.components))


    def __sub__(self, other : 'Vector') -> 'Vector':
        '''Subtract another vector element-wise from this vector.

        Args:
            other (Vector): Another Vector to subtract from this vector.

        Returns:
            Vector: A new Vector representing the element-wise difference of the two vectors.
        '''
        return Vector(*(self.components - other.components))

    def __eq__(self, other: "Vector") -> bool:
        '''Checks if two Vectors are equal
        
        Args:
            other (Vector) : Vector to check equality to

        Returns:    
            equal (bool) : True if vectors are equal, False otherwise
        '''
        return np.all(self.components == other.components)
        
    def calculate_cross_product(self, b : 'Vector') -> 'Vector':
        '''Calculates the cross product of 2 vectors
        
        Calculates the cross product of 2 vectors, which is the unit vector that is 
        perpendicular to both vectors.
        
        Args:
            (a1,a2,a3) (Vector) : x,y,z of the tip of first vector, assumes vector origin is (0,0,0)
            (b1,b2,b3) (Vector) : x,y,z of the tip of second vector, assumes vector origin is (0,0,0)
        Returns:
            c (Vector) : x,y,z of vector resulting from a x b
        '''
        c = np.cross(self.components, b.components)
        return Vector(*c)

    def magnitude(self):
        '''Calculate the magnitude (length) of the vector.

        Returns:
            float: The magnitude (length) of the vector.
        '''
        return np.linalg.norm(self.components)
        
    def calculate_projection(self, b : 'Vector') -> 'Vector':
        '''Calculates the projection of this vector onto vector b.
        
        Args: 
            (a1,a2,a3) (Vector) : x,y,z of the tip of first vector, assumes vector origin is (0,0,0)
            (b1,b2,b3) (Vector) : x,y,z of the tip of second vector, assumes vector origin is (0,0,0)
        Returns: 
            proj_vector (Vector) : x,y,z of vector resulting from proj_b(a)
        '''
        b_np = b.components
        projected_vector = (np.dot(self.components, b_np) / np.dot(b_np, b_np)) * b_np
        return Vector(*projected_vector)

    def __str__(self):
        '''Redefinition of printing for Vectors

        Redefines the output of print(Vector()) to display the x,y,z attributes
        '''
        x,y,z = self.components
        return "[ " + str(x) + "\n  " + str(y) + "\n  " + str(z) + " ]"
    
    def scale(self, a):
        '''Scale the self vector by a scalar a'''
        scaled_vector = a * self.components
        return Vector(*scaled_vector)
    
if __name__ == "__main__":
    assert (Vector(1,2,3) + Vector(3,2,1) == Vector(4,4,4))
    assert (Vector(1,2,3).y == 2)
    assert (Vector(1,2,3) - Vector(-1,0,4) == Vector(2,2,-1))
    assert (Vector(1,0,0).calculate_cross_product(Vector(0,1,0)) == Vector(0,0,1))
    assert (Vector(1,0,0).calculate_cross_product(Vector(0,0,0)) == Vector(0,0,0))
    assert (Vector(1,0,0).magnitude() == 1)
    assert (Vector(0,3,4).magnitude() == 5)
    assert (Vector(-1,-2,-2).magnitude() == 3)
    assert (Vector(3,1,0).calculate_projection(Vector(1,0,0)) == Vector(3,0,0))
    assert (Vector(1,2,3).calculate_projection(Vector(0,0,0)) == Vector(0,0,0))