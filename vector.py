import math

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

    def __add__(self, other : 'Vector') -> 'Vector':
        '''Add two vectors element-wise.

        Args:
            other (Vector): Another Vector to add to this vector.

        Returns:
            Vector: A new Vector representing the element-wise sum of the two vectors.
        '''
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)


    def __sub__(self, other : 'Vector') -> 'Vector':
        '''Subtract another vector element-wise from this vector.

        Args:
            other (Vector): Another Vector to subtract from this vector.

        Returns:
            Vector: A new Vector representing the element-wise difference of the two vectors.
        '''
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)


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
        c = Vector(self.y*b.z - self.z*b.y,
             self.z*b.x - self.x*b.z,
             self.x*b.y - self.y*b.x)
        return c

    def magnitude(self):
        '''Calculate the magnitude (length) of the vector.

        Returns:
            float: The magnitude (length) of the vector.
        '''
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        
    def calculate_projection(self, b : 'Vector') -> 'Vector':
        '''Calculates the projection of this vector onto vector b.
        
        Args: 
            (a1,a2,a3) (Vector) : x,y,z of the tip of first vector, assumes vector origin is (0,0,0)
            (b1,b2,b3) (Vector) : x,y,z of the tip of second vector, assumes vector origin is (0,0,0)
        Returns: 
            proj_vector (Vector) : x,y,z of vector resulting from proj_b(a)
        '''
        a_dot_product_b = self.x*b.x+self.y*b.y+self.z*b.z
        b_magnitude = b.magnitude()
        normalize_factor = a_dot_product_b / (b_magnitude**2)
        proj_vector = Vector(normalize_factor*b.x, normalize_factor*b.y, normalize_factor*b.z)
        return proj_vector
    
   
