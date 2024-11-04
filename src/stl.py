# Naive STL generation code 
# When run outputs STL ASCII for a hardcoded box with some simple braille-like dots to the standard output
import numpy as np

def solid(name, facets):
    return '\n'.join([
        f'solid {name}',
        *facets,
        'endsolid'
    ])

def facet(a, b, c, n):
    return '\n'.join([
        f'facet normal {n[0]} {n[1]} {n[2]}', 
        'outer loop', 
        f'vertex {a[0]} {a[1]} {a[2]}', 
        f'vertex {b[0]} {b[1]} {b[2]}', 
        f'vertex {c[0]} {c[1]} {c[2]}', 
        'endloop', 
        'endfacet', 
    ])


def create_box_triangles(width, length, height, center):
    corner_a = center - np.array([width/2, length/2, height/2])
    corner_b = corner_a + np.array([width, length, height])

    dx = np.array([width, 0, 0])
    dy = np.array([0, length, 0])
    dz = np.array([0, 0, height])

    return [
        # Bottom 
        [
            corner_a,
            corner_a + dx,
            corner_a + dy,
        ],
        [
            corner_a + dx,
            corner_a + dx + dy,
            corner_a + dy,
        ],
        # Left
        [
            corner_a + dz,
            corner_a,
            corner_a + dy,
        ],
        [
            corner_a + dz + dy,
            corner_a + dz,
            corner_a + dy,
        ],
        # Right
        [
            corner_b,
            corner_b - dz,
            corner_b - dy,
        ],
        [
            corner_b - dz,
            corner_b - dz - dy,
            corner_b - dy,
        ],
        # Front
        [
            corner_a,
            corner_a + dz,
            corner_a + dx,
        ],
        [
            corner_a + dz,
            corner_a + dz + dx,
            corner_a + dx,
        ],
        # Back
        [
            corner_b - dz,
            corner_b,
            corner_b - dx,
        ],
        [
            corner_b - dz - dx,
            corner_b - dz,
            corner_b - dx,
        ],
        # Top 
        [
            corner_b,
            corner_b - dy,
            corner_b - dx,
        ],
        [
            corner_b - dx - dy,
            corner_b - dx,
            corner_b - dy,
        ],
    ]

def rotate_point_around_axis(a, b, c, theta):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    a_translate = a - b
    axis_vector = c - b
    u = axis_vector / np.linalg.norm(axis_vector)

    # Rodrigues' rotation formula
    a_rotated = a_translate * cos_theta + np.cross(u, a_translate) * sin_theta + u * np.dot(u, a_translate) * (1 - cos_theta)
    
    return a_rotated + b

# Given points which form a circle, create bottom and top triangles of a strip with angle towards center 
# The height is before the rotation TODO: Make height be the height after rotation for ease of computing the overall height of the braille dot
def make_strip(points, height, angleRad):
    assert len(points) > 2

    bottom_triangles = []
    for i in range(len(points)):
        point_a = points[i-1]
        point_b = points[i]
        length = point_b - point_a
        point_c = point_a + length * 0.5 + np.array([0, 0, height])
        point_c = rotate_point_around_axis(point_c, point_a, point_b, angleRad) 
        bottom_triangles.append([point_a, point_b, point_c])

    top_triangles = []
    for i in range(0, len(bottom_triangles)):
        point_a = bottom_triangles[i-1][2]
        point_b = bottom_triangles[i][0]
        point_c = bottom_triangles[i][2]
        top_triangles.append([point_a, point_b, point_c])

    return (bottom_triangles, top_triangles)

# Creates the base circle of points for the overall structure
def create_circle_points(n, length):
    angle = 360 / n

    theta = np.radians(angle)
    # Rotates around the z-axis
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])

    # Initialize base points
    start_position = np.array([0, 0, 0])
    points = [start_position]
    
    current_position = start_position
    current_direction = np.array([1, 0, 0])
    for _ in range(n-1):
        next_position = current_position + current_direction * length
        points.append(next_position)

        current_position = next_position
        current_direction = current_direction.dot(rotation_matrix)
    return points

# Assumes points form a circle and that points are equally distant such that the middle of opposite points is the center
def compute_center(points):
    n = len(points)
    assert n % 2 == 0
    a = points[0]
    b = points[n//2]
    diff = b-a
    return a + 0.5 * diff 

def create_circle_triangles(number_triangles, diameter, height, center):
    assert number_triangles % 2 == 0  # The way we calculate the side lengths makes it necessary to have an even number of points to ensure there are two points on the circle which are opposite each other
    assert number_triangles > 2 # Need to form at least a triangle to have some notion of something circle-like
    length = calculate_side_length(number_triangles, diameter)
    points = create_circle_points(number_triangles, length)
    point_center = compute_center(points)

    # Make strips to give curve
    strip_height = length * 0.8660254037844386 # Magic number is sqrt(3)/2 
    triangles = []
    angles = [20, 25, 30, 35, 40,  75, 80, 85, ]
    # angles = [30]

    for angle in angles:
        (bottom_triangles, top_triangles) = make_strip(points, strip_height, np.radians(angle))
        # (bottom_triangles, top_triangles) = make_strip(points, height, np.radians(angle))
        points = get_top_points(bottom_triangles)
        triangles.extend(bottom_triangles)
        triangles.extend(top_triangles)

    # Make Top
    total_height = points[-1][2] # Based on y coordinate of the last point. This assumes the points 
    for i in range(len(points)):
        a = points[i-1]
        b = points[i]
        length = b - a
        triangles.append([a, b, point_center + np.array([0, 0, total_height])])

    triangles = list(map(lambda t: t - point_center + center, triangles))
    
    return triangles

# Assumes top points are stored as the 3rd point for bottom triangles.
def get_top_points(triangles):
    return list(map(lambda t: t[2], triangles))
    
# Calculates specific angles used to calculate the distance between points in the base circle to achieve a specific circle diameter
def _calculate_angles(points, n):
    a = points[0]
    b = points[n//2]
    c = points[(n//2)-1] # Probably gives same whether +1 or -1

    ab = b-a
    ac = c-a

    dot_product = np.dot(ab, ac)
    magnitude_ab = np.linalg.norm(ab)
    magnitude_ac = np.linalg.norm(ac)

    cos_theta = dot_product / (magnitude_ab * magnitude_ac)

    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    angle_rad = np.arccos(cos_theta)
    angle_deg = np.degrees(angle_rad)

    opposite_angle = (180 - (360 / n)) / 2
    return (angle_deg, opposite_angle, 180 - opposite_angle - angle_deg)

# Calculates distance between points such that n points forms a circle with the required diameter
def calculate_side_length(n, diameter):
    points = create_circle_points(n, 1)
    (angle_a, angle_b, angle_c) = _calculate_angles(points, n)
    A = np.radians(angle_a) # We want to know the side length of the opposite of this corner
    C = np.radians(angle_c) # We know the side length of the opposite of this
    return diameter * np.sin(A) / np.sin(C) 

# Represent braile as an 8 bit bitstring such that 10000000 is the braille pattern with only the top left dot and 00000001 is the pattern with only the bottom right dot
# Using this representation calculate the different offsets for placing the pattern given a start position for top left and spacing for row and column
def get_offsets_by_bitstring(start, braille_bit, spacing_row, spacing_column):
    result = []
    i = 0
    1 << i
    for i in range(8):
        if (braille_bit & (1 << (7 - i))) > 0:
            row = i // 2
            col = i % 2
            result.append(start +  spacing_row * row + spacing_column * col)
    return result


# Hard coded 3D description of a box with 3 braille patterns on top
def main():
    zero = np.array([0, 0, 0]) # TODO: Calculate normal for triangles

    box_triangles = create_box_triangles(width=5, length=3, height=0.5, center=np.array([0, 0, 0.25]))
    box_facets = list(map(lambda t: facet(*t, zero), box_triangles))

    start = np.array([-2, 1, 0.5])
    spacing = 0.4
    offsets = get_offsets_by_bitstring(start, 0b11111111, spacing_row=np.array([0, -spacing, 0]), spacing_column=np.array([spacing, 0, 0]))

    circle_facets = []
    for offset in offsets:
        circle_triangles = 48 
        circle_diameter = 0.2
        circle_triangles = create_circle_triangles(number_triangles=circle_triangles, diameter=circle_diameter, height=5, center=offset)
        circle_facets += list(map(lambda t: facet(*t, zero), circle_triangles))

    start = np.array([-0.5, 1, 0.5])
    spacing = 0.4
    offsets = get_offsets_by_bitstring(start, 0b10101011, spacing_row=np.array([0, -spacing, 0]), spacing_column=np.array([spacing, 0, 0]))
    for offset in offsets:
        circle_triangles = 48 
        circle_diameter = 0.2
        circle_triangles = create_circle_triangles(number_triangles=circle_triangles, diameter=circle_diameter, height=5, center=offset)
        circle_facets += list(map(lambda t: facet(*t, zero), circle_triangles))

    start = np.array([0.5, 1, 0.5])
    spacing = 0.4
    offsets = get_offsets_by_bitstring(start, 0b01010100, spacing_row=np.array([0, -spacing, 0]), spacing_column=np.array([spacing, 0, 0]))
    for offset in offsets:
        circle_triangles = 48 
        circle_diameter = 0.2
        circle_triangles = create_circle_triangles(number_triangles=circle_triangles, diameter=circle_diameter, height=5, center=offset)
        circle_facets += list(map(lambda t: facet(*t, zero), circle_triangles))

    print(solid(f'my_circle{circle_triangles}', box_facets + circle_facets))

if __name__ == "__main__":
    main()
