import numpy as np
from matplotlib import pyplot as plt

np.set_printoptions(suppress=True, edgeitems=10, linewidth=100000)

E = 3.7 * 1E6
v = 0.1

# Relates stress to strain
D = np.array([
    [1.0, v, 0.0],
    [v, 1.0, 0.0],
    [0.0, 0.0, 0.5 * (1 - v)]
])
D *= E / (1.0 - v * v)

coords = np.array([
    0, 0,
    1, 0,
    0, 1,
], dtype=np.float64)

triangle = [0, 1, 2, 3, 4, 5]
u = np.zeros_like(coords)

def flatten_2d_coords(coords):
    res = coords.flatten()
    return res.tolist() if len(res.shape) == 1 else res.tolist()[0]

def flat_to_2d_coords(coords):
    return np.reshape(coords, (-1, 2))

def generate_b(triangle):
    displacements = coords[triangle] + u[triangle]

    def x(i, j):
        i -= 1
        j -= 1
        return displacements[2 * i] - displacements[2 *j]
    def y(i, j):
        i -= 1
        j -= 1
        return displacements[2 * i + 1] - displacements[2 * j + 1]

    det_j = x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3)
    b = np.array([
        [y(2, 3), 0, y(3, 1), 0, y(1, 2), 0],
        [0, x(3, 2), 0, x(1, 3), 0, x(2, 1)],
        [x(3, 2), y(2, 3), x(1, 3), y(3, 1), x(2, 1), y(1, 2)],
    ])
    return b / det_j

def rotate(coords, theta):
    t = np.matrix([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    c = flat_to_2d_coords(coords)
    return flatten_2d_coords(t.dot(c.T).T)
def translate(coords, v):
    c = flat_to_2d_coords(coords)
    return flatten_2d_coords(c + v)

def draw_coords(points, *args, **kwargs):
    plt.fill((points[::2]), (points[1::2]), *args, **kwargs)

E0 = np.array([
    [coords[0], coords[2]],
    [coords[1], coords[3]]
])
E0_inv = np.linalg.inv(E0)

def compute_strain(u):
    d = coords[triangle] + u[triangle]
    E = np.array([
        [d[0], d[2]],
        [d[1], d[3]]
    ])

    F = E0_inv.dot(E)
    strain = 0.5 * (F.T.dot(F) - np.eye(2))
    print (strain[0][0], strain[0][1], strain[1][1]

u = translate(rotate(coords, np.radians(45)), [3, 1]) - coords
print (generate_b(triangle).dot(u))
compute_stress(u)
draw_coords(coords, fill=False, color="red")
draw_coords(rotate(coords, np.radians(1)), fill=False, color="blue")
# plt.show()
