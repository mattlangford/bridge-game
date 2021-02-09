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

def compute_stress(u):
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
    inv_det_j = 1 / det_j

    du_dx = inv_det_j * (y(2, 3) * (u[0] - u[4]) - y(1, 3) * (u[2] - u[4]))
    du_dy = inv_det_j * (-x(2, 3) * (u[0] - u[4]) + x(1, 3) * (u[2] - u[4]))
    dv_dx = inv_det_j * (y(2, 3) * (u[1] - u[5]) - y(1, 3) * (u[3] - u[5]))
    dv_dy = inv_det_j * (-x(2, 3) * (u[1] - u[5]) + x(1, 3) * (u[3] - u[5]))

    print ("Simple Strain:", np.array([
        du_dx,
        dv_dy,
        du_dy + dv_dx
    ]))

    print ("Green Strain:", np.array([
        du_dx + 0.5 * (du_dx ** 2 + dv_dx ** 2),
        dv_dy + 0.5 * (du_dx ** 2 + dv_dx ** 2),
        (du_dy + dv_dx) * 0.5 + 0.5 * (du_dx * du_dy + dv_dx * dv_dy)
    ]))

    print ("Almansi Strain:",  np.array([
        du_dx - 0.5 * (du_dx ** 2 + dv_dx ** 2),
        dv_dy - 0.5 * (du_dy ** 2 + dv_dy ** 2),
        0.5 * (du_dy + dv_dx - (du_dx * du_dy + dv_dx * dv_dy))
    ]))

u = translate(rotate(coords, np.radians(45)), [3, 1]) - coords
print (generate_b(triangle).dot(u))
compute_stress(u)
draw_coords(coords, fill=False, color="red")
draw_coords(rotate(coords, np.radians(1)), fill=False, color="blue")
# plt.show()
