import numpy as np
from matplotlib import pyplot as plt
import sys

np.set_printoptions(suppress=True, edgeitems=10, linewidth=100000)

E = 3.7 * 1E6
v = 0.1
m = 50 * 3.1
c = 0.0

# Relates stress to strain
D = np.array([
    [1.0, v, 0.0],
    [v, 1.0, 0.0],
    [0.0, 0.0, 0.5 * (1 - v)]
])
D *= E / (1.0 - v * v)

D = np.array([
    [2 * v + E, E, 0.0],
    [E, 2 * v + E, 0.0 ],
    [0.0, 0.0, v]
])

fps = 2000.
dt = 1 / fps
dt2 = dt * dt

d = 1
coords = np.array([
    0, 0,
    0, 1,
    1, 0,
    1, 1,
    2, 0,
    2, 1,
    3, 0,
    3, 1,
    # 1, 2,
    # 1, 3,
    # 2, 3,
    # 2, 2,
    # 2, 1,
    # 2, 0,
    # 3, 0,
    # 3, 1,
    # 3, 2,
    # 3, 3,
    # 4, 3,
    # 4, 2,
    # 4, 1,
    # 4, 0,
], dtype=np.float64)

C = c * np.eye(len(coords))

fixed = np.array([
    1, 1,
    1, 1,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
    # 0, 0,
])
assert len(fixed) == len(coords), f"{len(fixed)} != {len(coords)}"

def get_triangle(p0, p1, p2):
    return [
        2 * p0, 2 * p0 + 1,
        2 * p1, 2 * p1 + 1,
        2 * p2, 2 * p2 + 1
    ]

def get_triangles(from_i, to_i):
    l = []
    for i in range(from_i, to_i):
        l.append(get_triangle(i, i + 1, i + 2))
    return l

triangles = [
    get_triangle(0, 1, 2),
    get_triangle(3, 2, 1),
    get_triangle(2, 3, 4),
    get_triangle(5, 4, 3),
    get_triangle(4, 5, 6),
    get_triangle(7, 6, 5),
    # get_triangle(5, 8, 7),
    # get_triangle(4, 8, 9),
    # get_triangle(3, 9, 10),
    # get_triangle(10, 9, 12),
    # get_triangle(9, 8, 13),
    # get_triangle(8, 7, 14),
    # get_triangle(8, 14, 13),
    # get_triangle(9, 13, 12),
    # get_triangle(10, 12, 11),
    # get_triangle(11, 12, 17),
    # get_triangle(12, 13, 16),
    # get_triangle(13, 14, 15),
    # get_triangle(13, 15, 16),
    # get_triangle(12, 16, 17),
    # get_triangle(11, 12, 17),
    # get_triangle(11, 17, 18),
]
u = np.zeros_like(coords)
u_dot = np.zeros_like(u)
u_dot_dot = np.zeros_like(u)

def compute_M():
    triangle_volume = 1 * 1 * 1 / 2
    M = np.zeros((len(coords), len(coords)))
    for triangle in triangles:
        mass = triangle_volume * m / 6
        for pt in triangle:
            M[pt, pt] += mass
    return M
M = compute_M()

def compute_F(triangle):
    sign = 1 if triangles.index(triangle) % 2 == 0 else -1
    _u, _v = u[triangle][::2], u[triangle][1::2]

    u10 = _u[1] - _u[0]
    u20 = _u[2] - _u[0]
    v10 = _v[1] - _v[0]
    v20 = _v[2] - _v[0]
    F = np.array([
        [u10, u20],
        [v10, v20]])
    F /= d
    F += np.eye(2)
    F *= sign
    return F

def compute_F(triangle):
    sign = 1 if triangles.index(triangle) % 2 == 0 else -1
    c = coords[triangle]
    points = c + u[triangle]

    x = np.array([
        [c[0], c[1]],
        [c[2], c[3]],
        [c[4], c[5]],
    ])
    y = np.array([
        [points[0], points[1]],
        [points[2], points[3]],
        [points[4], points[5]],
    ])

    Y = np.array([y[1] - y[0], y[2] - y[0]])
    X = np.array([x[1] - x[0], x[2] - x[0]])
    F = Y.dot(np.linalg.inv(X))
    return F

def compute_normals(triangle):
    sign = 1 if triangles.index(triangle) % 2 == 0 else -1
    _u, _v = u[triangle][::2], u[triangle][1::2]

    _d = sign * d
    p10 = np.array([
        -(_d + _v[1] - _v[0]),
               _u[1] - _u[0]
    ])
    p20 = np.array([
               _v[2] - _v[0],
        -(_d + _u[2] - _u[0])
    ])
    p12 = np.array([
         _d + _v[1] - _v[2],
      -(-_d + _u[1] - _u[2]),
    ])

    return p10, p20, p12

def compute_forces(triangle):
    sign = 1 if triangles.index(triangle) % 2 == 0 else -1
    _u, _v = u[triangle][::2], u[triangle][1::2]

    F = compute_F(triangle)
    strain = 0.5 * (F.T.dot(F) - np.eye(2))
    stress = D.dot([strain[0, 0], strain[1, 1], strain[0, 1]])
    S = np.array([
        [stress[0], stress[2]],
        [stress[2], stress[1]],
    ])
    p10, p20, p12 = compute_normals(triangle)

    f0 = -0.5 * F.dot(S).dot(p20 + p10)
    f1 = -0.5 * F.dot(S).dot(p10 + p12)
    f2 = -0.5 * F.dot(S).dot(p12 + p20)
    return f0, f1, f2

# [2] Table 6.5 B
class UpdatedLagrangianMethod(object):
    def __init__(self):
        self.u = np.zeros_like(coords)
        self.u_v = np.zeros_like(coords)
        self.u_a = np.zeros_like(coords)

    def incremental_strains(self):
        _u = self.u[::2]
        _v = self.u[1::2]
        _x = coords[::2]
        _y = coords[1::2]

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        det_v = np.linalg.det(v)
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        du_dx = (1 / det_v) * (b0 * _u[0] + b1 * _u[1] + b2 * _u[2])
        du_dy = (1 / det_v) * (c0 * _u[0] + c1 * _u[1] + c2 * _u[2])
        dv_dx = (1 / det_v) * (b0 * _v[0] + b1 * _v[1] + b2 * _v[2])
        dv_dy = (1 / det_v) * (c0 * _v[0] + c1 * _v[1] + c2 * _v[2])

        return np.array([
            du_dx + 0.5 * (du_dx ** 2 + dv_dx ** 2),
            dv_dy + 0.5 * (du_dx ** 2 + dv_dx ** 2),
            0.5 * (du_dy + dv_dx) + 0.5 * (du_dx * du_dy + dv_dx * dv_dy)
        ])

    def linear_b(self):
        _x = coords[::2]
        _y = coords[1::2]

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        det_v = np.linalg.det(v)
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        return np.array([
            [b0,  0, b1,  0, b2,  0],
            [ 0, c0,  0, c1,  0, c2],
            [c0, b0, c1, b1, c2, b2]
        ])

    def nonlinear_b(self):
        _x = coords[::2]
        _y = coords[1::2]

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        det_v = np.linalg.det(v)
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        return np.array([
            [b0,  0, b1,  0, b2,  0],
            [c0,  0, c1,  0, c2,  0],
            [ 0, b0,  0, b1,  0, b2],
            [ 0, c0,  0, c1,  0, c2],
        ])

    def linear_k(self):
        b = self.linear_b()
        return 

def update():
    global u, u_dot, u_dot_dot
    forces = -9.8 * np.diag(M)

    for triangle in triangles:
        f0, f1, f2 = compute_forces(triangle)
        forces[triangle[0]] += f0[0]
        forces[triangle[1]] += f0[1]
        forces[triangle[2]] += f1[0]
        forces[triangle[3]] += f1[1]
        forces[triangle[4]] += f2[0]
        forces[triangle[5]] += f2[1]

    forces[fixed == 1] = 0
    #print (forces)

    u_dot_dot = forces / np.diag(M)
    u = u + dt * u_dot + 0.5 * dt2 * u_dot_dot
    u_dot = u_dot + dt * u_dot_dot
    print (u)

def draw_triangles(u):
    points = coords + u

    for i, triangle in enumerate(triangles):
        x0, y0, x1, y1, x2, y2 = points[triangle]
        plt.scatter((x0, x1, x2), (y0, y1, y2))

        f0, f1, f2 = compute_forces(triangle)
        f0 *= 1E-6
        f1 *= 1E-6
        f2 *= 1E-6

        _u = u[triangle][::2]
        _v = u[triangle][1::2]

        sign = 1 if i % 2 == 0 else -1
        p10, p20, p12 = compute_normals(triangle)

        plt.plot([x0, x1, x1, x2, x2, x0], [y0, y1, y1, y2, y2, y0])

        plt.plot([x0, x0 + f0[0]], [y0, y0 + f0[1]])
        plt.plot([x1, x1 + f1[0]], [y1, y1 + f1[1]])
        plt.plot([x2, x2 + f2[0]], [y2, y2 + f2[1]])

        # plt.plot([(x1 + x0) / 2.0, (x1 + x0) / 2.0 + p10[0]], [(y1 + y0) / 2.0, (y1 + y0) / 2.0 + p10[1]])
        # plt.plot([(x2 + x0) / 2.0, (x2 + x0) / 2.0 + p20[0]], [(y2 + y0) / 2.0, (y2 + y0) / 2.0 + p20[1]])
        # plt.plot([(x2 + x1) / 2.0, (x2 + x1) / 2.0 + p12[0]], [(y2 + y1) / 2.0, (y2 + y1) / 2.0 + p12[1]])


def draw(i, u):
    plt.clf()
    plt.gca().set_aspect('equal', adjustable='box')

    print (f"Saving frame {i}")
    plt.xlim(min(coords) - 2, max(coords) + 2)
    plt.ylim(min(coords) - 2, max(coords) + 2)
    draw_triangles(u)

    if isinstance(i, int):
        plt.savefig(f"/tmp/{i:07}.png")
    else:
        plt.savefig(f"/tmp/{i}.png")

if __name__ == "__main__":
    draw("00_test", u)

    max_i = 5000
    if len(sys.argv) == 2:
        max_i = int(sys.argv[1])

    for i in range(max_i):
        print (f"Generating frame {i}")
        update()

        if i % int(1 / 60 * fps) == 0:
            draw(i, u)

        i += 1
