import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy.optimize

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

fps = 30.
dt = 1 / fps
dt2 = dt * dt

thickness = 1

coords = np.array([
    0, 0,
    0, 1,
    1, 0,
    1, 1,
    2, 0,
    2, 1,
    3, 0,
    3, 1,
    4, 0,
    4, 1,
    5, 0,
    5, 1,
    6, 0,
    6, 1,
    7, 0,
    7, 1,
    8, 0,
    8, 1,
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
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
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

triangles = get_triangles(0, 16)
print (triangles)
assert max([max(triangle) for triangle in triangles]) < len(coords), f"Needs to be {int(len(coords) / 2 - 1)}"

def compute_M():
    triangle_volume = 1 * 1 * 1 / 2
    M = np.zeros((len(coords), len(coords)))
    for triangle in triangles:
        mass = triangle_volume * m / 6
        for pt in triangle:
            M[pt, pt] += mass
    return M
M = compute_M()


def area(triangle, u):
    displacements = coords[triangle] + u[triangle]

    def x(i, j):
        i -= 1
        j -= 1
        return displacements[2 * i] - displacements[2 * j]
    def y(i, j):
        i -= 1
        j -= 1
        return displacements[2 * i + 1] - displacements[2 * j + 1]

    det_j = x(1, 3) * y(2, 3) - y(1, 3) * x(2, 3)
    return 0.5 * abs(det_j)


def generate_nonfixed_matrix(matrix):
    coords_map = []
    for i, is_fixed in enumerate(fixed):
        if is_fixed:
            continue
        coords_map.append(i)

    result = np.zeros((len(coords_map), len(coords_map)))
    for row_to, row_from in enumerate(coords_map):
        for col_to, col_from in enumerate(coords_map):
            result[row_to, col_to] = matrix[row_from, col_from]

    return result


def generate_full_vector(nonfixed_vector, result_vector):
    coords_map = []
    for i, is_fixed in enumerate(fixed):
        if is_fixed:
            continue
        coords_map.append(i)

    for from_i, to_i in enumerate(coords_map):
        result_vector[to_i] = nonfixed_vector[from_i]

    return result_vector


def generate_nonfixed_vector(full_vector):
    coords_map = []
    for i, is_fixed in enumerate(fixed):
        if is_fixed:
            continue
        coords_map.append(i)

    result = np.zeros(len(coords_map))
    for to_i, from_i in enumerate(coords_map):
        result[to_i] = full_vector[from_i]
    return result

def gen_global_k(gen_local):
    K = np.zeros([len(coords), len(coords)])
    for i, triangle in enumerate(triangles):
        local_k = gen_local(triangle)
        for row in range(6):
            for col in range(6):
                K[triangle[col], triangle[row]] += local_k[row, col]
    return K

# [2] Table 6.5 B
class UpdatedLagrangianMethod(object):
    def __init__(self):
        self.u = np.zeros_like(coords)
        self.u_v = np.zeros_like(coords)
        self.u_a = np.zeros_like(coords)

    def nonlinear_strains(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2] + _u
        _y = coords[triangle][1::2] + _v

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

    def nonlinear_b(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2] + _u
        _y = coords[triangle][1::2] + _v

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        return np.array([
            [b0,  0, b1,  0, b2,  0],
            [c0,  0, c1,  0, c2,  0],
            [ 0, b0,  0, b1,  0, b2],
            [ 0, c0,  0, c1,  0, c2],
        ])

    def nonlinear_k(self, u):
        def impl(triangle):
            strains = self.nonlinear_strains(triangle, u)
            t11, t22, t12 = D.dot(strains)
            stress_matrix = np.array([
                [t11, t12,  0,   0],
                [t12, t22,  0,   0],
                [  0,  0, t11, t12],
                [  0,  0, t12, t22],
            ])
            b = self.nonlinear_b(triangle, u)
            return thickness * area(triangle, u) * b.T.dot(stress_matrix).dot(b)
        return gen_global_k(impl)

    def linear_strains(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2] + _u
        _y = coords[triangle][1::2] + _v

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
            du_dx,
            dv_dy,
            0.5 * (du_dy + dv_dx)
        ])

    def linear_b(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2] + _u
        _y = coords[triangle][1::2] + _v

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        return np.array([
            [b0,  0, b1,  0, b2,  0],
            [ 0, c0,  0, c1,  0, c2],
            [c0, b0, c1, b1, c2, b2]
        ])

    def linear_k(self, u):
        def impl(triangle):
            b = self.linear_b(triangle, u)
            return thickness * area(triangle, u) * b.T.dot(D).dot(b)
        return gen_global_k(impl)

    def stress_forces(self, u):
        k = self.linear_k(u) + self.nonlinear_k(u)
        return k.dot(u)

        forces = np.zeros_like(u)
        for triangle in triangles:
            strains = self.nonlinear_strains(triangle, u)
            stress = D.dot(strains)
            b = self.linear_b(triangle, u)
            force = b.T.dot(stress)

            for from_i, to_i in enumerate(triangle):
                forces[to_i] += force[from_i]

        return forces

    def iteration(self, new_u):
        gravity = np.zeros_like(self.u)
        gravity[1::2] = -9.8 * np.diag(M)[1::2]

        rhs = gravity - self.stress_forces(new_u)
        rhs -= M.dot((4 / dt2) * (new_u - self.u) - (4 / dt) * self.u_v - self.u_a)
        rhs -= C.dot((2 / dt) * (new_u - self.u) - self.u_v)

        print ("  Linear    :", self.linear_strains(triangles[-1], new_u))
        print ("  Non linear: ", self.nonlinear_strains(triangles[-1], new_u))
        # print (f"  Gravity: {generate_nonfixed_vector(gravity)}")
        # print (f"  Stress : {generate_nonfixed_vector(rhs - gravity)}")

        K = generate_nonfixed_matrix(self.linear_k(new_u) + self.nonlinear_k(new_u) + (4 / dt2) * M + (2 / dt) * C)
        du = np.linalg.inv(K).dot(generate_nonfixed_vector(rhs))
        new_u += generate_full_vector(du, np.zeros_like(self.u))

        # print (f"  Du     : {du}")

        # print (f" Residual: {np.linalg.norm(generate_nonfixed_vector(rhs))}")

        return new_u

    def update(self):
        new_u = np.copy(self.u)

        previous_u = np.zeros_like(new_u)
        for i in range(10):
            new_u = self.iteration(new_u)
            if (np.linalg.norm(new_u - previous_u) < 1E-6):
                print (f"Converged after {i} iterations")
                break
            previous_u = new_u
        else:
            print ("Failed to converge.")

        alpha = 0.5
        beta = 0.25 * (0.5 + alpha) ** 2.0

        new_u = generate_nonfixed_vector(new_u)
        u = generate_nonfixed_vector(self.u)
        u_v = generate_nonfixed_vector(self.u_v)
        u_a = generate_nonfixed_vector(self.u_a)

        n_1 = u_v + (1 - alpha) * dt * u_a
        n_2 = u + u_v * dt + (0.5 - beta) * dt * dt * u_a

        new_u_a = (new_u - n_2) / (beta * dt2)
        new_u_v = n_1 + alpha * new_u_a * dt

        self.u = generate_full_vector(new_u, np.zeros_like(self.u))
        self.u_v = generate_full_vector(new_u_v, np.zeros_like(self.u_v))
        self.u_a = generate_full_vector(new_u_a, np.zeros_like(self.u_a))


def draw_triangles(u):
    points = coords + u

    ul = UpdatedLagrangianMethod()
    ul.u = u

    for i, triangle in enumerate(triangles):
        x0, y0, x1, y1, x2, y2 = points[triangle]

        stresses = D.dot(ul.nonlinear_strains(triangle, u))
        stress = np.linalg.norm(stresses)

        def color(c):
            return max(min(c, 1.0), 0.0)

        max_stress = 50000
        red = color(stress / max_stress)
        green = color((max_stress - stress) / max_stress)

        plt.fill((x0, x1, x2), (y0, y1, y2), facecolor=(red, green, 0.0), edgecolor="black", linewidth=1, zorder=-10)

    plt.scatter(points[::2], points[1::2], color=["red" if i == 1 else "blue" for i in fixed[::2]])


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
    ul = UpdatedLagrangianMethod()
    draw("00_init", ul.u)

    max_i = 5000
    if len(sys.argv) == 2:
        max_i = int(sys.argv[1])

    for i in range(max_i):
        print (f"Generating frame {i}")
        ul.update()

        if i % int(1 / 30 * fps) == 0:
            draw(i, ul.u)

        i += 1
