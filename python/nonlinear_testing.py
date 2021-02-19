import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy.optimize

np.set_printoptions(suppress=True, edgeitems=10, linewidth=100000)

E = 3.7 * 1E5
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

fps = 60.
dt = 1 / fps
dt2 = dt * dt

thickness = 1

coords = np.array([
    0, 0,
    0, 1,
    1, 0,
    1, 1,
    1, 2,
    1, 3,
    2, 3,
    2, 2,
    2, 1,
    2, 0,
    3, 0,
    3, 1,
    3, 2,
    3, 3,
    4, 3,
    4, 2,
    4, 1,
    4, 0,
], dtype=np.float64)

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
    p0 = p0 - 1
    p1 = p1 - 1
    p2 = p2 - 1
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
    get_triangle(1, 3, 4),
    get_triangle(1, 2, 4),
    get_triangle(3, 4, 9),
    get_triangle(4, 5, 8),
    get_triangle(5, 6, 7),
    get_triangle(5, 8, 7),
    get_triangle(4, 8, 9),
    get_triangle(3, 9, 10),
    get_triangle(10, 9, 12),
    get_triangle(9, 8, 13),
    get_triangle(8, 7, 14),
    get_triangle(8, 14, 13),
    get_triangle(9, 13, 12),
    get_triangle(10, 12, 11),
    get_triangle(11, 12, 17),
    get_triangle(12, 13, 16),
    get_triangle(13, 14, 15),
    get_triangle(13, 15, 16),
    get_triangle(12, 16, 17),
    get_triangle(11, 12, 17),
    get_triangle(11, 17, 18),
]
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

# triangles = get_triangles(0, 16)
# print (triangles)
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


def compute_C():
    return c * np.eye(len(coords))
C = compute_C()


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

def area0(triangle):
    return area(triangle, np.zeros_like(coords))

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
class TotalLagrangianMethod(object):
    def __init__(self):
        self.u = np.zeros_like(coords)
        self.v = np.zeros_like(coords)
        self.a = np.zeros_like(coords)

    def nonlinear_strains(self, triangle, u):
        # Incremental u, v
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        # u, v at the previous timestamp
        _u_t = self.u[triangle][::2]
        _v_t = self.u[triangle][1::2]
        # Initial coordinates
        _x = coords[triangle][::2]
        _y = coords[triangle][1::2]

        # We're trying to solve for a, b, c where
        #   u0 = a + b*x0 + c*y0
        #   u1 = a + b*x1 + c*y1
        #   u2 = a + b*x2 + c*y2
        # or in matrix form:
        #   [u0, u1, u2] = Q(x, y) * [a, b, c]
        #   [a, b, c] = Q(x, y)^-1 * [u0, u1, u2]
        Q = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        inv_Q = np.linalg.inv(Q)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_Q.flatten()

        # So now we can write the displacement functions as:
        #   N0(x, y) = h0(x, y) = (a0 + b0 * x0 + c0 * y0)
        #   N1(x, y) = h1(x, y) = (a1 + b1 * x1 + c1 * y1)
        #   N2(x, y) = h2(x, y) = (a2 + b2 * x2 + c2 * y2)
        # where:
        #   u(x, y) = N0(x, y) * u0 + N1(x, y) * u1 + N2(x, y) * u2
        #   v(x, y) = N0(x, y) * v0 + N1(x, y) * v1 + N2(x, y) * v2

        # These are for the "initial displacement" terms, since they don't change for changes in the incremental u
        du_dx_t = b0 * _u_t[0] + b1 * _u_t[1] + b2 * _u_t[2]
        du_dy_t = c0 * _u_t[0] + c1 * _u_t[1] + c2 * _u_t[2]
        dv_dx_t = b0 * _v_t[0] + b1 * _v_t[1] + b2 * _v_t[2]
        dv_dy_t = c0 * _v_t[0] + c1 * _v_t[1] + c2 * _v_t[2]

        du_dx = b0 * _u[0] + b1 * _u[1] + b2 * _u[2]
        du_dy = c0 * _u[0] + c1 * _u[1] + c2 * _u[2]
        dv_dx = b0 * _v[0] + b1 * _v[1] + b2 * _v[2]
        dv_dy = c0 * _v[0] + c1 * _v[1] + c2 * _v[2]

        # Initial Displacement terms
        init11 = du_dx_t * du_dx + dv_dx_t * dv_dx
        init22 = du_dy_t * du_dy + dv_dy_t * dv_dy
        init12 = 0.5 * (du_dx_t * du_dy + dv_dx_t * dv_dy + du_dy_t * du_dx + dv_dy_t * dv_dx)

        # Linear terms
        e11 = du_dx + init11
        e22 = dv_dy + init22
        e12 = 0.5 * (du_dy + dv_dx) + init12

        # Nonlinear terms
        n11 = 0.5 * (du_dx ** 2 + dv_dx ** 2)
        n22 = 0.5 * (du_dy ** 2 + dv_dy ** 2)
        n12 = 0.5 * (du_dx * du_dy + dv_dx * dv_dy)

        return np.array([
            e11 + n11,
            e22 + n22,
            2 * (e12 + n12)
        ])

    def nonlinear_b(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2]
        _y = coords[triangle][1::2]

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
            return thickness * area0(triangle) * b.T.dot(stress_matrix).dot(b)
        return gen_global_k(impl)

    def linear_b(self, triangle, u):
        _u = u[triangle][::2]
        _v = u[triangle][1::2]
        _x = coords[triangle][::2]
        _y = coords[triangle][1::2]

        v = np.array([
            [1, _x[0], _y[0]],
            [1, _x[1], _y[1]],
            [1, _x[2], _y[2]],
        ])
        inv_v = np.linalg.inv(v)
        a0, a1, a2, b0, b1, b2, c0, c1, c2 = inv_v.flatten()

        B0 = np.array([
            [b0,  0, b1,  0, b2,  0],
            [ 0, c0,  0, c1,  0, c2],
            [c0, b0, c1, b1, c2, b2]
        ])

        _u_t = self.u[triangle][::2]
        _v_t = self.u[triangle][1::2]

        l11 = b0 * _u_t[0] + b1 * _u_t[1] + b2 * _u_t[2]
        l22 = c0 * _v_t[0] + c1 * _v_t[1] + c2 * _v_t[2]
        l21 = b0 * _v_t[0] + b1 * _v_t[1] + b2 * _v_t[2]
        l12 = c0 * _u_t[0] + c1 * _u_t[1] + c2 * _u_t[2]

        B1 = np.array([
            [           l11 * b0,            l21 * b0,            l11 * b1,            l21 * b1,            l11 * b2,            l21 * b2],
            [           l12 * c0,            l22 * c0,            l12 * c1,            l22 * c1,            l12 * c2,            l22 * c2],
            [l11 * c0 + l12 * b0, l21 * c0 + l22 * b0, l11 * c1 + l12 * b1, l21 * c1 + l22 * b1, l11 * c2 + l12 * b2, l21 * c2 + l22 * b2],
        ])

        return B0 + B1


    def linear_k(self, u):
        def impl(triangle):
            b = self.linear_b(triangle, u)
            return thickness * area0(triangle) * b.T.dot(D).dot(b)
        return gen_global_k(impl)

    def stress_forces(self, u):
        return (self.linear_k(u) + self.nonlinear_k(u)).dot(u)

        forces = np.zeros_like(u)
        for triangle in triangles:
            strains = self.nonlinear_strains(triangle, u)
            stress = D.dot(strains)
            b = self.linear_b(triangle, u)
            force = thickness * area(triangle, u) * b.T.dot(stress)

            print ("      U: ", u[triangle])
            print ("  Force: ", force)
            print ("Strains: ", strains)
            print (" Stress: ", stress)
            print ("")

            for from_i, to_i in enumerate(triangle):
                if fixed[to_i]:
                    continue
                forces[to_i] += force[from_i]

        return forces

# [2] Table 6.5 B
class UpdatedLagrangianMethod(object):
    def __init__(self):
        self.u = np.zeros_like(coords)
        self.v = np.zeros_like(coords)
        self.a = np.zeros_like(coords)

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

        du_dx = b0 * _u[0] + b1 * _u[1] + b2 * _u[2]
        du_dy = c0 * _u[0] + c1 * _u[1] + c2 * _u[2]
        dv_dx = b0 * _v[0] + b1 * _v[1] + b2 * _v[2]
        dv_dy = c0 * _v[0] + c1 * _v[1] + c2 * _v[2]

        return np.array([
            du_dx + 0.5 * (du_dx ** 2 + dv_dx ** 2),
            dv_dy + 0.5 * (du_dx ** 2 + dv_dx ** 2),
            2 * (0.5 * (du_dy + dv_dx) + 0.5 * (du_dx * du_dy + dv_dx * dv_dy))
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
        return (self.linear_k(u) + self.nonlinear_k(u)).dot(u)

        forces = np.zeros_like(u)
        for triangle in triangles:
            strains = self.nonlinear_strains(triangle, u)
            stress = D.dot(strains)
            b = self.linear_b(triangle, u)
            force = thickness * area(triangle, u) * b.T.dot(stress)

            for from_i, to_i in enumerate(triangle):
                if fixed[to_i]:
                    continue
                forces[to_i] += force[from_i]

        return forces


def update_newmark(interface):
    M = generate_nonfixed_matrix(compute_M())
    C = generate_nonfixed_matrix(compute_C())

    u = generate_nonfixed_vector(interface.u)
    v = generate_nonfixed_vector(interface.v)
    a = generate_nonfixed_vector(interface.a)

    def iteration(next_u):
        gravity = np.zeros_like(u)
        gravity[1::2] = -9.8 * np.diag(M)[1::2]

        next_full_u = generate_full_vector(next_u, np.zeros(len(coords)))

        external_forces = gravity
        internal_forces = generate_nonfixed_vector(interface.stress_forces(next_full_u))

        rhs = external_forces - internal_forces
        rhs -= M.dot((4 / dt2) * (next_u - u) - (4 / dt) * v - a)
        rhs -= C.dot((2 / dt) * (next_u - u) - v)

        K = generate_nonfixed_matrix(interface.linear_k(next_full_u) + interface.nonlinear_k(next_full_u)) + (4 / dt2) * M + (2 / dt) * C
        du = np.linalg.inv(K).dot(rhs)

        initial_next_u = np.copy(next_u)
        next_u += du

        next_a = (4 / dt2) * (next_u - u) - (4 / dt) * v - a
        internal_forces = generate_nonfixed_vector(interface.stress_forces(generate_full_vector(next_u, np.zeros(len(coords)))))
        error = np.linalg.norm(external_forces - internal_forces - M.dot(next_a))

        return next_u, error

    next_u = np.zeros_like(u)

    for i in range(10):
        next_u, error = iteration(next_u)
        print (f"Iteration {i} error: {error}")
        if (error < 1E-2):
            print (f"Converged after {i + 1} iteration(s)")
            break
    else:
        raise Exception("Failed to converge.")

    next_a = (4 / dt2) * (next_u - u) - (4 / dt) * v - a
    next_v = v + dt / 2 * (a + next_a)

    interface.u = generate_full_vector(next_u, np.zeros_like(interface.u))
    interface.v = generate_full_vector(next_v, np.zeros_like(interface.v))
    interface.a = generate_full_vector(next_a, np.zeros_like(interface.a))

def update_bathe(interface):
    M = generate_nonfixed_matrix(compute_M())
    C = generate_nonfixed_matrix(compute_C())

    u = generate_nonfixed_vector(interface.u)
    v = generate_nonfixed_vector(interface.v)
    a = generate_nonfixed_vector(interface.a)

    a0 = 16.0 / dt2
    a1 = 4.0 / dt
    a2 = 9.0 / dt2
    a3 = 3.0 / dt
    a4 = 2.0 * a1
    a5 = 12.0 / dt2
    a6 = -3.0 / dt2
    a7 = -1.0 / dt

    def iteration_first(half_u):
        gravity = np.zeros_like(u)
        gravity[1::2] = -9.8 * np.diag(M)[1::2]

        half_u_full = generate_full_vector(half_u, np.zeros(len(coords)))

        rhs = gravity - generate_nonfixed_vector(interface.stress_forces(half_u_full))
        rhs -= M.dot(a0 * (half_u - u) - a4 * v - a)
        rhs -= C.dot(a1 * (half_u - u) - v)

        K = generate_nonfixed_matrix(interface.linear_k(half_u_full) + interface.nonlinear_k(half_u_full)) + a0 * M + a1 * C
        du = np.linalg.inv(K).dot(rhs)
        half_u += du

        return half_u

    c1 = 0.5 / (0.5 * dt)
    c2 = -1 / ((1 - 0.5) * 0.5 * dt)
    c3 = (2 - 0.5) / ((1 - 0.5) * dt)
    def iteration_second(next_u, half_u, half_v):
        gravity = np.zeros_like(u)
        gravity[1::2] = -9.8 * np.diag(M)[1::2]

        next_u_full = generate_full_vector(next_u, np.zeros(len(coords)))

        rhs = gravity - generate_nonfixed_vector(interface.stress_forces(next_u_full))
        rhs -= M.dot(c3 * c3 * next_u + c3 * c2 * half_u + c3 * c1 * u + c2 * half_v + c1 * v)
        rhs -= C.dot(c1 * u + c2 * half_u + c3 * next_u)

        K = generate_nonfixed_matrix(interface.linear_k(next_u_full) + interface.nonlinear_k(next_u_full)) + c3 * c3 * M + c3 * C
        du = np.linalg.inv(K).dot(rhs)
        next_u += du

        return next_u

    half_u = np.zeros_like(u)
    next_u = np.zeros_like(u)

    previous_u = np.zeros_like(half_u)

    # First iteration
    for i in range(10):
        half_u = iteration_first(half_u)
        if (np.linalg.norm(half_u - previous_u) < 1E-8):
            print (f"Converged after {i + 1} iteration(s)")
            break
        previous_u = half_u
    else:
        print ("Failed to converge.")

    half_v = a1 * (half_u - u) - v

    # Second iteration
    previous_u = np.zeros_like(previous_u)
    for i in range(10):
        next_u = iteration_second(next_u, half_u, half_v)
        if (np.linalg.norm(next_u - previous_u) < 1E-8):
            print (f"Converged after {i + 1} iteration(s)")
            break
        previous_u = next_u
    else:
        print ("Failed to converge.")

    next_v = c1 * u + c2 * half_u + c3 * next_u
    next_a = c1 * v + c2 * half_v + c3 * next_v

    interface.u = generate_full_vector(next_u, np.zeros_like(interface.u))
    interface.v = generate_full_vector(next_v, np.zeros_like(interface.v))
    interface.a = generate_full_vector(next_a, np.zeros_like(interface.a))


def draw_triangles(interface):
    points = coords + interface.u

    for i, triangle in enumerate(triangles):
        x0, y0, x1, y1, x2, y2 = points[triangle]

        stresses = D.dot(interface.nonlinear_strains(triangle, interface.u))
        stress = np.linalg.norm(stresses)

        def color(c):
            return max(min(c, 1.0), 0.0)

        max_stress = 50000
        red = color(stress / max_stress)
        green = color((max_stress - stress) / max_stress)

        plt.fill((x0, x1, x2), (y0, y1, y2), facecolor=(red, green, 0.0), edgecolor="black", linewidth=1, zorder=-10)

    plt.scatter(points[::2], points[1::2], color=["red" if i == 1 else "blue" for i in fixed[::2]])

def compute_energy(interface):
    u = interface.u
    v = interface.v

    K = interface.linear_k(u) + interface.nonlinear_k(u)

    kinetic_energy = 0.5 * v.transpose().dot(M).dot(v)
    strain_potential_energy = 0.5 * u.transpose().dot(K).dot(u)
    gravitational_potential_energy = 9.8 * np.sum(M.dot(u)[1::2])
    return kinetic_energy + strain_potential_energy + gravitational_potential_energy

def draw(i, interface):
    plt.clf()
    plt.gca().set_aspect('equal', adjustable='box')

    print (f"Saving frame {i}")
    plt.xlim(min(coords) - 2, max(coords) + 2)
    plt.ylim(min(coords) - 2, max(coords) + 2)
    draw_triangles(interface)

    if isinstance(i, int):
        plt.savefig(f"/tmp/{i:07}.png")
    else:
        plt.savefig(f"/tmp/{i}.png")

if __name__ == "__main__":
    #interface = TotalLagrangianMethod()
    interface = UpdatedLagrangianMethod()
    draw("00_init", interface)

    max_i = 5000
    if len(sys.argv) == 2:
        max_i = int(sys.argv[1])

    initial_energy = None
    for i in range(max_i):
        print (f"Generating frame {i}")
        update_newmark(interface)

        if initial_energy is None:
            initial_energy = compute_energy(interface)
        print (f"Total Energy: initial: {initial_energy}, at {i}: {compute_energy(interface)} t={i / fps:.2f}s")

        if i % int(1 / 30 * fps) == 0:
            draw(i, interface)

        i += 1
