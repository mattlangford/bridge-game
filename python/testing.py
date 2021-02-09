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

fps = 250.
dt = 1 / fps
dt2 = dt * dt

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
u = np.zeros_like(coords)
u_dot = np.zeros_like(u)
u_dot_dot = np.zeros_like(u)
u_dot_dot[1::2] = -9.8

def M():
    triangle_volume = 1 * 1 * 1 / 2
    M = np.zeros((len(coords), len(coords)))
    for triangle in triangles:
        mass = triangle_volume * m / 6
        for pt in triangle:
            M[pt, pt] += mass
    return M

def area(triangle):
    displacements = coords[triangle]# + u[triangle]

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


def b(triangle):
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

def b_almansi(u, triangle):
    displacements = coords[triangle] + u[triangle]
    q = u[triangle]
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

    du_dx = inv_det_j * (y(2, 3) * (q[0] - q[4]) - y(1, 3) * (q[2] - q[4]))
    du_dy = inv_det_j * (-x(2, 3) * (q[0] - q[4]) + x(1, 3) * (q[2] - q[4]))
    dv_dx = inv_det_j * (y(2, 3) * (q[1] - q[5]) - y(1, 3) * (q[3] - q[5]))
    dv_dy = inv_det_j * (-x(2, 3) * (q[1] - q[5]) + x(1, 3) * (q[3] - q[5]))

    return np.array([
        du_dx - 0.5 * (du_dx ** 2 + dv_dx ** 2),
        dv_dy - 0.5 * (du_dy ** 2 + dv_dy ** 2),
        0.5 * (du_dy + dv_dx - (du_dx * du_dy + dv_dx * dv_dy))
    ])

def k(triangle):
    thickness = 1
    b_ = b(triangle)
    k = thickness * area(triangle) * b_.T.dot(D).dot(b_)
    return k

def gen_K():
    K = np.zeros([len(coords), len(coords)])
    for i, triangle in enumerate(triangles):
        local_k = k(triangle)
        for row in range(6):
            for col in range(6):
                K[triangle[col], triangle[row]] += local_k[row, col]
    return K

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

def compute_energy(u, u_dot):
    kinetic_energy = 0.5 * u_dot.transpose().dot(M).dot(u_dot)
    strain_potential_energy = 0.5 * u.transpose().dot(K).dot(u)
    gravitational_potential_energy = 9.8 * np.sum(M.dot(u)[1::2])
    return kinetic_energy + strain_potential_energy + gravitational_potential_energy

def update_newmark(iteration):
    global u, u_dot, u_dot_dot
    K_nonfixed = generate_nonfixed_matrix(K)
    M_nonfixed = generate_nonfixed_matrix(M)
    C_nonfixed = generate_nonfixed_matrix(C)
    u_nonfixed = generate_nonfixed_vector(u)
    u_dot_nonfixed = generate_nonfixed_vector(u_dot)
    u_dot_dot_nonfixed = generate_nonfixed_vector(u_dot_dot)

    alpha = 0.5
    beta = 0.25 * (0.5 + alpha) ** 2.0

    n_1 = u_dot_nonfixed + (1 - alpha) * dt * u_dot_dot_nonfixed
    n_2 = u_nonfixed + u_dot_nonfixed * dt + (0.5 - beta) * dt * dt * u_dot_dot_nonfixed

    gravity = np.array([0 if i % 2 == 0 else -9.8 * M_nonfixed[i, i] for i in range(len(u_nonfixed))])
    rhs = gravity + M_nonfixed.dot(n_2 / (beta * dt2)) - C_nonfixed.dot(n_1 - n_2 * alpha / beta)

    u_next = np.linalg.inv(M_nonfixed / (beta * dt2) + C_nonfixed * alpha / (beta * dt) + K_nonfixed).dot(rhs)
    u_dot_dot_next = (u_next - n_2) / (beta * dt2)
    u_dot_next = n_1 + alpha * u_dot_dot_next * dt

    u = generate_full_vector(u_next, np.zeros_like(u))
    u_dot = generate_full_vector(u_dot_next, np.zeros_like(u_dot))
    u_dot_dot = generate_full_vector(u_dot_dot_next, np.zeros_like(u_dot_dot))

def update_hht(iteration):
    global u, u_dot, u_dot_dot
    K_nonfixed = generate_nonfixed_matrix(gen_K())
    M_nonfixed = generate_nonfixed_matrix(M)
    C_nonfixed = generate_nonfixed_matrix(C)
    u_nonfixed = generate_nonfixed_vector(u)
    u_dot_nonfixed = generate_nonfixed_vector(u_dot)
    u_dot_dot_nonfixed = generate_nonfixed_vector(u_dot_dot)

    alpha = 0.0#-1 / 3.0
    gamma = 0.5 * (1 - 2 * alpha)
    beta = 0.25 * (1 - alpha) ** 2

    n_1 = u_nonfixed + u_dot_nonfixed * dt + (0.5 - beta) * dt2 * u_dot_dot_nonfixed
    n_2 = u_dot_nonfixed + (1 - gamma) * dt * u_dot_dot_nonfixed

    gravity = np.array([0 if i % 2 == 0 else -9.8 * M_nonfixed[i, i] for i in range(len(u_nonfixed))])
    rhs = gravity + M_nonfixed.dot(n_1 / (beta * dt2)) - (1 + alpha) * C_nonfixed.dot(n_2 - n_1 * gamma / (beta + dt)) + alpha * C_nonfixed.dot(u_dot_nonfixed) + alpha * K_nonfixed.dot(u_nonfixed)

    u_next = np.linalg.inv(M_nonfixed / (beta * dt2) + (1 + alpha) * C_nonfixed * gamma / (beta * dt) + (1 + alpha) * K_nonfixed).dot(rhs)
    u_dot_dot_next = (u_next - n_1) / (beta * dt2)
    u_dot_next = n_2 + gamma * u_dot_dot_next * dt

    u = generate_full_vector(u_next, np.zeros_like(u))
    u_dot = generate_full_vector(u_dot_next, np.zeros_like(u_dot))
    u_dot_dot = generate_full_vector(u_dot_dot_next, np.zeros_like(u_dot_dot))

K = gen_K()
M = M()

def draw_triangles(u):
    points = coords + u
    forces = M.dot(u_dot_dot)

    for triangle in triangles:
        x1, y1, x2, y2, x3, y3 = points[triangle]

        stresses = D.dot(b_almansi(u, triangle))
        stress = np.linalg.norm(stresses)

        def color(c):
            return max(min(c, 1.0), 0.0)
        max_stress = 10000
        red = color(stress / max_stress)
        green = color((max_stress - stress) / max_stress)

        plt.fill((x1, x2, x3), (y1, y2, y3), facecolor=(red, green, 0.0), edgecolor="black", linewidth=1, zorder=-10)

    plt.quiver(points[::2], points[1::2], forces[::2], forces[1::2], color="black")
    plt.scatter(points[::2], points[1::2], color=["red" if i == 1 else "blue" for i in fixed[::2]])

def draw_k_matrix():
    points = coords + u
    max_element = np.max(K)
    for from_i in range(len(coords))[::2]:
        for to_i in range(len(coords))[::2]:
            if from_i == to_i: continue
            if abs(K[from_i, to_i]) < 1E-3: continue

            from_x = points[from_i]
            from_y = points[from_i + 1]
            to_x = points[to_i]
            to_y = points[to_i + 1]

            # print (f"element: {K[from_i, to_i]}, at from: {from_i} ({from_x}, {from_y}, to: {to_i}, {to_x}, {to_y}")

            element = abs(K[from_i, to_i])
            color = (element / max_element, element / max_element, 0.5)
            plt.plot((from_x, to_x), (from_y, to_y), color=color, linewidth=2)

    plt.scatter(points[::2], points[1::2], colors)


def draw(i, u):
    plt.clf()
    print (f"Saving frame {i}")
    plt.xlim(min(coords) - 2, max(coords) + 2)
    plt.ylim(min(coords) - 2, max(coords) + 2)
    draw_triangles(u)
    #draw_k_matrix()

    if isinstance(i, int):
        plt.savefig(f"/tmp/{i:07}.png")
    else:
        plt.savefig(f"/tmp/{i}.png")

if __name__ == "__main__":
    draw("00_init", u)

    gravity = np.zeros_like(u)
    for i in range(len(gravity)):
        if i % 2 != 0:
            gravity[i] = -9.8 * M[i, i]

    ref_u = generate_full_vector(np.linalg.inv(generate_nonfixed_matrix(K)).dot(generate_nonfixed_vector(gravity)), np.zeros_like(u))
    draw("00_steady", ref_u)

    max_i = 5000
    if len(sys.argv) == 2:
        max_i = int(sys.argv[1])

    initial_energy = None

    for i in range(max_i):
        print (f"Generating frame {i}")
        update_newmark(i)
        if initial_energy is None:
            initial_energy = compute_energy(u, u_dot)
        print (f"Total Energy: initial: {initial_energy}, at {i}: {compute_energy(u, u_dot)}", )
        print (f"u: {u}")

        if i % int(1 / 30 * fps) == 0:
            draw(i, u)

        i += 1
