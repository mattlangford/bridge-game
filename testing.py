import numpy as np
from matplotlib import pyplot as plt

np.set_printoptions(edgeitems=10, linewidth=100000)

E = 30 * 10E6
v = 0.25
m = 100

# Relates stress to strain
D = np.array([
    [1.0, v, 0.0],
    [v, 1.0, 0.0],
    [0.0, 0.0, 0.5 * (1 - v)]
])
D *= E / (1.0 - v * v)

dt = 1 / 100.
dt2 = dt * dt

coords = np.array([
    0, 0,
    1, 0,
    0, 1,
    1, 1,
    1, 2,
    5, 1,
    5, 2,
])

M = m * np.eye(len(coords))

fixed = np.array([
    1, 1,
    1, 1,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
])

f = np.array([
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
    0, 0,
])
g = 9.8 * m
f[1::2] = -g
print (f)

triangles = [
    [0, 1, 2, 3, 4, 5],
    [0, 1, 4, 5, 6, 7],
    [4, 5, 6, 7, 8, 9],
    [8, 9, 10, 11, 12, 13],
]
u = np.zeros_like(coords, dtype=np.float64)
u_dot = np.zeros_like(u)
u_dot_dot = np.zeros_like(u)
u_dot_dot[1::2] = -9.8

def area(triangle):
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

def k(triangle):
    thickness = 0.5
    b_ = b(triangle)
    k = thickness * area(triangle) * b_.T.dot(D).dot(b_)
    # for i, is_fixed in enumerate(fixed[triangle]):
    #     if is_fixed:
    #         k[i, :] = 0
    #         k[:, i] = 0
    return k

def K():
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

def update(iteration):
    global u, u_dot, u_dot_dot
    K_nonfixed = generate_nonfixed_matrix(K())
    M_nonfixed = generate_nonfixed_matrix(M)
    f_nonfixed = generate_nonfixed_vector(f)
    u_nonfixed = generate_nonfixed_vector(u)
    u_dot_nonfixed = generate_nonfixed_vector(u_dot)
    u_dot_dot_nonfixed = generate_nonfixed_vector(u_dot_dot)

    alpha = 0.5
    beta = 0.25 * (0.5 + alpha) ** 2.0

    a_0 = 1.0 / (beta * dt2)
    a_1 = alpha / (beta * dt)
    a_2 = 1.0 / (beta * dt)
    a_3 = 1.0 / (2.0 * beta) - 1.0
    a_4 = alpha / beta - 1.0
    a_5 = (dt / 2.0) * (alpha / beta - 2.0)
    a_6 = dt * (1.0 - alpha)
    a_7 = alpha * dt

    K_nonfixed = K_nonfixed + a_0 * M_nonfixed

    u_next_nonfixed = f_nonfixed + M_nonfixed.dot(a_0 * u_nonfixed + a_2 * u_dot_nonfixed + a_3 * u_dot_dot_nonfixed)
    u_next_nonfixed = np.linalg.inv(K_nonfixed).dot(u_next_nonfixed)

    u_dot_dot_next_nonfixed = a_0 * (u_next_nonfixed - u_nonfixed) - a_2 * u_dot_nonfixed - a_3 * u_dot_dot_nonfixed
    u_dot_next_nonfixed = u_dot_nonfixed + a_6 * u_dot_dot_nonfixed + a_7 * u_dot_dot_next_nonfixed

    u = generate_full_vector(u_next_nonfixed, np.zeros_like(u))
    u_dot = generate_full_vector(u_dot_nonfixed, np.zeros_like(u_dot))
    u_dot_dot = generate_full_vector(u_dot_dot_nonfixed, np.zeros_like(u_dot_dot))


for i in range(5000):
    print (f"Generating frame {i}")
    update(i)

    points = coords + u
    print (f"u: {u}")
    print (f"vertex: {points[10]}, {points[11]}")

    print (f"Saving frame {i}")
    plt.clf()
    plt.scatter(points[::2], points[1::2])
    plt.savefig(f"/tmp/{i:03}.png")

# print (np.linalg.inv(generate_nonfixed_matrix(K())).dot(generate_nonfixed_vector(f)))
