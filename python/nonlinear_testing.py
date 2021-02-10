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

d = 1
coords = np.array([
    0, 0,
    0, 1,
    1, 0,
    1, 1,
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
    #get_triangle(3, 2, 1),
    # get_triangle(1, 2, 4),
    # get_triangle(3, 4, 9),
    # get_triangle(4, 5, 8),
    # get_triangle(5, 6, 7),
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
u_dot_dot[1::2] = -9.8

def draw_triangles(u):
    points = coords + u

    for i, triangle in enumerate(triangles):
        x0, y0, x1, y1, x2, y2 = points[triangle]
        plt.scatter((x0, x1, x2), (y0, y1, y2))

        _u = u[triangle][::2]
        _v = u[triangle][1::2]

        sign = 1 if i % 2 == 0 else -1
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
            d + _v[1] - _v[2],
          -(-d + _u[1] - _u[2]),
        ])
        print (_u[1] - _u[2], _v[1] - _v[2])

        print (_u, _v)
        print (p10, p20, p12)

        plt.plot(
            [x0, x1, x1, x2, x2, x0],
            [y0, y1, y1, y2, y2, y0])

        plt.plot(
            [(x1 + x0) / 2.0, (x1 + x0) / 2.0 + p10[0]],
            [(y1 + y0) / 2.0, (y1 + y0) / 2.0 + p10[1]])
        plt.plot(
            [(x2 + x0) / 2.0, (x2 + x0) / 2.0 + p20[0]],
            [(y2 + y0) / 2.0, (y2 + y0) / 2.0 + p20[1]])
        plt.plot(
            [(x2 + x1) / 2.0, (x2 + x1) / 2.0 + p12[0]],
            [(y2 + y1) / 2.0, (y2 + y1) / 2.0 + p12[1]])


def draw(i, u):
    plt.clf()
    print (f"Saving frame {i}")
    plt.xlim(min(coords) - 2, max(coords) + 2)
    plt.ylim(min(coords) - 2, max(coords) + 2)
    draw_triangles(u)

    plt.gca().set_aspect('equal', adjustable='box')

    if isinstance(i, int):
        plt.savefig(f"/tmp/{i:07}.png")
    else:
        plt.savefig(f"/tmp/{i}.png")

if __name__ == "__main__":
    u[2] -= 0.5
    u[4] += 3
    draw("00_test", u)
