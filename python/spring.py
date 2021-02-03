import numpy as np
import matplotlib.pyplot as plt

fps = 100.
dt = 1 / fps
dt2 = dt * dt

g = -9.8

m = 10
k = 1000

# With no mass, the spring would stretch to here
spring_equilibrium = 0

# Where there is no net force, it rests here
x_equilibrium = spring_equilibrium + m * g / k

def gravitational_potential_energy(x):
    x_ref = 0
    return -m * g * (x - x_ref)
def spring_potential_energy(x):
    return 0.5 * k * (x - spring_equilibrium) ** 2
def kinetic_energy(v):
    return 0.5 * m * (v ** 2)

w = np.sqrt(k / m)

u = m * g / k
u_v = 0
u_a = 0
equi = spring_potential_energy(u) + gravitational_potential_energy(u) + kinetic_energy(u_v)
print ("equi:", equi)
u = 0

def update():
    global u, u_v, u_a
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

    gravity = -9.8 * m

    u_next = np.copy(gravity)
    u_next += m * (a_0 * u + a_2 * u_v + a_3 * u_a)

    u_next = u_next / (k + a_0 * m)

    u_a_next = a_0 * (u_next - u) - a_2 * u_v - a_3 * u_a
    u_v_next = u_v + a_6 * u_a + a_7 * u_a_next

    u = u_next
    u_v = u_v_next
    u_a = u_a_next

us = []
energy = []
for i in range(500):
    update()
    us.append(u)
    energy.append(spring_potential_energy(u) + gravitational_potential_energy(u) + kinetic_energy(u_v))

t = np.arange(0, len(us))

fig, ax1 = plt.subplots()
ax1.set_xlabel('frame')
ax1.set_ylabel('u', color="red")
ax1.plot(t, us, color="red")
ax1.hlines(y=m * g / k, xmin=0, xmax=len(us), color="red")
ax1.tick_params(axis='y', labelcolor="red")

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('energy', color="green")
ax2.plot(t, energy, color="green")
ax2.hlines(y=equi, xmin=0, xmax=len(us), color="blue")
ax2.tick_params(axis='y', labelcolor="green")

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


# ts = np.arange(0, np.pi, 0.01)
# x = A * np.cos(w * ts) + x_equilibrium
# v = -A * w * np.sin(w * ts)
# 
# U_s = spring_potential_energy(x)
# U_g = gravitational_potential_energy(x)
# K = kinetic_energy(v)
# 
# plt.axhline(y=x_equilibrium)
# plt.plot(ts, x, c="yellow")
# 
# plt.plot(ts, U_s, c="black")
# plt.plot(ts, U_g, c="blue")
# plt.plot(ts, K, c="red")
# 
# total = K + U_g + U_s
# plt.plot(ts, total, c="green")
# plt.show()
