import matplotlib.pyplot as plt
from math import *
from numpy import *
from scipy.optimize import fsolve
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy import delete

try:
    import numpypy as np  # for compatibility with numpy in pypy
except:
    import numpy as np  # if using numpy in cpython


def TDMASolve(a, b, c, d):
    n = len(a)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    xc = []
    for j in range(2, n):
        if (bc[j - 1] == 0):
            ier = 1
            return
        ac[j] = ac[j] / bc[j - 1]
        bc[j] = bc[j] - ac[j] * cc[j - 1]
    if (b[n - 1] == 0):
        ier = 1
        return
    for j in range(2, n):
        dc[j] = dc[j] - ac[j] * dc[j - 1]
    dc[n - 1] = dc[n - 1] / bc[n - 1]
    for j in range(n - 2, -1, -1):
        dc[j] = (dc[j] - cc[j] * dc[j + 1]) / bc[j]
    return dc


C = 2.0
K = 0.01
R = 2.0
T_0 = 0.0
T_1 = 1.0
alpha = 0.01

SQRT_K = pow(K, 0.5)


def plotData(u):
    plt.imshow(u, aspect='auto', interpolation=None)
    # plt.axis([0, time])
    plt.gca().set_autoscale_on(False);
    plt.colorbar()
    plt.show(block='False')


beta = C / K  # here beta is inverted

N = 25  # number of steps by r
M = 15  # number of steps by t

# time = 75000.0
# time = 100.0
time = 200.0
h_t = time / M  # step by t
h_r = R / (N + 1)  # step by r

t = linspace(0, time, M)
r = linspace(0, R, N + 1)

solution = zeros([M, N + 1])
solution[:][0] = T_0

# print solution[0]

coef = csr_matrix((N + 1, N + 1), dtype=float64).todense()  # coefficients for matrix

gamma = h_t / (h_r * beta)

coef[0, 0] = 1 - 6*gamma/h_r
# coef[0, 0] = 6 * h_t / (beta * pow(h_r, 2)) + 1
coef[0, 1] = 6*gamma/h_r
# coef[0, 1] = - 6 * h_t / (beta * pow(h_r, 2))

coef[N, N - 1] = alpha
coef[N, N - 2] = -1 / (2 * h_r)
coef[N, N] = 1 / (2 * h_r)

for i in range(1, N):
    coef[i, i - 1] = -gamma / h_r + gamma / r[i]
    coef[i, i] = 2 * gamma / h_r + 1
    coef[i, i + 1] = -gamma / h_r - gamma / r[i]

# print coef

row = zeros(N+1)

for j in range (1, M):
    row = solution[j-1]
    row[N] = T_1*alpha
    solution[j] = spsolve(coef, row);

print "numeric solution found"

solution = delete(solution, obj=N, axis=1)
plotData(solution)


# --------------------
# now accurate decision is performed

def eigen_function(eta):
    v1 = SQRT_K * (1 / R - alpha)
    v2 = R / SQRT_K
    return v1 / eta - 1 / tan(eta * v2)


def get_eigens(count):
    a = 0
    period = pi * SQRT_K / R
    b = a + period

    cnt = 0

    result = zeros(count)

    while (cnt < count):
        result[cnt] = fsolve(eigen_function, a + 0.01)
        a += period
        b += period
        cnt += 1

    return result


def T(t, eta):
    return exp(-t * pow(eta, 2) / C)


def Ra(r, eta):
    if r < pow(10, -8):
        return eta / SQRT_K

    return sin(r * eta / SQRT_K) / r


def chi(eta):
    return R / 2 - (SQRT_K / (4 * eta)) * sin(2 * R * eta / SQRT_K)


def Ca(eta):
    return alpha * R * sin(R * eta / SQRT_K) * (T_1 - T_0) * K / (pow(eta, 2) * chi(eta))


def u(r, t, eigens):
    length = len(eigens)
    s = 0
    for i in range(length):
        s += (Ca(eigens[i]) * Ra(r, eigens[i]) * T(t, eigens[i]))

    return T_1 - s


eigens = get_eigens(200)

net_u = zeros([M, N])

for i in range(M):
    for j in range(N):
        net_u[i, j] = u(r[j],t[i], eigens)

print "exact solution found"
plotData(net_u)
