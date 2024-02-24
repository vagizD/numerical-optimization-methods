import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import pade


a_n = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0,
       1.0/40_320.0, 1.0/362_880.0, 1.0/3_628_800.0, 1.0/39_916_800.0, 1.0/479_001_600.0,
       1.0/6_227_020_800.0, 1.0/87_178_291_200.0, 1.0/1_307_674_368_000.0,
       1.0/20_922_789_888_000.0, 1.0/334_764_638_208_000.0, 1.0/5_690_998_849_536_000.0
       ]

tests = np.linspace(np.log(2)/2.0, np.log(2)/2.0, 10000)
answers = [np.exp(x) for x in tests]


best_loss = np.inf
best_pair = None
best_poly = None

for m in range(1, 20):
    for n in range(1, 20):

        if m + n >= len(a_n):
            continue

        print(f"------ Calc n={n}, m={m} ------")

        loss = -1.0
        try:
            p, q = pade(a_n, n, m)
        except:
            print(f"       Error: n={n}, m={m}        ")
            continue

        def f(x):
            return np.abs(p(x) / q(x) - np.exp(x))

        lower = tests[0]
        upper = tests[-1]
        loss = integrate.quad(f, lower, upper)[0]
        # for x, y_true in zip(tests, answers):
        #     y_pred = p(x) / q(x)
        #     loss = np.max([loss, np.abs(y_true - y_pred)])

        if loss < best_loss:
            best_loss = loss
            best_pair = (n, m)
            best_poly = (p, q)

p, q = best_poly
n, m = best_pair
print(best_loss)

print(f"m = {m}")
print('[')
for c in reversed(p.coef):
    print(c, ',')
# print(q)

print(']')

print(f"n = {n}")
print('[')
for c in reversed(q.coef):
    print(c, ',')
# print(p)

print(']')
