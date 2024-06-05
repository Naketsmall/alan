import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def calculate_dx(x1, T1, x2, T2):
    return (x1 * (T2 - 810) + x2 * (810 - T1)) / (T2 - T1)


def integrate(X, F):
    return sum([(F[i] + F[i - 1]) * (X[i] - X[i - 1]) / 2 for i in range(1, len(X))])

def integrate_var(X, F):
    I = [0]
    for i in range(1, len(X)):
        I.append(I[i-1] + (F[i-1] + F[i])/2 * (X[i]-X[i-1]))
    return I


def get_h(A, B, C, T):
    return A + B * T / 1000 + C * (T / 1000) ** 2


def get_n(p1, r1, p2, r2):
    with np.errstate(invalid='raise'):
        try:
            if abs(r2 - r1) <= 10 ** -10:
                return 0
            if (abs(p2) <= 10 ** -10):
                return p2
            print(np.log(p2 / p1) / np.log(r2 / r1))
            return np.log(p2 / p1) \
                / np.log(r2 / r1)
        except FloatingPointError:
            print('DBZ: ', p1, r1, p2, r2)


def get_xi1(n, T):
    return -n / (1 - n) * T

def get_pi1(h, R, m, xi):
    return h - R / m * xi


def get_AB(p1, t1, p2, t2):
    if t1 - t2 == 0:
        if p1 == p2:
            #print('A - problem1')
            return 0, p2
        if p1 != p2:
            #print('A - problem2')
            return 999999999999999, p2
        return 0, p2
    else:
        A = (p1 - p2) / (t1 - t2)
        return A, p1 - A * t1


def get_pi2(h, R, m, A, B, T):
    with np.errstate(invalid='raise'):
        try:
            if abs(A) > 10 ** -10 and abs(B) > 10 ** -10:
                if (A * T + B) / B < 0:
                    return 0
                return h - R / m * (T - B / A * np.log((A * T + B) / B))
            elif abs(A) <= 10 ** -10 and abs(B) > 10 ** -10:
                return h
            elif abs(A) > 10 ** -10 and abs(B) <= 10 ** -10:
                return h - R / m * T
            else:
                print("EXCEPTION: A, B are too small", A, B)
                return 0
        except FloatingPointError:
            print('wtf:', (A * T + B) / B)

class Converter:

    def get_xi2(self, i, xip):
        T1 = self.DS['T'][self.x0_index]
        A, B = get_AB(self.DS['p'][i-1], self.DS['T'][i-1], self.DS['p'][i], self.DS['T'][i])
        print('A, B', A, B)
        return T1 - B/A*np.log((A*T1 + B)/B) + self.vpt(i) + xip, self.vpt(i) + xip

    def vpt(self, i):
        if i <= self.x0_index: return 0
        v1, v2 = 1/self.DS['rho'][i-1], 1/self.DS['rho'][i]
        p1, p2 = self.DS['p'][i-1], self.DS['p'][i]
        T1, T2 = self.DS['T'][i-1], self.DS['T'][i]

        S = (((v1+v2) * (p2-p1)) / ((p1+p2) * (v2-v1)))
        print(v1, v2, p1, p2, T1, T2)
        return S/(1-S)*(T2-T1)

    def __init__(self, DS, M, Y_names, A, B, C, eqs):
        self.x0_index = None
        self.DS = DS
        self.M = M
        self.Y_names = Y_names
        self.A = A
        self.B = B
        self.C = C
        self.eqs = eqs
        self.r = len(eqs)
        self.k = len(A)

        self.shift_x()

    def shift_x(self):
        x0_index = self.DS[self.DS['T'] < int(810)].index[-1]
        dx = calculate_dx(self.DS['x'][x0_index], self.DS['T'][x0_index],
                          self.DS['x'][x0_index + 1], self.DS['T'][x0_index + 1])
        self.DS['x'] -= dx
        self.x0_index = x0_index

    def calculate_q1(self):
        q1 = []
        H0 = [get_h(self.A[k], self.B[k], self.C[k], 293.15) for k in range(self.k)]
        Y0 = list(self.DS[self.Y_names].iloc[self.x0_index])

        for i in range(len(self.DS)):
            q1.append(-sum([H0[k] * (self.DS[self.Y_names[k]][i] - Y0[k]) for k in range(self.k)]))

        return q1

    def calculate_Q1(self):
        Q1 = [[] for r in range(self.r)]

        for i in range(len(self.DS)):
            Hk = [get_h(self.A[k], self.B[k], self.C[k], self.DS['T'][i]) for k in range(self.k)]
            for r in range(self.r):
                Q1[r].append(sum([Hk[k] * self.M[k] * self.eqs[r][k] for k in range(self.k)]))

        return np.array(Q1).T

    def calculate_Q2(self):
        Q2 = [[0] for r in range(self.r)]  #TODO: fill start value

        for i in range(1, len(self.DS)):
            Hk = [get_h(self.A[k], self.B[k], self.C[k], self.DS['T'][i]) for k in range(self.k)]
            n = get_n(self.DS['p'][i - 1], self.DS['rho'][i - 1], self.DS['p'][i], self.DS['rho'][i])
            xi = get_xi1(n, self.DS['T'][i])
            pik = [get_pi1(Hk[k], 8.314, self.M[k], xi) for k in range(self.k)]

            for r in range(self.r):
                Q2[r].append(sum([pik[k] * self.M[k] * self.eqs[r][k] for k in range(self.k)]))

        return np.array(Q2).T

    def calculate_Q3(self):
        Q3 = [[0] for r in range(self.r)]

        for i in range(1, len(self.DS)):
            Hk = [get_h(self.A[k], self.B[k], self.C[k], self.DS['T'][i]) for k in range(self.k)]
            A, B = get_AB(self.DS['p'][i - 1], self.DS['T'][i - 1], self.DS['p'][i], self.DS['T'][i])
            pik = [get_pi2(Hk[k], 8.314, self.M[k], A, B, self.DS['T'][i]) for k in range(self.k)]

            for r in range(self.r):
                Q3[r].append(sum([pik[k] * self.M[k] * self.eqs[r][k] for k in range(self.k)]))

        return np.array(Q3).T

    def calculate_Q4(self):
        Q4 = [[0] for r in range(self.r)]

        xip = 0
        for i in range(1, len(self.DS)):
            Hk = [get_h(self.A[k], self.B[k], self.C[k], self.DS['T'][i]) for k in range(self.k)]
            xi, xip = self.get_xi2(i, xip)
            print(xi, xip)
            pik = [get_pi1(Hk[k], 8.314, self.M[k], xi) for k in range(self.k)]

            for r in range(self.r):
                Q4[r].append(sum([pik[k] * self.M[k] * self.eqs[r][k] for k in range(self.k)]))

        return np.array(Q4).T

    def calculate_dqQ(self, Q):
        dq = [0]

        for i in range(1, len(self.DS)):
            dq.append(sum([Q[i][r] * self.DS['Rate[' + str(r) + ']'][i] for r in range(self.r)])
                      / (self.DS['rho'][i] * self.DS['u'][i]))
        return dq
