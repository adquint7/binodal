import numpy as np

global n, chi, vs, mub, vb, mus
vb = np.array([0, 1], dtype='float64')


def fhelm(v):
    return v * np.log(v) / n + (1 - v) * np.log(1 - v) + chi * v * (1 - v)


def muhelm(v, neg=False):
    if neg:
        return -((1 + np.log(v)) / n - (1 + np.log(1 - v)) + chi * (1 - 2 * v))
    else:
        return (1 + np.log(v)) / n - (1 + np.log(1 - v)) + chi * (1 - 2 * v)


def poshelm(v, neg=False):
    pos = v * (1/n - 1) - np.log(1 - v) + chi * v**2
    if neg:
        return -((1 + np.log(v)) / n - (1 + np.log(1 - v)) + chi * (1 - 2 * v))
    else:
        return (1 + np.log(v)) / n - (1 + np.log(1 - v)) + chi * (1 - 2 * v)


def dmuhelm(v, neg=False):
    if neg:
        return -(1 / (n * v) + 1 / (1 - v) - 2 * chi)
    else:
        return 1 / (n * v) + 1 / (1 - v) - 2 * chi


def errf(x, y):
    return np.sqrt(np.sum((x - y) ** 2))


def adamopt(df, f, arg=None, tol=1e-6, a=0.001, b1=0.9, b2=0.99, ep=1e-8):
    if arg is None:
        arg = []
    x = f + 1
    m1 = np.zeros(shape=x.shape)
    m2 = np.zeros(shape=x.shape)
    t = 0
    while errf(f, x) > tol:
        t += 1
        x = f
        g = df(x, *arg)
        m1 = b1*m1 + (1-b1)*g
        m2 = b2*m2 + (1-b2)*g**2
        m1t = m1/(1-b1**t)
        m2t = m2/(1-b2**t)
        f = x - a * m1t/(np.sqrt(m2t) + ep)
    return f


def calvs():
    global vs, mus
    vs = np.empty(shape=2)
    vs[0] = adamopt(dmuhelm, np.array([0.001]), ['True'])
    vs[1] = adamopt(dmuhelm, np.array([1-0.001]))
    mus = muhelm(vs)


def muw(v, bmu, neg=False):
    if neg:
        return -(muhelm(v) - bmu)
    else:
        return muhelm(v) - bmu


def calvb(bmu):
    global vb
    vb[0] = adamopt(muw, vs[0], [bmu])
    vb[1] = adamopt(muw, vs[1], [bmu])


def dmub(bmu, neg=False):
    calvb(bmu)
    if neg:
        return np.diff(fhelm(vb))/np.diff(vb) - bmu
    else:
        return bmu - np.diff(fhelm(vb)) / np.diff(vb)


def binodal():
    global mub, vb, vs
    calvs()
    mub = adamopt(dmub, np.mean(muhelm(vs)))
    calvb(mub)


def floryhugg(x1, x2):
    global n, chi
    n = x1
    chi = x2
    binodal()


def main():
    return 0


if __name__ == '__main__':
    main()
