import numpy as np

v = np.array([846835985, 9834798552])
u = np.array([87502093, 123094980])

while True:
    if np.linalg.norm(u) < np.linalg.norm(v):
        v, u = u, v
    m = v.dot(u) // v.dot(v)
    if m == 0:
        print(u, v)
        break
    u = u - m*v
