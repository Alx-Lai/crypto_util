from functools import reduce
from sage.all import factor
# adapt from
# https://cocalc.com/share/public_paths/5df71434b58d27a77f417c5740c5fb2b66bfb0f1
# pow(g, ?, p) == h
# p=9695659496402058345675109753360179404032028613270359117817069167974113257013740920328993609480389467494990201488126732816853607851393176121811846442254337
# g=5
# h=1004729126034699826225228696663122643793015439367040380902848957420839629964379383356642006018196989283939275354018486038346278824830164179429161988655007

def xgcd(a,b):
    prevx, x = 1,0
    prevy, y = 0,1
    while b:
        q = a//b
        x, prevx = prevx - q*x, x
        y, prevy = prevy - q*y, y
        a, b = b, a%b
    return a, prevx, prevy

def pohlig_hellman(p, g, h):
    r = p-1
    factor_list = []
    tmp = factor(r)
    for i in tmp:
        factor_list.append((int(i[0]), int(i[1])))

    X = []
    a = []
    for q, e in factor_list:
        a.append(q**e)
        # print(q,e)
        A = pow(g,((p-1)//(q**e)),p)
        B = pow(h,((p-1)//(q**e)) ,p)
        xg = xgcd(A, p)
        # Find A^-1
        A_inv = pow(A, -1, p)
        # print(A ,B)
        # Solve A**x = B for x0,x1 ...
        x = []
        lhs = pow(A, q**(e-1), p)
        rhs = pow(B, q**(e-1), p)
        for x0 in range(p):
            if pow(lhs,x0,p) == rhs % p:
                x.append(x0)
                break
        for i in reversed(range(e - 1)):
            degree = sum([x_i*q**j for x_i, j in zip(x, range(e - i))])
            rhs = pow(B * pow(A_inv,degree, p), q**i, p)
            for xi in range(p):
                if pow(lhs, xi, p) == rhs % p:
                    x.append(xi)
                    break
        # print(x)
        X.append(sum([x_i*(q**j) for x_i, j in zip(x, range(e))]))

    # print(X, a)
    M = reduce(lambda x,y: x * y, a)
    M_i = [M//a_i for a_i in a]
    M_i_inv = []
    for m_i, a_i in zip(M_i, a):
        xg = xgcd(m_i, a_i)
        M_i_inv.append((xg[1] % a_i + a_i) % a_i)
    x_sol = sum([r_i*m_i*m_inv for r_i, m_i, m_inv in zip(*[X, M_i, M_i_inv])]) % M
    return x_sol


if __name__ == "__main__":
    p0=9695659496402058345675109753360179404032028613270359117817069167974113257013740920328993609480389467494990201488126732816853607851393176121811846442254337
    g0=5
    h0=1004729126034699826225228696663122643793015439367040380902848957420839629964379383356642006018196989283939275354018486038346278824830164179429161988655007
    ans = pohlig_hellman(p0, g0, h0)
    assert(pow(g0,ans, p0) == h0)
    print(f'x_sol = {ans}')

