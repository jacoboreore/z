from mpmath import mp

# ------------------ θ via Spouge(+reflection) ------------------
def _spouge_coeffs(a: int):
    c = [mp.sqrt(2*mp.pi)]
    for k in range(1, a):
        sign = -1 if (k % 2 == 0) else 1
        num  = sign * mp.power(a-k, k-mp.mpf('0.5')) * mp.e**(a-k)
        den  = mp.factorial(k-1)
        c.append(num/den)
    return c

def _loggamma_spouge_shift(z, a: int, coeffs=None):
    coeffs = coeffs or _spouge_coeffs(a)
    w = z + (a-1)
    s = coeffs[0]
    for k in range(1, a):
        s += coeffs[k] / (z + (k-1))
    return (z - mp.mpf('0.5'))*mp.log(w) - w + mp.log(s)

def _loggamma_spouge_reflect(z, a: int, coeffs=None):
    coeffs = coeffs or _spouge_coeffs(a)
    lg_1mz = _loggamma_spouge_shift(1 - z, a, coeffs)
    return mp.log(mp.pi) - mp.log(mp.sin(mp.pi*z)) - lg_1mz

def theta_spouge(t, a=45, coeffs=None):
    z = mp.mpf('0.25') + 0.5j * mp.mpf(t)
    lg = _loggamma_spouge_reflect(z, a, coeffs or _spouge_coeffs(a))
    return mp.im(lg) - (mp.mpf(t)/2)*mp.log(mp.pi)

# ------------------ Gabcke/Titchmarsh Ψ and corrections ------------------
def _psi(p):
    cden = mp.cos(2*mp.pi*p)
    if abs(cden) > mp.mpf('1e-30'):
        return mp.cos(2*mp.pi*(p*p - p - mp.mpf('1')/16)) / cden
    # removable singularities near p=1/4, 3/4 (first-order expansion)
    if abs(p-mp.mpf('0.25')) < mp.mpf('1e-8'):
        eps = p - mp.mpf('0.25')
        return mp.mpf('0.5') - eps
    if abs(p-mp.mpf('0.75')) < mp.mpf('1e-8'):
        eps = p - mp.mpf('0.75')
        return mp.mpf('0.5') + eps
    return mp.cos(2*mp.pi*(p*p - p - mp.mpf('1')/16)) / cden

def _psi_d(p, n):
    # bump precision if cos(2πp) small so mp.diff is stable
    extra = 30 if abs(mp.cos(2*mp.pi*p)) < mp.mpf('1e-12') else 10
    with mp.workdps(mp.dps + extra):
        return mp.diff(_psi, p, n)

def C0(p): return _psi(p)
def C1(p): return -(mp.mpf(1)/(96*mp.pi**2)) * _psi_d(p, 3)
def C2(p): return (mp.mpf(1)/(18432*mp.pi**4))*_psi_d(p,6) + (mp.mpf(1)/(64*mp.pi**2))*_psi_d(p,2)

# ------------------ RS parameters ------------------
def rs_a(t): return mp.sqrt(mp.mpf(t)/(2*mp.pi))
def rs_N(t): return int(mp.floor(rs_a(t))) or 1
def rs_p(t): a = rs_a(t); return a - mp.floor(a)

# ------------------ Main RS sum (NO half-weight) ------------------
def Z_rs_main(t, theta):
    N = rs_N(t)
    cphi, sphi = mp.cos(theta), mp.sin(theta)
    S = mp.mpf('0')
    for i in range(N):
        n = i + 1
        if i == 0:
            c, s = cphi, sphi
        else:
            dphi = t * mp.log1p(1/mp.mpf(i))
            sd, cd = mp.sin(dphi), mp.cos(dphi)
            c, s = c*cd + s*sd, s*cd - c*sd
        S += c/mp.sqrt(n)
    return 2*S

# ------------------ Gabcke correction R(t) ------------------
def R_gabcke(t, use_C1=True, use_C2='auto'):
    p = rs_p(t); N = rs_N(t)
    A0 = (2*mp.pi/mp.mpf(t))**(mp.mpf('1')/4)   # (t/2π)^(-1/4)
    sgn = -1 if (N % 2 == 0) else 1            # (-1)^(N-1)
    corr = C0(p)
    if use_C1:
        corr += C1(p) * (2*mp.pi/mp.mpf(t))**(mp.mpf('1')/2)
    # auto-enable C2 when needed (near p≈1/4 or 3/4 or small cos)
    if use_C2 == 'auto':
        use_C2 = (abs(mp.cos(2*mp.pi*p)) < mp.mpf('1e-6'))
    if use_C2:
        corr += C2(p) * (2*mp.pi/mp.mpf(t))
    return sgn * A0 * corr

# ------------------ Siegel Z via RS + Gabcke ------------------
def siegelz_rs(t, *, spouge_a=45, use_C1=True, use_C2='auto', dps=None):
    # optional adaptive precision
    if dps is not None:
        ctx = mp.workdps(dps)
    else:
        # a lightweight heuristic if not specified
        # scale = 70 + max(0, int(mp.log10(mp.mpf(t)) / 1))
        scale = 70 + max(0, int(mp.floor(mp.log10(t))))

        ctx = mp.workdps(min(scale, 140))
    with ctx:
        theta = theta_spouge(t, a=spouge_a)
        Zmain = Z_rs_main(t, theta)
        R1 = R_gabcke(t, use_C1=use_C1, use_C2=use_C2)
        return Zmain + R1, Zmain, R1, theta

# ------------------ reference projector (for testing) ------------------
def Z_ref_mpmath(t, theta):
    s = mp.mpf('0.5') + 1j*mp.mpf(t)
    z = mp.zeta(s)
    c, s_ = mp.cos(theta), mp.sin(theta)
    return mp.re(z)*c - mp.im(z)*s_

# ------------------ quick demo ------------------
if __name__ == "__main__":
    mp.dps = 120
    tests = [
        #mp.mpf("3047.70881764088986878223663882"),
        #mp.mpf("5045.94167550201313249597140915"),
        #mp.mpf("888877770.8075114138098565586614875565161299"),
        #mp.mpf("88888891.29471578310681783197410193108012386"),
        mp.mpf(" 4400961943972.13952128950874642275779953323243530890719745406680758506377745278236083160864758483734 "),
    ]
    for t in tests:
        Ztot, Zmain, R1, theta = siegelz_rs(t, use_C1=True, use_C2='auto', dps=80)
        Zref = Z_ref_mpmath(t, theta) if t < 5e4 else mp.mpf('0')  # avoid slow checks for huge t
        print("\n=== t =", mp.nstr(t, 30), "===")
        print("theta  =", mp.nstr(theta, 30))
        print("Z_main =", mp.nstr(Zmain, 30))
        print("R1     =", mp.nstr(R1, 30))
        print("Z_tot  =", mp.nstr(Ztot, 30))
        if Zref:
            print("Z_ref  =", mp.nstr(Zref, 30))
            print("|Δ|     =", mp.nstr(abs(Ztot - Zref), 12))
