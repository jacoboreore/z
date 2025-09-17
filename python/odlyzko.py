# odlyzko_block.py
# Prototype Odlyzko–Schönhage-style block evaluator + RS baseline
# Requires: numpy, mpmath

import numpy as np
import mpmath as mp

TWOPI = 2.0 * mp.pi

# ----------------------------
# Helpers: theta(t), RS baseline
# ----------------------------

def theta(t):
    # θ(t) = Im log Γ(1/4 + i t/2) − (t/2) log π
    z = mp.mpc(0.25, t/2.0)
    lg = mp.loggamma(z)
    return lg.imag - (t/2.0) * mp.log(mp.pi)

def Z_riemann_siegel(t, use_gabcke=True):
    """
    Baseline RS evaluator O(sqrt(t)).
    Good for checking correctness at moderate t.
    """
    if t <= 0:
        return mp.mpf('0.0')

    a = mp.sqrt(t / TWOPI)
    N = int(mp.floor(a))
    th = theta(t)

    # Stable rotation
    cphi = mp.cos(th)
    sphi = mp.sin(th)
    S = cphi
    for n in range(1, N):
        dphi = t * mp.log1p(1.0/n)
        cd = mp.cos(dphi)
        sd = mp.sin(dphi)
        new_c = cphi*cd + sphi*sd
        new_s = sphi*cd - cphi*sd
        cphi, sphi = new_c, new_s
        S += cphi / mp.sqrt(n+1)

    Zmain = 2*S

    if not use_gabcke:
        return Zmain

    # Minimal Gabcke-type correction (A0*C0 only)
    p = a - N
    # psi(p) = cos(2π(p^2 - p - 1/16))/cos(2πp) (simple version; no singular guards)
    num = mp.cos(TWOPI*(p*p - p - mp.mpf('0.0625')))
    den = mp.cos(TWOPI*p)
    C0  = num/den if den != 0 else mp.mpf('0.5')  # crude fallback
    A0  = (TWOPI/t)**(mp.mpf('0.25'))
    R1  = A0*C0
    if (N % 2) == 0:
        R1 = -R1
    return Zmain + R1


# -----------------------------------------------------------
# Odlyzko–Schönhage style block summation (quadratic expansion)
# -----------------------------------------------------------

def _amp_poly_coeffs(m, order=2):
    """
    Amplitude expansion for (m+j)^(-1/2) about j=0:
      (m+j)^(-1/2) = m^(-1/2) * [ 1 - j/(2m) + 3 j^2/(8 m^2) - 5 j^3/(16 m^3) + ... ]
    Return coefficients a0..ap where amp ≈ m^(-1/2) * sum_{r=0}^order a_r * (j^r / m^r)
    """
    a = [mp.mpf('1.0')]
    if order >= 1:
        a.append(mp.mpf('-0.5'))
    if order >= 2:
        a.append(mp.mpf('0.375'))  # 3/8
    if order >= 3:
        a.append(mp.mpf('-0.3125'))  # -5/16
    if order >= 4:
        a.append(mp.mpf('0.2734375'))  # 35/128
    return a[:order+1]

def _phase_coeffs(m, order=2):
    """
    Phase expansion for log(m+j) about j=0:
      log(m+j) = log m + j/m - j^2/(2 m^2) + j^3/(3 m^3) - ...
    Return coefficients b1..bp for log(m+j) - log m ≈ sum_{r=1}^order b_r * (j^r / m^r)
    """
    b = [None]  # 1-based indexing for clarity
    if order >= 1:
        b.append(mp.mpf('1.0'))           # j/m
    if order >= 2:
        b.append(mp.mpf('-0.5'))          # - j^2/(2 m^2)
    if order >= 3:
        b.append(mp.mpf('1.0')/3)         # + j^3/(3 m^3)
    if order >= 4:
        b.append(mp.mpf('-0.25'))         # - j^4/(4 m^4)
    return b

def Z_block_OS(T, K=64, h=mp.mpf('0.05'), dps=80, block_alpha=mp.mpf('0.66'),
               phase_order=2, amp_order=2, include_R1=True):
    """
    Evaluate Z(t) for t_k = T + (k - (K-1)/2) * h, k=0..K-1
    using OS-like block expansion:
      sum_{n<=N} n^{-1/2 - i t} ≈ sum over blocks m..m+B-1
           m^{-1/2 - i t} * sum_{j=0}^{B-1} P_amp(j/m) * exp( -i t * P_phase(j/m) )

    Parameters:
      T            : center height
      K, h         : number of points and spacing
      dps          : decimal precision (mpmath)
      block_alpha  : B ≈ m^alpha (typical alpha ~ 2/3)
      phase_order  : Taylor order for log expansion (2 is quadratic/chirp)
      amp_order    : Taylor order for amplitude
      include_R1   : add a light RS remainder (A0*C0)

    Returns:
      ts (mp.mpf array), Zs (mp.mpf array)
    """
    mp.mp.dps = int(dps)

    # RS parameters
    a  = mp.sqrt(T / TWOPI)
    N  = int(mp.floor(a))

    # sample points
    mid = (K-1)/2
    ts  = [T + (k - mid)*h for k in range(K)]

    # accumulate S_k = sum_{n<=N} n^{-1/2} e^{-i t_k log n}
    Sk = [mp.mpc(0) for _ in range(K)]

    # Partition n into blocks with length B ≈ m^alpha, growing with m
    m = 1
    while m <= N:
        B = max(8, int(mp.floor((m**block_alpha))))
        if m + B - 1 > N:
            B = N - m + 1
        # Precompute j powers for the block
        j = np.arange(B, dtype=np.float64)  # double is fine for indices
        # m and powers as mp for precision
        mm = mp.mpf(m)

        # amplitude coefficients
        acoef = _amp_poly_coeffs(mm, amp_order)
        # Build amplitude as polynomial in j/mm:
        # amp[j] ≈ mm^(-1/2) * sum_r a_r * (j/mm)^r
        inv_mm  = 1/mm
        inv_sqr = 1/mp.sqrt(mm)
        # Compute powers of (j/mm) in numpy (double), then cast to mp where needed
        # For stability at large m, j/mm is small.
        j_over_m = mp.mpf('1.0')  # placeholder for loop composition
        # We'll build a vector of mp for amp (costly), but block B is small.
        amp = [mp.mpf('0.0') for _ in range(B)]
        # precompute j^r as int powers via numpy then to mp
        jpow = [np.ones(B, dtype=np.float64)]
        for r in range(1, amp_order+1):
            jpow.append(jpow[-1]*j)  # j^r

        for idx in range(B):
            val = mp.mpf('0.0')
            for r in range(0, amp_order+1):
                coeff = acoef[r]
                if r == 0:
                    val += coeff
                else:
                    val += coeff * mp.mpf(jpow[r][idx]) * (inv_mm**r)
            amp[idx] = inv_sqr * val

        # phase coefficients (without log m)
        bcoef = _phase_coeffs(mm, phase_order)

        # For each t_k: factor out exp(-i t_k log m), then chirp over j
        logm = mp.log(mm)
        for k, t in enumerate(ts):
            # base factor
            base = mp.e**(-1j * t * logm)

            # Build phase: exp(-i t * sum_{r>=1} b_r*(j^r / m^r))
            # We'll do r=1,2 via vectorized numpy for speed, higher r with a loop.
            # a1 = t/m, a2 = -t/(2 m^2), etc.
            # Make numpy arrays of complex64 for the exponent (ok; the end cast is mp)
            # But to keep precision, compute with mp and then fill python list.
            # For r=1,2 we can cheap-evaluate via numpy and cast piecewise.

            # r=1 and r=2 (most impact)
            a1 = bcoef[1] * (t * inv_mm) if phase_order >= 1 else mp.mpf('0.0')
            a2 = bcoef[2] * (t * (inv_mm**2)) if phase_order >= 2 else mp.mpf('0.0')

            # exponent per j ≈ -i*(a1*j + a2*j^2 + ...)
            # Compute with python complex + mp for accuracy
            # (B is small: it’s fine.)
            sblk = mp.mpc(0)
            for idx in range(B):
                ph = mp.mpf('0.0')
                if phase_order >= 1:
                    ph += a1 * mp.mpf(jpow[1][idx])
                if phase_order >= 2:
                    ph += a2 * mp.mpf(jpow[2][idx])
                if phase_order >= 3:
                    a3 = bcoef[3]*(t*(inv_mm**3))
                    ph += a3 * mp.mpf(jpow[3][idx])
                if phase_order >= 4:
                    a4 = bcoef[4]*(t*(inv_mm**4))
                    ph += a4 * mp.mpf(jpow[4][idx])

                sblk += amp[idx] * mp.e**(-1j * ph)

            Sk[k] += base * sblk

        m += B

    # RS style: Z = 2 Re[ e^{i θ(t)} * S(t) ] + R1  (we’ve directly built S with cosine phase)
    # But our S is actually sum n^{-1/2 - i t}, so:
    # Z(t) = e^{iθ} S + e^{-iθ} \bar{S}  = 2 Re(e^{iθ} S)
    Zs = []
    for k, t in enumerate(ts):
        th = theta(t)
        val = 2 * (mp.e**(1j*th) * Sk[k]).real
        if include_R1:
            a = mp.sqrt(t / TWOPI)
            N = int(mp.floor(a))
            p = a - N
            num = mp.cos(TWOPI*(p*p - p - mp.mpf('0.0625')))
            den = mp.cos(TWOPI*p)
            C0  = num/den if den != 0 else mp.mpf('0.5')
            A0  = (TWOPI/t)**(mp.mpf('0.25'))
            R1  = A0*C0
            if (N % 2) == 0:
                R1 = -R1
            val += R1
        Zs.append(val)

    return ts, Zs


# ----------------------------
# Demo / quick check
# ----------------------------
if __name__ == "__main__":
    mp.mp.dps = 60
    T  = mp.mpf('440012345677')      # demo scale; raise with care
    K  = 8
    h  = mp.mpf('0.1')

    ts, Zs = Z_block_OS(T, K=K, h=h, dps=60, block_alpha=mp.mpf('0.66'),
                        phase_order=2, amp_order=2, include_R1=True)

    print("OS-block demo:")
    for t, z in zip(ts, Zs):
        print(f"t={t!s}  Z≈{z}")

    # Single-point cross-check via RS at same precision (only feasible for modest t)
    t0 = mp.mpf('5e6')
    print("\nRS baseline check at t=5e6:")
    print(Z_riemann_siegel(t0))
