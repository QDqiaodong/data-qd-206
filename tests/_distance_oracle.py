"""Independent Bayesian-distance oracle (adaptive Gauss-Legendre quadrature).

This module is deliberately independent of the production package under test:

* integrator: adaptive Gauss-Legendre panels (GL16 cross-checked against GL32),
  whereas production uses adaptive Simpson on a different variable;
* all CDFs, quantiles and pair integrals are evaluated directly from the GL
  panels (no interpolation grid), so their error is the quadrature error;
* Python standard library only, fully deterministic, no randomness.

It exists solely to produce fixed reference numbers for offline tests and is
never imported by the package under test. Run it directly to (re)generate the
saved oracle JSON::

    python3 -I oracle_gl.py
"""
import json
import math

# Deterministic Gauss-Legendre node/weight tables (Newton generated, embedded).
N16_X = [-0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
         -0.755404408355003, -0.6178762444026438, -0.45801677765722737,
         -0.2816035507792589, -0.09501250983763744, 0.09501250983763744,
         0.2816035507792589, 0.45801677765722737, 0.6178762444026438,
         0.755404408355003, 0.8656312023878318, 0.9445750230732326,
         0.9894009349916499]
N16_W = [0.027152459411754058, 0.062253523938647776, 0.0951585116824929,
         0.12462897125553395, 0.14959598881657682, 0.16915651939500256,
         0.1826034150449236, 0.18945061045506847, 0.18945061045506847,
         0.1826034150449236, 0.16915651939500256, 0.14959598881657682,
         0.12462897125553395, 0.0951585116824929, 0.062253523938647776,
         0.027152459411754058]
N32_X = [-0.9972638618494816, -0.9856115115452684, -0.9647622555875064,
         -0.9349060759377397, -0.8963211557660521, -0.84936761373257,
         -0.7944837959679424, -0.7321821187402897, -0.6630442669302152,
         -0.5877157572407623, -0.5068999089322294, -0.42135127613063533,
         -0.33186860228212767, -0.23928736225213706, -0.1444719615827965,
         -0.04830766568773831, 0.04830766568773831, 0.1444719615827965,
         0.23928736225213706, 0.33186860228212767, 0.42135127613063533,
         0.5068999089322294, 0.5877157572407623, 0.6630442669302152,
         0.7321821187402897, 0.7944837959679424, 0.84936761373257,
         0.8963211557660521, 0.9349060759377397, 0.9647622555875064,
         0.9856115115452684, 0.9972638618494816]
N32_W = [0.007018610009470136, 0.01627439473090571, 0.02539206530926214,
         0.03427386291302141, 0.042835898022226704, 0.050998059262376154,
         0.05868409347853558, 0.06582222277636195, 0.07234579410884862,
         0.07819389578707044, 0.08331192422694672, 0.08765209300440374,
         0.0911738786957639, 0.09384439908080454, 0.09563872007927485,
         0.09654008851472785, 0.09654008851472785, 0.09563872007927485,
         0.09384439908080454, 0.0911738786957639, 0.08765209300440374,
         0.08331192422694672, 0.07819389578707044, 0.07234579410884862,
         0.06582222277636195, 0.05868409347853558, 0.050998059262376154,
         0.042835898022226704, 0.03427386291302141, 0.02539206530926214,
         0.01627439473090571, 0.007018610009470136]

LOG2PI = math.log(2.0 * math.pi)
PANEL_TOL = 1e-12
MIN_PANEL = 1e-7


def gauss_logpdf(x, mu, sig):
    return -0.5 * ((x - mu) / sig) ** 2 - math.log(sig) - 0.5 * LOG2PI


class Source:
    __slots__ = ("y", "sig", "L")

    def __init__(self, y, sig, L):
        self.y, self.sig, self.L = y, sig, L

    def log_rel(self, s, z, peak):
        """log of unnormalized posterior density in s=ln r, minus ``peak``."""
        r = math.exp(s)
        val = (3.0 * s - r / self.L
               - 0.5 * ((self.y - z - 1000.0 / r) / self.sig) ** 2)
        return val - peak


def golden_max(fn, a, b):
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = fn(c), fn(d)
    for _ in range(120):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fn(d)
        else:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fn(c)
    m = 0.5 * (a + b)
    return m, fn(m)


def find_mode(fn, s_guess):
    """Deterministic uphill walk to bracket the single mode, then golden search."""
    step = 1.0
    a = s_guess
    fup, fdown = fn(a + step), fn(a - step)
    direction = 1.0 if fup >= fdown else -1.0
    b = a + direction * step
    fb = fup if direction > 0 else fdown
    for _ in range(500):
        c = b + direction * step
        fc = fn(c)
        if fc <= fb:
            return golden_max(fn, min(a, c), max(a, c))
        a, b, fb = b, c, fc
        step *= 1.05
    return golden_max(fn, b - step, b + step)


def _gl_scaled(src, z, peak, a, b, xs, ws):
    """GL integral of the relative density exp(log_rel) over panel [a, b]."""
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    total = 0.0
    for x, w in zip(xs, ws):
        s = mid + half * x
        total += w * math.exp(src.log_rel(s, z, peak))
    return half * total


def build_panels(src, z):
    """Adaptive GL panels covering the support; returns (panels, peak, logZ).

    Each panel is (a, b, integral). Panels split while GL32 disagrees with
    2xGL16 beyond ``PANEL_TOL``; underflowed deep-tail panels are accepted.
    """
    y_eff = src.y - z
    s_guess = math.log(1000.0 / y_eff) if y_eff > 0.0 else math.log(src.L)
    mode, peak = find_mode(lambda s: src.log_density_abs(s, z), s_guess)
    h = 1e-4 * (1.0 + abs(mode))
    g2 = (src.log_density_abs(mode + h, z) - 2.0 * peak
          + src.log_density_abs(mode - h, z)) / (h * h)
    width = 1.0 / math.sqrt(-g2) if g2 < 0 else 1.0
    lo, hi = mode - 8.0 * width, mode + 8.0 * width
    while src.log_density_abs(lo, z) > peak - 46.0 and mode - lo < 160.0:
        lo -= 2.0 * width
    while src.log_density_abs(hi, z) > peak - 46.0 and hi - mode < 160.0:
        hi += 2.0 * width

    stack = [(lo, hi)]
    panels = []
    while stack:
        a, b = stack.pop()
        i32 = _gl_scaled(src, z, peak, a, b, N32_X, N32_W)
        mid = 0.5 * (a + b)
        i16 = (_gl_scaled(src, z, peak, a, mid, N16_X, N16_W)
               + _gl_scaled(src, z, peak, mid, b, N16_X, N16_W))
        if i32 == 0.0 and i16 == 0.0:
            panels.append((a, b, 0.0))
        elif abs(i32 - i16) > PANEL_TOL * max(i32, 1e-300) and (b - a) > MIN_PANEL:
            stack.append((a, mid))
            stack.append((mid, b))
        else:
            panels.append((a, b, i32))
    panels.sort()
    total = sum(p[2] for p in panels)
    log_z = peak + math.log(total)
    return panels, total, log_z


def cdf_from_panels(src, z, peak, panels, total, s):
    """Normalized CDF at s via direct GL integration up to s (no grid)."""
    if s <= panels[0][0]:
        return 0.0
    if s >= panels[-1][1]:
        return 1.0
    acc = 0.0
    for a, b, integ in panels:
        if b <= s:
            acc += integ
        elif a < s < b:
            half = 0.5 * (s - a)
            mid = 0.5 * (a + s)
            part = 0.0
            for x, w in zip(N32_X, N32_W):
                ss = mid + half * x
                part += w * math.exp(src.log_rel(ss, z, peak))
            acc += half * part
            break
        else:
            break
    return min(1.0, max(0.0, acc / total))


def quantile(conditional, q):
    """Invert a (possibly mixture) CDF by deterministic bisection in s=ln r."""
    lo = min(tab["lo"] for tab in conditional)
    hi = max(tab["hi"] for tab in conditional)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        c = sum(item["w"] * cdf_from_panels(
            item["src"], item["z"], item["peak"], item["panels"],
            item["total"], mid) for item in conditional)
        if c < q:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def pair_probability(cond_near, cond_far):
    """GL evaluation of the ordering integral, panel by panel.

    P(r_near < r_far) = sum_z w_z int p_n(s|z) F_f(s|z) ds.
    """
    result = 0.0
    for n_item, f_item in zip(cond_near, cond_far):
        sub = 0.0
        for a, b, _integ in n_item["panels"]:
            half = 0.5 * (b - a)
            mid = 0.5 * (a + b)
            panel_sum = 0.0
            for x, w in zip(N32_X, N32_W):
                s = mid + half * x
                density = math.exp(
                    n_item["src"].log_rel(s, n_item["z"], n_item["peak"])) \
                    / n_item["total"]
                cdf = cdf_from_panels(
                    f_item["src"], f_item["z"], f_item["peak"],
                    f_item["panels"], f_item["total"], s)
                panel_sum += w * density * cdf
            sub += half * panel_sum
        result += n_item["w"] * sub
    return result


# attach absolute log density used by the mode finder and bracket logic
def _log_density_abs(self, s, z):
    r = math.exp(s)
    return (3.0 * s - r / self.L
            - 0.5 * ((self.y - z - 1000.0 / r) / self.sig) ** 2)


Source.log_density_abs = _log_density_abs


def run(data, z_n=401):
    L = float(data["scale_length_pc"])
    mu = float(data["zero_point"]["mean_mas"])
    tau = float(data["zero_point"]["sigma_mas"])
    srcs = [(s["id"], Source(float(s["parallax_mas"]),
                             float(s["sigma_mas"]), L))
            for s in data["sources"]]
    fixed = tau == 0.0

    def log_evidence(z):
        v = 0.0 if fixed else gauss_logpdf(z, mu, tau)
        for _id, sc in srcs:
            _p, _t, lz = build_panels(sc, z)
            v += lz
        return v

    if fixed:
        zgrid = [mu]
        weights = [1.0]
    else:
        zlo0, zhi0 = mu - 16.0 * tau, mu + 16.0 * tau
        nscan = 160
        zs = [zlo0 + (zhi0 - zlo0) * i / nscan for i in range(nscan + 1)]
        hs = [log_evidence(z) for z in zs]
        j = max(range(nscan + 1), key=lambda k: hs[k])
        zm, _hm = golden_max(
            log_evidence, zs[max(0, j - 1)], zs[min(nscan, j + 1)])
        h = 1e-5 * max(1.0, abs(zm))
        h2 = (log_evidence(zm + h) - 2.0 * log_evidence(zm)
              + log_evidence(zm - h)) / (h * h)
        wz = 1.0 / math.sqrt(-h2) if h2 < 0 else tau
        zlo = max(zlo0, zm - 12.0 * wz)
        zhi = min(zhi0, zm + 12.0 * wz)
        zgrid = [zlo + (zhi - zlo) * i / (z_n - 1) for i in range(z_n)]
        ev = [log_evidence(z) for z in zgrid]
        emax = max(ev)
        hv = [math.exp(e - emax) for e in ev]
        dz = (zhi - zlo) / (z_n - 1)
        wsum = (hv[0] * 0.5 + hv[-1] * 0.5 + sum(hv[1:-1])) * dz
        weights = [v * dz / wsum for v in hv]
        weights[0] *= 0.5
        weights[-1] *= 0.5

    # build every conditional once, sequentially (small memory footprint)
    conditionals = {}
    for sid, sc in srcs:
        items = []
        for w, z in zip(weights, zgrid):
            panels, total, _lz = build_panels(sc, z)
            y_eff = sc.y - z
            s_guess = math.log(1000.0 / y_eff) if y_eff > 0.0 else math.log(sc.L)
            _mode, peak = find_mode(lambda s: sc.log_density_abs(s, z), s_guess)
            items.append({
                "w": w, "z": z, "src": sc, "panels": panels,
                "total": total, "peak": peak,
                "lo": panels[0][0], "hi": panels[-1][1],
            })
        conditionals[sid] = items

    ez = sum(w * z for w, z in zip(weights, zgrid))
    ez2 = sum(w * z * z for w, z in zip(weights, zgrid))
    result = {
        "zero_point": {
            "mean_mas": ez,
            "sigma_mas": math.sqrt(max(0.0, ez2 - ez * ez)),
        },
        "sources": {},
        "pairs": {},
    }
    for sid, _sc in srcs:
        cond = conditionals[sid]
        result["sources"][sid] = {
            "median_pc": quantile(cond, 0.5),
            "ci_lower_pc": quantile(cond, 0.15865),
            "ci_upper_pc": quantile(cond, 0.84135),
        }
    for p in data["pairs"]:
        result["pairs"][p["id"]] = pair_probability(
            conditionals[p["nearer"]], conditionals[p["farther"]])
    return result


FIXTURE = {
    "scale_length_pc": "1200",
    "zero_point": {"mean_mas": "0.08", "sigma_mas": "0.15"},
    "sources": [
        {"id": "A", "parallax_mas": "-0.20", "sigma_mas": "0.50"},
        {"id": "B", "parallax_mas": "0.75", "sigma_mas": "0.35"},
    ],
    "pairs": [{"id": "A-before-B", "nearer": "A", "farther": "B"}],
}

if __name__ == "__main__":
    print(json.dumps(run(FIXTURE), indent=2, sort_keys=True))
