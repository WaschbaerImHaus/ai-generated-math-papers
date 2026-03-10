"""
@file pde.py
@brief PDE-Modul: Partielle Differentialgleichungen (Partial Differential Equations).
@description
    Enthält numerische und analytische Lösungsverfahren für PDEs:
    - Wärmeleitungsgleichung (explizit, implizit, Crank-Nicolson)
    - Wellengleichung (Leapfrog, d'Alembert)
    - Laplace/Poisson-Gleichung (SOR, sparse Direkt)
    - Schrödinger-Gleichung (stationär & zeitabhängig)
    - Finite-Elemente-Methode 1D (Galerkin, Hutfunktionen)
    - Methode der Charakteristiken
    - Viskose Burgers-Gleichung

    Numerische Schemata:
    - FTCS (Forward-Time Central-Space) für Wärme (explizit)
    - BTCS (Backward-Time Central-Space) für Wärme (implizit)
    - Crank-Nicolson θ=1/2 (O(Δt²+Δx²), unbedingt stabil)
    - Leapfrog für Wellengleichung (CFL-Bedingung)
    - SOR-Iteration für Laplace
    - Sparse-LGS für Poisson

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
from typing import Callable, Optional
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import eigh_tridiagonal


# =============================================================================
# Klassifikation von PDEs 2. Ordnung
# =============================================================================

def classify_pde(A: float, B: float, C: float) -> str:
    """
    @brief Klassifikation einer PDE 2. Ordnung in 2 Variablen.
    @description
        Klassifiziert die PDE Au_xx + Bu_xy + Cu_yy + ... = 0
        anhand der Diskriminante D = B² - 4AC:
        - D < 0: elliptisch (z.B. Laplace-Gleichung)
        - D = 0: parabolisch (z.B. Wärmeleitungsgleichung)
        - D > 0: hyperbolisch (z.B. Wellengleichung)

    @param A Koeffizient vor u_xx
    @param B Koeffizient vor u_xy
    @param C Koeffizient vor u_yy
    @return Klassifikationsstring ('elliptic', 'parabolic', 'hyperbolic')
    @lastModified 2026-03-10
    """
    # Diskriminante berechnen
    discriminant = B**2 - 4 * A * C
    # Toleranz für numerischen Vergleich mit 0
    tol = 1e-12
    if discriminant < -tol:
        return 'elliptic'
    elif discriminant > tol:
        return 'hyperbolic'
    else:
        return 'parabolic'


# =============================================================================
# Wärmeleitungsgleichung u_t = α² u_xx
# =============================================================================

def heat_equation_explicit(
    u0: list,
    L: float,
    T: float,
    alpha: float = 1.0,
    nx: int = 50,
    nt: int = 500
) -> np.ndarray:
    """
    @brief Explizites FTCS-Verfahren für die Wärmeleitungsgleichung.
    @description
        Löst u_t = α² u_xx mit Dirichlet-Randwerten u(0,t)=u(L,t)=0
        und Anfangsbedingung u(x,0)=u0(x_i).

        Schema: u_i^{n+1} = u_i^n + r(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
        mit r = α²Δt/Δx².
        Stabilitätsbedingung: r ≤ 1/2.

    @param u0   Anfangswerte (nx+1 Punkte oder interpolierbar)
    @param L    Länge des Intervalls [0, L]
    @param T    Endzeit
    @param alpha Wärmediffusivität α
    @param nx   Anzahl der Raumgitterpunkte (innere)
    @param nt   Anzahl der Zeitschritte
    @return     2D-Array der Form (nt+1, nx+1): Lösung u(x,t)
    @lastModified 2026-03-10
    """
    # Gitterschrittweiten
    dx = L / nx
    dt = T / nt
    # Courant-Zahl (Stabilitätsparameter)
    r = alpha**2 * dt / dx**2

    if r > 0.5:
        # Warnung ausgeben, aber trotzdem berechnen
        import warnings
        warnings.warn(
            f"Stabilitätsbedingung verletzt: r={r:.4f} > 0.5. "
            "Das Ergebnis kann instabil sein."
        )

    # Gitterpunkte x_0, x_1, ..., x_nx
    x = np.linspace(0, L, nx + 1)
    # Anfangsbedingung auf das Gitter interpolieren
    u0_arr = np.array(u0, dtype=float)
    if len(u0_arr) != nx + 1:
        # Lineare Interpolation auf das Gitter
        x0_arr = np.linspace(0, L, len(u0_arr))
        u0_arr = np.interp(x, x0_arr, u0_arr)

    # Lösungs-Array: Zeile i = Zeitschritt i
    u = np.zeros((nt + 1, nx + 1))
    u[0, :] = u0_arr
    # Randwerte für alle Zeiten = 0
    u[:, 0] = 0.0
    u[:, nx] = 0.0

    # Zeitintegration via explizitem Euler
    for n in range(nt):
        # Innere Punkte aktualisieren (Vektorisiert)
        u[n + 1, 1:-1] = (
            u[n, 1:-1]
            + r * (u[n, 2:] - 2 * u[n, 1:-1] + u[n, :-2])
        )

    return u


def heat_equation_implicit(
    u0: list,
    L: float,
    T: float,
    alpha: float = 1.0,
    nx: int = 50,
    nt: int = 100
) -> np.ndarray:
    """
    @brief Implizites BTCS-Verfahren für die Wärmeleitungsgleichung.
    @description
        Unbedingt stabiles Schema. In jedem Zeitschritt wird
        ein tridiagonales LGS gelöst:

        (I + r·K) u^{n+1} = u^n

        wobei K die zweite Differenzmatrix ist und r = α²Δt/Δx².

    @param u0   Anfangswerte
    @param L    Länge des Intervalls
    @param T    Endzeit
    @param alpha Wärmediffusivität
    @param nx   Anzahl der Raumgitterpunkte
    @param nt   Anzahl der Zeitschritte
    @return     2D-Array (nt+1, nx+1): Lösung u(x,t)
    @lastModified 2026-03-10
    """
    dx = L / nx
    dt = T / nt
    r = alpha**2 * dt / dx**2
    # Anzahl innerer Punkte
    m = nx - 1

    # Anfangsbedingung interpolieren
    x = np.linspace(0, L, nx + 1)
    u0_arr = np.array(u0, dtype=float)
    if len(u0_arr) != nx + 1:
        x0_arr = np.linspace(0, L, len(u0_arr))
        u0_arr = np.interp(x, x0_arr, u0_arr)

    # Lösungs-Array
    u = np.zeros((nt + 1, nx + 1))
    u[0, :] = u0_arr

    # Tridiagonale Systemmatrix (I + r·K) für innere Punkte
    # Diagonale: 1+2r, Nebendiagonalen: -r
    diag_main = np.full(m, 1.0 + 2.0 * r)
    diag_off = np.full(m - 1, -r)

    # Sparse-Tridiagonalmatrix aufbauen
    A_mat = sparse.diags(
        [diag_off, diag_main, diag_off],
        offsets=[-1, 0, 1],
        format='csc'
    )

    for n in range(nt):
        # Rechte Seite = innere Punkte des aktuellen Zeitschritts
        rhs = u[n, 1:-1].copy()
        # LGS lösen
        u_inner = spsolve(A_mat, rhs)
        u[n + 1, 1:-1] = u_inner
        # Randwerte bleiben 0

    return u


def heat_equation_crank_nicolson(
    u0: list,
    L: float,
    T: float,
    alpha: float = 1.0,
    nx: int = 50,
    nt: int = 100
) -> np.ndarray:
    """
    @brief Crank-Nicolson-Verfahren für die Wärmeleitungsgleichung.
    @description
        Zeitlich und räumlich zweiter Ordnung: O(Δt² + Δx²).
        Unbedingt stabil (θ = 1/2, Mittelung von explizit und implizit).

        Schema: (I + r/2·K) u^{n+1} = (I - r/2·K) u^n

    @param u0   Anfangswerte
    @param L    Länge des Intervalls
    @param T    Endzeit
    @param alpha Wärmediffusivität
    @param nx   Anzahl der Raumgitterpunkte
    @param nt   Anzahl der Zeitschritte
    @return     2D-Array (nt+1, nx+1): Lösung u(x,t)
    @lastModified 2026-03-10
    """
    dx = L / nx
    dt = T / nt
    r = alpha**2 * dt / dx**2
    m = nx - 1  # Anzahl innerer Punkte

    # Anfangsbedingung
    x = np.linspace(0, L, nx + 1)
    u0_arr = np.array(u0, dtype=float)
    if len(u0_arr) != nx + 1:
        x0_arr = np.linspace(0, L, len(u0_arr))
        u0_arr = np.interp(x, x0_arr, u0_arr)

    u = np.zeros((nt + 1, nx + 1))
    u[0, :] = u0_arr

    # Linke Seite: (I + r/2·K)
    a_diag = np.full(m, 1.0 + r)
    a_off = np.full(m - 1, -r / 2.0)
    A_mat = sparse.diags([a_off, a_diag, a_off], [-1, 0, 1], format='csc')

    # Rechte Seite: (I - r/2·K) als dichte Matrix (oder Anwendung als Funktion)
    # Für Effizienz: rechte Seite direkt berechnen
    for n in range(nt):
        u_in = u[n, 1:-1]
        # Rechte Seite: (I - r/2·K) u^n
        rhs = u_in.copy()
        rhs[1:] += (r / 2.0) * u_in[:-1]
        rhs[:-1] += (r / 2.0) * u_in[1:]
        rhs -= r * u_in
        # Korrekte Formel: rhs_i = (r/2)u_{i-1} + (1-r)u_i + (r/2)u_{i+1}
        rhs2 = np.zeros(m)
        rhs2[:] = (1.0 - r) * u_in
        rhs2[:-1] += (r / 2.0) * u_in[1:]
        rhs2[1:] += (r / 2.0) * u_in[:-1]

        u[n + 1, 1:-1] = spsolve(A_mat, rhs2)

    return u


def heat_equation_analytical(
    x: np.ndarray,
    t: float,
    L: float = np.pi,
    n_terms: int = 20
) -> np.ndarray:
    """
    @brief Analytische Lösung der Wärmeleitungsgleichung via Fourier-Reihe.
    @description
        Für u(x,0) = sin(x) auf [0, π] mit u(0,t)=u(π,t)=0:
        u(x,t) = sin(x) · e^{-t}

        Allgemeine Fourier-Reihe: u(x,t) = Σ bₙ sin(nπx/L) e^{-(nπ/L)²t}
        mit bₙ = (2/L) ∫₀ᴸ u(x,0) sin(nπx/L) dx

        Für Anfangsbedingung sin(x) auf [0,π]: nur b₁=1 != 0.

    @param x       Ortsvektor
    @param t       Zeitpunkt
    @param L       Länge des Intervalls (Standard: π)
    @param n_terms Anzahl der Fourier-Terme
    @return        Lösungsvektor u(x,t)
    @lastModified 2026-03-10
    """
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)

    # Fourier-Koeffizienten numerisch berechnen (Trapezregel)
    x_quad = np.linspace(0, L, 1000)
    u0_quad = np.sin(np.pi * x_quad / L)  # Beispiel-Anfangsbedingung

    for n in range(1, n_terms + 1):
        # Koeffizient: bₙ = (2/L) ∫ u0(x) sin(nπx/L) dx
        integrand = u0_quad * np.sin(n * np.pi * x_quad / L)
        bn = (2.0 / L) * np.trapezoid(integrand, x_quad)
        # Zeitabfall: e^{-(nπ/L)²t}
        lambda_n = (n * np.pi / L) ** 2
        result += bn * np.sin(n * np.pi * x / L) * np.exp(-lambda_n * t)

    return result


# =============================================================================
# Wellengleichung u_tt = c² u_xx
# =============================================================================

def wave_equation_explicit(
    u0: list,
    u0_t: list,
    L: float,
    T: float,
    c: float = 1.0,
    nx: int = 100,
    nt: int = 500
) -> np.ndarray:
    """
    @brief Leapfrog-Schema für die Wellengleichung.
    @description
        Löst u_tt = c² u_xx mit:
        - u(x,0) = u0 (Anfangsauslenkung)
        - u_t(x,0) = u0_t (Anfangsgeschwindigkeit)
        - Dirichlet-Rand: u(0,t)=u(L,t)=0

        Leapfrog-Schema (CFL-Bedingung: λ = cΔt/Δx ≤ 1):
        u_i^{n+1} = 2u_i^n - u_i^{n-1} + λ²(u_{i+1}^n - 2u_i^n + u_{i-1}^n)

    @param u0    Anfangsauslenkung (nx+1 Werte)
    @param u0_t  Anfangsgeschwindigkeit (nx+1 Werte)
    @param L     Länge des Intervalls
    @param T     Endzeit
    @param c     Wellengeschwindigkeit
    @param nx    Raumgitterpunkte
    @param nt    Zeitschritte
    @return      2D-Array (nt+1, nx+1): Lösung u(x,t)
    @lastModified 2026-03-10
    """
    dx = L / nx
    dt = T / nt
    # CFL-Zahl
    lam = c * dt / dx

    if lam > 1.0:
        import warnings
        warnings.warn(
            f"CFL-Bedingung verletzt: λ={lam:.4f} > 1. Instabilität möglich."
        )

    x = np.linspace(0, L, nx + 1)

    # Anfangsbedingungen auf Gitter
    u0_arr = np.array(u0, dtype=float)
    if len(u0_arr) != nx + 1:
        x0_arr = np.linspace(0, L, len(u0_arr))
        u0_arr = np.interp(x, x0_arr, u0_arr)

    u0t_arr = np.array(u0_t, dtype=float)
    if len(u0t_arr) != nx + 1:
        x0_arr = np.linspace(0, L, len(u0t_arr))
        u0t_arr = np.interp(x, x0_arr, u0t_arr)

    u = np.zeros((nt + 1, nx + 1))
    u[0, :] = u0_arr

    # Erster Zeitschritt mit Taylorentwicklung (halber Schritt für u_t(0))
    r2 = lam**2
    u[1, 1:-1] = (
        u[0, 1:-1]
        + dt * u0t_arr[1:-1]
        + 0.5 * r2 * (u[0, 2:] - 2 * u[0, 1:-1] + u[0, :-2])
    )
    u[1, 0] = 0.0
    u[1, nx] = 0.0

    # Leapfrog für n >= 1
    for n in range(1, nt):
        u[n + 1, 1:-1] = (
            2 * u[n, 1:-1]
            - u[n - 1, 1:-1]
            + r2 * (u[n, 2:] - 2 * u[n, 1:-1] + u[n, :-2])
        )
        u[n + 1, 0] = 0.0
        u[n + 1, nx] = 0.0

    return u


def wave_equation_analytical(
    x: np.ndarray,
    t: float,
    L: float = np.pi,
    c: float = 1.0,
    n_terms: int = 10
) -> np.ndarray:
    """
    @brief Analytische Fourier-Reihen-Lösung der Wellengleichung.
    @description
        Für u(x,0) = sin(πx/L), u_t(x,0) = 0, u(0,t)=u(L,t)=0:
        u(x,t) = Σ [aₙcos(nπct/L) + bₙsin(nπct/L)] sin(nπx/L)

        Nur a₁=1 ist von null verschieden für die gewählte Anfangsbedingung.
        u(x,t) = sin(πx/L) · cos(πct/L)

    @param x       Ortsvektor
    @param t       Zeitpunkt
    @param L       Intervalllänge
    @param c       Wellengeschwindigkeit
    @param n_terms Anzahl der Fourier-Terme
    @return        Lösungsvektor u(x,t)
    @lastModified 2026-03-10
    """
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)

    # Numerische Fourier-Koeffizienten (Anfangsbedingung sin(πx/L))
    x_quad = np.linspace(0, L, 1000)
    u0_quad = np.sin(np.pi * x_quad / L)
    u0t_quad = np.zeros_like(x_quad)  # Anfangsgeschwindigkeit = 0

    for n in range(1, n_terms + 1):
        omega_n = n * np.pi * c / L
        k_n = n * np.pi / L

        # aₙ = (2/L) ∫ u0(x) sin(nπx/L) dx
        an = (2.0 / L) * np.trapezoid(u0_quad * np.sin(k_n * x_quad), x_quad)
        # bₙ = (2/(L·ωₙ)) ∫ u0_t(x) sin(nπx/L) dx
        if abs(omega_n) > 1e-14:
            bn = (2.0 / (L * omega_n)) * np.trapezoid(
                u0t_quad * np.sin(k_n * x_quad), x_quad
            )
        else:
            bn = 0.0

        result += (an * np.cos(omega_n * t) + bn * np.sin(omega_n * t)) * np.sin(k_n * x)

    return result


def wave_dispersion_relation(k: float, c: float) -> float:
    """
    @brief Dispersionsrelation der Wellengleichung.
    @description
        Für u_tt = c²u_xx gilt ω = ck (nicht-dispersiv).
        Die Gruppengeschwindigkeit dω/dk = c = Phasengeschwindigkeit.

    @param k  Wellenzahl
    @param c  Wellengeschwindigkeit
    @return   Kreisfrequenz ω
    @lastModified 2026-03-10
    """
    return c * k


# =============================================================================
# Laplace/Poisson-Gleichung
# =============================================================================

def laplace_equation_2d(
    nx: int = 20,
    ny: int = 20,
    boundary: Optional[dict] = None
) -> np.ndarray:
    """
    @brief 2D Laplace-Gleichung via SOR (Successive Over-Relaxation).
    @description
        Löst Δu = 0 auf [0,1]×[0,1] mit vorgegebenen Dirichlet-Randwerten.
        SOR mit optimaler Relaxationsparameter ω_opt ≈ 2/(1+sin(π/N)).

        Finite-Differenzen-Schema (5-Punkt-Stern):
        u_{i,j} = (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}) / 4

    @param nx       Gitterpunkte in x-Richtung
    @param ny       Gitterpunkte in y-Richtung
    @param boundary dict mit 'top', 'bottom', 'left', 'right' (je als array)
    @return         2D-Array (ny+1, nx+1): Lösung u(x,y)
    @lastModified 2026-03-10
    """
    # Standardrandwerte
    if boundary is None:
        boundary = {
            'top': np.ones(nx + 1),
            'bottom': np.zeros(nx + 1),
            'left': np.zeros(ny + 1),
            'right': np.zeros(ny + 1)
        }

    u = np.zeros((ny + 1, nx + 1))

    # Randwerte setzen
    u[ny, :] = boundary.get('top', np.zeros(nx + 1))
    u[0, :] = boundary.get('bottom', np.zeros(nx + 1))
    u[:, 0] = boundary.get('left', np.zeros(ny + 1))
    u[:, nx] = boundary.get('right', np.zeros(ny + 1))

    # Optimaler SOR-Relaxationsparameter
    N = max(nx, ny)
    omega = 2.0 / (1.0 + np.sin(np.pi / N))

    # SOR-Iteration bis Konvergenz
    max_iter = 10000
    tol = 1e-8
    for _ in range(max_iter):
        u_old = u.copy()
        for j in range(1, ny):
            for i in range(1, nx):
                # Gauss-Seidel-Wert
                u_gs = 0.25 * (u[j, i+1] + u[j, i-1] + u[j+1, i] + u[j-1, i])
                # SOR-Update
                u[j, i] = (1.0 - omega) * u[j, i] + omega * u_gs

        # Konvergenzprüfung
        diff = np.max(np.abs(u - u_old))
        if diff < tol:
            break

    return u


def poisson_equation_2d(
    f: Callable,
    nx: int = 20,
    ny: int = 20
) -> np.ndarray:
    """
    @brief 2D Poisson-Gleichung via sparse Direktlösung.
    @description
        Löst -Δu = f(x,y) auf [0,1]×[0,1] mit u=0 am Rand.
        Diskretisierung mit 5-Punkt-Stern, direkte sparse LGS-Lösung.

        Gitterpunktindizes: innere Punkte (i,j) → globaler Index i*(ny-1)+j

    @param f   Quellterm f(x,y) als aufrufbare Funktion
    @param nx  Gitterpunkte in x-Richtung
    @param ny  Gitterpunkte in y-Richtung
    @return    2D-Array (ny+1, nx+1): Lösung u(x,y)
    @lastModified 2026-03-10
    """
    dx = 1.0 / nx
    dy = 1.0 / ny
    dx2 = dx**2
    dy2 = dy**2

    # Innere Gitterpunkte
    mx = nx - 1
    my = ny - 1
    N = mx * my  # Gesamtanzahl innerer Punkte

    # Sparse-Systemmatrix aufbauen
    rows, cols, vals = [], [], []

    def idx(i, j):
        """Globaler Index für inneren Punkt (i,j)."""
        return (i - 1) * my + (j - 1)

    for i in range(1, nx):
        for j in range(1, ny):
            k = idx(i, j)
            # Hauptdiagonal
            rows.append(k)
            cols.append(k)
            vals.append(2.0 / dx2 + 2.0 / dy2)
            # Nachbar links
            if i > 1:
                rows.append(k)
                cols.append(idx(i - 1, j))
                vals.append(-1.0 / dx2)
            # Nachbar rechts
            if i < nx - 1:
                rows.append(k)
                cols.append(idx(i + 1, j))
                vals.append(-1.0 / dx2)
            # Nachbar unten
            if j > 1:
                rows.append(k)
                cols.append(idx(i, j - 1))
                vals.append(-1.0 / dy2)
            # Nachbar oben
            if j < ny - 1:
                rows.append(k)
                cols.append(idx(i, j + 1))
                vals.append(-1.0 / dy2)

    A_mat = sparse.csc_matrix((vals, (rows, cols)), shape=(N, N))

    # Rechte Seite: Quellterm auswerten
    x_arr = np.linspace(0, 1, nx + 1)
    y_arr = np.linspace(0, 1, ny + 1)
    rhs = np.zeros(N)
    for i in range(1, nx):
        for j in range(1, ny):
            rhs[idx(i, j)] = f(x_arr[i], y_arr[j])

    # Lösen
    u_inner = spsolve(A_mat, rhs)

    # Lösung in 2D-Array eintragen
    u = np.zeros((ny + 1, nx + 1))
    for i in range(1, nx):
        for j in range(1, ny):
            u[j, i] = u_inner[idx(i, j)]

    return u


def green_function_laplace_2d(
    x0: float,
    y0: float,
    x: np.ndarray,
    y: np.ndarray
) -> np.ndarray:
    """
    @brief Fundamentallösung (Greens Funktion) des Laplace-Operators in 2D.
    @description
        G(x,y;x₀,y₀) = -1/(2π) · ln(r)
        mit r = √((x-x₀)² + (y-y₀)²).

        Singularität bei r=0 wird durch kleines ε regularisiert.

    @param x0  Quellpunkt x-Koordinate
    @param y0  Quellpunkt y-Koordinate
    @param x   Feldpunkte x (1D-Array)
    @param y   Feldpunkte y (1D-Array)
    @return    2D-Array G(x_i, y_j; x0, y0)
    @lastModified 2026-03-10
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    eps = 1e-14  # Regularisierung der Singularität
    # Abstandsmatrix
    X, Y = np.meshgrid(x, y)
    r = np.sqrt((X - x0)**2 + (Y - y0)**2) + eps
    return -1.0 / (2.0 * np.pi) * np.log(r)


def harmonic_function_properties(u: np.ndarray) -> dict:
    """
    @brief Überprüft numerisch Eigenschaften harmonischer Funktionen.
    @description
        Harmonische Funktionen (Δu = 0) erfüllen:
        1. Mittelwerteigenschaft: u(x₀) = Mittel auf jedem Kreis/Kugel
        2. Maximum-Prinzip: Maximum und Minimum werden am Rand angenommen

    @param u  2D-Array der Lösung (approximativ harmonisch)
    @return   dict mit 'max_interior', 'max_boundary', 'satisfies_max_principle',
                       'mean_value_error'
    @lastModified 2026-03-10
    """
    ny, nx = u.shape

    # Randwerte
    boundary_vals = np.concatenate([
        u[0, :], u[-1, :], u[1:-1, 0], u[1:-1, -1]
    ])
    interior_vals = u[1:-1, 1:-1].flatten()

    max_boundary = np.max(boundary_vals)
    min_boundary = np.min(boundary_vals)
    max_interior = np.max(interior_vals) if len(interior_vals) > 0 else float('nan')
    min_interior = np.min(interior_vals) if len(interior_vals) > 0 else float('nan')

    # Maximum-Prinzip: Inneres Maximum ≤ Rand-Maximum (mit Toleranz)
    tol = 1e-6
    satisfies_max = (max_interior <= max_boundary + tol)
    satisfies_min = (min_interior >= min_boundary - tol)

    # Mittelwerteigenschaft: Mittelpunkt vs. Mittel aller Punkte
    # Grober Test: Mittelpunkt ≈ Mittel über alle Gitterpunkte (nur für spezielle Fälle exakt)
    center_i = ny // 2
    center_j = nx // 2
    center_val = u[center_i, center_j]
    mean_val = np.mean(u)
    mean_value_error = abs(center_val - mean_val)

    return {
        'max_interior': float(max_interior),
        'max_boundary': float(max_boundary),
        'min_interior': float(min_interior),
        'min_boundary': float(min_boundary),
        'satisfies_max_principle': bool(satisfies_max),
        'satisfies_min_principle': bool(satisfies_min),
        'mean_value_error': float(mean_value_error),
    }


# =============================================================================
# Schrödinger-Gleichung
# =============================================================================

def schrodinger_stationary(
    V: Callable,
    x: np.ndarray,
    n_states: int = 5
) -> dict:
    """
    @brief Stationäre Schrödinger-Gleichung via Finite-Differenzen.
    @description
        Löst das Eigenwertproblem:
        -ℏ²/(2m) ψ'' + V(x)ψ = Eψ  (Einheiten: ℏ=m=1)

        Diskretisierung: -ψ_{i+1} + 2ψ_i - ψ_{i-1} = Δx²(E-V_i)ψ_i
        → Tridiagonales Eigenwertproblem (symmetrisch).
        Lösung via numpy.linalg.eigh (stabil, effizient).

    @param V        Potential V(x) als Funktion
    @param x        Ortsgitter (1D-Array)
    @param n_states Anzahl der gesuchten Eigenzustände
    @return dict mit 'energies' (E_n) und 'wavefunctions' (ψ_n)
    @lastModified 2026-03-10
    """
    n = len(x)
    dx = x[1] - x[0]
    dx2 = dx**2

    # Potential auf dem Gitter
    V_arr = np.array([V(xi) for xi in x], dtype=float)

    # Tridiagonale Hamilton-Matrix aufbauen (ℏ=m=1)
    # H_ii = 1/dx² + V_i, H_{i,i±1} = -1/(2dx²)
    diag = 1.0 / dx2 + V_arr
    off_diag = np.full(n - 1, -0.5 / dx2)

    # Eigenwertlösung (nur die ersten n_states)
    # eigh_tridiagonal nutzt LAPACK und ist effizient
    n_states_eff = min(n_states, n - 2)
    energies, wavefunctions = eigh_tridiagonal(
        diag, off_diag,
        select='i',
        select_range=(0, n_states_eff - 1)
    )

    # Normierung: ∫|ψ|²dx = 1 (Trapezregel)
    psi_list = []
    for k in range(n_states_eff):
        psi = wavefunctions[:, k]
        norm = np.sqrt(np.trapezoid(np.abs(psi)**2, x))
        if norm > 1e-14:
            psi = psi / norm
        psi_list.append(psi)

    return {
        'energies': energies,
        'wavefunctions': psi_list,
        'x': x,
    }


def schrodinger_harmonic_oscillator(n_states: int = 5) -> dict:
    """
    @brief Harmonischer Oszillator: numerisch vs. analytisch.
    @description
        Potential V(x) = x²/2, analytische Lösung:
        E_n = n + 1/2 (n = 0,1,2,...)
        ψ_n(x) = (2^n n! √π)^{-1/2} H_n(x) e^{-x²/2}

        H_n: Hermite-Polynome (physikalisch).

    @param n_states Anzahl der Zustände
    @return dict mit 'energies_numerical', 'energies_analytical', 'errors', 'wavefunctions'
    @lastModified 2026-03-10
    """
    from scipy.special import hermite, factorial

    # Numerisches Gitter (groß genug, damit Randwerte ≈ 0)
    x = np.linspace(-6, 6, 400)

    # Stationäre Schrödinger für V(x)=x²/2
    result = schrodinger_stationary(lambda xi: 0.5 * xi**2, x, n_states)
    e_num = result['energies']
    psi_num = result['wavefunctions']

    # Analytische Energieeigenwerte: E_n = n + 1/2
    e_ana = np.array([n + 0.5 for n in range(n_states)])

    # Fehler
    n_compare = min(len(e_num), len(e_ana))
    errors = np.abs(e_num[:n_compare] - e_ana[:n_compare])

    return {
        'energies_numerical': e_num,
        'energies_analytical': e_ana[:n_compare],
        'errors': errors,
        'wavefunctions': psi_num,
        'x': x,
    }


def schrodinger_time_dependent(
    psi0: np.ndarray,
    V: np.ndarray,
    x: np.ndarray,
    T: float,
    nt: int = 200
) -> np.ndarray:
    """
    @brief Zeitabhängige Schrödinger-Gleichung via Crank-Nicolson.
    @description
        Löst iℏ ∂ψ/∂t = Ĥψ = [-ℏ²/(2m)∂²/∂x² + V(x)]ψ (ℏ=m=1)
        Crank-Nicolson: unitäres Schema, erhält die Norm ‖ψ‖²=1.

        Schema: (I + iΔt/2·H) ψ^{n+1} = (I - iΔt/2·H) ψ^n

    @param psi0  Anfangszustand ψ(x,0) (normiert)
    @param V     Potential V(x_i) als Array
    @param x     Ortsgitter
    @param T     Gesamtzeit
    @param nt    Anzahl der Zeitschritte
    @return      2D-Array (nt+1, nx): Wellenfunktion ψ(x,t) (komplex)
    @lastModified 2026-03-10
    """
    n = len(x)
    dt = T / nt
    dx = x[1] - x[0]
    dx2 = dx**2

    # Diagonale des Hamilton-Operators
    h_diag = 1.0 / dx2 + np.asarray(V, dtype=float)
    h_off = np.full(n - 1, -0.5 / dx2)

    # Sparse-Matrizen für Crank-Nicolson
    # A = I + i*dt/2 * H (linke Seite)
    # B = I - i*dt/2 * H (rechte Seite, expliziter Teil)
    alpha = 0.5j * dt
    a_diag = 1.0 + alpha * h_diag
    a_off = alpha * h_off
    A_mat = sparse.diags(
        [a_off, a_diag, a_off],
        offsets=[-1, 0, 1],
        format='csc',
        dtype=complex
    )

    # Lösungs-Array
    psi = np.zeros((nt + 1, n), dtype=complex)
    psi[0, :] = np.asarray(psi0, dtype=complex)

    for step in range(nt):
        psi_cur = psi[step, :]
        # Rechte Seite: (I - i*dt/2 * H) ψ
        rhs = psi_cur.copy()
        rhs[:-1] -= alpha * h_off * psi_cur[1:]
        rhs[1:] -= alpha * h_off * psi_cur[:-1]
        rhs -= alpha * h_diag * psi_cur
        rhs += psi_cur  # Doppelt addiert, korrigieren:
        # Korrekte Berechnung:
        rhs2 = psi_cur - alpha * (
            h_diag * psi_cur
            + np.concatenate([[0], h_off * psi_cur[:-1]])
            + np.concatenate([h_off * psi_cur[1:], [0]])
        )
        psi[step + 1, :] = spsolve(A_mat, rhs2)

    return psi


# =============================================================================
# Finite-Elemente-Methode (FEM) 1D
# =============================================================================

def fem_1d_poisson(
    f: Callable,
    a: float,
    b: float,
    n: int = 20
) -> dict:
    """
    @brief 1D FEM für -u'' = f(x), u(a)=u(b)=0.
    @description
        Galerkin-Methode mit stückweise linearen Basisfunktionen (Hutfunktionen).

        Schwache Formulierung:
        ∫ u' v' dx = ∫ f v dx  ∀v im Test-Raum

        Steifigkeitsmatrix K_ij = ∫ φ_i' φ_j' dx
        Massenvektor F_i = ∫ f(x) φ_i(x) dx

        Element-Beitrag pro Intervall [x_k, x_{k+1}], h=Δx:
        K_loc = (1/h) [[1,-1],[-1,1]], F_loc = (h/6)[2f(x_k)+f(x_{k+1}), f(x_k)+2f(x_{k+1})]

    @param f  Quellterm f(x)
    @param a  Linker Randpunkt
    @param b  Rechter Randpunkt
    @param n  Anzahl der Elemente
    @return   dict mit 'x', 'u' (FEM-Lösung), 'nodes'
    @lastModified 2026-03-10
    """
    # Knotenpunkte (n+1 Punkte)
    x = np.linspace(a, b, n + 1)
    h = (b - a) / n  # Elementlänge (gleichmäßig)

    # Steifigkeitsmatrix (Tridiagonal)
    K = np.zeros((n + 1, n + 1))
    F = np.zeros(n + 1)

    for k in range(n):
        xk = x[k]
        xk1 = x[k + 1]
        hk = xk1 - xk

        # Lokale Steifigkeitsmatrix
        K_loc = np.array([[1.0, -1.0], [-1.0, 1.0]]) / hk

        # Lokaler Lastvektor (Gauss-Quadratur mit 3 Punkten)
        f_mid = f(0.5 * (xk + xk1))
        fk = f(xk)
        fk1 = f(xk1)
        # Trapezregel: F_loc_i = ∫ f φ_i dx ≈ hk/6 · [2fk+fk1, fk+2fk1]
        F_loc = hk / 6.0 * np.array([2 * fk + fk1, fk + 2 * fk1])

        # Assemblierung
        K[k:k+2, k:k+2] += K_loc
        F[k:k+2] += F_loc

    # Dirichlet-Randwerte u(a)=u(b)=0 eintragen
    # Erste und letzte Gleichung ersetzen
    K[0, :] = 0.0
    K[0, 0] = 1.0
    F[0] = 0.0
    K[-1, :] = 0.0
    K[-1, -1] = 1.0
    F[-1] = 0.0

    # LGS lösen
    u = np.linalg.solve(K, F)

    return {
        'x': x,
        'u': u,
        'nodes': n + 1,
        'h': h,
    }


def fem_convergence_demo() -> dict:
    """
    @brief Zeigt O(h²)-Konvergenz der FEM für -u''=f.
    @description
        Testet FEM für -u'' = π²sin(πx), u(0)=u(1)=0.
        Exakte Lösung: u(x) = sin(πx).
        Fehler in L²-Norm: ‖u_h - u‖ = O(h²).

    @return dict mit 'h_values', 'errors', 'order' (mittlere Konvergenzordnung)
    @lastModified 2026-03-10
    """
    # Quellterm und analytische Lösung
    f = lambda x: np.pi**2 * np.sin(np.pi * x)
    u_exact = lambda x: np.sin(np.pi * x)

    h_values = []
    errors = []

    for n in [4, 8, 16, 32, 64, 128]:
        result = fem_1d_poisson(f, 0.0, 1.0, n)
        x = result['x']
        u_h = result['u']
        h = result['h']

        # L²-Fehler (Trapezregel)
        err = np.sqrt(np.trapezoid((u_h - u_exact(x))**2, x))
        h_values.append(h)
        errors.append(err)

    h_arr = np.array(h_values)
    e_arr = np.array(errors)

    # Konvergenzordnung schätzen (log-log-Steigung)
    if len(h_arr) > 1:
        log_h = np.log(h_arr)
        log_e = np.log(e_arr + 1e-16)
        order = np.polyfit(log_h, log_e, 1)[0]
    else:
        order = float('nan')

    return {
        'h_values': h_arr.tolist(),
        'errors': e_arr.tolist(),
        'order': float(order),
    }


# =============================================================================
# Methode der Charakteristiken
# =============================================================================

def method_of_characteristics(
    a: float,
    b: float,
    c: float,
    u0: Callable,
    x_range: tuple,
    t_range: tuple,
    nx: int = 50,
    nt: int = 50
) -> dict:
    """
    @brief Methode der Charakteristiken für au_x + bu_t = c·u.
    @description
        Charakteristikenkurven: dx/ds = a, dt/ds = b.
        Entlang der Charakteristiken: du/ds = c·u → u = u0·e^{cs/b}.

        Für a=1, b=1 (Advektionsgleichung): u(x,t) = u0(x-at/b)·e^{ct/b}

    @param a       Koeffizient vor u_x
    @param b       Koeffizient vor u_t (Zeitkoeffizient)
    @param c       Koeffizient in der rechten Seite c·u
    @param u0      Anfangsbedingung u(x,0)
    @param x_range (x_min, x_max)
    @param t_range (t_min, t_max)
    @param nx      Gitterpunkte in x
    @param nt      Gitterpunkte in t
    @return dict mit 'x', 't', 'u' (2D-Lösung), 'characteristics'
    @lastModified 2026-03-10
    """
    x = np.linspace(x_range[0], x_range[1], nx)
    t_arr = np.linspace(t_range[0], t_range[1], nt)

    u = np.zeros((nt, nx))

    for j, t in enumerate(t_arr):
        for i, xi in enumerate(x):
            # Charakteristik rückwärts: x0 = xi - (a/b)*t
            if abs(b) > 1e-14:
                x0 = xi - (a / b) * t
                # Lösung: u = u0(x0) * e^{(c/b)*t}
                u[j, i] = u0(x0) * np.exp((c / b) * t)
            else:
                u[j, i] = u0(xi)

    # Beispiel-Charakteristiken (3 Linien)
    characteristics = []
    for x_start in np.linspace(x_range[0], x_range[1], 5):
        char_x = x_start + (a / b) * t_arr if abs(b) > 1e-14 else np.full_like(t_arr, x_start)
        characteristics.append({'t': t_arr.tolist(), 'x': char_x.tolist()})

    return {
        'x': x,
        't': t_arr,
        'u': u,
        'characteristics': characteristics,
    }


def burgers_equation(
    u0: list,
    L: float = 2 * np.pi,
    T: float = 0.5,
    nx: int = 100,
    nt: int = 200,
    nu: float = 0.01
) -> np.ndarray:
    """
    @brief Viskose Burgers-Gleichung u_t + u·u_x = ν·u_xx.
    @description
        Modellproblem für nichtlineare Wellenphänomene und Schockbildung.
        Kombination aus Advektion (nichtlinear) und Diffusion.

        Explizites Upwind-Schema:
        - Advektion: u_x ≈ (u_i - u_{i-1})/Δx für u>0 (upwind)
        - Diffusion: u_xx ≈ (u_{i+1} - 2u_i + u_{i-1})/Δx²

        Randbedingungen: periodisch (Fourier-Modell).

    @param u0  Anfangsbedingung (nx Werte, periodisch)
    @param L   Periodenlänge
    @param T   Endzeit
    @param nx  Raumgitterpunkte
    @param nt  Zeitschritte
    @param nu  kinematische Viskosität
    @return    2D-Array (nt+1, nx): Lösung u(x,t)
    @lastModified 2026-03-10
    """
    dx = L / nx
    dt = T / nt

    # Anfangsbedingung auf Gitter
    u0_arr = np.array(u0, dtype=float)
    if len(u0_arr) != nx:
        x_src = np.linspace(0, L, len(u0_arr))
        x_dst = np.linspace(0, L, nx)
        u0_arr = np.interp(x_dst, x_src, u0_arr)

    u = np.zeros((nt + 1, nx))
    u[0, :] = u0_arr

    for n in range(nt):
        un = u[n, :]

        # Nichtlinearer Advektionsterm (upwind mit periodischen Grenzen)
        u_adv = np.where(
            un >= 0,
            (un - np.roll(un, 1)) / dx,   # Vorwärtsdifferenz (upwind für u>0)
            (np.roll(un, -1) - un) / dx   # Rückwärtsdifferenz (upwind für u<0)
        )

        # Diffusionsterm (zentrale Differenz)
        u_diff = (np.roll(un, -1) - 2 * un + np.roll(un, 1)) / dx**2

        # Euler-Schritt
        u[n + 1, :] = un + dt * (-un * u_adv + nu * u_diff)

    return u
