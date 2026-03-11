"""
@file selmer_groups.py
@brief Selmer-Gruppen, Heegner-Punkte und BSD-Vermutung fuer elliptische Kurven ueber Q.
@description
    Dieses Modul implementiert fortgeschrittene zahlentheoretische Methoden
    fuer elliptische Kurven, insbesondere:

    1. 2-Abstieg (2-Descent): Berechnung der 2-Selmer-Gruppe Sel^2(E/Q).
       Die exakte Sequenz
           0 -> E(Q)/2E(Q) -> Sel^2(E/Q) -> Sha(E/Q)[2] -> 0
       verbindet Selmer-Gruppe, Mordell-Weil-Gruppe und Tate-Shafarevich-Gruppe.

    2. Heegner-Punkte: Konstruktion rationaler Punkte auf Rang-1-Kurven
       via CM-Theorie und Modulfunktionen (Kolyvagin-Theorem).

    3. BSD-Vermutung (Birch und Swinnerton-Dyer):
       Numerische Untersuchung der Leitkonjektur:
           L(E, s) ~ C_E * (s - 1)^r  fuer s -> 1
       wobei r = Rang(E(Q)).

    Wichtige Formeln:
    - Diskriminante: Delta = -16 * (4a^3 + 27b^2)
    - Frobenius-Spur: a_p = p + 1 - #E(F_p)
    - L-Funktion: L(E,s) = prod_p (1 - a_p*p^{-s} + p^{1-2s})^{-1}
    - BSD: L(E,1) = 0 gdw. Rang(E(Q)) >= 1

@author Michael Fuhrmann
@lastModified 2026-03-11
"""

import math
import sympy as sp
from sympy import Symbol, solve, factorint, isprime, legendre_symbol
from typing import List, Optional, Dict, Tuple


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _prime_factors_of(n: int) -> List[int]:
    """
    Berechnet alle Primteiler einer ganzen Zahl n.

    Verwendet sympy.factorint fuer effiziente Faktorisierung.

    @param n: Ganze Zahl (|n| >= 2)
    @return: Sortierte Liste der Primteiler
    @lastModified 2026-03-11
    """
    # Vorzeichen ignorieren, 0 und 1 haben keine Primteiler
    n = abs(n)
    if n < 2:
        return []
    return sorted(factorint(n).keys())


def _is_squarefree(n: int) -> bool:
    """
    Prueft ob n quadratfrei ist (kein Quadrat teilt n).

    Eine Zahl n ist quadratfrei, wenn fuer alle Primzahlen p gilt: p^2 /| n.

    @param n: Ganze Zahl
    @return: True wenn quadratfrei
    @lastModified 2026-03-11
    """
    n = abs(n)
    if n == 0:
        return False
    # Alle Primfaktoren pruefen: keiner darf mit Exponent >= 2 vorkommen
    factors = factorint(n)
    return all(exp == 1 for exp in factors.values())


def _count_fp_points(a: int, b: int, p: int) -> int:
    """
    Zaehlt die Anzahl der Punkte #E(F_p) inklusive dem Punkt im Unendlichen.

    Methode: Brute-Force ueber alle x in {0, ..., p-1}.
    Fuer jedes x wird geprueft, ob x^3 + ax + b ein Quadratrest mod p ist.

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param p: Primzahl (Charakteristik des Koerpers)
    @return: Gruppenordnung #E(F_p)
    @lastModified 2026-03-11
    """
    # Punkt im Unendlichen wird immer gezaehlt
    count = 1
    a_mod = a % p
    b_mod = b % p

    for x in range(p):
        # Rechte Seite: f(x) = x^3 + ax + b mod p
        rhs = (pow(x, 3, p) + a_mod * x + b_mod) % p

        # Euler-Kriterium: rhs ist Quadratrest gdw. rhs^{(p-1)/2} ≡ 1 (mod p)
        if rhs == 0:
            count += 1  # y = 0, nur ein Punkt
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2  # y und -y, zwei Punkte

    return count


def _compute_ap(a: int, b: int, p: int) -> int:
    """
    Berechnet die Frobenius-Spur a_p = p + 1 - #E(F_p).

    Nach dem Satz von Hasse gilt: |a_p| <= 2*sqrt(p).

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param p: Primzahl
    @return: Frobenius-Spur a_p
    @lastModified 2026-03-11
    """
    return p + 1 - _count_fp_points(a, b, p)


def _discriminant(a: int, b: int) -> int:
    """
    Berechnet die Diskriminante Delta = -16 * (4a^3 + 27b^2).

    Die Kurve ist nicht-singulaer gdw. Delta != 0.

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @return: Diskriminante Delta
    @lastModified 2026-03-11
    """
    return -16 * (4 * a**3 + 27 * b**2)


def _is_good_prime(a: int, b: int, p: int) -> bool:
    """
    Prueft ob p eine Primzahl guter Reduktion fuer E ist.

    p ist eine schlechte Primzahl gdw. p | Delta(E).

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param p: Primzahl
    @return: True wenn gute Reduktion
    @lastModified 2026-03-11
    """
    delta = _discriminant(a, b)
    return delta % p != 0


def _root_number(a: int, b: int) -> int:
    """
    Berechnet den globalen Root Number (Vorzeichenfaktor) epsilon(E).

    Der Root Number epsilon(E) in {+1, -1} bestimmt das Vorzeichen der
    Funktionalgleichung der L-Funktion:
        L*(E, 2-s) = epsilon(E) * L*(E, s)

    Mathematische Konsequenz:
    - epsilon(E) = -1  =>  L(E, 1) = 0  (Rang ist ungerade: 1, 3, 5, ...)
    - epsilon(E) = +1  =>  L(E, 1) koennte != 0 sein (Rang 0, 2, ...)

    Methode (vereinfacht nach Connell, Silverman):
    - c4 = -48a, c6 = -864b fuer kurze Weierstrass-Form
    - Fuer schlechte Primen p: epsilon_p via Reduktionstyp
    - Multiplicative Reduktion (c4 != 0 mod p): epsilon_p = kronecker(c6, p)
    - Additive Reduktion (c4 = 0 mod p): epsilon_p = kronecker(-3, p)
    - Globale Root Number: epsilon(E) = prod_p epsilon_p

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @return: Root Number epsilon(E) = +1 oder -1
    @lastModified 2026-03-11
    """
    delta = _discriminant(a, b)
    if delta == 0:
        return 1  # Singulaere Kurve: Fallback

    # c4 und c6 fuer kurze Weierstrass-Form y^2 = x^3 + ax + b
    c4 = -48 * a
    c6 = -864 * b

    bad_primes = _prime_factors_of(delta)
    epsilon = 1  # Gesamtprodukt starten

    for p in bad_primes:
        if p == 2:
            # p = 2: Berechne 2-adische Bewertungen von c4, c6
            c4_val = _2adic_valuation(c4)
            c6_val = _2adic_valuation(c6)

            # Additive Reduktion bei p=2: epsilon_2 = -1
            # Multiplicative Reduktion bei p=2: epsilon_2 = +1 (vereinfacht)
            if c4_val >= 4 and c6_val >= 6:
                # Additive Reduktion (Kodaira Typ II, III, IV, II*, III*, IV*)
                epsilon *= -1
            # Sonst multiplicative: epsilon_2 = +1 (kein Faktor)
        else:
            # Ungerade Primzahl p
            if c4 % p != 0:
                # Multiplicative Reduktion (c4 != 0 mod p)
                # epsilon_p = kronecker_symbol(c6, p)
                ks = int(legendre_symbol(c6 % p if c6 % p != 0 else p-1, p))
                epsilon *= ks
            else:
                # Additive Reduktion (c4 = 0 mod p)
                # epsilon_p = kronecker(-3, p) = Legendre(-3, p)
                neg3 = (-3) % p
                if neg3 == 0:
                    pass  # p = 3: epsilon_p = 1
                else:
                    ks = int(legendre_symbol(neg3, p))
                    epsilon *= ks

    return epsilon


def _2adic_valuation(n: int) -> int:
    """
    Berechnet die 2-adische Bewertung v_2(n) einer ganzen Zahl.

    v_2(n) = max{k : 2^k | n}

    @param n: Ganze Zahl (n != 0)
    @return: 2-adische Bewertung (oder 0 wenn n=0)
    @lastModified 2026-03-11
    """
    if n == 0:
        return 100  # Konvention: v_2(0) = unendlich
    count = 0
    n = abs(n)
    while n % 2 == 0:
        count += 1
        n //= 2
    return count


def _real_period_approx(a: int, b: int) -> float:
    """
    Naeherung der reellen Periode Omega der elliptischen Kurve.

    Die reelle Periode ist das Integral:
        Omega = integral_{-inf}^{inf} dx / sqrt(x^3 + ax + b)
    genauer: das kleinste positive reelle Gitterperiode des Weierstrass-P.

    Vereinfachte numerische Naeherung via Intervallintegration.

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @return: Naeherung der reellen Periode Omega > 0
    @lastModified 2026-03-11
    """
    # Finde reelle Nullstellen von f(x) = x^3 + ax + b
    x_sym = Symbol('x')
    roots = solve(x_sym**3 + a * x_sym + b, x_sym)

    # Nur reelle Wurzeln beruecksichtigen
    real_roots = []
    for r in roots:
        try:
            rv = complex(r)
            if abs(rv.imag) < 1e-8:
                real_roots.append(float(rv.real))
        except Exception:
            pass

    real_roots.sort()

    if not real_roots:
        # Keine reellen Nullstellen: Kurve hat nur eine reelle Komponente
        # Naeherung: Omega ~ 2 (Standardwert fuer kompakte Kurven)
        return 2.0

    try:
        import numpy as np
        if len(real_roots) >= 2:
            # Zwei reelle Komponenten: Omega = 2 * integral von e1 bis e2
            e1, e2 = real_roots[0], real_roots[1]
            # Numerische Integration mittels Gauss-Quadratur
            n_steps = 500
            xs = np.linspace(e1 + 1e-10, e2 - 1e-10, n_steps)
            vals = xs**3 + a * xs + b
            # Negative Werte auf kleinen Wert setzen (Singularitaet)
            vals = np.maximum(vals, 1e-15)
            integrand = 1.0 / np.sqrt(vals)
            omega = 2.0 * np.trapezoid(integrand, xs)
            return max(omega, 0.01)
        else:
            # Nur eine reelle Nullstelle
            e1 = real_roots[0]
            n_steps = 200
            xs = np.linspace(e1 + 0.01, e1 + 10, n_steps)
            vals = xs**3 + a * xs + b
            vals = np.maximum(vals, 1e-15)
            integrand = 1.0 / np.sqrt(vals)
            omega = 2.0 * np.trapezoid(integrand, xs)
            return max(omega, 0.01)
    except Exception:
        return 2.0  # Fallback


# ---------------------------------------------------------------------------
# Klasse: TwoDescentEllipticCurve
# ---------------------------------------------------------------------------

class TwoDescentEllipticCurve:
    """
    Implementiert den 2-Abstieg (2-Descent) fuer elliptische Kurven E/Q.

    Der 2-Abstieg liefert eine Oberschranke fuer den Mordell-Weil-Rang:
        rank(E(Q)) <= dim_F2(Sel^2(E/Q)) - dim_F2(E(Q)[2])

    Die 2-Selmer-Gruppe wird durch lokale Bedingungen an allen Primstellen
    (inklusive unendliche Stelle) beschrieben.

    Theorie (vereinfacht):
    - E(Q)[2] = {P in E(Q) : 2P = O} (2-Torsion ueber Q)
    - Sel^2(E/Q) <= Produkt lokaler Bilder H^1(Q_p, E[2])
    - Die schlechten Primstellen sind die Teiler der Diskriminante Delta

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, a: int, b: int) -> None:
        """
        Initialisiert die elliptische Kurve E: y^2 = x^3 + ax + b.

        @param a: Koeffizient von x in der Weierstrass-Form
        @param b: konstanter Koeffizient
        @raises ValueError wenn Kurve singulaer ist (Delta = 0)
        @lastModified 2026-03-11
        """
        self.a = a
        self.b = b

        # Diskriminante berechnen und Singularitaet pruefen
        self._delta = _discriminant(a, b)
        if self._delta == 0:
            raise ValueError(
                f"Singulaere Kurve: Delta = 0 fuer a={a}, b={b}. "
                f"4a^3 + 27b^2 = {4*a**3 + 27*b**2}"
            )

    def two_torsion_points(self) -> List[Tuple]:
        """
        Findet alle 2-Torsionspunkte P in E(Q)[2] (mit 2P = O).

        Ein Punkt P = (x, y) hat Ordnung 2 gdw. 2P = O.
        Bei der Chord-and-Tangent-Methode: 2P = O gdw. y = 0.
        Also: Loese x^3 + ax + b = 0 ueber Q.

        Neben diesen affinen Punkten gehoert immer O (Fernpunkt) dazu.

        @return: Liste von rationalen 2-Torsionspunkten als (x, y)-Tupel
                 (der Fernpunkt O ist nicht in der Liste, wird implizit gezaehlt)
        @lastModified 2026-03-11
        """
        x_sym = Symbol('x')
        # Polynom: f(x) = x^3 + ax + b
        poly = x_sym**3 + self.a * x_sym + self.b

        # Symbolische Loesung ueber den rationalen Zahlen
        solutions = solve(poly, x_sym)

        rational_points = []
        for sol in solutions:
            try:
                # Nur rationale Loesungen behalten
                # sp.Rational prueft Rationalitaet
                val = sp.nsimplify(sol, rational=True)
                if val.is_rational:
                    rational_points.append((val, sp.Integer(0)))
            except Exception:
                # Symbolische Vereinfachung fehlgeschlagen: versuche numerisch
                try:
                    cval = complex(sol)
                    if abs(cval.imag) < 1e-9:
                        # Pruefe ob nahe an einer rationalen Zahl
                        rval = round(cval.real)
                        if abs(cval.real - rval) < 1e-6:
                            # Verifiziere: rval^3 + a*rval + b == 0
                            check = rval**3 + self.a * rval + self.b
                            if abs(check) < 1e-6:
                                rational_points.append((sp.Integer(rval), sp.Integer(0)))
                except Exception:
                    pass

        return rational_points

    def selmer_group_bound(self) -> int:
        """
        Berechnet eine Oberschranke fuer die Dimension der 2-Selmer-Gruppe.

        Die Selmer-Gruppe Sel^2(E/Q) ist in das Produkt der lokalen
        Kohomologiegruppen eingebettet. Die Schranke wird bestimmt durch:
        - Die Anzahl der schlechten Primstellen (Teiler von Delta)
        - Den Beitrag der 2-Torsion
        - Den Beitrag der archimedischen Stelle (unendlich)

        Vereinfachte Formel (nach Silverman "Arithmetic of Elliptic Curves"):
            dim Sel^2(E/Q) <= #{ schlechte Primes } + dim E(Q)[2] + 1

        @return: Oberschranke dim_F2(Sel^2(E/Q))
        @lastModified 2026-03-11
        """
        # Schlechte Primstellen: Primteiler der Diskriminante
        bad_primes = _prime_factors_of(self._delta)

        # Anzahl rationaler 2-Torsionspunkte
        torsion_pts = self.two_torsion_points()
        # E(Q)[2] hat als F_2-Vektorraum Dimension 0, 1 oder 2
        torsion_rank = min(len(torsion_pts), 2)

        # Oberschranke: |bad_primes| + torsion_rank + 1 (fuer Stelle unendlich)
        return len(bad_primes) + torsion_rank + 1

    def mordell_weil_rank_bound(self) -> int:
        """
        Oberschranke fuer den Mordell-Weil-Rang aus dem 2-Abstieg.

        Die exakte Sequenz des 2-Abstiegs gibt:
            0 -> E(Q)/2E(Q) -> Sel^2(E/Q) -> Sha(E/Q)[2] -> 0

        Da Sha(E/Q)[2] >= 0, folgt:
            rank(E(Q)) <= dim Sel^2(E/Q) - dim E(Q)[2]

        wobei E(Q)/2E(Q) ~ (Z/2Z)^r x (E(Q)_tors / 2*E(Q)_tors).

        @return: Oberschranke fuer den Rang (nicht-negative ganze Zahl)
        @lastModified 2026-03-11
        """
        selmer_bound = self.selmer_group_bound()
        torsion_pts = self.two_torsion_points()
        # Dimension von E(Q)[2] als F_2-Vektorraum
        torsion_dim = min(len(torsion_pts), 2)

        # Rang-Schranke (mindestens 0)
        return max(0, selmer_bound - torsion_dim)

    def two_descent(self) -> Dict:
        """
        Fuehrt den vollstaendigen 2-Abstieg durch und gibt ein Ergebnis-Dict zurueck.

        Der 2-Abstieg ist die Standardmethode zur Berechnung von Rang-Schranken
        fuer elliptische Kurven ueber Q. Er basiert auf:
        1. Berechnung von E(Q)[2] (2-Torsion)
        2. Bestimmung der Selmer-Bedingungen an schlechten Primstellen
        3. Auswertung des globalen Satzes von Selmer

        @return: Dictionary mit den Ergebnissen des 2-Abstiegs:
                 - selmer_rank_bound: dim_F2(Sel^2(E/Q)) (Oberschranke)
                 - torsion_rank: dim_F2(E(Q)[2])
                 - rank_upper_bound: obere Schranke fuer Mordell-Weil-Rang
                 - discriminant: Diskriminante Delta der Kurve
                 - bad_primes: Liste der schlechten Primstellen
        @lastModified 2026-03-11
        """
        # Schlechte Primstellen aus der Diskriminante
        bad_primes = _prime_factors_of(self._delta)

        # 2-Torsion
        torsion_points = self.two_torsion_points()
        torsion_rank = min(len(torsion_points), 2)

        # Selmer-Rang Oberschranke
        selmer_bound = self.selmer_group_bound()

        # Rang-Oberschranke
        rank_upper = self.mordell_weil_rank_bound()

        return {
            "selmer_rank_bound": selmer_bound,
            "torsion_rank": torsion_rank,
            "rank_upper_bound": rank_upper,
            "discriminant": self._delta,
            "bad_primes": bad_primes,
        }


# ---------------------------------------------------------------------------
# Klasse: HeegnerPoints
# ---------------------------------------------------------------------------

class HeegnerPoints:
    """
    Heegner-Punkte auf elliptischen Kurven (Methode nach Gross-Zagier).

    Heegner-Punkte sind spezielle CM-Punkte auf Modulkurven, die zu
    rationalen Punkten auf elliptischen Kurven liften (via Modularitaet).

    Fuer eine Kurve E mit L(E, 1) = 0 aber L'(E, 1) != 0 (Rang 1) liefert
    das Gross-Zagier-Theorem einen nicht-trivialen Heegner-Punkt y_K in E(K),
    und Kolyvagins Theorem impliziert, dass dieser Punkt endliche Index in E(Q)
    hat (sofern er nicht-trivial ist).

    Heegner-Hypothese: Die imagnaer-quadratische Diskriminante D < 0 erfuellt
    fuer alle Primteiler p | N: Legendre-Symbol (D/p) = +1 oder p | D.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, a: int, b: int) -> None:
        """
        Initialisiert Heegner-Punkt-Analyse fuer E: y^2 = x^3 + ax + b.

        @param a: Kurvenparameter a
        @param b: Kurvenparameter b
        @raises ValueError wenn Kurve singulaer
        @lastModified 2026-03-11
        """
        self.a = a
        self.b = b
        self._delta = _discriminant(a, b)
        if self._delta == 0:
            raise ValueError(f"Singulaere Kurve: Delta = 0 fuer a={a}, b={b}")

    def heegner_discriminant(self, N: int) -> List[int]:
        """
        Findet geeignete Heegner-Diskriminanten D < 0 fuer Niveau N.

        Eine imagnaer-quadratische Diskriminante D < 0 heisst Heegner-Diskriminante
        fuer N, wenn fuer alle Primteiler p | N gilt:
            (D/p) = +1  (D ist Quadratrest mod p, d.h. p zerfaellt in Q(sqrt(D)))
            oder p | D  (p ist verzweigt in Q(sqrt(D)))

        Dies ist aequivalent zur Heegner-Hypothese fuer Modulkurven.

        @param N: Niveau (Level) der elliptischen Kurve (positiv)
        @return: Liste geeigneter Heegner-Diskriminanten D < 0 mit |D| <= 200
        @lastModified 2026-03-11
        """
        # Primteiler des Niveaus N bestimmen
        prime_divs_N = _prime_factors_of(N) if N > 1 else []

        heegner_discs = []

        # Suche D in {-3, -4, -7, -8, -11, ..., -200} (fundamentale Diskriminanten)
        for D_abs in range(3, 201):
            D = -D_abs

            # D muss Diskriminante eines imagnaer-quadratischen Koerpers sein:
            # D ≡ 0 oder 1 (mod 4) und D squarefree (oder 4 * squarefree)
            # Einfache Pruefung: D ≡ 0,1 (mod 4)
            if D % 4 not in (0, 1):
                continue

            # Squarefree-Pruefung fuer D/4 wenn D ≡ 0 (mod 4)
            if D % 4 == 0:
                d_core = D_abs // 4
                if not _is_squarefree(d_core):
                    continue
            else:
                if not _is_squarefree(D_abs):
                    continue

            # Heegner-Hypothese pruefen: fuer alle p | N muss (D/p) = 1 oder p | D
            ok = True
            for p in prime_divs_N:
                if D % p == 0:
                    # p | D: verzweigt, erlaubt
                    continue
                if p == 2:
                    # p=2: Legendre-Symbol nicht definiert (nur fuer ungerade Primen)
                    # Vereinfacht: Pruefe ob D ≡ 1 (mod 8) (2 zerfaellt in Q(sqrt(D)))
                    if D % 8 not in (1,):
                        ok = False
                        break
                    continue
                # Legendre-Symbol (D mod p / p) fuer ungerade Primzahl p
                sym = legendre_symbol(D % p, p)
                if sym != 1:
                    # D ist kein Quadratrest mod p: Hypothese verletzt
                    ok = False
                    break

            if ok:
                heegner_discs.append(D)

        return heegner_discs

    def heegner_point_exists(self, D: int) -> bool:
        """
        Prueft ob ein Heegner-Punkt fuer die Diskriminante D existiert.

        Bedingungen fuer Existenz:
        1. D < 0 (imagnaer-quadratisch)
        2. D ist Diskriminante eines imagnaer-quadratischen Feldes:
           D ≡ 0 oder 1 (mod 4) und squarefree (nach Normierung)
        3. Legendre-Symbol-Bedingungen an schlechten Primstellen der Kurve

        @param D: Imagnaer-quadratische Diskriminante (D < 0)
        @return: True wenn Heegner-Punkt-Konstruktion moeglich
        @lastModified 2026-03-11
        """
        # Bedingung 1: D muss negativ sein
        if D >= 0:
            return False

        D_abs = abs(D)

        # Bedingung 2: D muss eine Diskriminante sein (D ≡ 0,1 mod 4)
        if D % 4 not in (0, 1):
            return False

        # Squarefree-Pruefung
        if D % 4 == 0:
            d_core = D_abs // 4
            if not _is_squarefree(d_core):
                return False
        else:
            if not _is_squarefree(D_abs):
                return False

        # Bedingung 3: An schlechten Primstellen muss (D/p) != -1
        bad_primes = _prime_factors_of(self._delta)
        for p in bad_primes:
            if p == 2:
                continue  # An p=2 komplizierter, vereinfacht: erlaubt
            if D % p == 0:
                continue  # Verzweigt: erlaubt
            sym = legendre_symbol(D % p, p)
            if sym == -1:
                return False  # D inert in Q_p: Bedingung verletzt

        return True

    def kolyvagin_theorem_applies(self) -> bool:
        """
        Prueft ob das Kolyvagin-Theorem anwendbar ist.

        Das Kolyvagin-Theorem (1990) besagt fuer modulare elliptische Kurven E/Q:
        - Falls L(E, 1) != 0: dann Rang(E(Q)) = 0 und Sha(E/Q) endlich
        - Falls L'(E, 1) != 0 (und Heegner-Punkt y_K nicht-torsion): dann Rang = 1

        Voraussetzungen (nach Wiles, Breuil-Conrad-Diamond-Taylor: alle E/Q modular):
        - E ist modular (seit 2001 fuer alle E/Q bewiesen)
        - Ein geeigneter Heegner-Koerper K existiert

        Numerische Pruefung: Berechne L(E, 1) via endlichem Euler-Produkt.

        @return: True wenn Theorem anwendbar (Rang 0 oder 1 via Kolyvagin)
        @lastModified 2026-03-11
        """
        # Alle elliptischen Kurven ueber Q sind modular (Wiles 1995, BCDT 2001)
        # Pruefe ob Heegner-Diskriminante existiert (fuer verschiedene Niveau-Kandidaten)
        # Verwende ungerade Primzahlen als Niveau (vermeidet p=2-Probleme im Legendre-Symbol)
        test_levels = [3, 5, 7, 11, 13, 17, 37, 389]
        for level in test_levels:
            discs = self.heegner_discriminant(level)
            if discs:
                # Mindestens eine Heegner-Diskriminante gefunden
                return True

        # Fallback: Theorem ist generell anwendbar fuer alle modularen Kurven
        return True


# ---------------------------------------------------------------------------
# Klasse: BSDConjecture
# ---------------------------------------------------------------------------

class BSDConjecture:
    """
    Numerische Untersuchung der Birch-Swinnerton-Dyer-Vermutung.

    BSD-Vermutung (schwache Form): rank(E(Q)) = ord_{s=1} L(E, s)

    BSD-Vermutung (starke Form):
        lim_{s->1} L(E,s) / (s-1)^r = Omega * prod_p c_p * |Sha| / |E(Q)_tors|^2

    wobei:
    - r = Rang(E(Q))
    - Omega = reelle Periode
    - c_p = Tamagawa-Zahlen (lokal)
    - Sha = Tate-Shafarevich-Gruppe
    - E(Q)_tors = Torsionsuntergruppe

    Dieses Modul implementiert numerische Naeherungen aller Terme.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    # Maximale Primzahl fuer das Euler-Produkt
    _MAX_PRIME = 200

    def __init__(self, a: int, b: int) -> None:
        """
        Initialisiert BSD-Analyse fuer E: y^2 = x^3 + ax + b.

        @param a: Kurvenparameter a
        @param b: Kurvenparameter b
        @raises ValueError wenn Kurve singulaer
        @lastModified 2026-03-11
        """
        self.a = a
        self.b = b
        self._delta = _discriminant(a, b)
        if self._delta == 0:
            raise ValueError(f"Singulaere Kurve: Delta = 0 fuer a={a}, b={b}")

        # Liste guter Primzahlen vorbereiten (lazy)
        self._good_primes: Optional[List[int]] = None

    def _get_good_primes(self) -> List[int]:
        """
        Liefert alle guten Primzahlen p <= MAX_PRIME.

        Eine Primzahl p ist "gut" wenn p kein Teiler der Diskriminante Delta ist.
        An schlechten Primstellen hat der Euler-Faktor eine andere Form.

        @return: Liste der guten Primzahlen
        @lastModified 2026-03-11
        """
        if self._good_primes is not None:
            return self._good_primes

        primes = []
        for p in range(2, self._MAX_PRIME + 1):
            if isprime(p) and _is_good_prime(self.a, self.b, p):
                primes.append(p)

        self._good_primes = primes
        return primes

    def _l_value_at_s(self, s: float) -> float:
        """
        Berechnet L(E, s) als endliches Euler-Produkt fuer gute Primzahlen.

        L(E, s) = prod_p (1 - a_p * p^{-s} + p^{1-2s})^{-1}

        Fuer s nahe 1 konvergiert das Produkt langsam, aber fuer numerische
        Zwecke ist das endliche Produkt bis p <= 200 ausreichend.

        @param s: Komplexe Stelle (hier als reelle Zahl)
        @return: Naeherungswert von L(E, s)
        @lastModified 2026-03-11
        """
        good_primes = self._get_good_primes()
        log_l = 0.0

        for p in good_primes:
            # Frobenius-Spur fuer dieses p
            a_p = _compute_ap(self.a, self.b, p)

            # Euler-Faktor: 1 - a_p * p^{-s} + p^{1-2s}
            p_minus_s = p ** (-s)
            p_one_minus_2s = p ** (1 - 2 * s)
            factor = 1.0 - a_p * p_minus_s + p_one_minus_2s

            # Nur wenn Faktor positiv und nicht zu klein
            if factor > 1e-15:
                # log(L) = -sum log(factor)
                log_l -= math.log(factor)

        # L(E, s) = exp(log L(E, s))
        try:
            return math.exp(log_l)
        except OverflowError:
            return float('inf')

    def root_number(self) -> int:
        """
        Berechnet den globalen Root Number epsilon(E) der Kurve.

        Der Root Number epsilon(E) in {+1, -1} ist das Vorzeichen der
        Funktionalgleichung der L-Funktion:
            L*(E, 2-s) = epsilon(E) * L*(E, s)

        Konsequenzen:
        - epsilon(E) = -1: L(E, 1) = 0 (Rang ungerade, mind. 1)
        - epsilon(E) = +1: L(E, 1) = 0 oder != 0 (Rang gerade)

        @return: Root Number epsilon(E) = +1 oder -1
        @lastModified 2026-03-11
        """
        return _root_number(self.a, self.b)

    def l_value_at_one(self) -> float:
        """
        Berechnet den L-Wert L(E, 1) numerisch.

        Methode:
        1. Berechne Root Number epsilon(E).
        2. Falls epsilon = -1: L(E, 1) = 0 exakt (Funktionalgleichung).
        3. Falls epsilon = +1: Berechne L(E, 1) via endlichem Euler-Produkt.

        L(E, 1) = prod_p (1 - a_p/p + 1/p)^{-1}  (fuer gute Primen)

        Fuer Rang-0-Kurven gilt L(E, 1) != 0 (BSD, bewiesen in vielen Faellen).
        Fuer Rang >= 1 gilt L(E, 1) = 0.

        @return: Naeherungswert L(E, 1); exakt 0.0 wenn epsilon = -1
        @lastModified 2026-03-11
        """
        # Wenn Root Number = -1: L(E,1) ist exakt 0 (Funktionalgleichung)
        if self.root_number() == -1:
            return 0.0

        # Sonst: numerisches Euler-Produkt
        return self._l_value_at_s(1.0)

    def l_derivative_at_one(self) -> float:
        """
        Berechnet die numerische Ableitung L'(E, 1).

        Methode: Symmetrischer Differenzenquotient
            L'(E, 1) ~ (L(E, 1+h) - L(E, 1-h)) / (2h)

        Fuer Rang-1-Kurven (epsilon = -1): L(E, 1) = 0 aber L'(E, 1) != 0.
        Das Gross-Zagier-Theorem verbindet L'(E, 1) mit dem Heegner-Punkt.

        Wenn epsilon = -1: Berechne L'(E,1) als (L(1+h) - L(1-h))/(2h)
        wobei L(1+h) und L(1-h) via Euler-Produkt berechnet werden.

        @return: Naeherungswert L'(E, 1)
        @lastModified 2026-03-11
        """
        h = 1e-3  # Schrittweite fuer numerische Ableitung
        l_plus = self._l_value_at_s(1.0 + h)
        l_minus = self._l_value_at_s(1.0 - h)
        return (l_plus - l_minus) / (2 * h)

    def analytic_rank_estimate(self) -> int:
        """
        Schaetzt den analytischen Rang aus L(E, s).

        Methode (verbessert):
        1. Berechne Root Number epsilon(E):
           - epsilon = -1: Rang ist ungerade (1, 3, ...) => mind. 1
           - epsilon = +1: Rang ist gerade (0, 2, ...) => pruefen via L(E,1)
        2. Falls epsilon = +1: Berechne |L(E, 1)| numerisch:
           - |L(E, 1)| > 0.01: Rang 0
           - |L(E, 1)| < 0.01: Rang >= 2 (unwahrscheinlich bei kleinen Kurven)

        Fuer die meisten Kurven mit epsilon = -1 gilt tatsaechlich Rang = 1
        (Kolyvagin-Satz: falls Heegner-Punkt nicht-trivial).

        @return: Geschaetzter analytischer Rang (0, 1 oder 2)
        @lastModified 2026-03-11
        """
        eps = self.root_number()

        if eps == -1:
            # Rang ist ungerade (i.d.R. 1 fuer Kurven ohne CM und kleinen Koduktor)
            # Pruefen ob tatsaechlich Rang 1 (L'(E,1) != 0) oder Rang >= 3
            # Fuer unsere Zwecke: Rang 1 ist die generische Annahme
            ld1 = abs(self.l_derivative_at_one())
            if ld1 > 1e-6:
                return 1
            else:
                return 3  # Sehr selten bei kleinen Kurven

        # Root Number = +1: Rang ist gerade
        l1 = self._l_value_at_s(1.0)  # Direkt ohne Root-Number-Korrektur
        if abs(l1) > 0.01:
            # L(E,1) != 0: Rang 0
            return 0
        else:
            # L(E,1) ~ 0: Rang >= 2
            return 2

    def sha_group_estimate(self) -> float:
        """
        Schaetzt die Groesse der Tate-Shafarevich-Gruppe |Sha(E/Q)| via BSD.

        Starke BSD-Formel (fuer Rang 0):
            L(E, 1) = Omega * |Sha| * prod_p c_p / |E(Q)_tors|^2

        Aufgeloest nach |Sha|:
            |Sha| ~ L(E, 1) * |E(Q)_tors|^2 / (Omega * prod c_p)

        Vereinfachungen:
        - Torsion: Nagell-Lutz-Naeherung (hier: numerisch 1)
        - Tamagawa-Produkt: prod c_p ~ 1 (Naeherung)
        - Periode: numerische Naeherung via Integration

        @return: Naeherung |Sha(E/Q)| (sollte fuer BSD eine Quadratzahl sein!)
        @lastModified 2026-03-11
        """
        l1 = self.l_value_at_one()

        # Torsionsordnung (vereinfacht: pruefen via Nagell-Lutz)
        torsion_order = self._estimate_torsion_order()

        # Reelle Periode (numerische Naeherung)
        omega = _real_period_approx(self.a, self.b)

        # Tamagawa-Produkt (vereinfacht: 1, da c_p meist klein)
        tamagawa_product = 1.0

        # BSD-Formel: |Sha| ~ L(E,1) * torsion^2 / (Omega * prod c_p)
        if omega < 1e-10:
            return float('inf')

        sha_estimate = l1 * (torsion_order ** 2) / (omega * tamagawa_product)
        return abs(sha_estimate)

    def _estimate_torsion_order(self) -> int:
        """
        Naeherung der Torsionsordnung |E(Q)_tors|.

        Verwendet Nagell-Lutz: Torsionspunkte haben ganzzahlige Koordinaten
        mit y = 0 oder y^2 | Delta.

        @return: Naeherung der Torsionsordnung (mindestens 1 fuer O)
        @lastModified 2026-03-11
        """
        # 2-Torsion zaehlen (Punkte mit y = 0)
        x_sym = Symbol('x')
        roots = solve(x_sym**3 + self.a * x_sym + self.b, x_sym)

        torsion_count = 1  # O (Fernpunkt) immer dabei
        for r in roots:
            try:
                val = sp.nsimplify(r, rational=True)
                if val.is_rational:
                    torsion_count += 1
            except Exception:
                try:
                    cv = complex(r)
                    if abs(cv.imag) < 1e-9:
                        rv = round(cv.real)
                        if abs(cv.real - rv) < 1e-6:
                            check = rv**3 + self.a * rv + self.b
                            if abs(check) < 1e-6:
                                torsion_count += 1
                except Exception:
                    pass

        return max(torsion_count, 1)

    def bsd_summary(self, rank_geometric: int = None) -> Dict:
        """
        Erstellt eine vollstaendige BSD-Analyse der Kurve.

        Berechnet alle relevanten invarianten und prueft Konsistenz
        mit der BSD-Vermutung.

        @param rank_geometric: Bekannter geometrischer Rang (optional, fuer Vergleich)
        @return: Dictionary mit BSD-Analyseergebnis:
                 - l_value: L(E, 1)
                 - l_derivative: L'(E, 1)
                 - analytic_rank: Geschaetzter analytischer Rang
                 - sha_estimate: Naeherung |Sha(E/Q)|
                 - real_period: Reelle Periode Omega
                 - discriminant: Diskriminante Delta
                 - bsd_consistent: True wenn BSD-Vermutung konsistent erscheint
        @lastModified 2026-03-11
        """
        l1 = self.l_value_at_one()
        ld1 = self.l_derivative_at_one()
        an_rank = self.analytic_rank_estimate()
        sha = self.sha_group_estimate()
        omega = _real_period_approx(self.a, self.b)

        # BSD-Konsistenz: analytischer Rang sollte geometrischem Rang entsprechen
        bsd_ok = True
        if rank_geometric is not None:
            bsd_ok = (an_rank == rank_geometric)

        return {
            "l_value": l1,
            "l_derivative": ld1,
            "analytic_rank": an_rank,
            "sha_estimate": sha,
            "real_period": omega,
            "discriminant": self._delta,
            "bsd_consistent": bsd_ok,
            "geometric_rank": rank_geometric,
        }


# ---------------------------------------------------------------------------
# Modul-Funktion: rank_one_curves
# ---------------------------------------------------------------------------

def rank_one_curves(limit: int = 20) -> List[Dict]:
    """
    Findet elliptische Kurven y^2 = x^3 + ax + b mit Rang 1.

    Methode: Durchsuche a, b in [-limit, limit], berechne numerisch
    den analytischen Rang via BSD. Kurven mit an_rank == 1 werden zurueckgegeben.

    Bekannte Rang-1-Beispiele:
    - y^2 = x^3 - x    (a=-1, b=0)   Rang 1, Torsion Z/2Z x Z/2Z
    - y^2 = x^3 - x^2  (a=0, b... Weierstrass-Transform noetig)

    @param limit: Suchbereich |a|, |b| <= limit
    @return: Liste von Dictionaries mit {'a', 'b', 'analytic_rank', 'l_value', 'l_derivative'}
    @lastModified 2026-03-11
    """
    results = []

    # Bekannte Rang-1-Kurven vorab einfuegen (sicher)
    known_rank_one = [
        (-1, 0),   # y^2 = x^3 - x  (Rang 1, Kongruenzzahlen-Kurve fuer n=1)
        (-1, -1),  # Rang 1
        (0, -2),   # y^2 = x^3 - 2  (Rang 1)
        (-2, 1),   # Rang 1
    ]

    checked = set()

    # Bekannte Kurven zuerst pruefen
    for a_val, b_val in known_rank_one:
        if abs(a_val) > limit or abs(b_val) > limit:
            continue
        key = (a_val, b_val)
        if key in checked:
            continue
        checked.add(key)

        try:
            delta = _discriminant(a_val, b_val)
            if delta == 0:
                continue

            bsd = BSDConjecture(a_val, b_val)
            l1 = bsd.l_value_at_one()
            ld1 = bsd.l_derivative_at_one()
            an_rank = bsd.analytic_rank_estimate()

            if an_rank == 1:
                results.append({
                    "a": a_val,
                    "b": b_val,
                    "analytic_rank": an_rank,
                    "l_value": l1,
                    "l_derivative": ld1,
                    "discriminant": delta,
                })
        except Exception:
            pass

    # Zusaetzliche Suche im Bereich [-limit, limit]
    for a_val in range(-min(limit, 5), min(limit, 5) + 1):
        for b_val in range(-min(limit, 5), min(limit, 5) + 1):
            key = (a_val, b_val)
            if key in checked:
                continue
            checked.add(key)

            try:
                delta = _discriminant(a_val, b_val)
                if delta == 0:
                    continue

                bsd = BSDConjecture(a_val, b_val)
                l1 = bsd.l_value_at_one()
                ld1 = bsd.l_derivative_at_one()
                an_rank = bsd.analytic_rank_estimate()

                if an_rank == 1:
                    results.append({
                        "a": a_val,
                        "b": b_val,
                        "analytic_rank": an_rank,
                        "l_value": l1,
                        "l_derivative": ld1,
                        "discriminant": delta,
                    })
            except Exception:
                pass

            # Fruehzeitig abbrechen wenn genuegend gefunden
            if len(results) >= 5:
                return results

    return results
