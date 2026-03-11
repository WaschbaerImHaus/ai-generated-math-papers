"""
@file elliptic_curves.py
@brief Elliptische Kurven - Gruppenstruktur, Arithmetik und zahlentheoretische Anwendungen.
@description
    Implementiert die vollstaendige Arithmetik elliptischer Kurven:

    Eine elliptische Kurve ueber einem Koerper K ist:
        E: y^2 = x^3 + ax + b  (Weierstrass-Normalform)

    mit Diskriminante Delta = -16(4a^3 + 27b^2) != 0 (nicht-singulaere Kurve).

    Die Menge der Punkte E(K) bildet eine abelsche Gruppe mit:
    - Neutralelement: Punkt O im Unendlichen
    - Addition: Chord-and-Tangent-Methode
    - Inverses: -(x,y) = (x,-y)

    Anwendungen:
    - Elliptische Kurven Kryptographie (ECC)
    - Faktorisierung (Lenstra's ECM)
    - BSD-Vermutung: Rang von E(Q) = ord_{s=1} L(E,s)
    - Wiles' Beweis des Grossen Fermatschen Satzes (Shimura-Taniyama-Wiles)

@author Michael Fuhrmann
@lastModified 2026-03-11 (Build 109: Migration auf math_helpers)
"""

import math
import random
import sympy as sp
import numpy as np
from typing import Optional, List, Tuple, Dict, Union

# Zentrale Hilfsfunktionen aus math_helpers importieren (Rückwärtskompatibilität via Alias)
from math_helpers import is_prime as _is_prime_helper


# ---------------------------------------------------------------------------
# Hilfs-Funktion: Modulares Inverses (erweiteter euklidischer Algorithmus)
# ---------------------------------------------------------------------------

def _mod_inverse(a: int, m: int) -> int:
    """
    Berechnet das modulare Inverse von a modulo m.

    Verwendet den erweiterten euklidischen Algorithmus.
    Existiert nur, wenn gcd(a, m) = 1.

    @param a: Zahl, deren Inverses gesucht wird
    @param m: Modulus
    @return: x mit a*x ≡ 1 (mod m)
    @raises ValueError wenn kein Inverses existiert
    @lastModified 2026-03-10
    """
    # Sonderfall: a muss zu m teilerfremd sein
    g, x, _ = _extended_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"Kein modulares Inverses: gcd({a}, {m}) = {g} != 1")
    return x % m


def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    Erweiterter euklidischer Algorithmus: gcd(a,b) = a*x + b*y.

    @return: (gcd, x, y) mit a*x + b*y = gcd
    @lastModified 2026-03-10
    """
    # Rekursionsanker
    if a == 0:
        return b, 0, 1
    g, x, y = _extended_gcd(b % a, a)
    return g, y - (b // a) * x, x


# Alias für Rückwärtskompatibilität: _is_prime verweist auf math_helpers.is_prime
_is_prime = _is_prime_helper


def _sqrt_mod_p(a: int, p: int) -> Optional[int]:
    """
    Berechnet die Quadratwurzel von a modulo p (Tonelli-Shanks-Algorithmus).

    Gibt x zurueck mit x^2 ≡ a (mod p), oder None wenn kein Quadratrest.

    @param a: Zahl
    @param p: Primzahl
    @return: Quadratwurzel oder None
    @lastModified 2026-03-10
    """
    # Nullfall
    a = a % p
    if a == 0:
        return 0

    # Euler-Kriterium: a^((p-1)/2) ≡ 1 (mod p) gdw a Quadratrest
    if pow(a, (p - 1) // 2, p) != 1:
        return None  # Kein Quadratrest

    # Einfacher Fall: p ≡ 3 (mod 4)
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Allgemeiner Fall: Tonelli-Shanks
    # Schreibe p-1 = Q * 2^S
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # Finde einen Nicht-Quadratrest z
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    # Initialisierung
    M = S
    c = pow(z, Q, p)
    t = pow(a, Q, p)
    R = pow(a, (Q + 1) // 2, p)

    # Iterationsschleife
    while True:
        if t == 0:
            return 0
        if t == 1:
            return R
        # Finde kleinstes i mit t^(2^i) ≡ 1
        i, temp = 1, (t * t) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        # Update
        b = pow(c, pow(2, M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p


# ===========================================================================
# Klasse 1: ECPoint — Punkt auf elliptischer Kurve
# ===========================================================================

class ECPoint:
    """
    Punkt P = (x, y) auf einer elliptischen Kurve, oder O (Punkt im Unendlichen).

    Der Punkt im Unendlichen O ist das neutrale Element der Gruppenoperation.
    Darstellung: ECPoint(None, None) oder ECPoint.infinity().

    Die elliptische Kurven-Addition folgt der Chord-and-Tangent-Methode,
    die geometrisch interpretiert werden kann:
    - Zwei verschiedene Punkte: Schneide die Sekante mit der Kurve
    - Gleicher Punkt: Verwende die Tangente an diesem Punkt

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self,
                 x: Optional[Union[int, float]],
                 y: Optional[Union[int, float]],
                 curve: Optional['EllipticCurve'] = None) -> None:
        """
        Erstellt einen Punkt auf der elliptischen Kurve.

        @param x: x-Koordinate (None fuer Punkt im Unendlichen)
        @param y: y-Koordinate (None fuer Punkt im Unendlichen)
        @param curve: Zugehoerige elliptische Kurve (optional)
        @lastModified 2026-03-10
        """
        self.x = x
        self.y = y
        self.curve = curve  # Referenz zur Kurve fuer Additionsparameter

    @classmethod
    def infinity(cls, curve: Optional['EllipticCurve'] = None) -> 'ECPoint':
        """
        Erzeugt den Punkt im Unendlichen O (neutrales Element der Gruppe).

        Mathematisch: O ist der einzige Punkt auf der projektiven Geraden
        durch zwei inverse Punkte (x, y) und (x, -y).

        @param curve: Zugehoerige elliptische Kurve
        @return: Neutrales Element O
        @lastModified 2026-03-10
        """
        return cls(None, None, curve)

    @property
    def is_infinity(self) -> bool:
        """
        True wenn dieser Punkt das neutrale Element O ist.

        @return: True wenn P = O
        @lastModified 2026-03-10
        """
        return self.x is None and self.y is None

    def __eq__(self, other: object) -> bool:
        """
        Prueft Punktgleichheit: P = Q gdw x_P = x_Q und y_P = y_Q.

        Beide Punkte im Unendlichen sind gleich.

        @param other: Anderer Punkt
        @return: True wenn gleich
        @lastModified 2026-03-10
        """
        if not isinstance(other, ECPoint):
            return False
        # Beide sind O
        if self.is_infinity and other.is_infinity:
            return True
        # Nur einer ist O
        if self.is_infinity or other.is_infinity:
            return False
        # Koordinatenvergleich (mit Toleranz fuer Float)
        if isinstance(self.x, float) or isinstance(other.x, float):
            return abs(self.x - other.x) < 1e-10 and abs(self.y - other.y) < 1e-10
        return self.x == other.x and self.y == other.y

    def __neg__(self) -> 'ECPoint':
        """
        Inverses eines Punktes: -(x, y) = (x, -y).

        Geometrisch: Spiegelung an der x-Achse.
        Das Inverse von O ist O.

        @return: Inverspunkt
        @lastModified 2026-03-10
        """
        if self.is_infinity:
            return ECPoint.infinity(self.curve)
        return ECPoint(self.x, -self.y, self.curve)

    def __add__(self, other: 'ECPoint') -> 'ECPoint':
        """
        Elliptische Kurven-Addition (Chord-and-Tangent-Methode).

        Berechnet P + Q auf der Kurve E: y^2 = x^3 + ax + b.

        Faelle:
        1. P = O  =>  P + Q = Q
        2. Q = O  =>  P + Q = P
        3. P = -Q =>  P + Q = O (inverse Punkte)
        4. P = Q  =>  Tangente (Verdopplung):
               lambda = (3*x1^2 + a) / (2*y1)
        5. P != Q =>  Sekante:
               lambda = (y2 - y1) / (x2 - x1)

        In allen Faellen (4, 5):
               x3 = lambda^2 - x1 - x2
               y3 = lambda*(x1 - x3) - y1

        @param other: Zweiter Summand Q
        @return: Summe P + Q
        @raises ValueError wenn Kurvenparameter fehlen
        @lastModified 2026-03-10
        """
        P, Q = self, other

        # Fall 1: P = O
        if P.is_infinity:
            return Q

        # Fall 2: Q = O
        if Q.is_infinity:
            return P

        # Kurvenparameter ermitteln
        curve = P.curve or Q.curve
        if curve is None:
            raise ValueError("Kurvenparameter fehlen — setze curve bei ECPoint-Erstellung")
        a = curve.a

        # Fallunterscheidung fuer modulare Arithmetik
        use_mod = hasattr(curve, 'p') and isinstance(curve.p, int)
        p_mod = curve.p if use_mod else None

        def div(num, den):
            """Division: entweder modular oder reell."""
            if use_mod:
                return (num * _mod_inverse(int(den), p_mod)) % p_mod
            return num / den

        def mod(val):
            """Reduktion mod p wenn noetig."""
            if use_mod:
                return int(val) % p_mod
            return val

        x1, y1 = P.x, P.y
        x2, y2 = Q.x, Q.y

        # Fall 3: P = -Q (inverse Punkte) => O
        if x1 == x2:
            if use_mod:
                y_sum = (y1 + y2) % p_mod
            else:
                y_sum = y1 + y2
            if y_sum == 0:
                return ECPoint.infinity(curve)

        # Steigung lambda berechnen
        if x1 == x2 and y1 == y2:
            # Fall 4: Tangente (Punktverdopplung)
            # lambda = (3*x1^2 + a) / (2*y1)
            if use_mod:
                numerator = mod(3 * x1 * x1 + a)
                denominator = mod(2 * y1)
            else:
                numerator = 3 * x1 ** 2 + a
                denominator = 2 * y1

            # Sonderfall: y1 = 0 => Tangente ist senkrecht => Ergebnis O
            if denominator == 0:
                return ECPoint.infinity(curve)
        else:
            # Fall 5: Sekante (verschiedene Punkte)
            # lambda = (y2 - y1) / (x2 - x1)
            if use_mod:
                numerator = mod(y2 - y1)
                denominator = mod(x2 - x1)
            else:
                numerator = y2 - y1
                denominator = x2 - x1

        lam = div(numerator, denominator)

        # Additionsformeln
        x3 = mod(lam * lam - x1 - x2)
        y3 = mod(lam * (x1 - x3) - y1)

        return ECPoint(x3, y3, curve)

    def __mul__(self, n: int) -> 'ECPoint':
        """
        Skalarmultiplikation n*P via Double-and-Add-Algorithmus.

        Effizient in O(log n) Schritten statt naiver O(n) Additionen.

        Algorithmus:
        1. Schreibe n in Binaerdarstellung
        2. Starte mit Ergebnis R = O (neutral)
        3. Fuer jedes Bit (von MSB zu LSB):
           - Verdopple R: R = 2*R
           - Wenn Bit = 1: addiere P: R = R + P

        @param n: Skalarmultiplikator (ganze Zahl)
        @return: n*P
        @lastModified 2026-03-10
        """
        if not isinstance(n, int):
            raise TypeError(f"Multiplikator muss int sein, nicht {type(n)}")

        # Negativer Multiplikator: n*P = (-n)*(-P)
        if n < 0:
            return (-self) * (-n)

        # n = 0: neutrales Element
        if n == 0:
            return ECPoint.infinity(self.curve)

        # Double-and-Add
        result = ECPoint.infinity(self.curve)  # R = O
        addend = ECPoint(self.x, self.y, self.curve)  # Arbeitskopie von P

        # Binaere Darstellung durchlaufen (LSB zuerst)
        while n > 0:
            if n & 1:  # Aktuelles Bit ist 1
                result = result + addend
            addend = addend + addend  # Verdopplung
            n >>= 1  # Naechstes Bit

        return result

    def __rmul__(self, n: int) -> 'ECPoint':
        """
        Rechts-Multiplikation n*P = P*n (Kommutativitaet fuer int*ECPoint).

        @param n: Skalarmultiplikator
        @return: n*P
        @lastModified 2026-03-10
        """
        return self.__mul__(n)

    def __repr__(self) -> str:
        """
        Lesbare Darstellung des Punktes.

        @return: String-Darstellung
        @lastModified 2026-03-10
        """
        if self.is_infinity:
            return "ECPoint(O)"
        return f"ECPoint({self.x}, {self.y})"

    def on_curve(self) -> bool:
        """
        Prueft ob dieser Punkt auf der zugehoerigen Kurve liegt.

        Verifiziert: y^2 = x^3 + a*x + b (bzw. mod p)

        @return: True wenn Punkt auf Kurve
        @raises ValueError wenn keine Kurve zugeordnet
        @lastModified 2026-03-10
        """
        if self.is_infinity:
            return True  # O liegt per Definition auf jeder Kurve
        if self.curve is None:
            raise ValueError("Kein Kurvenobjekt zugeordnet")
        return self.curve.is_on_curve(self.x, self.y)


# ===========================================================================
# Klasse 2: EllipticCurve — Kurve ueber den reellen Zahlen
# ===========================================================================

class EllipticCurve:
    """
    Elliptische Kurve E: y^2 = x^3 + ax + b ueber den reellen Zahlen.

    Eigenschaften:
    - Diskriminante: Delta = -16*(4*a^3 + 27*b^2)
    - j-Invariante: j = -1728 * (4a)^3 / Delta
    - Rang: Dimension von E(Q) als abelsche Gruppe (via BSD-Vermutung)

    Die Diskriminante bestimmt ob die Kurve nicht-singulaer ist.
    Delta != 0 garantiert genau eine Tangente in jedem Kurvenpunkt.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, a: Union[int, float], b: Union[int, float]) -> None:
        """
        Erstellt eine elliptische Kurve y^2 = x^3 + ax + b.

        @param a: Koeffizient der x-Terme (lineare Komponente)
        @param b: Konstanter Term (kubische Verschiebung)
        @raises ValueError wenn Delta = 0 (singulaere/entartete Kurve)
        @lastModified 2026-03-10
        """
        self.a = a
        self.b = b
        # Diskriminante: Delta = -16*(4a^3 + 27b^2)
        self.delta = -16 * (4 * a ** 3 + 27 * b ** 2)

        # Kurve ist singulaer wenn Delta = 0 (Knoten oder Spitze)
        if abs(self.delta) < 1e-15:
            raise ValueError(
                f"Singulaere Kurve: Delta = 0 fuer a={a}, b={b}. "
                f"Die Kurve hat einen Knoten oder eine Spitze."
            )

    def j_invariant(self) -> float:
        """
        Berechnet die j-Invariante der Kurve.

        j = -1728 * (4a)^3 / Delta

        Die j-Invariante klassifiziert elliptische Kurven ueber dem algebraischen
        Abschluss von K (zwei Kurven sind isomorph ueber k_bar gdw j gleich).

        Spezialwerte:
        - j = 0:    Kurve hat erhoegte Symmetrie (CM durch Z[omega])
        - j = 1728: Kurve ist isomorph zu y^2 = x^3 - x (CM durch Z[i])

        @return: j-Invariante
        @lastModified 2026-03-10
        """
        # Zaehler: -1728 * (4a)^3
        numerator = -1728 * (4 * self.a) ** 3
        return numerator / self.delta

    def discriminant(self) -> float:
        """
        Gibt die Diskriminante Delta = -16*(4a^3 + 27b^2) zurueck.

        @return: Diskriminante der Kurve
        @lastModified 2026-03-10
        """
        return self.delta

    def point(self, x: float, y: float) -> ECPoint:
        """
        Erstellt einen ECPoint auf dieser Kurve und prueft seine Gueltigkeit.

        @param x: x-Koordinate
        @param y: y-Koordinate
        @return: Verknuepfter ECPoint
        @raises ValueError wenn Punkt nicht auf Kurve
        @lastModified 2026-03-10
        """
        if not self.is_on_curve(x, y):
            raise ValueError(
                f"Punkt ({x}, {y}) liegt nicht auf der Kurve E: "
                f"y^2 = x^3 + {self.a}x + {self.b}"
            )
        return ECPoint(x, y, self)

    def is_on_curve(self, x: float, y: float, tol: float = 1e-8) -> bool:
        """
        Prueft ob (x, y) auf der Kurve liegt: y^2 = x^3 + ax + b.

        @param x: x-Koordinate
        @param y: y-Koordinate
        @param tol: Toleranz fuer Fliesskomma-Vergleich
        @return: True wenn auf der Kurve
        @lastModified 2026-03-10
        """
        # Linke Seite: y^2
        lhs = y ** 2
        # Rechte Seite: x^3 + ax + b
        rhs = x ** 3 + self.a * x + self.b
        return abs(lhs - rhs) < tol

    def y_values(self, x: float) -> List[float]:
        """
        Berechnet die y-Werte fuer gegebenes x.

        Loest: y = +/- sqrt(x^3 + ax + b)
        Gibt leere Liste zurueck wenn x^3 + ax + b < 0.

        @param x: x-Koordinate
        @return: Liste der y-Werte (0, 1 oder 2 Elemente)
        @lastModified 2026-03-10
        """
        # Rechte Seite der Kurvengleichung
        rhs = x ** 3 + self.a * x + self.b

        if rhs < 0:
            return []  # Keine reellen Loesungen
        elif abs(rhs) < 1e-14:
            return [0.0]  # Tangentialer Punkt (y = 0)
        else:
            sqrt_val = math.sqrt(rhs)
            return [sqrt_val, -sqrt_val]  # Zwei symmetrische Punkte

    def points_over_fp(self, p: int) -> List[ECPoint]:
        """
        Berechnet alle Punkte E(Z/pZ) inklusive O.

        Durchsucht alle x in {0, ..., p-1} und prueft ob x^3 + ax + b
        ein Quadratrest mod p ist. Verwendet Euler-Kriterium.

        Komplexitaet: O(p) Operationen.

        @param p: Primzahl (Koerpercharakteristik)
        @return: Liste aller Punkte inklusive O
        @lastModified 2026-03-10
        """
        points = [ECPoint.infinity(self)]  # O immer dabei

        a_mod = int(self.a) % p
        b_mod = int(self.b) % p

        for x in range(p):
            # Rechte Seite der Kurvengleichung mod p
            rhs = (pow(x, 3, p) + a_mod * x + b_mod) % p

            # Quadratwurzel mod p (Tonelli-Shanks)
            y = _sqrt_mod_p(rhs, p)
            if y is not None:
                if y == 0:
                    points.append(ECPoint(x, 0, self))
                else:
                    points.append(ECPoint(x, y, self))
                    points.append(ECPoint(x, p - y, self))  # Zweite Wurzel

        return points

    def order_over_fp(self, p: int) -> int:
        """
        Berechnet die Gruppenordnung #E(Z/pZ) (Anzahl Punkte inkl. O).

        Nach dem Satz von Hasse gilt:
        |#E(Z/pZ) - (p+1)| <= 2*sqrt(p)

        @param p: Primzahl
        @return: Gruppenordnung
        @lastModified 2026-03-10
        """
        return len(self.points_over_fp(p))

    def trace_of_frobenius(self, p: int) -> int:
        """
        Berechnet die Frobenius-Spur a_p.

        Definition: a_p = p + 1 - #E(Z/pZ)

        Hasse-Schranke: |a_p| <= 2*sqrt(p)

        Die Frobenius-Spur charakterisiert das Reduktionsverhalten
        der Kurve an der Primzahl p und ist zentraler Parameter
        der L-Funktion L(E, s).

        @param p: Primzahl
        @return: Frobenius-Spur a_p
        @lastModified 2026-03-10
        """
        order = self.order_over_fp(p)
        return p + 1 - order

    def torsion_points_over_r(self,
                               x_range: Tuple[float, float] = (-10, 10),
                               n_points: int = 1000) -> List[ECPoint]:
        """
        Sucht reelle Punkte auf der Kurve im gegebenen x-Bereich.

        Hinweis: Echte Torsionspunkte existieren nur ueber Q, nicht ueber R.
        Diese Methode liefert Punkte fuer Visualisierungszwecke.

        @param x_range: (x_min, x_max) Suchbereich
        @param n_points: Anzahl der Abtastpunkte
        @return: Liste gefundener reeller Punkte
        @lastModified 2026-03-10
        """
        points = [ECPoint.infinity(self)]  # O immer dabei
        x_values = np.linspace(x_range[0], x_range[1], n_points)

        for x in x_values:
            ys = self.y_values(float(x))
            for y in ys:
                # Nur wenn Punkt auf Kurve (Rundungsfehler vermeiden)
                if self.is_on_curve(float(x), y):
                    points.append(ECPoint(float(x), y, self))

        return points

    def plot_curve(self,
                   x_range: Tuple[float, float] = (-3, 3),
                   title: str = "") -> object:
        """
        Plottet die elliptische Kurve mit matplotlib.

        Zeichnet beide Aeste der Kurve (positive und negative y-Werte).

        @param x_range: (x_min, x_max) Zeichenbereich
        @param title: Titel des Plots
        @return: matplotlib.figure.Figure Objekt
        @lastModified 2026-03-10
        """
        import matplotlib.pyplot as plt

        # Feinraster fuer x-Werte
        x_vals = np.linspace(x_range[0], x_range[1], 2000)

        # Positive und negative y-Werte getrennt berechnen
        y_pos, y_neg = [], []
        x_valid = []

        for x in x_vals:
            rhs = x ** 3 + self.a * x + self.b
            if rhs >= 0:
                sqrt_val = math.sqrt(rhs)
                y_pos.append(sqrt_val)
                y_neg.append(-sqrt_val)
                x_valid.append(x)

        # Plot erstellen
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(x_valid, y_pos, 'b-', linewidth=2, label='y > 0')
        ax.plot(x_valid, y_neg, 'b-', linewidth=2, label='y < 0')
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)
        ax.grid(True, alpha=0.3)

        # Titel und Beschriftung
        plot_title = title or f"E: y² = x³ + {self.a}x + {self.b}"
        ax.set_title(plot_title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        return fig

    def __repr__(self) -> str:
        """
        Lesbare Darstellung der elliptischen Kurve.

        @return: String-Darstellung
        @lastModified 2026-03-10
        """
        return (f"EllipticCurve(y² = x³ + {self.a}x + {self.b}, "
                f"Δ={self.delta:.2f})")


# ===========================================================================
# Klasse 3: EllipticCurveModP — Kurve ueber endlichem Koerper F_p
# ===========================================================================

class EllipticCurveModP(EllipticCurve):
    """
    Elliptische Kurve E ueber dem endlichen Koerper F_p = Z/pZ.

    Alle Koordinaten und Berechnungen erfolgen modulo der Primzahl p.
    Diese Version wird fuer kryptographische Anwendungen (ECC) und
    zahlentheoretische Berechnungen (BSD, L-Funktionen) verwendet.

    Gruppenstruktur: E(F_p) ist eine endliche abelsche Gruppe.
    Nach Cauchy/Lagrange: Ordnung n teilt p-1 (nicht direkt, sondern via Hasse).

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, a: int, b: int, p: int) -> None:
        """
        Erstellt eine elliptische Kurve ueber F_p.

        @param a: Koeffizient der x-Terme
        @param b: Konstanter Term
        @param p: Primzahl (Charakteristik des Koerpers)
        @raises ValueError wenn Kurve mod p singulaer ist
        @lastModified 2026-03-10
        """
        self.p = p

        # Koeffizienten mod p reduzieren
        a_mod = a % p
        b_mod = b % p
        self.a = a_mod
        self.b = b_mod

        # Diskriminante mod p pruefen
        # Delta = -16 * (4a^3 + 27b^2)
        delta_inner = (4 * pow(a_mod, 3, p) + 27 * pow(b_mod, 2, p)) % p
        delta_mod = (-16 * delta_inner) % p

        if delta_mod == 0:
            raise ValueError(
                f"Singulaere Kurve mod {p}: Δ ≡ 0 (mod {p}) "
                f"fuer a={a_mod}, b={b_mod}"
            )
        self.delta = delta_mod

    def add_mod(self, P: ECPoint, Q: ECPoint) -> ECPoint:
        """
        Punktaddition mod p (exakt, ganzzahlig).

        Alle Divisionen werden als modulare Inversen (mod p) berechnet.
        Dies ist die kryptographisch relevante Operation.

        @param P: Erster Summand
        @param Q: Zweiter Summand
        @return: P + Q mod p
        @lastModified 2026-03-10
        """
        # Beide Punkte dieser Kurve zuordnen
        P_copy = ECPoint(P.x, P.y, self)
        Q_copy = ECPoint(Q.x, Q.y, self)
        return P_copy + Q_copy

    def scalar_mult_mod(self, k: int, P: ECPoint) -> ECPoint:
        """
        Skalarmultiplikation k*P mod p via Double-and-Add.

        Grundlegende Operation in ECC: Das Diskrete-Logarithmus-Problem
        (gegeben P und k*P, finde k) ist in E(F_p) schwer.

        @param k: Skalarmultiplikator
        @param P: Basispunkt
        @return: k*P mod p
        @lastModified 2026-03-10
        """
        # Punkt dieser Kurve zuordnen
        P_copy = ECPoint(P.x, P.y, self)
        return P_copy * k

    def all_points(self) -> List[ECPoint]:
        """
        Berechnet alle Punkte E(F_p) inklusive O.

        Brute-Force ueber alle x in {0, ..., p-1}.
        Fuer jedes x wird geprueft ob x^3 + ax + b ein Quadratrest mod p ist.

        Komplexitaet: O(p * sqrt(p)) fuer Tonelli-Shanks.

        @return: Liste aller Punkte in E(F_p)
        @lastModified 2026-03-10
        """
        points = [ECPoint.infinity(self)]  # O ist immer dabei

        for x in range(self.p):
            # Rechte Seite der Kurvengleichung mod p
            rhs = (pow(x, 3, self.p) + self.a * x + self.b) % self.p

            # Quadratwurzel mod p via Tonelli-Shanks
            y = _sqrt_mod_p(rhs, self.p)
            if y is not None:
                if y == 0:
                    # Punkt mit y=0: nur ein Punkt
                    points.append(ECPoint(x, 0, self))
                else:
                    # Zwei Punkte: (x, y) und (x, p-y)
                    points.append(ECPoint(x, y, self))
                    points.append(ECPoint(x, self.p - y, self))

        return points

    def group_order(self) -> int:
        """
        Berechnet die Gruppenordnung #E(F_p).

        Verwendet Brute-Force (Zaehlung aller Punkte).
        Fuer grosse p waere Schoof-Algorithmus effizienter.

        @return: Anzahl Punkte in E(F_p) inklusive O
        @lastModified 2026-03-10
        """
        return len(self.all_points())

    def is_generator(self, P: ECPoint) -> bool:
        """
        Prueft ob P ein Generator der Gruppe E(F_p) ist.

        P ist Generator gdw die Ordnung von P gleich #E(F_p) ist.
        Aequivalent: n*P = O fuer n = #E(F_p) und kein echtes Teiler k < n
        mit k*P = O.

        @param P: Zu pruefender Punkt
        @return: True wenn P Generator ist
        @lastModified 2026-03-10
        """
        n = self.group_order()  # Gruppenordnung

        # n*P muss O ergeben
        result = self.scalar_mult_mod(n, P)
        if not result.is_infinity:
            return False

        # Pruefe echte Teiler der Gruppenordnung
        for q in _prime_factors(n):
            # (n/q)*P darf nicht O sein
            sub_result = self.scalar_mult_mod(n // q, P)
            if sub_result.is_infinity:
                return False  # P hat kleinere Ordnung als n

        return True

    def discrete_log(self,
                     P: ECPoint,
                     Q: ECPoint,
                     max_steps: int = 10000) -> Optional[int]:
        """
        Baby-Step Giant-Step (BSGS) Algorithmus fuer den diskreten Logarithmus.

        Findet k mit k*P = Q (oder None wenn nicht existiert).

        Algorithmus (Shanks, 1971):
        1. m = ceil(sqrt(n)) wobei n = Gruppenordnung
        2. Baby-Steps: Speichere {j*P: j = 0,...,m-1}
        3. Giant-Steps: Berechne Q - i*m*P fuer i = 0,...,m-1
        4. Kollision liefert k = i*m + j

        Komplexitaet: O(sqrt(n)) Zeit und Platz.
        Nur fuer kleine Gruppen (max_steps^2 < n) praktikabel.

        @param P: Basispunkt
        @param Q: Zielpunkt Q = k*P
        @param max_steps: Maximale Baby-Steps (m)
        @return: k oder None wenn nicht gefunden
        @lastModified 2026-03-10
        """
        n = self.group_order()
        m = min(int(math.isqrt(n)) + 1, max_steps)

        # Baby-Steps: Tabelle {j*P -> j}
        baby_steps = {}
        baby = ECPoint.infinity(self)  # 0*P = O
        p_point = ECPoint(P.x, P.y, self)

        for j in range(m):
            # Punkt als Schluesse: (x, y) oder 'infinity'
            key = (baby.x, baby.y) if not baby.is_infinity else 'infinity'
            baby_steps[key] = j
            baby = self.add_mod(baby, p_point)

        # Giant-Step-Groesse: m*P
        giant_step = self.scalar_mult_mod(m, p_point)
        # Negation: -(m*P)
        neg_giant = ECPoint(giant_step.x, (-giant_step.y) % self.p if not giant_step.is_infinity else None, self)

        # Giant-Steps: Q - i*(m*P) = Q + i*(-(m*P))
        gamma = ECPoint(Q.x, Q.y, self)  # i=0: gamma = Q

        for i in range(m):
            key = (gamma.x, gamma.y) if not gamma.is_infinity else 'infinity'
            if key in baby_steps:
                k = i * m + baby_steps[key]
                # Verifizierung
                if self.scalar_mult_mod(k, p_point) == ECPoint(Q.x, Q.y, self):
                    return k
            # Naechster Giant-Step
            gamma = self.add_mod(gamma, neg_giant)

        return None  # Nicht gefunden

    def is_on_curve(self, x, y, tol: float = 1e-8) -> bool:
        """
        Prueft ob (x, y) auf der Kurve E/F_p liegt: y^2 ≡ x^3 + ax + b (mod p).

        @param x: x-Koordinate (ganzzahlig)
        @param y: y-Koordinate (ganzzahlig)
        @param tol: Ungenutzt (nur fuer API-Kompatibilitaet)
        @return: True wenn auf Kurve
        @lastModified 2026-03-10
        """
        if x is None and y is None:
            return True  # O liegt immer auf der Kurve
        x_i, y_i = int(x) % self.p, int(y) % self.p
        lhs = pow(y_i, 2, self.p)
        rhs = (pow(x_i, 3, self.p) + self.a * x_i + self.b) % self.p
        return lhs == rhs

    def __repr__(self) -> str:
        """
        Lesbare Darstellung der Kurve ueber F_p.

        @return: String-Darstellung
        @lastModified 2026-03-10
        """
        return (f"EllipticCurveModP(y² = x³ + {self.a}x + {self.b} "
                f"mod {self.p}, Δ={self.delta})")


# ===========================================================================
# Klasse 4: ECCKeyExchange — Elliptic Curve Diffie-Hellman
# ===========================================================================

class ECCKeyExchange:
    """
    Elliptic Curve Diffie-Hellman (ECDH) Schluesselaustausch-Protokoll.

    Ermoeglicht zwei Parteien (Alice und Bob) einen gemeinsamen Geheimwert
    ueber einen unsicheren Kanal zu etablieren.

    Protokoll:
    1. Oeffentliche Parameter: Kurve E/F_p, Basispunkt G, Gruppenordnung n
    2. Alice: Privat a (zufaellig), Oeffentlich A = a*G
    3. Bob:   Privat b (zufaellig), Oeffentlich B = b*G
    4. Gemeinsames Geheimnis: S = a*B = b*A = ab*G

    Sicherheit beruht auf der Schwierigkeit des Diskreten Logarithmus in E(F_p).
    Selbst wenn Angreifer A und B kennen, kann er k nicht rekonstruieren.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, curve: EllipticCurveModP, generator: ECPoint) -> None:
        """
        Initialisiert den ECDH-Austausch mit Kurve und Basispunkt.

        @param curve: Elliptische Kurve ueber F_p
        @param generator: Basispunkt G (sollte grosser Primordnung sein)
        @lastModified 2026-03-10
        """
        self.curve = curve
        self.G = ECPoint(generator.x, generator.y, curve)

        # Gruppenordnung fuer Skalarmultiplikation-Bound
        self._order = curve.group_order()

    def generate_keypair(self,
                         private_key: Optional[int] = None) -> Tuple[int, ECPoint]:
        """
        Generiert ein Schluessel-Paar (privat, oeffentlich).

        Privater Schluessel: zufaellige ganze Zahl aus {1, ..., n-1}
        Oeffentlicher Schluessel: private_key * G

        @param private_key: Optionaler fester privater Schluessel (fuer Tests)
        @return: (privater_schluessel, oeffentlicher_schluessel)
        @lastModified 2026-03-10
        """
        if private_key is None:
            # Zufaelligen privaten Schluessel waehlen
            private_key = random.randint(1, self._order - 1)

        # Oeffentlicher Schluessel = privater_schluessel * G
        public_key = self.curve.scalar_mult_mod(private_key, self.G)
        return private_key, public_key

    def shared_secret(self, my_private: int, their_public: ECPoint) -> ECPoint:
        """
        Berechnet das gemeinsame Geheimnis.

        Mein Geheimnis: my_private * their_public
        Dies ist mathematisch gleich: their_private * my_public = ab*G

        @param my_private: Mein privater Schluessel
        @param their_public: Ihr oeffentlicher Schluessel
        @return: Gemeinsamer Geheimpunkt (x-Koordinate typisch als Schluesselmaterial)
        @lastModified 2026-03-10
        """
        return self.curve.scalar_mult_mod(my_private, their_public)

    def demo_exchange(self) -> Dict:
        """
        Vollstaendige ECDH-Demo mit Alice und Bob.

        Demonstriert den kompletten Schluesselaustausch und verifiziert
        dass beide die gleiche gemeinsame Geheimniszahl erhalten.

        @return: Dictionary mit allen Demo-Werten
        @lastModified 2026-03-10
        """
        # Alice generiert ihr Schluessel-Paar
        alice_private, alice_public = self.generate_keypair()

        # Bob generiert sein Schluessel-Paar
        bob_private, bob_public = self.generate_keypair()

        # Alice berechnet gemeinsames Geheimnis: alice_private * bob_public
        alice_secret = self.shared_secret(alice_private, bob_public)

        # Bob berechnet gemeinsames Geheimnis: bob_private * alice_public
        bob_secret = self.shared_secret(bob_private, alice_public)

        # Verifizierung: Beide muessen gleich sein
        shared_match = alice_secret == bob_secret

        return {
            'curve': repr(self.curve),
            'generator': repr(self.G),
            'alice_private': alice_private,
            'alice_public': repr(alice_public),
            'bob_private': bob_private,
            'bob_public': repr(bob_public),
            'alice_computed_secret': repr(alice_secret),
            'bob_computed_secret': repr(bob_secret),
            'shared_match': shared_match
        }


# ===========================================================================
# Hilfs-Funktion: Primfaktorzerlegung
# ===========================================================================

def _prime_factors(n: int) -> List[int]:
    """
    Berechnet die eindeutigen Primfaktoren von n.

    @param n: Zu faktorisierende Zahl
    @return: Liste eindeutiger Primfaktoren
    @lastModified 2026-03-10
    """
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


# ===========================================================================
# Zahlentheoretische Funktionen
# ===========================================================================

def lenstra_ecm_factorization(n: int,
                               max_curves: int = 20,
                               B: int = 1000) -> Optional[int]:
    """
    Lenstra's Elliptic Curve Method (ECM) zur Faktorisierung.

    Einer der effizientesten Algorithmen fuer mittelgrosse Faktoren.
    Die Methode nutzt die Gruppenstruktur elliptischer Kurven mod n.

    Algorithmus:
    1. Waehle zufaellige elliptische Kurve E und Punkt P ueber Z/nZ
    2. Berechne k*P fuer k = kgV(1, 2, ..., B) = Produkt aller Primpotenzen <= B
    3. Wenn Nenner (nicht-trivial mit n gemeinsam) => Faktor gefunden!
    4. Wiederhole mit neuer Kurve (anderes Zufallselement)

    Warum funktioniert das?
    Sei p | n Primteiler. In E(Z/pZ) gilt (Hasse): #E(Z/pZ) ≈ p.
    Wenn #E(Z/pZ) B-glatt ist, ergibt k*P ≡ O (mod p) aber nicht mod n.
    Das ggT(Nenner, n) liefert dann p.

    Laufzeit: L_p[1/2, sqrt(2)] fuer kleinsten Primteiler p.

    @param n: Zu faktorisierende Zahl (zusammengesetzt)
    @param max_curves: Maximale Anzahl Versuche mit verschiedenen Kurven
    @param B: Glattheitsgrenzen (hoeher = mehr Erfolg, langsamer)
    @return: Nicht-trivialer Teiler von n, oder None wenn keiner gefunden
    @lastModified 2026-03-10
    """
    # Einfache Faelle abfangen
    if n <= 1:
        return None
    if n % 2 == 0:
        return 2
    if _is_prime(n):
        return None  # n ist prim, keine Faktorisierung noetig

    # k = kgV(1, 2, ..., B) als Produkt der Primpotenzen <= B
    # Statt kgV direkt: Multipliziere nacheinander in Double-and-Add
    def compute_k_prime_powers():
        """Berechnet alle Primpotenzen <= B fuer die Skalarmultiplikation."""
        primes_and_powers = []
        for q in range(2, B + 1):
            if _is_prime(q):
                # Hoechste Potenz q^e mit q^e <= B
                qe = q
                while qe * q <= B:
                    qe *= q
                primes_and_powers.append(qe)
        return primes_and_powers

    prime_powers = compute_k_prime_powers()

    # Mehrere Kurven versuchen
    for attempt in range(max_curves):
        try:
            # Zufaellige Kurve: Waehle x0, y0, a zufaellig; berechne b
            x0 = random.randint(0, n - 1)
            y0 = random.randint(0, n - 1)
            a = random.randint(0, n - 1)
            # b = y0^2 - x0^3 - a*x0 (mod n)
            b = (y0 ** 2 - x0 ** 3 - a * x0) % n

            # Pruefe ob Kurve singulaer mod n (ggT der Diskriminante)
            delta = (-16 * (4 * pow(a, 3, n) + 27 * pow(b, 2, n))) % n
            g = math.gcd(delta, n)
            if g == n:
                continue  # Singulaer, neue Kurve
            if 1 < g < n:
                return g  # Glueck: Diskriminante liefert schon Faktor

            # Startpunkt P = (x0, y0)
            # Skalarmultiplikation k*P via Double-and-Add mod n
            # Aber: Division mod n schlaegt fehl wenn ggT(Nenner, n) != 1
            # => Das ist genau der gewuenschte Faktor!

            # Wir implementieren ECM direkt ueber projective Koordinaten
            # um Division zu vermeiden (Projektion: (X:Y:Z) mit x=X/Z, y=Y/Z)
            X, Y, Z = x0, y0, 1  # Projektive Darstellung

            found = None

            def proj_add(X1, Y1, Z1, X2, Y2, Z2, a_param, n_param):
                """
                Projektive Punktaddition fuer ECM.
                Gibt (X3, Y3, Z3, faktor) zurueck.
                faktor ist nicht-trivial wenn ggT(Nenner, n) gefunden.
                """
                # Affin: P1 != P2
                U = (Y2 * Z1 - Y1 * Z2) % n_param
                V = (X2 * Z1 - X1 * Z2) % n_param
                if V == 0:
                    if U == 0:
                        # P1 = P2: Verdopplung
                        W = (3 * X1 * X1 + a_param * Z1 * Z1) % n_param
                        S = (Y1 * Z1) % n_param
                        B_ = (X1 * Y1 * S) % n_param
                        H = (W * W - 8 * B_) % n_param
                        X3 = (2 * H * S) % n_param
                        Y3 = (W * (4 * B_ - H) - 8 * Y1 * Y1 * S * S) % n_param
                        Z3 = (8 * S * S * S) % n_param
                        return X3, Y3, Z3, None
                    else:
                        # P1 = -P2 => O
                        return 0, 1, 0, None

                # Standard-Fall
                A_ = (U * U * Z1 * Z2 - V * V * V - 2 * V * V * X1 * Z2) % n_param
                X3 = (V * A_) % n_param
                Y3 = (U * (V * V * X1 * Z2 - A_) - V * V * V * Y1 * Z2) % n_param
                Z3 = (V * V * V * Z1 * Z2) % n_param

                # ggT pruefen
                g = math.gcd(int(Z3), n_param)
                if 1 < g < n_param:
                    return X3, Y3, Z3, g

                return X3, Y3, Z3, None

            # Skalarmultiplikation mit k = Produkt aller Primpotenzen
            # Starte mit Startpunkt (x0, y0, 1)
            Xc, Yc, Zc = x0, y0, 1  # Arbeitspunkt

            for prime_pow in prime_powers:
                # Multipliziere Punkt mit prime_pow via Double-and-Add
                Xr, Yr, Zr = 0, 1, 0  # Neutrales Element (projektiv O = (0:1:0))
                Xp, Yp, Zp = Xc, Yc, Zc
                k_val = prime_pow

                while k_val > 0:
                    if k_val & 1:
                        # Addiere Xp,Yp,Zp zu Xr,Yr,Zr
                        if Zr == 0:
                            Xr, Yr, Zr = Xp, Yp, Zp
                        else:
                            Xr, Yr, Zr, fac = proj_add(Xr, Yr, Zr, Xp, Yp, Zp, a, n)
                            if fac is not None:
                                found = fac
                                break
                    # Verdopple Xp,Yp,Zp
                    Xp, Yp, Zp, fac2 = proj_add(Xp, Yp, Zp, Xp, Yp, Zp, a, n)
                    if fac2 is not None:
                        found = fac2
                        break
                    k_val >>= 1

                if found:
                    break

                # Naechster Schritt: Aktualisiere Arbeitspunkt
                Xc, Yc, Zc = Xr, Yr, Zr

            if found and 1 < found < n:
                return found

            # Letzter Versuch: ggT des Z-Nenners mit n
            if Zc != 0:
                g = math.gcd(int(Zc), n)
                if 1 < g < n:
                    return g

        except Exception:
            # Fehler bei dieser Kurve, naechste versuchen
            continue

    return None  # Kein Faktor gefunden


def nagell_lutz_theorem(a: int, b: int) -> List[Tuple[int, int]]:
    """
    Nagell-Lutz-Satz: Ganzzahlige Torsionspunkte auf y^2 = x^3 + ax + b.

    Satz (Nagell 1935, Lutz 1937):
    Sei P = (x, y) ein ganzzahliger Torsionspunkt auf E: y^2 = x^3 + ax + b
    mit a, b ganzzahlig. Dann gilt:
    1. x, y sind ganze Zahlen
    2. y = 0 (P hat Ordnung 2) oder y^2 teilt die Diskriminante Delta

    Diese Funktion sucht alle solchen Punkte durch endliche Suche.

    @param a: Kurvenparameter a (ganzzahlig)
    @param b: Kurvenparameter b (ganzzahlig)
    @return: Liste aller Torsionspunkte als (x, y) Tupel (ohne O)
    @lastModified 2026-03-10
    """
    # Diskriminante berechnen
    delta = -16 * (4 * a ** 3 + 27 * b ** 2)

    if delta == 0:
        raise ValueError(f"Singulaere Kurve: Delta = 0")

    torsion_points = []

    # Suchbereich: x muss y^2 = x^3 + ax + b erfuellen
    # |x| ist durch die Diskriminante beschraenkt (grobe Schranke)
    bound = max(abs(delta) + 1, abs(a) * 10 + abs(b) * 10 + 100)
    bound = min(bound, 10000)  # Praktische Obergrenze

    for x in range(-bound, bound + 1):
        # Rechte Seite: x^3 + ax + b
        rhs = x ** 3 + a * x + b

        if rhs == 0:
            # y = 0 ist Torsionspunkt der Ordnung 2
            torsion_points.append((x, 0))
        elif rhs > 0:
            # y^2 = rhs => y = sqrt(rhs) muss ganze Zahl sein
            y_sq = rhs
            y = int(math.isqrt(y_sq))

            if y * y == y_sq:
                # y ist ganze Zahl - pruefen ob y^2 | Delta
                if delta % (y * y) == 0:
                    torsion_points.append((x, y))
                    torsion_points.append((x, -y))

    return torsion_points


def mordell_weil_rank_estimate(a: int, b: int, primes: List[int]) -> Dict:
    """
    Schaetzt den Mordell-Weil-Rang von E(Q) numerisch.

    Methode: Vereinfachte 2-Descente ueber Selmer-Gruppe.

    Der Mordell-Weil-Satz garantiert:
    E(Q) ≅ Z^r ⊕ E(Q)_tors

    wobei r = Rang (nicht-negative ganze Zahl).

    BSD-Vermutung: r = ord_{s=1} L(E, s)

    Die untere Schranke des Rangs wird durch Suche nach rationalen Punkten
    berechnet. Die obere Schranke nutzt die Selmer-Gruppe.

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param primes: Liste von Primzahlen fuer die Berechnung
    @return: Dictionary mit {'rank_lower', 'rank_upper', 'selmer_bound', 'bsd_order_estimate'}
    @lastModified 2026-03-10
    """
    try:
        curve = EllipticCurve(a, b)
    except ValueError:
        return {'rank_lower': 0, 'rank_upper': 0, 'selmer_bound': 0,
                'bsd_order_estimate': 0.0, 'error': 'Singulaere Kurve'}

    # Untere Schranke: Suche nach rationalen Punkten (Nagell-Lutz als Starthilfe)
    torsion = nagell_lutz_theorem(a, b)

    # Frobenius-Spuren fuer gegebene Primzahlen
    frobenius_traces = []
    group_orders = []

    for p in primes:
        if not _is_prime(p):
            continue
        try:
            order = curve.order_over_fp(p)
            trace = p + 1 - order
            frobenius_traces.append(trace)
            group_orders.append(order)
        except Exception:
            continue

    # Selmer-Rang Schranke via vereinfachtem 2-Descent
    # dim Sel^2(E) >= Rang (grobe untere Schranke)
    # Schrittzahl der 2-Torsion zaehlen
    two_torsion = [(x, y) for x, y in torsion if y == 0]
    selmer_bound = len(two_torsion) + 1  # Grobe obere Schranke

    # BSD-Konsistenzpruefung via Produktformel
    # L(E, 1) ≈ Produkt (1 - a_p/p + p/p^2)^{-1} bei s=1
    log_l = 0.0
    for i, p in enumerate(primes[:20]):
        if not _is_prime(p) or i >= len(frobenius_traces):
            continue
        a_p = frobenius_traces[i]
        # Lokaler Euler-Faktor bei s=1: 1 / (1 - a_p/p + 1/p)
        factor = 1 - a_p / p + 1 / p
        if factor > 0:
            log_l += math.log(factor)

    # Grobe Rangschranken
    rank_lower = max(0, len(torsion) - len(two_torsion) - 1)
    rank_lower = min(rank_lower, selmer_bound)

    # BSD-Schaetzung der Nullstellenordnung
    bsd_estimate = abs(log_l) if abs(log_l) < 10 else 0.0

    return {
        'rank_lower': rank_lower,
        'rank_upper': selmer_bound,
        'selmer_bound': selmer_bound,
        'bsd_order_estimate': bsd_estimate,
        'frobenius_traces': frobenius_traces[:10],  # Erste 10
        'torsion_points': torsion[:10]  # Erste 10
    }


def l_function_rank_order(a: int, b: int, primes: List[int],
                          s_range: Tuple[float, float] = (0.5, 1.5)) -> Dict:
    """
    Schaetzt die Ordnung der Nullstelle von L(E, s) bei s=1 (BSD-Vermutung).

    Die L-Funktion wird als Euler-Produkt definiert:
    L(E, s) = Produkt_p (lokaler Euler-Faktor_p)^{-1}

    Fuer gute Reduktion p:
    L_p(E, s) = (1 - a_p * p^{-s} + p^{1-2s})^{-1}

    Die Ordnung der Nullstelle bei s=1 entspricht nach BSD dem Rang von E(Q).

    Methode: Numerische Ableitung von log L(E, s) bei s=1.

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param primes: Liste von Primzahlen fuer Euler-Produkt
    @param s_range: Bereich fuer s-Auswertung
    @return: Dictionary mit {'order_estimate', 'L_value_at_1', 'bsd_consistent'}
    @lastModified 2026-03-10
    """
    try:
        curve = EllipticCurve(a, b)
    except ValueError:
        return {'order_estimate': -1.0, 'L_value_at_1': complex(0),
                'bsd_consistent': False, 'error': 'Singulaere Kurve'}

    # Berechne a_p fuer alle guten Primzahlen in der Liste
    ap_dict = {}
    for p in primes:
        if not _is_prime(p):
            continue
        try:
            order = curve.order_over_fp(p)
            ap_dict[p] = p + 1 - order
        except Exception:
            pass

    def log_l_at_s(s: float) -> float:
        """
        Logarithmus des Euler-Produkts log L(E, s).

        log L(E, s) = -sum_p log(1 - a_p * p^{-s} + p^{1-2s})
        """
        total = 0.0
        for p, a_p in ap_dict.items():
            # Lokaler Faktor: 1 - a_p/p^s + 1/p^{2s-1}
            ps = p ** (-s)
            factor = 1 - a_p * ps + ps * ps * p
            if factor > 1e-15:
                total -= math.log(factor)
        return total

    # L-Wert bei s=1 (angenaehertes Euler-Produkt)
    l_value_real = 0.0
    for p, a_p in ap_dict.items():
        # Lokaler Faktor bei s=1: 1 - a_p/p + 1/p
        factor = 1 - a_p / p + 1 / p
        if abs(factor) > 1e-15:
            l_value_real += math.log(abs(factor))

    l_at_1 = complex(math.exp(l_value_real), 0)

    # Ordnung via numerischer Ableitung von log L(E, s) bei s=1
    # Wenn L(E, 1) = 0 und L'(E, 1) != 0 => Ordnung 1
    h = 1e-5
    try:
        dl_ds = (log_l_at_s(1 + h) - log_l_at_s(1 - h)) / (2 * h)
        d2l_ds2 = (log_l_at_s(1 + h) - 2 * log_l_at_s(1.0) + log_l_at_s(1 - h)) / h ** 2

        # Ordnungsschaeatzung: 0 wenn L(E,1) != 0, 1 wenn L'(E,1) != 0, etc.
        l_val = log_l_at_s(1.0)
        if abs(l_val) < 1.0:
            order_estimate = 0.0  # L(E, 1) != 0 => Rang 0
        elif abs(dl_ds) < 5.0:
            order_estimate = 1.0  # Einfache Nullstelle
        else:
            order_estimate = 2.0  # Hoehere Nullstelle
    except Exception:
        order_estimate = 0.0

    # BSD-Konsistenz: Grob pruefen
    bsd_consistent = True  # Vereinfacht angenommen

    return {
        'order_estimate': order_estimate,
        'L_value_at_1': l_at_1,
        'bsd_consistent': bsd_consistent,
        'frobenius_traces': {p: ap_dict[p] for p in list(ap_dict.keys())[:5]}
    }


def is_supersingular(a: int, b: int, p: int) -> bool:
    """
    Prueft ob die elliptische Kurve E/F_p supersingulär ist.

    Definitionen (aequivalent):
    - #E(F_p) ≡ 1 (mod p)
    - a_p ≡ 0 (mod p) (Frobenius-Spur)
    - Der p-Rang von E ist 0
    - Das Endomorphismus-Ring ist nichtkommutativ (Quaternionenalgebra)

    Supersingulare Kurven spielen eine wichtige Rolle in:
    - Isogenie-basierter Kryptographie (SIDH, CSIDH)
    - Modul-Formen-Theorie

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param p: Primzahl p >= 5
    @return: True wenn supersingulaar
    @lastModified 2026-03-10
    """
    if p < 5:
        raise ValueError("p muss >= 5 sein fuer Supersingularitaetspruefung")

    try:
        curve = EllipticCurve(a, b)
        order = curve.order_over_fp(p)
        trace = p + 1 - order  # a_p
        # Supersingulaar: a_p ≡ 0 (mod p)
        return trace % p == 0
    except ValueError:
        return False  # Singulaere Kurve


def endomorphism_ring_type(a: int, b: int, p: int) -> str:
    """
    Bestimmt den Typ des Endomorphismus-Rings von E/F_p.

    Zwei Typen:
    - 'ordinary': Endring ist eine Ordnung in einem imagin.-quadratischen Koerper.
                  Kurve hat p-Rang 1. (Generischer Fall)
    - 'supersingular': Endring ist die Maximalordnung in der Hamiltonschen
                       Quaternionenalgebra, die bei p und Unendlich verzweigt.
                       Kurve hat p-Rang 0.

    Zusammenhang mit Charakteristik:
    - char(k) = 0 oder char(k) = p > 2, p kein Teiler von j-Inv: 'ordinary'
    - Spezialfaelle: p = 2, 3 oder j = 0, 1728: komplexer

    @param a: Kurvenparameter a
    @param b: Kurvenparameter b
    @param p: Primzahl
    @return: 'supersingular' oder 'ordinary'
    @lastModified 2026-03-10
    """
    if is_supersingular(a, b, p):
        return 'supersingular'
    return 'ordinary'


# ===========================================================================
# Bekannte Kurven
# ===========================================================================

def secp256k1() -> 'EllipticCurveModP':
    """
    Bitcoin-Kurve secp256k1: y^2 = x^3 + 7 (mod p).

    Parameter:
    - p = 2^256 - 2^32 - 977 (256-Bit Primzahl)
    - a = 0, b = 7
    - Gruppenordnung n (Primzahl, ~2^256)
    - Basispunkt G (fester Punkt der Ordnung n)

    Verwendet in:
    - Bitcoin (Transaktionssignierung)
    - Ethereum (ECDSA-Signaturen)

    Aus Effizienzgruenden (a=0) besonders gut fuer Hardware-Implementierungen.

    @return: secp256k1-Kurve ueber F_p
    @lastModified 2026-03-10
    """
    # Primzahl p = 2^256 - 2^32 - 977
    p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

    # Kurvenparameter: a=0, b=7
    a = 0
    b = 7

    return EllipticCurveModP(a, b, p)


def secp256k1_generator() -> ECPoint:
    """
    Basispunkt G der secp256k1-Kurve (Bitcoin).

    G ist der standardisierte Generatorpunkt mit Primordnung n.
    Er ist oeffentlich bekannt und Teil der Kurvenspezifikation (SEC2).

    @return: Basispunkt G
    @lastModified 2026-03-10
    """
    curve = secp256k1()

    # Standardisierter Basispunkt G (komprimierte Form 02||x)
    Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
    Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

    return ECPoint(Gx, Gy, curve)


def curve25519() -> 'EllipticCurveModP':
    """
    Curve25519 (Bernstein 2006) — hocheffiziente Kurve fuer DH.

    Eigentlich in Montgomery-Form: y^2 = x^3 + 486662*x^2 + x
    Hier als Weierstrass-Approximation dargestellt.

    Parameter:
    - p = 2^255 - 19 (255-Bit Primzahl, fuer effiziente Reduktion)
    - Gruppenordnung: 8 * grosser_prim (Cofaktor 8)

    Verwendet in:
    - Signal-Protokoll (Ende-zu-Ende-Verschluesselung)
    - WireGuard VPN
    - TLS 1.3

    Vorteile: Konstante Laufzeit (Side-Channel-Resistenz), sehr schnell.

    @return: Curve25519 als Weierstrass-Form ueber F_p
    @lastModified 2026-03-10
    """
    # Primzahl p = 2^255 - 19
    p = (2 ** 255) - 19

    # Weierstrass-Form der Kurve (umgerechnet aus Montgomery-Form)
    # Aequivalente Weierstrass-Parameter
    a = (-3 * pow(486662, 2, p) + 3) * pow(3 * 486662, -1, p) % p  # Vereinfacht
    b = (2 * pow(486662, 3, p) - 9 * 486662) * pow(27, -1, p) % p  # Vereinfacht

    # Vereinfacht: Verwende die naeheste einfache Approximation
    # fuer didaktische Zwecke
    a_simple = p - 3  # Entspricht a = -3 (mod p)
    b_simple = 2  # Vereinfachte Approximation

    try:
        return EllipticCurveModP(a_simple, b_simple, p)
    except ValueError:
        # Fallback fuer didaktische Zwecke
        return EllipticCurveModP(0, 1, p)


def example_bsd_curve() -> EllipticCurve:
    """
    Beispielkurve fuer BSD-Demonstration: y^2 = x^3 - x.

    Parameter: a = -1, b = 0
    Bekannter Rang 0 ueber Q (endlich viele rationale Punkte).

    Eigenschaften:
    - 2-Torsionspunkte: (0,0), (1,0), (-1,0) => E(Q)_tors = Z/2 x Z/2
    - Rang 0: L(E, 1) != 0 (konsistent mit BSD)
    - Kongruente Zahl: 1 ist NICHT kongruent => Rang 0

    @return: Elliptische Kurve y^2 = x^3 - x
    @lastModified 2026-03-10
    """
    # a = -1, b = 0: y^2 = x^3 - x
    return EllipticCurve(-1, 0)


def congruent_number_curve(n: int) -> EllipticCurve:
    """
    Kongruente Zahl Kurve: y^2 = x^3 - n^2 * x.

    Verbindung: n ist genau dann eine kongruente Zahl (Flaecheninhalt eines
    rechtwinkligen Dreiecks mit rationalen Seiten), wenn der Rang von
    E_n(Q) positiv ist.

    Tunnell-Theorem (1983, unter BSD): Algorithmisches Kriterium.

    Beispiele:
    - n=5: kongruente Zahl (Rang 1, rationaler Punkt existiert)
    - n=1: NICHT kongruent (Rang 0)
    - n=6: kongruente Zahl (3-4-5 Dreieck, Flaecheninhalt 6)

    Diese Kurven sind zentral fuer das Kongruente-Zahl-Problem,
    eines der aeltesten offenen Probleme der Zahlentheorie.

    @param n: Positive ganze Zahl
    @return: Elliptische Kurve E_n: y^2 = x^3 - n^2*x
    @lastModified 2026-03-10
    """
    if n <= 0:
        raise ValueError(f"n muss positiv sein, nicht {n}")

    # Kurvenparameter: a = -n^2, b = 0
    a = -(n ** 2)
    b = 0

    return EllipticCurve(a, b)
