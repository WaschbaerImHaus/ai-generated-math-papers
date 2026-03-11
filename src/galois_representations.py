"""
@file galois_representations.py
@brief Galois-Darstellungen: ρ: Gal(Q̄/Q) → GL_n(ℤ_p), Langlands-Programm, Tate-Modul.
@description
    Implementiert die zentralen Konzepte der Theorie der Galois-Darstellungen,
    welche eine der tiefsten Verbindungen zwischen Zahlentheorie, algebraischer
    Geometrie und automorphen Formen herstellt.

    Eine **ℓ-adische Galois-Darstellung** ist ein stetiger Gruppenhomomorphismus

        ρ: Gal(Q̄/Q) → GL_n(ℤ_ℓ)

    wobei Gal(Q̄/Q) die absolute Galois-Gruppe von ℚ ist, versehen mit der
    profiniten Topologie (projektiver Limes endlicher Galois-Gruppen), und
    GL_n(ℤ_ℓ) die invertierbare n×n-Matrizen über den ℓ-adischen ganzen Zahlen.

    **Frobenius-Elemente**: Für eine bei p unverzweigte Darstellung existiert
    ein wohlbestimmtes konjugiertes Frobenius-Element Frob_p ∈ Gal(Q̄/Q) mit

        Frob_p(ζ) = ζ^p  für alle ℓ^n-ten Einheitswurzeln ζ.

    Die **Spur** tr(ρ(Frob_p)) und die **charakteristischen Polynome** det(1 - ρ(Frob_p)X)
    codieren wichtige arithmetische Information.

    **Tate-Modul**: Für eine elliptische Kurve E/ℚ ist der ℓ-adische Tate-Modul

        T_ℓ(E) = lim_{←} E[ℓ^n]  ≅  ℤ_ℓ²

    als ℤ_ℓ-Modul, und liefert eine 2-dimensionale Galois-Darstellung.
    Die Frobenius-Spur ist a_p = p + 1 - #E(𝔽_p) (bei guter Reduktion).

    **Langlands-Programm**: Die Vermutung (teilweise bewiesen durch Wiles, Taylor,
    Harris-Taylor, Ngô u.a.) besagt, dass jede "motivische" Galois-Darstellung
    einer automorphen Darstellung einer algebraischen Gruppe entspricht.

    **Ramanujan-Vermutung**: Für die Tate-Modul-Darstellung einer elliptischen Kurve
    gilt |a_p| ≤ 2√p (bewiesen durch Hasse für elliptische Kurven, allgemein durch
    Deligne als Teil von Weil II).

    Dieses Modul implementiert:
    1. GaloisRepresentation      — Abstrakte Basisklasse
    2. TateModuleRepresentation  — Tate-Modul einer elliptischen Kurve
    3. CyclotomicCharacter       — Zyklotomischer Charakter χ_ℓ
    4. DirichletCharacterRepresentation — 1-dim. Darstellung via Dirichlet-Charakter
    5. SymmetricPowerRepresentation     — Sym^k einer 2-dim. Darstellung
    6. LanglandsCorrespondence   — Numerische Evidenz für Langlands
    7. Hilfsfunktionen           — p_adic_matrix, galois_group_order, artin_conductor, ...

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple
from functools import lru_cache

import numpy as np
from sympy import (
    symbols, Poly, factor_list, ZZ, isprime as sympy_isprime,
    galois_group, Rational, sqrt as sym_sqrt, I
)
from sympy.abc import x as sym_x


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _is_prime(n: int) -> bool:
    """
    @brief Primzahltest via Trial Division (für kleine n) und sympy (für große n).
    @param n Zu testende Zahl
    @return True wenn n eine Primzahl ist
    @lastModified 2026-03-11
    """
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _primes_up_to(n: int) -> List[int]:
    """
    @brief Erzeugt alle Primzahlen bis n via Sieb des Eratosthenes.
    @param n Obere Schranke
    @return Liste aller Primzahlen ≤ n
    @lastModified 2026-03-11
    """
    if n < 2:
        return []
    # Sieb initialisieren
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            # Alle Vielfachen von i markieren
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return [i for i in range(2, n + 1) if sieve[i]]


def _sym_power_trace(alpha: complex, beta: complex, k: int) -> complex:
    """
    @brief Berechnet die Spur der k-ten symmetrischen Potenz Sym^k(ρ) bei einem
           2-dimensionalen Modul mit Eigenvalues α, β via der Formel:
               tr(Sym^k(ρ)(Frob_p)) = Σ_{j=0}^{k} α^j · β^{k-j}
    @param alpha Erster Eigenwert des Frobenius
    @param beta  Zweiter Eigenwert des Frobenius
    @param k     Potenz der symmetrischen Potenz
    @return Spur des Frobenius auf Sym^k
    @lastModified 2026-03-11
    """
    # Summe α^j * β^(k-j) für j = 0, ..., k
    result = sum(alpha**j * beta**(k - j) for j in range(k + 1))
    return result


# ===========================================================================
# 1. ABSTRAKTE BASISKLASSE
# ===========================================================================

class GaloisRepresentation(ABC):
    """
    @brief Abstrakte Basisklasse für ℓ-adische Galois-Darstellungen.
    @description
        Modelliert einen stetigen Gruppenhomomorphismus
            ρ: Gal(Q̄/Q) → GL_n(ℤ_ℓ)
        über seine lokalen Daten (Frobenius-Spuren) bei guten Primzahlen.
        Konkrete Unterklassen implementieren trace_frobenius() für spezifische
        mathematische Objekte (elliptische Kurven, Dirichlet-Charaktere etc.).
    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, ell: int, dimension: int, name: str) -> None:
        """
        @brief Initialisiert eine Galois-Darstellung.
        @param ell       Primzahl ℓ (charakteristisch für die p-adischen Koeffizienten)
        @param dimension Dimension n der Darstellung (Größe der GL_n-Matrizen)
        @param name      Bezeichner der Darstellung
        @raises ValueError wenn ell keine Primzahl oder dimension < 1
        @lastModified 2026-03-11
        """
        if not _is_prime(ell):
            raise ValueError(f"ell muss eine Primzahl sein, erhalten: {ell}")
        if dimension < 1:
            raise ValueError(f"Dimension muss ≥ 1 sein, erhalten: {dimension}")
        # Primzahl ℓ der Darstellung (Koeffizientenring ℤ_ℓ)
        self.ell = ell
        # Dimension der Darstellung (= Rang des GL_n)
        self.dim = dimension
        # Name/Bezeichner der Darstellung
        self.name = name

    @abstractmethod
    def trace_frobenius(self, p: int) -> complex:
        """
        @brief Berechnet die Frobenius-Spur tr(ρ(Frob_p)) an der Primzahl p.
        @description
            Für fast alle Primzahlen p (alle p außer endlich vielen Verzweigungsprimzahlen)
            ist Frob_p ∈ Gal(Q̄/Q) ein wohldefiniertes konjugiertes Element.
            Die Spur tr(ρ(Frob_p)) ∈ ℤ_ℓ ist eine fundamentale arithmetische Invariante.
        @param p Primzahl p ≠ ℓ
        @return Frobenius-Spur (komplex für allgemeine Darstellungen)
        @lastModified 2026-03-11
        """
        pass

    def is_unramified(self, p: int) -> bool:
        """
        @brief Prüft, ob die Darstellung bei p unverzweigt ist.
        @description
            Eine Darstellung ρ ist bei p unverzweigt, wenn die Trägheitsgruppe
            I_p ⊂ Gal(Q̄/Q) trivial auf dem Darstellungsraum wirkt, d.h.
            ρ|_{I_p} = id. In diesem Fall ist ρ(Frob_p) wohldefiniert.
            Standardmäßig: unverzweigt bei p ≠ ℓ. Unterklassen können überschreiben.
        @param p Zu prüfende Primzahl
        @return True wenn bei p unverzweigt
        @lastModified 2026-03-11
        """
        # Standard: unverzweigt für alle p ≠ ℓ
        return p != self.ell

    def l_function_euler_factor(self, p: int, s: complex) -> complex:
        """
        @brief Berechnet den lokalen L-Faktor bei p: L_p(ρ, s) = det(1 - Frob_p · p^{-s})^{-1}.
        @description
            Für eine n-dimensionale Darstellung mit Frobenius-Eigenwerten α_1,...,α_n ist
                L_p(ρ, s) = ∏_{i=1}^{n} (1 - α_i · p^{-s})^{-1}
                           = det(I_n - ρ(Frob_p) · p^{-s})^{-1}
            Für dim=2 mit charakteristischem Polynom x² - a_p x + p:
                L_p(ρ, s) = (1 - a_p·p^{-s} + p·p^{-2s})^{-1}
            Dies ist der Baustein des Euler-Produkts L(ρ, s) = ∏_p L_p(ρ, s).
        @param p Primzahl
        @param s Komplexes Argument der L-Funktion
        @return Wert des lokalen Euler-Faktors L_p(ρ, s) ∈ ℂ
        @lastModified 2026-03-11
        """
        if not self.is_unramified(p):
            # Verzweigte Primzahlen: vereinfachter lokaler Faktor
            return 1.0 + 0j

        # p^{-s} berechnen
        p_to_minus_s = p ** (-s)

        if self.dim == 1:
            # 1-dimensionaler Fall: L_p = (1 - α·p^{-s})^{-1}
            alpha = self.trace_frobenius(p)
            denom = 1.0 - alpha * p_to_minus_s
            if abs(denom) < 1e-15:
                return float('inf')
            return 1.0 / denom

        elif self.dim == 2:
            # 2-dimensionaler Fall: L_p = (1 - a_p·p^{-s} + p·p^{-2s})^{-1}
            a_p = self.trace_frobenius(p)
            denom = 1.0 - a_p * p_to_minus_s + p * (p_to_minus_s ** 2)
            if abs(denom) < 1e-15:
                return float('inf')
            return 1.0 / denom

        else:
            # Allgemeiner Fall: Näherung via Spur (nur für dim ≤ 4 sinnvoll)
            a_p = self.trace_frobenius(p)
            # Newton-Identitäten Näherung: det(1-AX) ≈ 1 - tr(A)X + ...
            denom = 1.0 - a_p * p_to_minus_s
            if abs(denom) < 1e-15:
                return float('inf')
            return 1.0 / denom

    def __repr__(self) -> str:
        """@brief String-Repräsentation der Darstellung. @lastModified 2026-03-11"""
        return f"GaloisRepresentation(name='{self.name}', ell={self.ell}, dim={self.dim})"


# ===========================================================================
# 2. TATE-MODUL EINER ELLIPTISCHEN KURVE
# ===========================================================================

class TateModuleRepresentation(GaloisRepresentation):
    """
    @brief Tate-Modul T_ℓ(E) einer elliptischen Kurve E/ℚ als Galois-Darstellung.
    @description
        Für eine elliptische Kurve E/ℚ ist der ℓ-adische Tate-Modul definiert als

            T_ℓ(E) = lim_{←n} E[ℓ^n]

        der projektive Limes der ℓ^n-Torsionspunkte. Als abelsche Gruppe gilt

            T_ℓ(E) ≅ ℤ_ℓ²  (2-dimensionaler ℤ_ℓ-Modul)

        und Gal(Q̄/Q) wirkt auf T_ℓ(E), was eine Darstellung

            ρ_{E,ℓ}: Gal(Q̄/Q) → GL_2(ℤ_ℓ)

        liefert. Bei guter Reduktion an p ≠ ℓ gilt für das charakteristische Polynom
        des Frobenius:

            det(1 - Frob_p · X | T_ℓ(E)) = 1 - a_p·X + p·X²

        mit a_p = p + 1 - #E(𝔽_p) (Hasse-Weil). Die Hasse-Schranke lautet |a_p| ≤ 2√p.

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, a_values: Dict[int, int], conductor: int, ell: int) -> None:
        """
        @brief Initialisiert den Tate-Modul einer elliptischen Kurve.
        @param a_values  Dict p → a_p mit den Frobenius-Spuren für gute Primzahlen p
                         (a_p = p + 1 - #E(𝔽_p))
        @param conductor Führer N der elliptischen Kurve (Produkt der Verzweigungsprimzahlen
                         mit zugehörigen Exponenten)
        @param ell       Primzahl ℓ für den Tate-Modul (ℓ ≠ charakteristik der Reduktion)
        @raises ValueError wenn conductor < 1
        @lastModified 2026-03-11
        """
        super().__init__(ell, 2, f"TateModule(N={conductor}, ell={ell})")
        if conductor < 1:
            raise ValueError(f"Führer muss ≥ 1 sein, erhalten: {conductor}")
        # a_p-Werte: p → a_p = p + 1 - #E(𝔽_p)
        self.a_values: Dict[int, int] = dict(a_values)
        # Führer der elliptischen Kurve
        self.conductor: int = conductor

    def trace_frobenius(self, p: int) -> complex:
        """
        @brief Gibt a_p = tr(Frob_p | T_ℓ(E)) zurück.
        @description
            Bei guter Reduktion (p ∤ N) ist a_p = p + 1 - #E(𝔽_p) ∈ ℤ.
            Bei schlechter Reduktion (p | N) gilt:
            - Additive Reduktion:     a_p = 0
            - Multiplikative Reduktion split:     a_p = +1
            - Multiplikative Reduktion non-split: a_p = -1
            Hier wird 0 zurückgegeben, wenn p nicht im a_values-Dict vorhanden.
        @param p Primzahl
        @return a_p ∈ ℤ ⊂ ℂ
        @lastModified 2026-03-11
        """
        # Aus gespeichertem Dict lesen, sonst 0 (schlechte Reduktion)
        return complex(self.a_values.get(p, 0))

    def characteristic_polynomial_frobenius(self, p: int) -> Tuple[float, float, float]:
        """
        @brief Gibt die Koeffizienten des charakteristischen Polynoms x² - a_p·x + p zurück.
        @description
            Bei guter Reduktion (p ∤ N) hat Frob_p das charakteristische Polynom

                P_p(X) = X² - a_p·X + p  ∈  ℤ[X]

            Dies folgt aus dem Weil-Theorem und der Tate-Konjektur (bewiesen durch Faltings).
            Die Wurzeln α, β erfüllen α·β = p und α + β = a_p, sowie |α| = |β| = √p
            (Hasse-Deligne-Riemann-Hypothese für elliptische Kurven).
        @param p Primzahl bei guter Reduktion
        @return Tuple (1, -a_p, p) als Koeffizienten von x² - a_p·x + p
        @lastModified 2026-03-11
        """
        a_p = self.a_values.get(p, 0)
        # Rückgabe: (Leitkoeffizient, linearer Term, konstanter Term)
        return (1.0, float(-a_p), float(p))

    def is_unramified(self, p: int) -> bool:
        """
        @brief Prüft ob E gute Reduktion bei p hat (äquivalent zu p ∤ N und p ≠ ℓ).
        @param p Zu prüfende Primzahl
        @return True wenn gute Reduktion und p ≠ ℓ
        @lastModified 2026-03-11
        """
        return (self.conductor % p != 0) and (p != self.ell)

    def frobenius_eigenvalues(self, p: int) -> Tuple[complex, complex]:
        """
        @brief Berechnet die Frobenius-Eigenwerte α, β mit α+β=a_p und α·β=p.
        @description
            Löst x² - a_p·x + p = 0 via Mitternachtsformel.
            Nach Hasse gilt |α| = |β| = √p (sofern die Hasse-Schranke erfüllt ist).
        @param p Primzahl bei guter Reduktion
        @return Paar (α, β) der komplexen Frobenius-Eigenwerte
        @lastModified 2026-03-11
        """
        a_p = self.a_values.get(p, 0)
        # Diskriminante des charakteristischen Polynoms: a_p² - 4p
        discriminant = a_p * a_p - 4 * p
        sqrt_disc = cmath.sqrt(discriminant)
        alpha = (a_p + sqrt_disc) / 2
        beta = (a_p - sqrt_disc) / 2
        return (alpha, beta)


# ===========================================================================
# 3. ZYKLOTOMISCHER CHARAKTER
# ===========================================================================

class CyclotomicCharacter:
    """
    @brief Zyklotomischer Charakter χ_ℓ: Gal(Q̄/Q) → ℤ_ℓ*.
    @description
        Der zyklotomische Charakter χ_ℓ ist definiert durch die Wirkung von
        Gal(Q̄/Q) auf den ℓ-Potenz-Einheitswurzeln: Für σ ∈ Gal(Q̄/Q) und
        eine primitive ℓ^n-te Einheitswurzel ζ gilt

            σ(ζ) = ζ^{χ_ℓ(σ)}

        Das Frobenius-Element Frob_p (für p ≠ ℓ) wirkt durch

            χ_ℓ(Frob_p) = p  (in ℤ_ℓ*)

        d.h. Frob_p(ζ) = ζ^p, was genau der Frobenius-Endomorphismus x ↦ x^p
        auf 𝔽_p ist. Der zyklotomische Charakter ist die "tautologische"
        1-dimensionale Galois-Darstellung, der Prototyp eines Hecke-Charakters.

        **Parität**: Die komplexe Konjugation c ∈ Gal(ℂ/ℝ) wirkt durch
        ζ ↦ ζ^{-1}, also χ_ℓ(c) = -1. Daher heißt χ_ℓ "ungerade" (odd).

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, ell: int) -> None:
        """
        @brief Initialisiert den zyklotomischen Charakter mod ℓ.
        @param ell Primzahl ℓ
        @raises ValueError wenn ell keine Primzahl ist
        @lastModified 2026-03-11
        """
        if not _is_prime(ell):
            raise ValueError(f"ell muss eine Primzahl sein, erhalten: {ell}")
        # Primzahl ℓ des zyklotomischen Charakters
        self.ell: int = ell

    def evaluate(self, p: int) -> int:
        """
        @brief Wertet den zyklotomischen Charakter am Frobenius Frob_p aus.
        @description
            Gibt χ_ℓ(Frob_p) = p mod ℓ zurück, da Frob_p auf ℓ-ten Einheitswurzeln
            durch den Frobenius-Endomorphismus ζ ↦ ζ^p wirkt. In ℤ_ℓ* entspricht
            dies dem Bild von p unter dem Reduktionshomomorphismus ℤ → ℤ/ℓℤ.
        @param p Primzahl (≠ ℓ)
        @return p mod ℓ ∈ {1, ..., ℓ-1} (p ≠ ℓ garantiert ≠ 0)
        @lastModified 2026-03-11
        """
        return p % self.ell

    def is_odd(self) -> bool:
        """
        @brief Gibt zurück, ob der Charakter ungerade (odd) ist.
        @description
            Ein Charakter χ: Gal(Q̄/Q) → ℂ* heißt ungerade, wenn
            χ(komplexe Konjugation) = -1. Der zyklotomische Charakter χ_ℓ
            ist stets ungerade, da komplexe Konjugation auf Einheitswurzeln
            durch ζ ↦ ζ^{-1} wirkt, also χ_ℓ(c) = -1 in ℤ_ℓ*.
        @return True (zyklotomischer Charakter ist immer ungerade)
        @lastModified 2026-03-11
        """
        # Zyklotomischer Charakter ist stets ungerade (komplexe Konjugation → -1)
        return True

    def tensor_power(self, n: int) -> "CyclotomicCharacterPower":
        """
        @brief Gibt χ_ℓ^n, die n-te Tensorpotenz des zyklotomischen Charakters, zurück.
        @param n Exponent (kann negativ sein für Twist)
        @return Objekt, das χ_ℓ^n repräsentiert
        @lastModified 2026-03-11
        """
        return CyclotomicCharacterPower(self.ell, n)

    def __repr__(self) -> str:
        """@brief String-Repräsentation. @lastModified 2026-03-11"""
        return f"CyclotomicCharacter(ell={self.ell})"


class CyclotomicCharacterPower:
    """
    @brief Tensorpotenz χ_ℓ^n des zyklotomischen Charakters.
    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, ell: int, power: int) -> None:
        """
        @brief Initialisiert χ_ℓ^n.
        @param ell   Primzahl ℓ
        @param power Exponent n
        @lastModified 2026-03-11
        """
        self.ell = ell
        self.power = power

    def evaluate(self, p: int) -> int:
        """
        @brief Gibt χ_ℓ^n(Frob_p) = p^n mod ℓ zurück.
        @param p Primzahl
        @return p^n mod ℓ
        @lastModified 2026-03-11
        """
        return pow(p, self.power, self.ell)


# ===========================================================================
# 4. DIRICHLET-CHARAKTER ALS GALOIS-DARSTELLUNG
# ===========================================================================

class DirichletCharacterRepresentation(GaloisRepresentation):
    """
    @brief 1-dimensionale Galois-Darstellung via Dirichlet-Charakter.
    @description
        Nach dem Kronecker-Weber-Theorem ist jede abelsche Erweiterung K/ℚ
        in einem zyklotomischen Körper ℚ(ζ_n) enthalten. Die zugehörigen
        1-dimensionalen Galois-Darstellungen

            ρ_χ: Gal(ℚ(ζ_q)/ℚ) ≅ (ℤ/qℤ)* → ℂ*

        entsprechen genau den Dirichlet-Charakteren mod q via Class-Field-Theory.

        Die L-Funktion ist:
            L(ρ_χ, s) = L(s, χ) = ∑_{n=1}^∞ χ(n)/n^s

        Das Frobenius-Element Frob_p (für p ∤ q) entspricht n ↦ p mod q,
        also ρ_χ(Frob_p) = χ(p).

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, chi_values: Dict[int, complex], conductor: int, ell: int) -> None:
        """
        @brief Initialisiert die Dirichlet-Charakter-Darstellung.
        @param chi_values Dict n → χ(n) mit den Charakterwerten (n von 1 bis conductor)
        @param conductor  Führer q des Dirichlet-Charakters
        @param ell        Primzahl ℓ (für Koeffizientenring)
        @raises ValueError wenn conductor < 1
        @lastModified 2026-03-11
        """
        super().__init__(ell, 1, f"DirichletChar(q={conductor}, ell={ell})")
        if conductor < 1:
            raise ValueError(f"Führer muss ≥ 1 sein, erhalten: {conductor}")
        # Charakterwerte χ(n) für n = 1, ..., q
        self.chi_values: Dict[int, complex] = dict(chi_values)
        # Führer des Dirichlet-Charakters
        self.conductor: int = conductor

    def trace_frobenius(self, p: int) -> complex:
        """
        @brief Gibt χ(p) = ρ_χ(Frob_p) zurück.
        @description
            Das Frobenius-Element bei p ∤ q entspricht p mod q in (ℤ/qℤ)*.
            Der Charakter nimmt dort den Wert χ(p mod q) ∈ ℂ* an.
            Bei p | q (Verzweigung) ist χ(p) = 0 (der Charakter verschwindet).
        @param p Primzahl
        @return χ(p) wenn p ∤ q, sonst 0
        @lastModified 2026-03-11
        """
        if self.conductor % p == 0:
            # Verzweigte Primzahl: Charakter verschwindet
            return 0j
        # Periodiziät: χ(p) = χ(p mod q)
        residue = p % self.conductor
        return complex(self.chi_values.get(residue, 0))

    def is_unramified(self, p: int) -> bool:
        """
        @brief Dirichlet-Charakter-Darstellung ist bei p unverzweigt ⟺ p ∤ q.
        @param p Primzahl
        @return True wenn p ∤ conductor und p ≠ ℓ
        @lastModified 2026-03-11
        """
        return (self.conductor % p != 0) and (p != self.ell)


# ===========================================================================
# 5. SYMMETRISCHE POTENZ EINER 2-DIM. DARSTELLUNG
# ===========================================================================

class SymmetricPowerRepresentation(GaloisRepresentation):
    """
    @brief Sym^k einer 2-dimensionalen Galois-Darstellung.
    @description
        Für eine 2-dimensionale Darstellung ρ: Gal(Q̄/Q) → GL_2(ℤ_ℓ) mit
        Frobenius-Eigenwerten α, β (bei guter Primzahl p) gilt

            Sym^k(ρ): Gal(Q̄/Q) → GL_{k+1}(ℤ_ℓ)

        mit Frobenius-Spur

            tr(Sym^k(ρ)(Frob_p)) = Σ_{j=0}^{k} α^j · β^{k-j}

        Dies ist ein Spezialfall der **Langlands-Funktorialität** (Sym^k-Lifting):
        Die Vermutung, dass Sym^k(π_E) eine automorphe Form auf GL_{k+1} ist,
        ist für k ≤ 4 bewiesen (Kim-Shahidi) und allgemein offen.

        **Ramanujan-Implikation**: Aus |a_p| ≤ 2√p folgt für die Eigenwerte
        α = √p·e^{iθ}, β = √p·e^{-iθ}, und
            tr(Sym^k(Frob_p)) = p^{k/2} · sin((k+1)θ)/sin(θ)

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, base_rep: GaloisRepresentation, k: int) -> None:
        """
        @brief Initialisiert Sym^k der Basis-Darstellung.
        @param base_rep Basis-Darstellung (muss 2-dimensional sein)
        @param k        Grad der symmetrischen Potenz (k ≥ 0)
        @raises ValueError wenn base_rep nicht 2-dimensional oder k < 0
        @lastModified 2026-03-11
        """
        if base_rep.dim != 2:
            raise ValueError(
                f"Basis-Darstellung muss 2-dimensional sein, "
                f"erhalten: dim={base_rep.dim}"
            )
        if k < 0:
            raise ValueError(f"k muss ≥ 0 sein, erhalten: {k}")
        # Sym^k hat Dimension k+1
        super().__init__(
            base_rep.ell,
            k + 1,
            f"Sym^{k}({base_rep.name})"
        )
        # Basis-Darstellung (2-dimensional)
        self.base_rep = base_rep
        # Grad der symmetrischen Potenz
        self.k = k

    def dimension(self) -> int:
        """
        @brief Gibt die Dimension von Sym^k zurück: k+1.
        @description
            Sym^k(V) für einen 2-dimensionalen Vektorraum V hat Dimension k+1,
            da eine Basis durch {v₁^j ⊗ v₂^{k-j} : j=0,...,k} gegeben ist.
        @return k+1
        @lastModified 2026-03-11
        """
        return self.k + 1

    def trace_frobenius(self, p: int) -> complex:
        """
        @brief Berechnet tr(Sym^k(ρ)(Frob_p)) via Newton-Power-Sum-Formel.
        @description
            Seien α, β die Frobenius-Eigenwerte der Basis-Darstellung bei p.
            Die Spur von Sym^k(Frob_p) ist:

                tr = Σ_{j=0}^{k} α^j · β^{k-j}

            Diese "geometrische Summe" in α/β lässt sich auch schreiben als:
            - Falls α = β: (k+1)·α^k
            - Sonst: (α^{k+1} - β^{k+1}) / (α - β)

            Für den Tate-Modul einer elliptischen Kurve mit a_p = α + β und
            det = α·β = p (bei guter Reduktion) berechnen wir α, β explizit.
        @param p Primzahl bei guter Reduktion
        @return Spur des Frobenius auf Sym^k
        @lastModified 2026-03-11
        """
        if not self.base_rep.is_unramified(p):
            return 0j

        # Basis-Darstellung muss TateModuleRepresentation sein für genaue Eigenwerte
        if isinstance(self.base_rep, TateModuleRepresentation):
            alpha, beta = self.base_rep.frobenius_eigenvalues(p)
        else:
            # Allgemeinfall: a_p bekannt, Determinante unbekannt → Näherung
            a_p = self.base_rep.trace_frobenius(p)
            # Für unbekannte Determinante: nehme p (Standard für Gewicht-2-Formen)
            det_p = complex(p)
            disc = a_p * a_p - 4 * det_p
            sqrt_disc = cmath.sqrt(disc)
            alpha = (a_p + sqrt_disc) / 2
            beta = (a_p - sqrt_disc) / 2

        return _sym_power_trace(alpha, beta, self.k)


# ===========================================================================
# 6. LANGLANDS-KORRESPONDENZ (NUMERISCHE EVIDENZ)
# ===========================================================================

class LanglandsCorrespondence:
    """
    @brief Numerische Evidenz für das Langlands-Programm.
    @description
        Das **Langlands-Programm** ist eine weitreichende Reihe von Vermutungen und
        Theoremen, die "motivische" Galois-Darstellungen mit "automorphen" Darstellungen
        algebraischer Gruppen verknüpfen.

        Zentrale Instanz: Für eine elliptische Kurve E/ℚ ist die zugehörige 2-dimensionale
        Galois-Darstellung ρ_{E,ℓ} einer Modulform f_E von Gewicht 2 zugeordnet
        (Modularity Theorem, bewiesen durch Wiles 1995 für semistabile E, allgemein
        durch Breuil-Conrad-Diamond-Taylor 2001):

            a_p(f_E) = a_p(E) = p + 1 - #E(𝔽_p)  für fast alle p

        Die L-Funktionen stimmen überein: L(ρ_{E,ℓ}, s) = L(f_E, s).

        Diese Klasse implementiert numerische Checks:
        1. L(ρ, s) als Euler-Produkt ∏_p L_p(ρ, s)
        2. Funktionalgleichung Λ(ρ, s) = ε · Λ(ρ̌, 1-s)
        3. Ramanujan-Vermutung |a_p| ≤ 2√p

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self) -> None:
        """
        @brief Initialisiert die Langlands-Korrespondenz-Klasse.
        @lastModified 2026-03-11
        """
        pass

    def automorphic_l_function(
        self,
        rep: GaloisRepresentation,
        s: complex,
        num_primes: int = 50
    ) -> complex:
        """
        @brief Berechnet L(ρ, s) = ∏_p L_p(ρ, s) als Euler-Produkt.
        @description
            Das Euler-Produkt L(ρ, s) = ∏_p det(1 - Frob_p · p^{-s})^{-1}
            wird über die ersten num_primes Primzahlen approximiert.
            Für Re(s) >> 1 konvergiert das Produkt absolut.

            Die Verbindung zum Langlands-Programm: L(ρ, s) soll mit einer
            automorphen L-Funktion L(π, s) übereinstimmen, wobei π eine
            automorphe Darstellung von GL_n(𝔸_ℚ) ist.
        @param rep        Galois-Darstellung
        @param s          Komplexes Argument (Re(s) sollte groß sein für gute Konvergenz)
        @param num_primes Anzahl der Primzahlen im Euler-Produkt
        @return Näherungswert für L(ρ, s)
        @lastModified 2026-03-11
        """
        # Liste der ersten num_primes Primzahlen
        primes = _primes_up_to(500)[:num_primes]
        result = complex(1.0)
        for p in primes:
            # Lokalen Euler-Faktor multiplizieren
            factor = rep.l_function_euler_factor(p, s)
            if not cmath.isfinite(factor):
                continue
            result *= factor
        return result

    def check_functional_equation(
        self,
        rep: GaloisRepresentation,
        s: complex
    ) -> Dict[str, complex]:
        """
        @brief Prüft numerisch die Funktionalgleichung Λ(ρ, s) ≈ ε · Λ(ρ̌, 1-s).
        @description
            Die vollständige L-Funktion eines Tate-Moduls ist:

                Λ(E, s) = (√N / 2π)^s · Γ(s) · L(E, s)

            und erfüllt die Funktionalgleichung

                Λ(E, s) = ε_E · Λ(E, 2-s)

            mit ε_E ∈ {+1, -1} (Vorzeichen der Funktionalgleichung).

            Diese Methode berechnet Λ(ρ, s) und Λ(ρ, 1-s) numerisch und
            gibt beide Werte sowie den geschätzten ε-Faktor zurück.
        @param rep Galois-Darstellung (muss TateModuleRepresentation sein für genaue Ergebnisse)
        @param s   Komplexes Argument
        @return Dict mit 'lambda_s', 'lambda_1ms', 'epsilon_estimate'
        @lastModified 2026-03-11
        """
        # Vollständige L-Funktion: Γ-Faktor × L(ρ, s)
        def gamma_factor(z: complex) -> complex:
            """Gamma-Funktion via Stirling-Näherung für |z| > 5."""
            if abs(z) < 0.1:
                return float('inf') + 0j
            # Stirling-Näherung: Γ(z) ≈ √(2π/z) · (z/e)^z
            return cmath.sqrt(2 * cmath.pi / z) * (z / cmath.e) ** z

        # L-Werte berechnen
        l_at_s = self.automorphic_l_function(rep, s, num_primes=30)
        l_at_1ms = self.automorphic_l_function(rep, 1 - s, num_primes=30)

        # Gamma-Faktoren
        gamma_s = gamma_factor(s)
        gamma_1ms = gamma_factor(1 - s)

        # Vollständige L-Funktionen (vereinfacht, ohne Leiter-Normierung)
        lambda_s = gamma_s * l_at_s
        lambda_1ms = gamma_1ms * l_at_1ms

        # Epsilon-Schätzung: ε ≈ Λ(s) / Λ(1-s) (nur sinnvoll wenn Λ(1-s) ≠ 0)
        if abs(lambda_1ms) > 1e-10:
            epsilon_estimate = lambda_s / lambda_1ms
        else:
            epsilon_estimate = complex(float('nan'))

        return {
            "lambda_s": lambda_s,
            "lambda_1ms": lambda_1ms,
            "epsilon_estimate": epsilon_estimate,
        }

    def check_ramanujan_conjecture(
        self,
        rep: GaloisRepresentation,
        prime_list: List[int]
    ) -> Dict[str, object]:
        """
        @brief Prüft die Ramanujan-Vermutung |a_p| ≤ 2√p für alle p in prime_list.
        @description
            Die **Ramanujan-Vermutung** (für elliptische Kurven bewiesen durch Hasse,
            allgemein durch Deligne als "Weil II"):
            Für den Tate-Modul einer elliptischen Kurve E/ℚ gilt

                |a_p| ≤ 2√p  für alle p bei guter Reduktion.

            Äquivalent: Die Frobenius-Eigenwerte α, β erfüllen |α| = |β| = √p,
            d.h. sie liegen auf dem Kreis mit Radius √p.

            Für allgemeine Galois-Darstellungen aus Modulformen der Gewicht k
            lautet die Vermutung |a_p| ≤ 2·p^{(k-1)/2}.
        @param rep        Galois-Darstellung (2-dimensional)
        @param prime_list Liste von zu prüfenden Primzahlen
        @return Dict mit 'all_satisfied', 'violations', 'max_ratio'
        @lastModified 2026-03-11
        """
        violations = []
        ratios = []

        for p in prime_list:
            if not rep.is_unramified(p):
                continue
            a_p = rep.trace_frobenius(p)
            abs_ap = abs(a_p)
            bound = 2 * math.sqrt(p)
            ratio = abs_ap / bound if bound > 0 else float('inf')
            ratios.append(ratio)

            if abs_ap > bound + 1e-9:
                violations.append({
                    "p": p,
                    "a_p": a_p,
                    "|a_p|": abs_ap,
                    "2*sqrt(p)": bound,
                })

        return {
            "all_satisfied": len(violations) == 0,
            "violations": violations,
            "max_ratio": max(ratios) if ratios else 0.0,
            "num_checked": len(ratios),
        }


# ===========================================================================
# 7. p-ADISCHE DARSTELLUNGSMATRIX
# ===========================================================================

def p_adic_representation_matrix(
    n: int,
    p: int,
    precision: int = 3
) -> np.ndarray:
    """
    @brief Erzeugt eine zufällige GL_n(ℤ/p^k)-Matrix als Demo einer p-adischen Darstellung.
    @description
        Eine p-adische Darstellung ρ: G → GL_n(ℤ_p) kann via Reduktion modulo p^k
        auf GL_n(ℤ/p^k ℤ) reduziert werden. Diese Funktion erzeugt eine invertierbare
        Matrix über ℤ/p^k als Beispiel einer solchen reduzierten Darstellung.

        Die Matrix wird so konstruiert, dass sie modulo p invertierbar ist (det ≢ 0 mod p),
        was die Bedingung für Elemente von GL_n(ℤ/p^k) ist.

        Anwendung: In der Iwasawa-Theorie werden Familien von mod-p^k-Darstellungen
        studiert, die einen projekiven Limes bilden.
    @param n         Dimension der Matrix (n×n)
    @param p         Primzahl p
    @param precision Genauigkeit k (Koeffizienten in ℤ/p^k)
    @return numpy-Matrix A ∈ GL_n(ℤ/p^k), invertierbar modulo p
    @raises ValueError wenn n < 1, p keine Primzahl, oder precision < 1
    @lastModified 2026-03-11
    """
    if n < 1:
        raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")
    if not _is_prime(p):
        raise ValueError(f"p muss eine Primzahl sein, erhalten: {p}")
    if precision < 1:
        raise ValueError(f"precision muss ≥ 1 sein, erhalten: {precision}")

    # Modulus: p^k
    modulus = p ** precision

    # Konstruiere invertierbare Matrix: starte mit Einheitsmatrix + kleine Perturbationen
    rng = np.random.default_rng(seed=42)

    max_attempts = 100
    for _ in range(max_attempts):
        # Zufällige Matrix mit Einträgen in {0, ..., p^k - 1}
        A = rng.integers(0, modulus, size=(n, n))
        # Diagonale auf 1 setzen für bessere Invertierbarkeit
        for i in range(n):
            A[i, i] = 1

        # Prüfe ob det(A) mod p invertierbar ist (also ≢ 0 mod p)
        A_mod_p = A % p
        det_mod_p = int(round(np.linalg.det(A_mod_p.astype(float)))) % p
        if det_mod_p != 0:
            return A % modulus

    # Fallback: Einheitsmatrix (immer invertierbar)
    return np.eye(n, dtype=int)


# ===========================================================================
# 8. GALOIS-GRUPPEN-ORDNUNG
# ===========================================================================

def galois_group_order(polynomial: str) -> int:
    """
    @brief Berechnet die Ordnung der Galois-Gruppe eines irreduziblen Polynoms über ℚ.
    @description
        Die Galois-Gruppe Gal(K/ℚ) eines irreduziblen Polynoms f ∈ ℚ[x]
        ist eine transitive Untergruppe der symmetrischen Gruppe S_n (n = deg f).

        Bekannte Fälle:
        - deg = 1: |Gal| = 1
        - deg = 2: |Gal| = 2 (immer ℤ/2ℤ)
        - deg = 3 irreduzibel: |Gal| = 3 (zyklisch) oder 6 (S₃)
        - deg = 4 irreduzibel: |Gal| ∈ {4, 8, 12, 24}
        - x^n - p (p prim): Gal = ℤ/nℤ ⋊ (ℤ/nℤ)* (Kreisteilungstheorie)

        Diese Funktion nutzt SymPy's galois_group()-Funktion für exakte Berechnung.
        Fallback: Schätzung über Grad des Zerfällungskörpers.

        **Beispiele** (für Polynome über ℚ):
        - x² - 2:  Gal ≅ ℤ/2ℤ,   |Gal| = 2
        - x³ - 2:  Gal ≅ S₃,      |Gal| = 6
        - x⁴ - 2:  Gal ≅ D₄,      |Gal| = 8
        - x⁵ - x - 1: Gal ≅ S₅ (nach Tschirnhaus), |Gal| = 120

    @param polynomial String-Repräsentation des Polynoms in x, z.B. "x**4 - 2"
    @return Ordnung der Galois-Gruppe
    @raises ValueError wenn Polynom nicht geparst werden kann
    @lastModified 2026-03-11
    """
    from sympy import parse_expr, Symbol
    from sympy.polys.numberfields import field_isomorphism
    from sympy import minimal_polynomial, rootof, Poly

    x = symbols('x')

    try:
        # Polynom parsen
        expr = parse_expr(polynomial, local_dict={'x': x})
        poly_obj = Poly(expr, x, domain='QQ')

        # SymPy galois_group() für exakte Berechnung verwenden
        try:
            G, is_alt = galois_group(poly_obj)
            return G.order()
        except (NotImplementedError, AttributeError, TypeError):
            pass

        # Fallback: Schätzung für bekannte Muster
        deg = poly_obj.degree()
        coeffs = poly_obj.all_coeffs()

        if deg == 1:
            return 1
        elif deg == 2:
            # x^2 + bx + c: Galois-Gruppe ist ℤ/2ℤ (irreduzibel über ℚ)
            return 2
        elif deg == 3:
            # x^3 - p: Gal = S₃ wenn nicht reine Kubikwurzel eines Quadrats
            # Diskriminante Δ: Gal = A₃ = ℤ/3ℤ wenn Δ Quadrat, sonst S₃
            a, b, c, d = coeffs[0], coeffs[1], coeffs[2], coeffs[3]
            discriminant = (18*a*b*c*d - 4*b**3*d + b**2*c**2
                           - 4*a*c**3 - 27*a**2*d**2)
            # Ist Δ ein perfektes Quadrat? Falls ja: ℤ/3ℤ (|Gal|=3), sonst S₃ (|Gal|=6)
            from sympy import sqrt as sym_sqrt2, Rational
            disc_val = float(discriminant)
            sqrt_disc = disc_val**0.5 if disc_val >= 0 else None
            if sqrt_disc is not None and abs(sqrt_disc - round(sqrt_disc)) < 1e-9:
                return 3
            return 6
        elif deg == 4:
            # x^4 - n: generisch D₄ mit Ordnung 8
            # Genauere Analyse nötig, aber 8 ist häufigster Fall für x^n-p
            return 8
        else:
            # Allgemeine Schätzung: S_n hat Ordnung n!
            return math.factorial(deg)

    except Exception as e:
        raise ValueError(f"Konnte Polynom nicht analysieren: {polynomial}. Fehler: {e}")


# ===========================================================================
# 9. ARTIN-LEITER
# ===========================================================================

def artin_conductor(
    rep: GaloisRepresentation,
    prime_list: List[int]
) -> int:
    """
    @brief Schätzt den Artin-Leiter einer Galois-Darstellung.
    @description
        Der **Artin-Leiter** (conductor) N(ρ) einer Galois-Darstellung
        ρ: Gal(Q̄/Q) → GL_n(ℤ_ℓ) ist definiert als

            N(ρ) = ∏_p p^{f(ρ, p)}

        wobei der Exponent f(ρ, p) die Verzweigung bei p misst:

            f(ρ, p) = dim(ρ^{I_p}) + Swan(ρ, p) ≥ 0

        mit I_p = Trägheitsgruppe, ρ^{I_p} = unverzweigter Anteil,
        Swan(ρ, p) = wilder Verzweigungsterm.

        Für eine unverzweigte Darstellung bei p ist f(ρ, p) = 0.
        Für zahme Verzweigung (p ≠ ℓ): f(ρ, p) ≥ 1.

        Diese Funktion schätzt den Leiter via: Prüfe ob is_unramified(p)
        für jedes p in prime_list, und multipliziere verzweigte p.

        Für TateModuleRepresentation wird der tatsächliche Führer verwendet.
    @param rep        Galois-Darstellung
    @param prime_list Liste der zu prüfenden Primzahlen
    @return Geschätzte Ordnung des Artin-Leiters (Produkt der verzweigten Primzahlen)
    @lastModified 2026-03-11
    """
    # Für TateModule: direkter Zugriff auf den Führer
    if isinstance(rep, TateModuleRepresentation):
        return rep.conductor

    # Für DirichletCharacter: Führer ist der Modulus
    if isinstance(rep, DirichletCharacterRepresentation):
        return rep.conductor

    # Allgemeiner Fall: Produkt der verzweigten Primzahlen
    conductor_estimate = 1
    for p in prime_list:
        if not rep.is_unramified(p):
            # Verzweigte Primzahl: zum Leiter multiplizieren
            conductor_estimate *= p

    return conductor_estimate if conductor_estimate > 1 else 1


# ===========================================================================
# 10. WEIL-DELIGNE-DARSTELLUNG (DEMO)
# ===========================================================================

def weil_deligne_representation_demo(p: int) -> Dict[str, object]:
    """
    @brief Demo der Weil-Deligne-Darstellung für eine p-adische elliptische Kurve.
    @description
        Eine **Weil-Deligne-Darstellung** ist ein Paar (r, N) bestehend aus:
        - r: W_p → GL_n(ℂ) — Darstellung der Weil-Gruppe W_p von ℚ_p
        - N: nilpotente Matrix (Monodromie-Operator) mit r(σ)·N·r(σ)^{-1} = ||σ||·N

        Hier ist ||σ|| = p^{-v(σ)} die p-adische Norm (Weil-Gruppe-Charakter).

        Für eine elliptische Kurve E/ℚ_p gibt es drei Typen von lokalen Darstellungen
        (Klassifikation nach Kodaira-Néron-Typ):

        1. **Gute Reduktion** (E hat gute Reduktion bei p):
           - r(Frob_p) = Frobenius-Matrix mit Eigenwerten α, β (|α|=|β|=√p)
           - N = 0 (keine Monodromie)

        2. **Multiplikative Reduktion** (Tate-Kurve E_q: y²+xy=x³-...):
           - r = Frobenius-Twist des zyklotomischen Charakters
           - N = [[0,1],[0,0]] (unipotenter Monodromieoperator)

        3. **Additive Reduktion** (schlimmste Verzweigung):
           - Komplizierter; wird hier nicht vollständig implementiert.

        Literatur: J.-P. Serre, "Local Fields" (1979);
                   B. Conrad, "Ramification and finiteness" (Seminaire Bourbaki).
    @param p Primzahl (Charakteristik der Residuenkörper-Reduktion)
    @return Dict mit Feldern 'prime', 'frobenius_matrix', 'monodromy_matrix', 'type', 'description'
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"p muss eine Primzahl sein, erhalten: {p}")

    # Beispiel: Kurve y² = x³ - x (Conductor 32, Mordell-Weil-Rang 0)
    # a_p-Werte für kleine Primzahlen:
    a_p_table = {2: 0, 3: 0, 5: -2, 7: 0, 11: 0, 13: 6, 17: 2, 19: 0, 23: 0, 29: 6}

    if p == 2:
        # Additive Reduktion bei p=2 (Führer 32 = 2^5)
        frobenius = np.array([[0, -1], [1, 0]], dtype=complex)  # Demo-Frobenius
        monodromy = np.array([[0, 1], [0, 0]], dtype=complex)   # Nilpotent
        reduction_type = "additive"
        description = (
            f"p={p}: Additive Reduktion. Weil-Deligne-Darstellung hat "
            f"nicht-triviale wilde Verzweigung. N ≠ 0."
        )
    elif a_p_table.get(p, None) is not None:
        a_p = a_p_table[p]
        disc = complex(a_p**2 - 4*p)
        sqrt_disc = cmath.sqrt(disc)
        alpha = (a_p + sqrt_disc) / 2
        beta = (a_p - sqrt_disc) / 2

        # Frobenius-Matrix in diagonaler Form
        frobenius = np.array([[alpha, 0], [0, beta]], dtype=complex)
        monodromy = np.zeros((2, 2), dtype=complex)  # Keine Monodromie bei guter Red.
        reduction_type = "good"
        description = (
            f"p={p}: Gute Reduktion. a_p={a_p}. "
            f"Frobenius-Eigenwerte: α≈{alpha:.4f}, β≈{beta:.4f}. "
            f"N=0 (keine Monodromie)."
        )
    else:
        # Unbekannte Primzahl: generische multiplikative Reduktion
        frobenius = np.array([[p, 0], [0, 1]], dtype=complex)
        monodromy = np.array([[0, 1], [0, 0]], dtype=complex)
        reduction_type = "multiplicative"
        description = (
            f"p={p}: Angenommene multiplikative Reduktion. "
            f"Tate-Kurve-Parametrisierung. N ≠ 0 (Monodromie-Operator)."
        )

    return {
        "prime": p,
        "frobenius_matrix": frobenius,
        "monodromy_matrix": monodromy,
        "type": reduction_type,
        "description": description,
        "nilpotent_check": bool(np.allclose(monodromy @ monodromy, np.zeros((2, 2)))),
    }


# ===========================================================================
# HILFSFUNKTION: a_p-Werte für y² = x³ - x berechnen
# ===========================================================================

def compute_ap_for_y2_x3_minus_x(p: int) -> int:
    """
    @brief Berechnet a_p = p + 1 - #E(𝔽_p) für E: y² = x³ - x.
    @description
        Zählt die affinen Punkte auf E(𝔽_p) + den Punkt im Unendlichen:

            #E(𝔽_p) = 1 + Σ_{x=0}^{p-1} (1 + Legendre(x³-x, p))

        Nach dem Legendre-Symbol-Kalkül gilt für p ≡ 3 (mod 4): a_p = 0.
        Für p ≡ 1 (mod 4): a_p = 2·Re(π) wo π Gaußsche Primzahl mit π·π̄ = p.

        Diese Funktion zählt brute-force für kleine p.
    @param p Primzahl p > 3 (bei p=2 oder p=3 hat E schlechte Reduktion)
    @return a_p = p + 1 - #E(𝔽_p)
    @lastModified 2026-03-11
    """
    if p <= 3:
        return 0  # Schlechte Reduktion oder zu kleines p

    # Quadratische Reste mod p berechnen
    quadratic_residues = set()
    for i in range(p):
        quadratic_residues.add((i * i) % p)

    # Punkte auf E: y² = x³ - x zählen
    count = 1  # Punkt im Unendlichen
    for x in range(p):
        rhs = (pow(x, 3, p) - x) % p
        if rhs == 0:
            count += 1  # Nur der Punkt y=0
        elif rhs in quadratic_residues:
            count += 2  # Zwei Punkte y und -y

    return p + 1 - count
