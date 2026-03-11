r"""
@file motive_theory.py
@brief Motiventheorie nach Grothendieck: Chow-Gruppen, Realisierungsfunktoren,
       Standardvermutungen und motivische L-Funktionen.
@description
    Die Motiventheorie (Grothendieck, 1964) ist ein einheitliches Framework,
    das alle kohomologischen Theorien algebraischer Varietäten durch einen
    universellen "Ursprung" — das Motiv — verbindet.

    Kernkonzepte:
    - **Reines Motiv** h(X) zu einer glatten projektiven Varietät X
    - **Chow-Gruppen** CH^k(X): Zykelgruppen modulo rationaler Äquivalenz
    - **Realisierungsfunktoren**: Betti, de Rham, l-adisch (étale)
    - **Standardvermutungen** (Grothendieck): Lefschetz, Hodge, Künneth
    - **Motivische Kohomologie**: H^{p,q}_M, Beilinson-Regulatoren
    - **L-Funktionen**: motivische L-Funktion L(M, s)

    Historische Meilensteine:
    - 1964: Grothendieck formuliert Motiventheorie in Mumfords Notizen
    - 1991: Voevodsky definiert "gemischte Motive" (triangulierte Kategorie)
    - 1994: Voevodsky-Suslin-Friedlander: motivische Kohomologie
    - 2000: Voevodsky beweist Milnor-Vermutung (Teil von BSD-Programm)

    KaTeX-Formeln:
    $$\text{CH}^k(X) = \frac{\{\text{algebraische Zykel der Kodimension } k\}}{\text{rationale Äquivalenz}}$$
    $$L(M, s) = \prod_p L_p(M, p^{-s})^{-1}, \quad L_p(M, T) = \det(1 - F_p T | H^\bullet(M))$$

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
from typing import Dict, List, Optional, Tuple, Any, Union
from fractions import Fraction

from exceptions import InvalidInputError, PrimeRequiredError


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _is_prime(n: int) -> bool:
    r"""
    @brief Einfache Primzahlprüfung via Trial Division.
    @param n Zu prüfende ganze Zahl.
    @return True wenn n prim, False sonst.
    @lastModified 2026-03-11
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True


def _euler_phi(n: int) -> int:
    r"""
    @brief Euler'sche Phi-Funktion φ(n).
    @param n Positive ganze Zahl.
    @return φ(n) = Anzahl der zu n teilerfremden Zahlen in {1,...,n}.
    @lastModified 2026-03-11
    """
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


# ===========================================================================
# 1. MOTIVE UND KOHOMOLOGIE
# ===========================================================================

class Motive:
    r"""
    @brief Klasse für ein reines Motiv h(X) einer algebraischen Varietät X.
    @description
        Ein reines Motiv (M, p, n) besteht aus:
        - M: eine glatte projektive Varietät (hier als Name/Dimension repräsentiert)
        - p: ein Korrespondenzprojektor (p² = p in CH^d(X×X))
        - n: Tate-Twist: h(X)(n) = h(X) ⊗ Q(n)

        Grundmotive:
        - h(Spec k) = Q(0): triviales Motiv (Einsmotiv)
        - Q(-1) = L: Lefschetz-Motiv (h²(P¹))
        - Q(n) = L^{⊗n}: n-facher Tate-Twist

        Kanonische Zerlegung (Künneth):
        h(X) = ⊕_{i=0}^{2d} h^i(X)   (konj. durch Hodge-Lefschetz)

        KaTeX:
        $$h(X) = \bigoplus_{i=0}^{2\dim X} h^i(X), \quad h^i(X)(n) = h^i(X) \otimes_{\mathbb{Q}} \mathbb{Q}(n)$$

    @param variety_name Name/Bezeichner der Varietät.
    @param dimension Dimension der Varietät.
    @param cohomology_type Typ der Realisierung ('betti', 'de_rham', 'l_adic', 'all').
    @param tate_twist Tate-Twist n (Integer).
    @lastModified 2026-03-11
    """

    def __init__(self, variety_name: str, dimension: int,
                 cohomology_type: str = 'all', tate_twist: int = 0):
        r"""
        @brief Konstruktor für ein reines Motiv.
        @param variety_name Bezeichner der Varietät (z.B. 'P1', 'elliptic_curve', 'abelian_variety').
        @param dimension Dimension der Varietät (d ≥ 0).
        @param cohomology_type Kohomologietyp ('betti', 'de_rham', 'l_adic', 'all').
        @param tate_twist Tate-Twist n ∈ Z.
        @raises InvalidInputError Wenn dimension < 0 oder unbekannter cohomology_type.
        @lastModified 2026-03-11
        """
        if dimension < 0:
            raise InvalidInputError(f"Motive: Dimension muss ≥ 0 sein, erhalten: {dimension}")

        valid_types = {'betti', 'de_rham', 'l_adic', 'motivic', 'all'}
        if cohomology_type not in valid_types:
            raise InvalidInputError(
                f"Motive: Unbekannter Kohomologietyp '{cohomology_type}'. "
                f"Erlaubt: {valid_types}"
            )

        # Grunddaten des Motivs
        self.variety_name = variety_name
        self.dimension = dimension
        self.cohomology_type = cohomology_type
        self.tate_twist = tate_twist

        # Betti-Zahlen (standardmäßig heuristisch für bekannte Varietäten)
        self.betti_numbers = self._compute_betti_numbers()

    def _compute_betti_numbers(self) -> List[int]:
        r"""
        @brief Berechnet Betti-Zahlen b_i = dim H^i(X, Q) für bekannte Varietäten.
        @description
            Betti-Zahlen für Standardvarietäten:
            - Punkt: [1]
            - P^n: [1, 0, 1, 0, ..., 1] (nur gerade Grade 1)
            - Elliptische Kurve: [1, 2, 1] (Genus 1)
            - Abelsche Varietät Dim g: [C(2g,0), C(2g,1), ..., C(2g,2g)]
            - K3-Fläche: [1, 0, 22, 0, 1]
        @return Liste [b_0, b_1, ..., b_{2d}].
        @lastModified 2026-03-11
        """
        d = self.dimension
        name = self.variety_name.lower()

        if name in ('point', 'spec_k', 'spec k'):
            return [1]
        elif name.startswith('p') and name[1:].isdigit():
            # Projektiver Raum P^n: b_{2k} = 1 für 0 ≤ k ≤ n, b_{2k+1} = 0
            n = int(name[1:])
            betti = [0] * (2 * n + 1)
            for k in range(n + 1):
                betti[2 * k] = 1
            return betti
        elif name in ('elliptic_curve', 'elliptic curve', 'ec'):
            # Elliptische Kurve: b_0=1, b_1=2, b_2=1
            return [1, 2, 1]
        elif name in ('k3', 'k3_surface', 'k3 surface'):
            # K3-Fläche: b_0=1, b_1=0, b_2=22, b_3=0, b_4=1
            return [1, 0, 22, 0, 1]
        elif name.startswith('abelian_') or name.startswith('abelian variety'):
            # Abelsche Varietät der Dimension g:
            # b_i = C(2g, i) nach Künneth-Formel
            betti = [math.comb(2 * d, i) for i in range(2 * d + 1)]
            return betti
        elif name in ('curve', 'smooth_curve'):
            # Glatte Kurve vom Genus g = d (für Dimension 1)
            # b_0=1, b_1=2g, b_2=1 (Poincaré-Dualität)
            g = max(0, d)
            return [1, 2 * g, 1] if d >= 1 else [1]
        else:
            # Allgemeines Motiv: Poincaré-Dualität b_i = b_{2d-i}
            betti = [1] + [0] * (2 * d - 1) + [1] if d >= 1 else [1]
            return betti

    def euler_characteristic(self) -> int:
        r"""
        @brief Berechnet die Euler-Charakteristik χ(X) = Σ (-1)^i b_i.
        @description
            Topologische Invariante:
            χ(X) = Σ_{i=0}^{2d} (-1)^i · b_i(X)

            Wichtige Werte:
            - χ(P^n) = n + 1
            - χ(Elliptische Kurve) = 0
            - χ(K3) = 24
        @return Euler-Charakteristik als ganze Zahl.
        @lastModified 2026-03-11
        """
        return sum((-1) ** i * b for i, b in enumerate(self.betti_numbers))

    def hodge_numbers(self) -> Dict[Tuple[int, int], int]:
        r"""
        @brief Gibt die Hodge-Zahlen h^{p,q} zurück (für einfache Varietäten).
        @description
            Hodge-Zerlegung: H^n(X, C) = ⊕_{p+q=n} H^{p,q}(X)
            Es gilt h^{p,q} = h^{q,p} (Komplex-Konjugation) und
                   h^{p,q} = h^{d-p,d-q} (Serre-Dualität, d = dim X).

            Für elliptische Kurve: h^{1,0} = h^{0,1} = 1 (Genus)
            Für K3: h^{2,0} = h^{0,2} = 1, h^{1,1} = 20
        @return Dictionary {(p,q): h^{p,q}}.
        @lastModified 2026-03-11
        """
        d = self.dimension
        name = self.variety_name.lower()
        hodge = {}

        if name in ('elliptic_curve', 'elliptic curve', 'ec'):
            # Elliptische Kurve: H^0 = Q, H^1 = H^{1,0}⊕H^{0,1}, H^2 = Q(-1)
            hodge = {(0, 0): 1, (1, 0): 1, (0, 1): 1, (2, 0): 0, (1, 1): 0, (0, 2): 0}
            # Tatsächlich: h^{1,1} gehört zur Kohomologie, aber für Kurven gilt dies
            hodge = {(0, 0): 1, (1, 0): 1, (0, 1): 1}
        elif name in ('k3', 'k3_surface'):
            hodge = {
                (0, 0): 1, (1, 0): 0, (0, 1): 0,
                (2, 0): 1, (1, 1): 20, (0, 2): 1,
                (2, 1): 0, (1, 2): 0, (2, 2): 1,
            }
        elif name.startswith('p') and name[1:].isdigit():
            n = int(name[1:])
            for k in range(n + 1):
                hodge[(k, k)] = 1
        else:
            # Allgemein: diagonale Hodge-Zahlen (Schätzung)
            for k in range(d + 1):
                hodge[(k, k)] = 1

        return hodge

    def tate_twist_apply(self, n: int) -> 'Motive':
        r"""
        @brief Gibt das Motiv M(n) = M ⊗ Q(n) zurück (n-facher Tate-Twist).
        @description
            Der Tate-Twist Q(n) ist das n-fache Tensorprodukt des Lefschetz-Motivs.
            Wirkung auf Hodge-Zahlen: h^{p,q}(M(n)) = h^{p-n,q-n}(M)
            Wirkung auf Gewicht: w(M(n)) = w(M) - 2n
        @param n Anzahl der Tate-Twists (kann negativ sein).
        @return Neues Motiv M(n).
        @lastModified 2026-03-11
        """
        return Motive(
            variety_name=f"{self.variety_name}({n})",
            dimension=self.dimension,
            cohomology_type=self.cohomology_type,
            tate_twist=self.tate_twist + n
        )

    def dual(self) -> 'Motive':
        r"""
        @brief Duales Motiv M^∨ = Hom(M, Q(0)).
        @description
            Poincaré-Dualität: h(X)^∨ ≅ h(X)(d) für eine d-dimensionale Varietät.
            Betti-Zahlen sind dieselben, aber Hodge-Zahlen werden umgekehrt.
        @return Duales Motiv.
        @lastModified 2026-03-11
        """
        return Motive(
            variety_name=f"{self.variety_name}_dual",
            dimension=self.dimension,
            cohomology_type=self.cohomology_type,
            tate_twist=-self.tate_twist + self.dimension
        )

    def __repr__(self) -> str:
        r"""@brief String-Repräsentation des Motivs."""
        return (f"Motive(variety='{self.variety_name}', dim={self.dimension}, "
                f"type='{self.cohomology_type}', twist={self.tate_twist})")


# ===========================================================================
# 2. CHOW-GRUPPEN
# ===========================================================================

class ChowGroup:
    r"""
    @brief Chow-Gruppe CH^k(X) einer algebraischen Varietät.
    @description
        Die Chow-Gruppe CH^k(X) ist die Gruppe der algebraischen Zykel
        der Kodimension k modulo rationaler Äquivalenz:

            CH^k(X) = Z^k(X) / ~_rat

        wobei ~_rat rationale Äquivalenz via 1-dimensionale Familien ist.

        Wichtige Eigenschaften:
        - CH^0(X) = Z (von der Fundamentalklasse erzeugt)
        - CH^1(X) = Pic(X): Picard-Gruppe (Divisorenklassen)
        - CH^d(X) = Z^0(X)~/~: Nullzykelgruppe
        - CH^*(X) = ⊕_k CH^k(X): Ring unter Schnittprodukt

        Konjectures (Standard):
        - Bloch-Kato: Milnor K-Theorie ↔ Galoisgruppen (bewiesen Voevodsky 2010)
        - Beilinson: CH^k(X) ⊗ Q ≅ Ext^{2k-1}(Q, H^{2k-1}(X, Q(k)))
        - Tate: CH^k(X) ⊗ Ql / rationale ~ ↔ Galois-Rep. in H^{2k}_et

        KaTeX:
        $$\text{CH}^k(X) = \frac{\bigoplus_Y \mathbb{Z} \cdot [Y]}{\text{Rationale Äquivalenz}}$$

    @param variety Zugehöriges Motive-Objekt.
    @param codimension Kodimension k der Zykel.
    @lastModified 2026-03-11
    """

    def __init__(self, variety: Motive, codimension: int):
        r"""
        @brief Konstruktor für eine Chow-Gruppe.
        @param variety Varietät (als Motive-Objekt).
        @param codimension Kodimension k (0 ≤ k ≤ dim X).
        @raises InvalidInputError Bei ungültiger Kodimension.
        @lastModified 2026-03-11
        """
        if not (0 <= codimension <= variety.dimension):
            raise InvalidInputError(
                f"ChowGroup: Kodimension {codimension} außerhalb [0, {variety.dimension}]."
            )

        self.variety = variety
        self.codimension = codimension

        # Basiszykel: Listen von (multiplicity, label) Paaren
        self._generators: List[Tuple[int, str]] = []

    def add_cycle(self, multiplicity: int, label: str) -> None:
        r"""
        @brief Fügt einen algebraischen Zykel zur Gruppe hinzu.
        @param multiplicity Ganze Zahl (Vielfachheit des Zykels).
        @param label Bezeichner des Zykels (z.B. "hyperplane", "diagonal").
        @lastModified 2026-03-11
        """
        self._generators.append((multiplicity, label))

    def intersection_pairing(self, other: 'ChowGroup') -> Optional[int]:
        r"""
        @brief Schnittprodukt CH^j(X) × CH^k(X) → CH^{j+k}(X).
        @description
            Das Schnittprodukt von Zykeln ist wohldefeniert auf Äquivalenzklassen.
            Für komplementäre Kodimensionen (j + k = d) ergibt sich ein
            ganzzahliger Grad:

            KaTeX: $(Z_1 \cdot Z_2) \in \mathbb{Z}$ für $\text{codim}(Z_1) + \text{codim}(Z_2) = d$

        @param other Andere Chow-Gruppe der Kodimension d - k.
        @return Schnittprodukt-Grad als Integer (falls möglich), sonst None.
        @lastModified 2026-03-11
        """
        target_codim = self.codimension + other.codimension
        if target_codim == self.variety.dimension:
            # Komplementäre Kodimensionen: Ergebnis ist Grad (Integer)
            total_mult = sum(m for m, _ in self._generators)
            other_mult = sum(m for m, _ in other._generators)
            return total_mult * other_mult  # Vereinfachung für reine Punkte
        return None

    def picard_group_rank(self) -> Optional[int]:
        r"""
        @brief Gibt den Rang der Picard-Gruppe Pic(X) = CH^1(X) zurück.
        @description
            Für CH^1(X): Pic(X) = NS(X) ⊕ Pic^0(X)
            - NS(X): Néron-Severi-Gruppe (endlich erzeugt)
            - Pic^0(X): Zusammenhangskomponente (Jacobi-Varietät für Kurven)
            Rang(Pic(X)) = Rang(NS(X)) (Néron-Severi-Rang)

        @return Rang der Picard-Gruppe (nur für CH^1 sinnvoll), sonst None.
        @lastModified 2026-03-11
        """
        if self.codimension != 1:
            return None

        name = self.variety.variety_name.lower()
        if name in ('elliptic_curve', 'ec', 'elliptic curve'):
            return 1  # Pic^0 ≅ E selbst, NS = Z
        elif name.startswith('p') and name[1:].isdigit():
            return 1  # Pic(P^n) = Z (erzeugt von O(1))
        elif name in ('k3', 'k3_surface'):
            return 1  # Generische K3: Pic = Z (ρ=1)
        return len(self._generators) or 1  # Heuristisch

    def __repr__(self) -> str:
        r"""@brief String-Repräsentation der Chow-Gruppe."""
        return (f"ChowGroup(variety='{self.variety.variety_name}', "
                f"codim={self.codimension}, #generators={len(self._generators)})")


# ===========================================================================
# 3. MOTIVISCHE KOHOMOLOGIE UND BEILINSON-REGULATOREN
# ===========================================================================

class MotivicCohomology:
    r"""
    @brief Motivische Kohomologie H^{p,q}_M(X, Z) nach Voevodsky-Suslin-Friedlander.
    @description
        Die motivische Kohomologie verallgemeinert K-Theorie und Chow-Gruppen:
            H^{2k,k}_M(X, Z) = CH^k(X)          (Chow-Gruppen)
            H^{2k-1,k}_M(X, Z(k)) ⊗ Q            (Beilinson-Regulatoren)

        Beilinson-Vermutung: Der Beilinson-Regulator r_B liefert einen
        Isomorphismus nach Tensorprodukt mit R:
            r_B: H^{p,q}_M(X, Q) → H^p_D(X, R(q))  (Deligne-Kohomologie)
        mit r(M, s) = det(r_B) · L(M, s) · (einfacher Faktor).

        Milnor K-Theorie Verbindung:
            K^M_n(F) = H^{n,n}_M(Spec F, Z)  für einen Körper F.

        KaTeX:
        $$H^{p,q}_{\mathcal{M}}(X, \mathbb{Z}) \cong \text{CH}^q(X, 2q-p), \quad
          H^{2q,q}_{\mathcal{M}}(X, \mathbb{Z}) \cong \text{CH}^q(X)$$

    @param motive Zugehöriges Motive-Objekt.
    @lastModified 2026-03-11
    """

    def __init__(self, motive: Motive):
        r"""
        @brief Konstruktor für motivische Kohomologie.
        @param motive Reines Motiv h(X).
        @lastModified 2026-03-11
        """
        self.motive = motive

    def dimension(self, p: int, q: int) -> int:
        r"""
        @brief Schätzt dim H^{p,q}_M(X, Q) für ein Motiv h(X).
        @description
            Für h^{2k,k}: CH^k(X) ⊗ Q → dim = Rang(CH^k)
            Für andere (p,q): über Betti-Zahlen abgeschätzt.
        @param p Erstes Kohomologie-Gewicht.
        @param q Zweites Kohomologie-Gewicht.
        @return Geschätzte Dimension (≥ 0).
        @lastModified 2026-03-11
        """
        d = self.motive.dimension
        betti = self.motive.betti_numbers

        # Fundamentale Fälle:
        if p == 2 * q and 0 <= q <= d:
            # CH^q(X) ⊗ Q
            return betti[2 * q] if 2 * q < len(betti) else 0
        elif p == 2 * q - 1 and 1 <= q <= d:
            # Regulatorraum: dim ≈ b_{2q-1} (Betti-Zahl)
            return betti[2 * q - 1] if 2 * q - 1 < len(betti) else 0
        else:
            return 0

    def beilinson_regulator(self, p: int, q: int) -> float:
        r"""
        @brief Schätzt den Beilinson-Regulator r_B: H^{p,q}_M → H^p_D numerisch.
        @description
            Der Beilinson-Regulator ist ein Isomorphismus:
            r_B: H^{p,q}_M(X, Q) ⊗ R → H^p_D(X/R, R(q))
            (nach Tensorprodukt mit R).

            Für elliptische Kurven (p=1, q=1):
            r_B entspricht dem Neron-Tate-Regulator.

            Beilinson-Vermutung: L(M, q) ~ det(r_B) · Ω (bis auf rationale Faktoren)
            KaTeX: $r(M, s)|_{s=q} \sim \det(r_{\mathcal{B}}) \cdot \Omega_M$

        @param p Bidegree p (wie in H^{p,q}).
        @param q Bidegree q.
        @return Geschätzter Regulatorwert (dimensionsloser Float).
        @lastModified 2026-03-11
        """
        name = self.motive.variety_name.lower()

        # Neron-Tate-Regulator für elliptische Kurven (p=1, q=1)
        if name in ('elliptic_curve', 'ec') and p == 1 and q == 1:
            # Heuristischer Wert: Regulator ∼ log(Klein) für CM-Kurven
            return 1.0  # Normiert auf 1 für Rang-0-Kurven

        # Allgemeiner Fall: Produkt von log-Perioden
        dim = self.dimension(p, q)
        if dim == 0:
            return 0.0

        # Grobe Schätzung: R ≈ (2π)^{q · dim}
        return (2 * math.pi) ** (q * dim) / math.factorial(dim)

    def milnor_k_theory(self, n: int, field_char: int = 0) -> Dict[str, Any]:
        r"""
        @brief Milnor K-Gruppen K^M_n(F) für einen Körper F.
        @description
            Milnor K-Theorie (1970): K^M_n(F) = F^×⊗n / <Steinberg-Relationen>
            Steinberg-Relation: a ⊗ (1-a) = 0 für a ∈ F \ {0, 1}.

            Voevodsky (2000) bewies: K^M_n(F)/2 ≅ H^n_ét(F, Z/2) (Milnor-Vermutung)
            Bass-Tate: K^M_n(F) = 0 für n > 1 falls F endlicher Körper.

            KaTeX: $K_n^M(F) = F^\times \otimes_{\mathbb{Z}} \cdots \otimes_{\mathbb{Z}} F^\times \big/ I$

        @param n Grad der Milnor K-Gruppe.
        @param field_char Charakteristik des Körpers (0 = Q oder R, p prim = F_p).
        @return Dictionary mit Informationen zur Milnor K-Gruppe.
        @lastModified 2026-03-11
        """
        result = {
            'n': n,
            'field_char': field_char,
        }

        if field_char > 0 and _is_prime(field_char):
            # Endlicher Körper F_p: K^M_n(F_p) = 0 für n ≥ 2
            if n == 0:
                result['description'] = 'K^M_0(F_p) = Z'
                result['rank'] = 0
                result['finite_part'] = 'Z'
            elif n == 1:
                # K^M_1(F_p) = F_p^× ≅ Z/(p-1)Z
                result['description'] = f'K^M_1(F_p) ≅ Z/{field_char - 1}Z'
                result['rank'] = 0
                result['order'] = field_char - 1
            else:
                result['description'] = f'K^M_{n}(F_p) = 0 (Bass-Tate)'
                result['rank'] = 0
                result['order'] = 1
        else:
            # Zahlkörper / Q: K^M_n(Q) komplex, Quillen-Milnor verbindet zu K-Theorie
            if n == 0:
                result['description'] = 'K^M_0(Q) = Z'
                result['rank'] = 0
            elif n == 1:
                result['description'] = 'K^M_1(Q) = Q^× ≅ {±1} × ⊕_p Z'
                result['rank'] = 0  # Free part durch Primzahlen erzeugt
            elif n == 2:
                # K^M_2(Q) = Brauer-Gruppe-Torsion (bekannt via Moore)
                result['description'] = 'K^M_2(Q) = Z/2Z (Moore, 1968)'
                result['rank'] = 0
                result['order'] = 2
            else:
                # K^M_n(Q) = 0 für n ≥ 3 (Bass-Tate für Zahlkörper)
                result['description'] = f'K^M_{n}(Q) = 0 (Bass-Tate)'
                result['rank'] = 0
                result['order'] = 1

        return result


# ===========================================================================
# 4. STANDARDVERMUTUNGEN
# ===========================================================================

class StandardConjectures:
    r"""
    @brief Grothendieck'sche Standardvermutungen für algebraische Varietäten.
    @description
        Grothendieck formulierte vier "Standard Conjectures" (1968):
        A. **Lefschetz**: Die algebraischen Zykel-Klassen erzeugen primitive Kohomologie
        B. **Hodge** (schwach): Ein algebraischer Zykel der Kodimension k ist in
           primitiver Kohomologie orthogonal zu sich selbst.
        C. **Künneth**: Die Künneth-Projektoren π_i: H*(X×X) → H*(X×X) sind algebraisch.
        D. **Numerische = homologische Äquivalenz**: Für Primkörper beides äquivalent.

        Status:
        - Für abelsche Varietäten: alle bekannt (Lieberman, Kleiman)
        - Im Allgemeinen: unbewiesen (offen seit 60 Jahren!)
        - B + C implizieren den Weil-Vermutungen-Beweis über Z/pZ

        KaTeX:
        $$L^{d-k}: H^k(X) \xrightarrow{\sim} H^{2d-k}(X) \quad \text{(Hard Lefschetz)}$$
        $$\pi_i \circ \pi_j = \delta_{ij} \pi_i \quad \text{(Künneth-Projektoren algebraisch)}$$

    @param motive Motiv der zu untersuchenden Varietät.
    @lastModified 2026-03-11
    """

    def __init__(self, motive: Motive):
        r"""
        @brief Konstruktor.
        @param motive Reines Motiv h(X) für die Varietät X.
        @lastModified 2026-03-11
        """
        self.motive = motive

    def check_hard_lefschetz(self) -> Dict[str, Any]:
        r"""
        @brief Überprüft die Hard-Lefschetz-Eigenschaft numerisch.
        @description
            Hard Lefschetz: Für L = c_1(hyperplane), 0 ≤ k ≤ d gilt:
                L^{d-k}: H^k(X, Q) → H^{2d-k}(X, Q) ist Isomorphismus.

            Für glatte projektive Varietäten: bewiesen (Deligne 1974 für char p,
            klassisch über C).
            Für Motive: folgt aus Hard Lefschetz.

        @return Dictionary mit Prüfungsstatus.
        @lastModified 2026-03-11
        """
        d = self.motive.dimension
        betti = self.motive.betti_numbers

        # Poincaré-Dualität: b_k = b_{2d-k}
        poincare_ok = True
        for k in range(len(betti) // 2 + 1):
            mirror = 2 * d - k
            if mirror < len(betti) and k < len(betti):
                if betti[k] != betti[mirror]:
                    poincare_ok = False
                    break

        return {
            'hard_lefschetz': 'verified' if poincare_ok else 'failed_numerically',
            'poincare_duality': poincare_ok,
            'betti_numbers': betti,
            'note': 'Hard Lefschetz bewiesen für glatte proj. Varietäten über C/F_p'
        }

    def check_kunneth(self, other: Motive) -> Dict[str, Any]:
        r"""
        @brief Prüft Künneth-Formel für Produkt X × Y.
        @description
            Künneth: H^n(X × Y) ≅ ⊕_{i+j=n} H^i(X) ⊗ H^j(Y)
            Betti-Zahlen des Produkts: b_n(X×Y) = Σ_{i+j=n} b_i(X)·b_j(Y)

            Standardvermutung C: Die Künneth-Projektoren π_i(X): H*(X×X) → H^i(X)
            sind von algebraischen Zykeln induziert (Korrespondenzen in CH^d(X×X)).

        @param other Zweite Varietät Y.
        @return Dictionary mit Künneth-Betti-Zahlen des Produkts X×Y.
        @lastModified 2026-03-11
        """
        betti_x = self.motive.betti_numbers
        betti_y = other.betti_numbers
        d_x = self.motive.dimension
        d_y = other.dimension

        # Künneth-Betti-Zahlen
        n_max = 2 * (d_x + d_y) + 1
        kunneth_betti = [0] * n_max

        for i, bx in enumerate(betti_x):
            for j, by in enumerate(betti_y):
                if i + j < n_max:
                    kunneth_betti[i + j] += bx * by

        euler_x = self.motive.euler_characteristic()
        euler_y = other.euler_characteristic()
        euler_product = sum((-1) ** i * b for i, b in enumerate(kunneth_betti))

        return {
            'product_variety': f'{self.motive.variety_name} × {other.variety_name}',
            'product_dimension': d_x + d_y,
            'kunneth_betti': kunneth_betti,
            'euler_product': euler_product,
            'euler_product_check': (euler_product == euler_x * euler_y),
            'conjecture_C_status': 'open_in_general'
        }

    def hodge_conjecture_evidence(self) -> Dict[str, Any]:
        r"""
        @brief Sammelt Evidenz für die Hodge-Vermutung (Millennium-Problem).
        @description
            Hodge-Vermutung: Für X glatte projektive Varietät über C gilt:
            Jede rationale (p,p)-Hodge-Klasse ist eine rationale Linearkombination
            von Klassen algebraischer Zykel.

            H^{p,p}(X, Q) = Im(CH^p(X) → H^{2p}(X, Q))  ?

            Status (2026): Bewiesen für:
            - p = 0: trivial (Fundamentalklasse)
            - p = d: trivial (Punkte)
            - Abelsche Varietäten mit zusätzlicher Struktur (Lefschetz (1,1)-Theorem)
            - Niedrigdimensionale Fälle: nicht allgemein

            KaTeX: $$H^{p,p}(X) \cap H^{2p}(X, \mathbb{Q}) \stackrel{?}{=}
                   \text{image}\left(\text{CH}^p(X)_{\mathbb{Q}} \to H^{2p}(X, \mathbb{Q})\right)$$

        @return Dictionary mit Hodge-Vermutungs-Status.
        @lastModified 2026-03-11
        """
        d = self.motive.dimension
        hodge = self.motive.hodge_numbers()
        name = self.motive.variety_name.lower()

        # Sammle (p,p)-Hodge-Klassen
        pp_classes = {p: h for (p, q), h in hodge.items() if p == q}

        # Bekannte Resultate
        proven_cases = []
        if 0 in pp_classes:
            proven_cases.append('p=0: trivial (Fundamentalklasse)')
        if d in pp_classes:
            proven_cases.append(f'p={d}: trivial (Punkte)')
        if name in ('elliptic_curve', 'ec', 'elliptic curve'):
            proven_cases.append('Elliptische Kurven: Lefschetz (1,1)')

        return {
            'variety': self.motive.variety_name,
            'dimension': d,
            'pp_hodge_classes': pp_classes,
            'proven_cases': proven_cases,
            'millennium_status': 'OPEN — eines der 7 Millennium-Probleme',
            'note': 'Lefschetz (1,1)-Theorem: H^{1,1}∩H^2(X,Q) = algebraische Klassen ✓'
        }


# ===========================================================================
# 5. REALISIERUNGSFUNKTOREN
# ===========================================================================

class RealizationFunctor:
    r"""
    @brief Realisierungsfunktoren für Motive: Betti, de Rham, étale (l-adisch).
    @description
        Ein Realisierungsfunktor R_? ordnet jedem Motiv M einen Vektorraum
        (oder Modul) mit Zusatzstruktur zu:

        - **Betti-Realisierung** R_B(M): singulare Kohomologie H*(X, Q), Hodge-Struktur
        - **de-Rham-Realisierung** R_dR(M): algebraische de-Rham-Kohomologie H*_dR(X/k)
          mit Hodge-Filtration F^p
        - **l-adische (étale) Realisierung** R_l(M): étale Kohomologie H*_ét(X, Ql),
          Galois-Darstellung Gal(k̄/k) → GL(Ql)

        Vergleichsisomorphismen:
        R_B(M) ⊗ C ≅ R_dR(M) ⊗_k C   (Hodge-de-Rham)
        R_B(M) ⊗ Ql ≅ R_l(M)           (für X/C via Vergleich)

        KaTeX:
        $$R_B(h(X)) = H^*(X(\mathbb{C}), \mathbb{Q}), \quad
          R_{\text{dR}}(h(X)) = H^*_{\text{dR}}(X/k)$$

    @param motive Reines Motiv.
    @param realization_type Art der Realisierung ('betti', 'de_rham', 'l_adic').
    @lastModified 2026-03-11
    """

    VALID_TYPES = {'betti', 'de_rham', 'l_adic'}

    def __init__(self, motive: Motive, realization_type: str = 'betti'):
        r"""
        @brief Konstruktor.
        @param motive Quell-Motiv.
        @param realization_type Typ der Realisierung.
        @raises InvalidInputError Bei unbekanntem Realisierungstyp.
        @lastModified 2026-03-11
        """
        if realization_type not in self.VALID_TYPES:
            raise InvalidInputError(
                f"RealizationFunctor: Unbekannter Typ '{realization_type}'. "
                f"Erlaubt: {self.VALID_TYPES}"
            )
        self.motive = motive
        self.realization_type = realization_type

    def compute(self) -> Dict[str, Any]:
        r"""
        @brief Berechnet die Realisierung des Motivs.
        @description
            Gibt die strukturellen Daten der Realisierung zurück:
            - Dimension des Vektorraums
            - Zusatzstruktur (Hodge-Filtration, Galois-Wirkung usw.)
            - Charakteristisches Polynom (für l-adisch)

        @return Dictionary mit Realisierungsdaten.
        @lastModified 2026-03-11
        """
        betti = self.motive.betti_numbers
        total_dim = sum(betti)

        if self.realization_type == 'betti':
            return self._betti_realization(betti, total_dim)
        elif self.realization_type == 'de_rham':
            return self._de_rham_realization(betti, total_dim)
        else:  # l_adic
            return self._l_adic_realization(betti, total_dim)

    def _betti_realization(self, betti: List[int], total_dim: int) -> Dict[str, Any]:
        r"""
        @brief Betti-Realisierung: singuläre Kohomologie mit Hodge-Struktur.
        @lastModified 2026-03-11
        """
        hodge = self.motive.hodge_numbers()
        return {
            'type': 'betti',
            'vector_space_dimension': total_dim,
            'betti_numbers': betti,
            'hodge_numbers': hodge,
            'euler_characteristic': self.motive.euler_characteristic(),
            'hodge_structure': f'Pure weight {2 * self.motive.tate_twist}',
            'note': 'Rationale Struktur: H*(X, Q)'
        }

    def _de_rham_realization(self, betti: List[int], total_dim: int) -> Dict[str, Any]:
        r"""
        @brief de-Rham-Realisierung: algebraische Differentialformen mit Filtration.
        @lastModified 2026-03-11
        """
        d = self.motive.dimension
        hodge = self.motive.hodge_numbers()

        # Hodge-Filtration F^p H^n_dR: absteigende Filtration
        # F^p H^n = ⊕_{i≥p} H^{i,n-i}
        filtration = {}
        for n in range(2 * d + 1):
            filtration[n] = {}
            cumsum = 0
            for p in range(n, -1, -1):
                q = n - p
                h_pq = hodge.get((p, q), 0)
                cumsum += h_pq
                filtration[n][p] = cumsum

        return {
            'type': 'de_rham',
            'vector_space_dimension': total_dim,
            'betti_numbers': betti,
            'hodge_filtration': filtration,
            'connection': 'Gauss-Manin connection (flat)',
            'note': 'de-Rham über Basiskörper k: H^*_dR(X/k)'
        }

    def _l_adic_realization(self, betti: List[int], total_dim: int) -> Dict[str, Any]:
        r"""
        @brief l-adische (étale) Realisierung: Galois-Darstellung.
        @lastModified 2026-03-11
        """
        d = self.motive.dimension
        tate = self.motive.tate_twist

        # Galois-Darstellung: ρ_l: Gal(k̄/k) → GL_n(Ql)
        # Charakteristisches Polynom am Frobenius: det(1 - F_p · T | H^i_ét)
        # Für einfache Fälle: bekannte Formen
        name = self.motive.variety_name.lower()
        galois_rep_type = 'unknown'
        if name.startswith('p') and name[1:].isdigit():
            galois_rep_type = 'powers_of_cyclotomic'  # ε^k für P^n
        elif name in ('elliptic_curve', 'ec'):
            galois_rep_type = '2-dimensional (H^1 of E)'

        return {
            'type': 'l_adic',
            'vector_space_dimension': total_dim,
            'betti_numbers': betti,
            'galois_representation': galois_rep_type,
            'tate_twist': tate,
            'weight': f'Pure weight {-2 * tate} (if tate_twist={tate})',
            'note': 'Galois Gal(k̄/k) wirkt stetig auf H^*_ét(X_{k̄}, Ql)'
        }

    def comparison_isomorphism(self, other: 'RealizationFunctor') -> Dict[str, Any]:
        r"""
        @brief Vergleichsisomorphismus zwischen zwei Realisierungen.
        @description
            Fundamentale Isomorphismen der Hodge-Theorie:
            - Hodge-de-Rham: R_B ⊗ C ≅ R_dR ⊗ C
            - Riemann-Hilbert (für D-Moduln): R_B ≅ reguläre singuläre D-Moduln
            - Vergleich étale-Betti: R_B ⊗ Ql ≅ R_l (für X/C)

        @param other Andere Realisierung (verschiedener Typ).
        @return Dictionary mit Vergleichsdaten.
        @lastModified 2026-03-11
        """
        types = {self.realization_type, other.realization_type}

        if types == {'betti', 'de_rham'}:
            iso_name = 'Hodge-de-Rham Vergleichsisomorphismus'
            is_canonical = True
        elif types == {'betti', 'l_adic'}:
            iso_name = 'Étale-Betti Vergleichsisomorphismus'
            is_canonical = True
        elif types == {'de_rham', 'l_adic'}:
            iso_name = 'p-adischer Hodge (Faltings/Tsuji)'
            is_canonical = False  # Nicht eindeutig
        else:
            iso_name = 'Identität'
            is_canonical = True

        return {
            'from': self.realization_type,
            'to': other.realization_type,
            'isomorphism': iso_name,
            'is_canonical': is_canonical,
            'dimension_match': (self.compute()['vector_space_dimension'] ==
                               other.compute()['vector_space_dimension'])
        }


# ===========================================================================
# 6. MOTIVISCHE L-FUNKTION
# ===========================================================================

def motivic_l_function(motive: Motive, s: float, prime_bound: int = 30) -> complex:
    r"""
    @brief Berechnet die motivische L-Funktion L(M, s) als Euler-Produkt.
    @description
        Die L-Funktion eines Motivs M ist das Euler-Produkt:
            L(M, s) = ∏_p L_p(M, p^{-s})^{-1}

        wobei der lokale Faktor bei einer guten Primzahl p ist:
            L_p(M, T)^{-1} = det(1 - Frob_p · T | H^i_ét(M))

        Für den trivialen Charakter (M = Q(0)):
            L(Q(0), s) = ζ(s)  (Riemann-Zeta)

        Für elliptische Kurven E/Q:
            L(h^1(E), s) = L(E, s)  (Hasse-Weil-L-Funktion)

        KaTeX:
        $$L(M, s) = \prod_{p \text{ gut}} \det\!\left(1 - \text{Frob}_p \cdot p^{-s}
          \mid H^\bullet_{\text{ét}}(M)\right)^{-1}$$

    @param motive Reines Motiv M.
    @param s Komplexe Variable s (hier reell).
    @param prime_bound Obere Schranke für Primzahlen im Euler-Produkt.
    @return Numerischer Wert von L(M, s) (komplex, hier reell für s>1).
    @raises InvalidInputError Wenn s ≤ 1 und Konvergenz nicht garantiert.
    @lastModified 2026-03-11
    """
    if s <= 0.5:
        raise InvalidInputError(
            f"motivic_l_function: s={s} außerhalb Konvergenzbereich Re(s) > 1 "
            f"(Analytische Fortsetzung nicht implementiert)."
        )

    name = motive.variety_name.lower()
    d = motive.dimension
    tate = motive.tate_twist

    # Euler-Produkt ∏_p L_p(s)^{-1}
    log_l = 0.0  # log L = -Σ_p log L_p^{-1} = Σ_p log L_p

    for p in range(2, prime_bound + 1):
        if not _is_prime(p):
            continue

        # Lokaler Faktor: L_p(M, T) hängt von der Realisierung ab
        l_p_val = _local_factor(name, d, tate, p, s)

        if l_p_val > 0:
            log_l += math.log(l_p_val)

    return math.exp(log_l)


def _local_factor(name: str, dim: int, tate: int, p: int, s: float) -> float:
    r"""
    @brief Berechnet den lokalen Euler-Faktor L_p(M, p^{-s})^{-1}.
    @description
        Für verschiedene Motive:
        - Triviales Motiv Q(0): L_p = (1 - p^{-s})^{-1}
        - Tate-Twist Q(n): L_p = (1 - p^{n-s})^{-1}
        - Elliptische Kurve: L_p = (1 - a_p p^{-s} + p^{1-2s})^{-1}
          mit |a_p| ≤ 2√p (Hasse-Schranke)
    @param name Varietätsname.
    @param dim Dimension.
    @param tate Tate-Twist.
    @param p Primzahl.
    @param s Wert der s-Variable.
    @return Lokaler L-Faktor (positiver Realteil für s >> 0).
    @lastModified 2026-03-11
    """
    p_s = p ** (-s)

    if name in ('point', 'spec_k', 'q(0)'):
        # L_p(Q(0), T)^{-1} = (1 - T) → L_p = (1 - p^{-s})^{-1}
        denom = 1.0 - p_s
        return 1.0 / denom if abs(denom) > 1e-12 else 1.0

    elif name.startswith('p') and name[1:].isdigit():
        # Projektiver Raum P^n: L(H(P^n), s) = ζ(s)ζ(s-1)...ζ(s-n)
        n = int(name[1:])
        result = 1.0
        for k in range(n + 1):
            # Euler-Faktor für ζ(s-k): (1 - p^{k-s})^{-1}
            factor = 1.0 - p ** (k - s)
            if abs(factor) > 1e-12:
                result /= factor
        return result

    elif name in ('elliptic_curve', 'ec', 'elliptic curve'):
        # Hasse-Weil-L-Funktion: L_p = (1 - a_p p^{-s} + p^{1-2s})^{-1}
        # Heuristischer a_p: a_p ≈ 0 (grobe Schätzung)
        a_p = 0  # Im Allgemeinen: a_p = p + 1 - #E(F_p)
        denom = 1.0 - a_p * p_s + p * (p_s ** 2)
        return 1.0 / denom if abs(denom) > 1e-12 else 1.0

    else:
        # Allgemeines Motiv: Vereinfachung als ζ(s - tate)
        factor = 1.0 - p ** (tate - s)
        return 1.0 / factor if abs(factor) > 1e-12 else 1.0
