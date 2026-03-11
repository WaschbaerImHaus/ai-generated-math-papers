"""
@file multilinear_algebra.py
@brief Multilineare Algebra: Tensorprodukte, äußere Algebra, symmetrische Algebra.
@description
    Dieses Modul implementiert grundlegende Konzepte der multilinearen Algebra:

    - Tensorprodukt: u⊗v (Rang-1-Tensor), Kronecker-Produkt A⊗B
    - Allgemeine (p,q)-Tensoren mit Kontraktion, Symmetrisierung, Antisymmetrisierung
    - Äußere Algebra Λ(V): Wedge-Produkt u∧v, k-Formen, Hodge-Dual
    - Gram-Matrix und Gram-Determinante
    - Symmetrische Algebra S(V): Polarisierung, Dimension
    - Multilineare Formen: Symmetrietests, Rang, Sylvester-Signatur

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import numpy as np
import math
from itertools import combinations, permutations
from typing import Optional


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _sign_permutation(perm: list) -> int:
    """
    Berechnet das Vorzeichen einer Permutation (+1 oder -1).
    Verwendet Zählverfahren für Inversionen.

    @param perm  Liste mit Permutation der Indizes 0..n-1
    @return      +1 für gerade, -1 für ungerade Permutation
    """
    n = len(perm)
    visited = [False] * n
    sign = 1
    # Zähle Zyklen → Vorzeichen = (-1)^(n - Anzahl_Zyklen)
    for i in range(n):
        if not visited[i]:
            cycle_len = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = perm[j]
                cycle_len += 1
            if cycle_len % 2 == 0:
                sign *= -1
    return sign


# ---------------------------------------------------------------------------
# Tensorprodukt
# ---------------------------------------------------------------------------

def tensor_product_vectors(u: list, v: list) -> list:
    """
    Berechnet das Tensorprodukt (äußeres Produkt) zweier Vektoren: u⊗v.

    Das Ergebnis ist ein Rang-1-Tensor (Matrix), dessen Eintrag an Position (i,j)
    gleich u[i] * v[j] ist:

        (u⊗v)_{ij} = uᵢ · vⱼ

    @param u     Erster Vektor (Liste mit n Elementen)
    @param v     Zweiter Vektor (Liste mit m Elementen)
    @return      n×m-Matrix als verschachtelte Liste
    @lastModified 2026-03-10
    """
    u_arr = np.array(u, dtype=float)
    v_arr = np.array(v, dtype=float)
    # np.outer berechnet u_i * v_j für alle i,j
    result = np.outer(u_arr, v_arr)
    return result.tolist()


def tensor_product_matrices(A: list, B: list) -> list:
    """
    Berechnet das Kronecker-Produkt (Tensorprodukt) zweier Matrizen: A⊗B.

    Wenn A eine m×n-Matrix und B eine p×q-Matrix ist, ist A⊗B eine
    (mp)×(nq)-Matrix mit Blockstruktur:

        A⊗B = [[a₁₁·B, a₁₂·B, ...],
                [a₂₁·B, a₂₂·B, ...], ...]

    @param A     Erste Matrix als verschachtelte Liste (m×n)
    @param B     Zweite Matrix als verschachtelte Liste (p×q)
    @return      (mp)×(nq)-Matrix als verschachtelte Liste
    @lastModified 2026-03-10
    """
    A_arr = np.array(A, dtype=float)
    B_arr = np.array(B, dtype=float)
    # np.kron berechnet das Kronecker-Produkt
    result = np.kron(A_arr, B_arr)
    return result.tolist()


class Tensor:
    """
    Allgemeiner (p,q)-Tensor: p kovariante, q kontravariante Indizes.
    Intern als numpy-Array der Form [n₁,...,nₚ,m₁,...,mq] gespeichert.

    Konvention in diesem Modul:
    - kovariante Indizes (untere Indizes) kommen zuerst
    - kontravariante Indizes (obere Indizes) kommen danach

    Beispiele:
    - (2,0)-Tensor: Bilinearform (Metrik g_{ij})
    - (0,2)-Tensor: Bivektorfeld
    - (1,1)-Tensor: Lineare Abbildung (Endomorphismus)
    """

    def __init__(self, components: np.ndarray, covariant: int, contravariant: int):
        """
        Initialisiert einen (covariant, contravariant)-Tensor.

        @param components      numpy-Array der Komponenten
        @param covariant       Anzahl kovarianter Indizes (untere Indizes)
        @param contravariant   Anzahl kontravarianter Indizes (obere Indizes)
        @lastModified 2026-03-10
        """
        self.components = np.array(components, dtype=float)
        self.covariant = covariant
        self.contravariant = contravariant
        # Rang des Tensors: Gesamtzahl der Indizes
        self.rank = covariant + contravariant
        # Überprüfe, dass die Dimensionen der Komponenten stimmen
        if len(self.components.shape) != self.rank and self.rank > 0:
            raise ValueError(
                f"Tensor Rang {self.rank} erwartet Array mit {self.rank} Dimensionen, "
                f"aber {len(self.components.shape)} wurden übergeben."
            )

    def contract(self, i: int, j: int) -> 'Tensor':
        """
        Kontraktion (Spurbildung) des Tensors über die Indizes i und j.

        Die Kontraktion summiert über einen kovarianten und einen kontravarianten
        Index: T^{...a...}_{...a...} (Einstein-Summenkonvention).

        Die resultierenden Tensor hat Rang (p-1, q-1).

        @param i   Erster Indexposition (0-basiert)
        @param j   Zweite Indexposition (0-basiert), muss verschieden von i sein
        @return    Kontrahierter Tensor mit zwei weniger Indizes
        @lastModified 2026-03-10
        """
        if i == j:
            raise ValueError("Die beiden Kontraktionsindizes müssen verschieden sein.")
        if i >= self.rank or j >= self.rank:
            raise ValueError(f"Indexposition außerhalb des Rangs {self.rank}")

        # np.trace führt Spur über zwei Achsen durch
        result = np.trace(self.components, axis1=i, axis2=j)

        # Neuer Typ: beide kontrahierten Indizes entfernt
        # Wir müssen wissen, welche Indizes kovariant/kontravariant waren
        new_cov = self.covariant - 1
        new_contra = self.contravariant - 1
        if new_cov < 0 or new_contra < 0:
            # Allgemeine Kontraktion zweier gleichartiger Indizes
            new_cov = max(0, self.covariant - 1)
            new_contra = max(0, self.contravariant - 1)

        return Tensor(result, new_cov, new_contra)

    def symmetrize(self) -> 'Tensor':
        """
        Symmetrisierung des Tensors über alle Indizes.

        T_{(i₁...iₙ)} = (1/n!) Σ_{σ∈Sₙ} T_{σ(i₁)...σ(iₙ)}

        @return  Vollständig symmetrisierter Tensor gleichen Rangs
        @lastModified 2026-03-10
        """
        n = self.rank
        if n == 0:
            return Tensor(self.components.copy(), self.covariant, self.contravariant)

        shape = self.components.shape
        result = np.zeros(shape)
        indices = list(range(n))
        count = 0

        # Summiere über alle Permutationen der Indizes
        for perm in permutations(indices):
            # Transponiere den Tensor gemäß der Permutation
            result += np.transpose(self.components, perm)
            count += 1

        # Normiere durch Anzahl der Permutationen (= n!)
        result /= count
        return Tensor(result, self.covariant, self.contravariant)

    def antisymmetrize(self) -> 'Tensor':
        """
        Antisymmetrisierung (Alternierung) des Tensors über alle Indizes.

        T_{[i₁...iₙ]} = (1/n!) Σ_{σ∈Sₙ} sgn(σ) T_{σ(i₁)...σ(iₙ)}

        @return  Vollständig antisymmetrisierter Tensor gleichen Rangs
        @lastModified 2026-03-10
        """
        n = self.rank
        if n == 0:
            return Tensor(self.components.copy(), self.covariant, self.contravariant)

        shape = self.components.shape
        result = np.zeros(shape)
        indices = list(range(n))
        count = 0

        # Summiere über alle Permutationen mit Vorzeichen
        for perm in permutations(indices):
            sgn = _sign_permutation(list(perm))
            result += sgn * np.transpose(self.components, perm)
            count += 1

        result /= count
        return Tensor(result, self.covariant, self.contravariant)

    def __add__(self, other: 'Tensor') -> 'Tensor':
        """
        Addition zweier Tensoren gleichen Typs und gleicher Form.

        @param other  Tensor gleichen Typs
        @return       Summentensor
        @lastModified 2026-03-10
        """
        if self.covariant != other.covariant or self.contravariant != other.contravariant:
            raise ValueError("Tensoren müssen denselben Typ (p,q) haben.")
        if self.components.shape != other.components.shape:
            raise ValueError("Tensoren müssen dieselbe Form haben.")
        return Tensor(
            self.components + other.components,
            self.covariant,
            self.contravariant
        )

    def __mul__(self, scalar) -> 'Tensor':
        """
        Skalarmultiplikation eines Tensors.

        @param scalar  Skalar (int oder float)
        @return        Skalierter Tensor
        @lastModified 2026-03-10
        """
        return Tensor(
            scalar * self.components,
            self.covariant,
            self.contravariant
        )

    def __repr__(self) -> str:
        return f"Tensor(Typ=({self.covariant},{self.contravariant}), Form={self.components.shape})"


def tensor_contraction(T: np.ndarray, i: int, j: int) -> np.ndarray:
    """
    Spurbildung (Kontraktion) eines numpy-Arrays über zwei Achsen i und j.

    Summiert über den gemeinsamen Index: T^{i...j...} mit gleichen Positionen i=j.
    Das Ergebnis hat zwei Indizes weniger als T.

        (Spur über Achsen i,j) = Σₐ T[..., a, ..., a, ...]

    @param T   numpy-Array (Tensor beliebigen Rangs)
    @param i   Erste Achse (0-basiert)
    @param j   Zweite Achse (0-basiert), muss i ≠ j gelten
    @return    Array mit zwei Dimensionen weniger
    @lastModified 2026-03-10
    """
    if i == j:
        raise ValueError("Die Kontraktionsachsen i und j müssen verschieden sein.")
    # np.trace führt Spur über zwei Achsen aus
    return np.trace(T, axis1=i, axis2=j)


def outer_product(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Äußeres Produkt zweier Arrays: (a⊗b)_{ij} = aᵢ · bⱼ.

    Für Vektoren a ∈ ℝⁿ und b ∈ ℝᵐ ergibt sich eine n×m-Matrix.

    @param a   Erster Vektor/Array
    @param b   Zweiter Vektor/Array
    @return    Äußeres Produkt als numpy-Array
    @lastModified 2026-03-10
    """
    return np.outer(a, b)


# ---------------------------------------------------------------------------
# Äußere (Graßmann-)Algebra
# ---------------------------------------------------------------------------

def wedge_product(u: list, v: list) -> list:
    """
    Äußeres (Wedge-)Produkt u∧v zweier Vektoren in ℝⁿ.

    Für Vektoren u,v ∈ ℝⁿ ist u∧v ein Bivector mit Komponenten:
        (u∧v)_{ij} = uᵢvⱼ - uⱼvᵢ  für i < j

    Die Dimension des Raums Λ²(ℝⁿ) ist C(n,2).

    Gibt den Koeffizientenvektor in der Standardbasis
    {e_i∧e_j | i < j} zurück.

    @param u   Erster Vektor (n-dimensional)
    @param v   Zweiter Vektor (n-dimensional)
    @return    Koeffizientenliste der Länge C(n,2)
    @lastModified 2026-03-10
    """
    n = len(u)
    if len(v) != n:
        raise ValueError(f"Vektoren müssen gleiche Dimension haben: {n} ≠ {len(v)}")

    u_arr = np.array(u, dtype=float)
    v_arr = np.array(v, dtype=float)

    # Berechne alle Koeffizienten (u∧v)_{ij} = uᵢvⱼ - uⱼvᵢ für i < j
    result = []
    for i, j in combinations(range(n), 2):
        result.append(u_arr[i] * v_arr[j] - u_arr[j] * v_arr[i])

    return result


def exterior_power(vectors: list, k: int) -> list:
    """
    k-te äußere Potenz: u₁∧u₂∧...∧uₖ für k Vektoren in ℝⁿ.

    Das Ergebnis ist ein k-Vektor mit Komponenten:
        (u₁∧...∧uₖ)_{i₁...iₖ} = det[uⱼ[iₗ]_{j,l=1}^k]

    Die Koeffizienten werden über alle geordneten k-Teilmengen von {0,...,n-1}
    berechnet.

    Spezialfall: k = n liefert die Determinante als einzigen Koeffizienten.

    @param vectors   Liste von k Vektoren, jeder n-dimensional
    @param k         Grad der äußeren Potenz
    @return          Koeffizientenliste der Länge C(n,k)
    @lastModified 2026-03-10
    """
    if len(vectors) != k:
        raise ValueError(f"{k} Vektoren erwartet, aber {len(vectors)} erhalten.")
    if k == 0:
        return [1.0]

    n = len(vectors[0])
    if any(len(v) != n for v in vectors):
        raise ValueError("Alle Vektoren müssen dieselbe Dimension haben.")

    mat = np.array(vectors, dtype=float)  # k×n-Matrix

    result = []
    # Iteriere über alle aufsteigenden k-Teilmengen der Indizes
    for indices in combinations(range(n), k):
        # Berechne die Unterdeterminante für diese Indexauswahl
        submatrix = mat[:, list(indices)]  # k×k-Matrix
        det_val = np.linalg.det(submatrix)
        result.append(det_val)

    return result


class ExteriorAlgebra:
    """
    Äußere Algebra Λ(V) = ⊕ₖ Λᵏ(V) für n-dimensionalen Vektorraum V = ℝⁿ.

    Λᵏ(V) hat Dimension C(n,k).
    Gesamtdimension: 2ⁿ (direkte Summe aller Grade).

    Eine k-Form wird als Dictionary dargestellt:
        {(i₁,i₂,...,iₖ): Koeffizient}
    wobei i₁ < i₂ < ... < iₖ (standardisierte aufsteigende Basis).
    """

    def __init__(self, dimension: int):
        """
        Initialisiert die äußere Algebra für ℝ^dimension.

        @param dimension   Dimension n des Vektorraums
        @lastModified 2026-03-10
        """
        self.n = dimension

    def basis_k_forms(self, k: int) -> list:
        """
        Gibt die Basisvektoren von Λᵏ(V) als Liste von Tupeln zurück.

        Jedes Tupel (i₁,...,iₖ) mit i₁ < i₂ < ... < iₖ repräsentiert
        die Basisform eᵢ₁∧eᵢ₂∧...∧eᵢₖ.

        @param k   Grad der k-Formen
        @return    Liste aller Basis-Tupel, Länge = C(n,k)
        @lastModified 2026-03-10
        """
        return list(combinations(range(self.n), k))

    def wedge(self, form1: dict, form2: dict) -> dict:
        """
        Keilprodukt (Wedge-Produkt) zweier Differentialformen.

        Wenn form1 eine p-Form und form2 eine q-Form ist, ist form1∧form2 eine
        (p+q)-Form:

            (α∧β)(v₁,...,vₚ₊ᵩ) = Σ_{σ} sgn(σ) α(v_{σ(1)},...) β(...)

        Die Formen sind als Dictionaries {Tupel: Koeffizient} kodiert.
        Tupel mit i₁ < i₂ < ... < iₖ sind kanonische Basis-Indizes.

        @param form1   p-Form als Dictionary
        @param form2   q-Form als Dictionary
        @return        (p+q)-Form als Dictionary
        @lastModified 2026-03-10
        """
        result = {}

        for idx1, coeff1 in form1.items():
            for idx2, coeff2 in form2.items():
                # Prüfe auf gemeinsame Indizes (äußeres Produkt = 0 bei Wiederholung)
                combined = list(idx1) + list(idx2)
                if len(set(combined)) < len(combined):
                    # Wiederholter Index → Beitrag = 0
                    continue

                # Sortiere die Indizes und berechne das Vorzeichen der Permutation
                sorted_combined = sorted(combined)
                # Berechne Vorzeichen: Anzahl der Vertauschungen beim Sortieren
                perm = [combined.index(x) for x in sorted_combined]
                # Entferne bereits verwendete Positionen für Duplikate
                perm_unique = []
                used = []
                for x in sorted_combined:
                    positions = [i for i, v in enumerate(combined) if v == x and i not in used]
                    perm_unique.append(positions[0])
                    used.append(positions[0])
                sgn = _sign_permutation(perm_unique)

                key = tuple(sorted_combined)
                result[key] = result.get(key, 0.0) + sgn * coeff1 * coeff2

        # Entferne Nulleinträge
        return {k: v for k, v in result.items() if abs(v) > 1e-15}

    def hodge_dual(self, form: dict, metric: Optional[np.ndarray] = None) -> dict:
        """
        Hodge-Dual ⋆α einer k-Form in ℝⁿ (oder mit Metrik g).

        Für eine k-Form α ist ⋆α eine (n-k)-Form:
            ⋆(eᵢ₁∧...∧eᵢₖ) = ε^{i₁...iₖ}_{j₁...j_{n-k}} eʲ¹∧...∧eʲ^{n-k}

        Für die Standardmetrik: ⋆(eᵢ₁∧...∧eᵢₖ) = sgn(σ) eⱼ₁∧...∧eⱼ_{n-k}
        wobei {i₁,...,iₖ,j₁,...,j_{n-k}} eine Permutation von {0,...,n-1} ist.

        @param form    k-Form als Dictionary
        @param metric  Optionaler Metriktensor (n×n-Matrix); None = Euklidisch
        @return        Hodge-Dual als (n-k)-Form
        @lastModified 2026-03-10
        """
        n = self.n

        # Metrik-Determinante (Standard: 1 für Euklidische Metrik)
        if metric is None:
            metric_sqrt_det = 1.0
            metric_inv = np.eye(n)
        else:
            det_g = np.linalg.det(metric)
            metric_sqrt_det = math.sqrt(abs(det_g))
            metric_inv = np.linalg.inv(metric)

        result = {}
        all_indices = set(range(n))

        for idx_tuple, coeff in form.items():
            k = len(idx_tuple)
            # Komplementäre Indizes (für die Dual-Form)
            complement = sorted(all_indices - set(idx_tuple))

            if len(complement) != n - k:
                continue

            # Berechne das Vorzeichen der Permutation (idx_tuple + complement) → (0,1,...,n-1)
            full_perm = list(idx_tuple) + complement
            sgn = _sign_permutation(full_perm)

            key = tuple(complement)
            result[key] = result.get(key, 0.0) + sgn * metric_sqrt_det * coeff

        return {k: v for k, v in result.items() if abs(v) > 1e-15}


def gram_matrix(vectors: list) -> list:
    """
    Berechnet die Gram-Matrix G der Vektoren.

    G_{ij} = ⟨vᵢ, vⱼ⟩ (inneres Produkt / Skalarprodukt)

    Die Gram-Matrix ist positiv semidefinit für beliebige Vektoren
    und positiv definit genau dann, wenn die Vektoren linear unabhängig sind.

    @param vectors   Liste von Vektoren (alle gleicher Dimension)
    @return          k×k Gram-Matrix als verschachtelte Liste
    @lastModified 2026-03-10
    """
    mat = np.array(vectors, dtype=float)
    # G = V · Vᵀ (Skalarprodukte aller Paare)
    G = mat @ mat.T
    return G.tolist()


def gram_determinant(vectors: list) -> float:
    """
    Berechnet die Gram-Determinante det(G) = |v₁∧...∧vₖ|².

    Die Gram-Determinante ist das Quadrat des Volumens des Parallelotops,
    das von den Vektoren aufgespannt wird:

        det(G) = |v₁∧...∧vₖ|² ≥ 0

    Spezialfall k=n: det(G) = det(V)²

    @param vectors   Liste von k Vektoren (k ≤ n)
    @return          Gram-Determinante (nicht-negativ)
    @lastModified 2026-03-10
    """
    G = np.array(gram_matrix(vectors))
    return float(np.linalg.det(G))


# ---------------------------------------------------------------------------
# Symmetrische Algebra
# ---------------------------------------------------------------------------

class SymmetricAlgebra:
    """
    Symmetrische Algebra S(V) = ⊕ₖ Sᵏ(V) für V = ℝⁿ.

    Sᵏ(V) = Raum der vollständig symmetrischen k-Tensoren.
    Isomorph zum Polynomring ℝ[x₁,...,xₙ] in Grad k.

    dim(Sᵏ(V)) = C(n+k-1, k) = "n+k-1 über k"
    (Anzahl Monomials vom Grad k in n Variablen)
    """

    def __init__(self, dimension: int):
        """
        Initialisiert die symmetrische Algebra für ℝ^dimension.

        @param dimension   Dimension n des Vektorraums
        @lastModified 2026-03-10
        """
        self.n = dimension

    def symmetric_power_dim(self, k: int) -> int:
        """
        Berechnet die Dimension von Sᵏ(V).

        dim(Sᵏ(V)) = C(n+k-1, k) = (n+k-1)! / (k! · (n-1)!)

        @param k   Grad der symmetrischen Potenz
        @return    Dimension des Raums Sᵏ(V)
        @lastModified 2026-03-10
        """
        return math.comb(self.n + k - 1, k)

    def polarization(self, poly_coeffs: dict, k: int) -> dict:
        """
        Polarisierung eines homogenen Polynoms vom Grad k.

        Die Polarisierung eines Polynoms p(x) vom Grad k liefert eine
        symmetrische k-lineare Form P(v₁,...,vₖ), so dass P(v,...,v) = p(v).

        poly_coeffs: {(e₁,...,eₙ): Koeffizient} mit |e₁+...+eₙ| = k
        (Multiindex-Darstellung der Monomials)

        Formel: P(v₁,...,vₖ) = (1/k!) Σ_{S⊆{1,...,k}} (-1)^{k-|S|} p(Σᵢ∈S vᵢ)

        @param poly_coeffs   Polynomkoeffizienten als Multiindex-Dictionary
        @param k             Grad des Polynoms
        @return              Symmetrisierte Form (Multiindex-Dictionary)
        @lastModified 2026-03-10
        """
        # Normiere die Koeffizienten durch den Symmetrisierungsfaktor
        # Jeder Multiindex (e₁,...,eₙ) bekommt Faktor k! / (e₁!·...·eₙ!)
        result = {}
        for multiindex, coeff in poly_coeffs.items():
            # Berechne den Symmetrisierungsfaktor
            denom = 1
            for e in multiindex:
                denom *= math.factorial(e)
            factor = math.factorial(k) / denom
            result[multiindex] = coeff / factor

        return result


# ---------------------------------------------------------------------------
# Multilineare Formen
# ---------------------------------------------------------------------------

def is_alternating(T: np.ndarray) -> bool:
    """
    Prüft ob ein Tensor alternierend (antisymmetrisch) ist.

    Ein Tensor T ist alternierend, wenn für jede Transposition zweier Indizes
    das Vorzeichen wechselt:
        T[...,i,...,j,...] = -T[...,j,...,i,...]

    Äquivalent: T = 0, sobald zwei Indizes gleich sind.

    @param T   numpy-Array (Tensor)
    @return    True wenn alternierend, False sonst
    @lastModified 2026-03-10
    """
    ndim = T.ndim
    if ndim < 2:
        return True  # Skalare und Vektoren sind trivial alternierend

    # Prüfe alle Paare von Achsen
    for i in range(ndim):
        for j in range(i + 1, ndim):
            # Transponiere Achsen i und j
            T_transposed = np.swapaxes(T, i, j)
            # Prüfe ob T = -T_transposed (mit numerischer Toleranz)
            if not np.allclose(T, -T_transposed, atol=1e-10):
                return False
    return True


def is_symmetric_tensor(T: np.ndarray) -> bool:
    """
    Prüft ob ein Tensor vollständig symmetrisch ist.

    Ein Tensor T ist symmetrisch, wenn die Komponenten unter allen
    Permutationen der Indizes invariant bleiben:
        T[...,i,...,j,...] = T[...,j,...,i,...]

    @param T   numpy-Array (Tensor)
    @return    True wenn symmetrisch, False sonst
    @lastModified 2026-03-10
    """
    ndim = T.ndim
    if ndim < 2:
        return True

    # Prüfe alle Paare von Achsen
    for i in range(ndim):
        for j in range(i + 1, ndim):
            T_transposed = np.swapaxes(T, i, j)
            if not np.allclose(T, T_transposed, atol=1e-10):
                return False
    return True


def multilinear_form_rank(M: np.ndarray) -> int:
    """
    Berechnet den Rang eines (0,2)-Tensors (Bilinearform) als Matrix-Rang.

    Für eine symmetrische Bilinearform B: V×V → ℝ ist der Rang gleich
    der Anzahl der nicht-verschwindenden Eigenwerte der darstellenden Matrix.

    @param M   2D-numpy-Array (quadratische oder rechteckige Matrix)
    @return    Rang der Matrix (Anzahl linear unabhängiger Zeilen/Spalten)
    @lastModified 2026-03-10
    """
    M_arr = np.array(M, dtype=float)
    return int(np.linalg.matrix_rank(M_arr, tol=1e-10))


def signature_bilinear_form(M: np.ndarray) -> tuple:
    """
    Berechnet die Signatur (p, q) einer symmetrischen Bilinearform.

    Nach dem Sylvester'schen Trägheitssatz ist die Signatur eine Invariante
    der Bilinearform (unabhängig von der Basiswahl):
    - p = Anzahl positiver Eigenwerte
    - q = Anzahl negativer Eigenwerte
    - Nullstellen (Entartung) werden ignoriert

    Beispiele:
    - Euklidisches Skalarprodukt: Signatur (n, 0)
    - Minkowski-Metrik: Signatur (1, 3) oder (3, 1)
    - Indefinite Form: p > 0 und q > 0

    @param M   Symmetrische n×n-Matrix der Bilinearform
    @return    Tupel (p, q) mit Anzahl positiver und negativer Eigenwerte
    @lastModified 2026-03-10
    """
    M_arr = np.array(M, dtype=float)

    # Symmetrisiere die Matrix (für numerische Stabilität)
    M_sym = (M_arr + M_arr.T) / 2

    # Berechne Eigenwerte (nur reell für symmetrische Matrizen)
    eigenvalues = np.linalg.eigvalsh(M_sym)

    # Zähle positive und negative Eigenwerte
    tol = 1e-10
    p = int(np.sum(eigenvalues > tol))   # Anzahl positiver Eigenwerte
    q = int(np.sum(eigenvalues < -tol))  # Anzahl negativer Eigenwerte

    return (p, q)
