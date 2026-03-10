"""
@file geometric_number_theory.py
@brief Geometrische Zahlentheorie (Geometry of Numbers): Gitter, LLL-Reduktion,
       Minkowski-Sätze, Quadratische Formen, Gauß-Komposition, Kugelpackungen.
@description
    Dieses Modul implementiert die Kernkonzepte der Geometrischen Zahlentheorie:
    - Gitter im ℝⁿ (Basis, Determinante, duales Gitter, kürzester Vektor)
    - LLL-Algorithmus (Lenstra-Lenstra-Lovász) zur Gitterbasisreduktion
    - Minkowski-Sätze (konvexer Körper, Linearformen, Minkowski-Schranke)
    - Binäre quadratische Formen (Diskriminante, Reduktion, Äquivalenz)
    - Gauß'sche Komposition von quadratischen Formen (Klassenzahl, Klassengruppe)
    - Kugelpackungstheorie (Kusszahlen, Packungsdichten, Leech-Gitter)

    Mathematische Grundlagen:
    - Ein Gitter Λ = {Σ aᵢbᵢ | aᵢ ∈ ℤ} mit Basisvektoren b₁,...,bₙ
    - det(Λ) = |det(B)| ist das Volumen des Fundamentalbereichs
    - Minkowski: Vol(K) > 2ⁿ·det(Λ) → K enthält einen von 0 verschiedenen Gitterpunkt

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import numpy as np
from typing import List, Optional, Tuple
import math
import fractions


# ============================================================
# Klasse: Lattice (Gitter im ℝⁿ)
# ============================================================

class Lattice:
    """
    @class Lattice
    @brief Repräsentiert ein Gitter im ℝⁿ.
    @description
        Ein Gitter Λ ⊂ ℝⁿ ist die Menge aller ganzzahligen Linearkombinationen
        von linear unabhängigen Basisvektoren b₁, ..., bₙ:
            Λ = { Σᵢ aᵢ bᵢ | aᵢ ∈ ℤ }
        Die Basismatrix B hat die Basisvektoren als Zeilen.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, basis_matrix: np.ndarray):
        """
        @brief Initialisiert das Gitter mit einer Basismatrix.
        @description
            Jede Zeile der Basismatrix entspricht einem Basisvektor des Gitters.
            Die Matrix muss quadratisch sein (n×n) oder kann m×n sein (m ≤ n).

        @param basis_matrix numpy-Array der Form (m, n) - jede Zeile ein Basisvektor
        @raises ValueError wenn die Basis singulär ist (Determinante ≈ 0)
        @lastModified 2026-03-10
        """
        # Sicherstellen, dass die Basis als Float gespeichert wird
        self.basis = np.array(basis_matrix, dtype=float)

        # Dimensionen speichern: m Basisvektoren in ℝⁿ
        if self.basis.ndim == 1:
            # Ein einzelner Vektor → 1D-Gitter
            self.basis = self.basis.reshape(1, -1)

        self.m, self.n = self.basis.shape  # m Zeilen, n Spalten

    def determinant(self) -> float:
        """
        @brief Berechnet die Gitterdeterminante det(B).
        @description
            Für eine quadratische Basis (n×n) ist det(Λ) = det(B).
            Für m < n gilt det(Λ) = √(det(BᵀB)) (Gram-Determinante).
            Die Determinante gibt das n-dimensionale Volumen des Fundamentalbereichs an.

        @return Gitterdeterminante (vorzeichenbehaftet für quadratische Basen)
        @lastModified 2026-03-10
        """
        if self.m == self.n:
            # Quadratische Basis: direkte Determinante
            return np.linalg.det(self.basis)
        else:
            # Nicht-quadratische Basis: Gram-Determinante √(det(BᵀB))
            gram = self.gram_matrix()
            return np.sqrt(max(0.0, np.linalg.det(gram)))

    def fundamental_domain_volume(self) -> float:
        """
        @brief Berechnet das Volumen des Fundamentalbereichs |det(B)|.
        @description
            Das Volumen des Fundamentalparallelotops ist |det(B)|.
            Dies entspricht der Dichte des Gitters: Je kleiner det(Λ),
            desto dichter liegt das Gitter im Raum.

        @return Volumen des Fundamentalbereichs (immer nicht-negativ)
        @lastModified 2026-03-10
        """
        return abs(self.determinant())

    def contains_point(self, point: np.ndarray, tol: float = 1e-10) -> bool:
        """
        @brief Prüft, ob ein Punkt zum Gitter gehört.
        @description
            Ein Punkt p liegt im Gitter Λ genau dann, wenn p = B·a für ein a ∈ ℤⁿ.
            Das heißt, die Koordinaten in der Basis müssen ganzzahlig sein.
            Berechnet a = B⁻¹·p und prüft ob alle Einträge ganzzahlig sind.

        @param point numpy-Vektor im ℝⁿ
        @param tol Toleranz für Ganzzahligkeitsprüfung (Standard: 1e-10)
        @return True wenn der Punkt im Gitter liegt, sonst False
        @lastModified 2026-03-10
        """
        p = np.array(point, dtype=float)

        if self.m == self.n:
            # Quadratische Basis: Löse B·a = p nach a
            det = np.linalg.det(self.basis)
            if abs(det) < 1e-14:
                # Singuläre Basis: Sonderfall
                return False
            try:
                coords = np.linalg.solve(self.basis.T, p)
            except np.linalg.LinAlgError:
                return False
        else:
            # Überbestimmtes System: Kleinste-Quadrate-Lösung
            coords, residuals, _, _ = np.linalg.lstsq(self.basis.T, p, rcond=None)
            # Überprüfe ob die Lösung exakt ist (Residuum ≈ 0)
            reconstructed = self.basis.T @ coords
            if np.linalg.norm(reconstructed - p) > tol:
                return False

        # Alle Koordinaten müssen ganzzahlig sein
        rounded = np.round(coords)
        return np.all(np.abs(coords - rounded) < tol)

    def shortest_vector_approx(self) -> np.ndarray:
        """
        @brief Näherung des kürzesten Gittervektors via LLL-Reduktion.
        @description
            Das Problem des kürzesten Vektors (SVP) ist NP-hart.
            LLL liefert eine Approximation mit Garantie:
            ||b₁|| ≤ 2^{(n-1)/4} · λ₁(Λ)
            wobei λ₁(Λ) das erste Sukzessivminimum ist.

        @return Näherung des kürzesten Gittervektors als numpy-Array
        @lastModified 2026-03-10
        """
        # LLL-Reduktion anwenden
        lll = LLLReduction()
        reduced_basis = lll.reduce(self.basis)

        # Der erste Vektor der LLL-reduzierten Basis ist die beste Näherung
        # Wähle den kürzesten Vektor aus allen reduzierten Basisvektoren
        norms = [np.linalg.norm(reduced_basis[i]) for i in range(reduced_basis.shape[0])]
        idx = np.argmin(norms)
        return reduced_basis[idx]

    def successive_minima_approx(self, k: int) -> float:
        """
        @brief Näherung des k-ten Sukzessivminimums λₖ(Λ).
        @description
            Das k-te Sukzessivminimum λₖ ist der kleinste Radius r, sodass
            k linear unabhängige Gittervektoren in der Kugel B(0, r) liegen.
            Via LLL werden die k kürzesten Basisvektoren als Näherung verwendet.

        @param k Ordnung des Sukzessivminimums (1 ≤ k ≤ n)
        @return Näherungswert für λₖ(Λ)
        @raises ValueError wenn k außerhalb des gültigen Bereichs liegt
        @lastModified 2026-03-10
        """
        if k < 1 or k > self.m:
            raise ValueError(f"k muss zwischen 1 und {self.m} liegen, aber k={k}")

        # LLL-Reduktion
        lll = LLLReduction()
        reduced_basis = lll.reduce(self.basis)

        # Sortiere Basisvektoren nach Länge
        norms = sorted([np.linalg.norm(reduced_basis[i]) for i in range(reduced_basis.shape[0])])

        # k-tes Sukzessivminimum ≈ k-te kürzeste LLL-Basisvektorlänge
        return norms[k - 1]

    def gram_matrix(self) -> np.ndarray:
        """
        @brief Berechnet die Gram-Matrix G = B·Bᵀ.
        @description
            Die Gram-Matrix enthält die inneren Produkte aller Basisvektoren:
            G[i,j] = <bᵢ, bⱼ>
            Sie ist positiv semidefinit und charakterisiert die Geometrie des Gitters.

        @return Gram-Matrix als (m×m) numpy-Array
        @lastModified 2026-03-10
        """
        # G[i,j] = bᵢ · bⱼ = Skalarprodukt der i-ten und j-ten Basisvektoren
        return self.basis @ self.basis.T

    def dual_lattice_basis(self) -> np.ndarray:
        """
        @brief Berechnet die Basis des dualen Gitters Λ* = (B⁻¹)ᵀ.
        @description
            Das duale Gitter Λ* = {x ∈ ℝⁿ | <x, λ> ∈ ℤ für alle λ ∈ Λ}.
            Für eine quadratische Basis B gilt: B* = (B⁻¹)ᵀ = (Bᵀ)⁻¹.
            Es gilt det(Λ*) = 1/det(Λ).

        @return Basismatrix des dualen Gitters als (n×n) numpy-Array
        @raises ValueError wenn die Basis singulär ist
        @lastModified 2026-03-10
        """
        if self.m != self.n:
            raise ValueError("Duale Gitterbasis nur für quadratische Basis definiert")

        det = np.linalg.det(self.basis)
        if abs(det) < 1e-14:
            raise ValueError("Singuläre Basis: Duales Gitter nicht definiert")

        # (B⁻¹)ᵀ = (Bᵀ)⁻¹
        return np.linalg.inv(self.basis).T


# ============================================================
# Klasse: LLLReduction (Lenstra-Lenstra-Lovász Algorithmus)
# ============================================================

class LLLReduction:
    """
    @class LLLReduction
    @brief Implementiert den LLL-Gitterbasisreduktionsalgorithmus.
    @description
        Der LLL-Algorithmus (1982, Lenstra-Lenstra-Lovász) findet in polynomialer Zeit
        eine "fast kurze" und "fast orthogonale" Gitterbasis.

        Eigenschaften der LLL-reduzierten Basis (b₁,...,bₙ):
        - Größenkondition: |μᵢⱼ| ≤ 1/2 für alle i > j
        - Lovász-Bedingung: δ·||b*ₖ₋₁||² ≤ ||b*ₖ + μₖ,ₖ₋₁·b*ₖ₋₁||²

        Anwendungen: Kryptanalyse, Faktorisierung, diophantische Approximation.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, delta: float = 0.75):
        """
        @brief Initialisiert LLL mit dem Lovász-Parameter δ.
        @description
            Der Parameter δ steuert den Trade-off zwischen Qualität und Laufzeit.
            Typisch: δ = 3/4 = 0.75 (originale LLL-Wahl).
            Gültig: 1/4 < δ < 1. Größeres δ → bessere Reduktion, mehr Iterationen.

        @param delta Lovász-Parameter (Standard: 0.75)
        @raises ValueError wenn δ außerhalb (1/4, 1) liegt
        @lastModified 2026-03-10
        """
        if not (0.25 < delta < 1.0):
            raise ValueError(f"delta muss in (0.25, 1.0) liegen, aber delta={delta}")
        self.delta = delta

    def gram_schmidt_orthogonalization(self, basis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Gram-Schmidt-Orthogonalisierung ohne Normierung.
        @description
            Berechnet die Gram-Schmidt-Koeffizienten μᵢⱼ und die orthogonalen Vektoren b*ᵢ:
                b*₁ = b₁
                b*ᵢ = bᵢ - Σⱼ<ᵢ μᵢⱼ · b*ⱼ
                μᵢⱼ = <bᵢ, b*ⱼ> / <b*ⱼ, b*ⱼ>

        @param basis numpy-Array (m×n) der Basisvektoren (Zeilen)
        @return Tupel (B_star, mu): orthogonalisierte Basis und Gram-Schmidt-Koeffizientenmatrix
        @lastModified 2026-03-10
        """
        m, n = basis.shape
        # b_star: orthogonalisierte Basisvektoren (ohne Normierung)
        b_star = np.zeros_like(basis, dtype=float)
        # mu: Gram-Schmidt-Koeffizienten μᵢⱼ (untere Dreiecksmatrix)
        mu = np.zeros((m, m), dtype=float)

        for i in range(m):
            b_star[i] = basis[i].copy()
            for j in range(i):
                # μᵢⱼ = <bᵢ, b*ⱼ> / <b*ⱼ, b*ⱼ>
                norm_sq = np.dot(b_star[j], b_star[j])
                if norm_sq > 1e-14:
                    mu[i, j] = np.dot(basis[i], b_star[j]) / norm_sq
                else:
                    mu[i, j] = 0.0
                # Projektionsanteil abziehen
                b_star[i] -= mu[i, j] * b_star[j]

        return b_star, mu

    def reduce(self, basis_matrix: np.ndarray) -> np.ndarray:
        """
        @brief Führt die LLL-Reduktion durch und gibt die reduzierte Basis zurück.
        @description
            LLL-Algorithmus (Lenstra-Lenstra-Lovász, 1982):
            1. Gram-Schmidt-Orthogonalisierung berechnen
            2. Größenreduktion: μᵢⱼ ← round(μᵢⱼ), bᵢ ← bᵢ - round(μᵢⱼ)·bⱼ
            3. Lovász-Bedingung prüfen: δ·||b*ₖ₋₁||² ≤ ||b*ₖ + μₖ,ₖ₋₁·b*ₖ₋₁||²
            4. Falls verletzt: bₖ und bₖ₋₁ tauschen, zurück zu Schritt 1

        @param basis_matrix numpy-Array (m×n) der Gitterbasis (Zeilen = Basisvektoren)
        @return LLL-reduzierte Basis als numpy-Array
        @lastModified 2026-03-10
        """
        # Arbeitskopie der Basis als Float
        B = np.array(basis_matrix, dtype=float)
        m, n = B.shape

        if m == 1:
            # Trivialfall: Ein Vektor ist immer reduziert
            return B.copy()

        # Maximale Iterationen zur Sicherheit (verhindert Endlosschleifen)
        max_iter = 1000 * m * m
        iteration = 0

        k = 1  # Aktueller Index im LLL-Algorithmus

        while k < m and iteration < max_iter:
            iteration += 1

            # Gram-Schmidt neu berechnen
            b_star, mu = self.gram_schmidt_orthogonalization(B)

            # Größenreduktion für B[k] bezüglich B[k-1], ..., B[0]
            for j in range(k - 1, -1, -1):
                # Runde μₖⱼ auf die nächste ganze Zahl
                mu_round = round(mu[k, j])
                if mu_round != 0:
                    # bₖ ← bₖ - round(μₖⱼ) · bⱼ
                    B[k] -= mu_round * B[j]
                    # Gram-Schmidt nach Änderung neu berechnen
                    b_star, mu = self.gram_schmidt_orthogonalization(B)

            # Lovász-Bedingung prüfen:
            # δ·||b*_{k-1}||² ≤ ||b*_k + μ_{k,k-1}·b*_{k-1}||²
            # Vereinfacht: δ·||b*_{k-1}||² ≤ (δ - μ_{k,k-1}²)·||b*_{k-1}||² + ||b*_k||²
            norm_k_minus_1_sq = np.dot(b_star[k - 1], b_star[k - 1])
            norm_k_sq = np.dot(b_star[k], b_star[k])

            lovasz_lhs = self.delta * norm_k_minus_1_sq
            lovasz_rhs = norm_k_sq + mu[k, k - 1] ** 2 * norm_k_minus_1_sq

            if lovasz_lhs > lovasz_rhs + 1e-10:
                # Lovász-Bedingung verletzt: Tausche B[k] und B[k-1]
                B[[k, k - 1]] = B[[k - 1, k]]
                # Einen Schritt zurückgehen
                k = max(k - 1, 1)
            else:
                # Weiter zum nächsten Vektor
                k += 1

        return B

    def is_lll_reduced(self, basis_matrix: np.ndarray) -> bool:
        """
        @brief Prüft, ob eine Basis LLL-reduziert ist.
        @description
            Eine Basis ist LLL-reduziert wenn:
            1. Größenkondition: |μᵢⱼ| ≤ 1/2 für alle i > j
            2. Lovász-Bedingung: δ·||b*ₖ₋₁||² ≤ (μₖ,ₖ₋₁² + 1)·||b*ₖ₋₁||² korrekt

        @param basis_matrix numpy-Array der zu prüfenden Basis
        @return True wenn LLL-reduziert, sonst False
        @lastModified 2026-03-10
        """
        B = np.array(basis_matrix, dtype=float)
        m = B.shape[0]

        if m == 1:
            return True

        b_star, mu = self.gram_schmidt_orthogonalization(B)

        for i in range(m):
            for j in range(i):
                # Größenkondition: |μᵢⱼ| ≤ 1/2
                if abs(mu[i, j]) > 0.5 + 1e-10:
                    return False

        for k in range(1, m):
            # Lovász-Bedingung prüfen
            norm_k_minus_1_sq = np.dot(b_star[k - 1], b_star[k - 1])
            norm_k_sq = np.dot(b_star[k], b_star[k])

            lovasz_lhs = self.delta * norm_k_minus_1_sq
            lovasz_rhs = norm_k_sq + mu[k, k - 1] ** 2 * norm_k_minus_1_sq

            if lovasz_lhs > lovasz_rhs + 1e-10:
                return False

        return True


# ============================================================
# Minkowski-Sätze (Standalone-Funktionen)
# ============================================================

def minkowski_convex_body_theorem(volume: float, det_lattice: float, dim: int) -> bool:
    """
    @brief Prüft, ob Minkowskis Satz über konvexe Körper einen Gitterpunkt garantiert.
    @description
        Minkowskis Fundamentalsatz (1889):
        Sei K ⊂ ℝⁿ ein konvexer, symmetrischer Körper (K = -K) mit
            Vol(K) > 2ⁿ · det(Λ)
        Dann enthält K mindestens einen von 0 verschiedenen Gitterpunkt.

        Der Satz ist scharf: Das Volumen 2ⁿ·det(Λ) kann ohne Gitterpunkt erreicht werden
        (z.B. halboffene Würfel).

    @param volume Volumen des konvexen Körpers K
    @param det_lattice Determinante des Gitters |det(Λ)|
    @param dim Dimension n des Raumes
    @return True wenn Minkowski's Satz einen Gitterpunkt garantiert
    @lastModified 2026-03-10
    """
    # Schwellenwert: 2ⁿ · det(Λ)
    threshold = (2 ** dim) * abs(det_lattice)
    # Satz gilt wenn Vol(K) > Schwellenwert (strikt)
    return volume > threshold


def minkowski_linear_forms(lambda_coefficients: np.ndarray, bounds: np.ndarray) -> bool:
    """
    @brief Minkowskis Satz über simultane Linearformen.
    @description
        Minkowskis Satz über simultane Diophantische Approximation:
        Für n reelle Zahlen λ₁,...,λₙ und positive c₁,...,cₙ mit c₁·...·cₙ ≥ 1
        existieren ganze Zahlen q₁,...,qₙ (nicht alle 0) mit
            |λ₁q₁ + ... + λₙqₙ| < c₁
        und |qᵢ| ≤ cᵢ für i=2,...,n.

        Diese vereinfachte Version prüft: Π bounds[i] ≥ 1

    @param lambda_coefficients Koeffizienten der Linearformen (numpy-Array)
    @param bounds Schranken c₁,...,cₙ (numpy-Array, alle > 0)
    @return True wenn Minkowskis Satz anwendbar ist (Produkt der Schranken ≥ 1)
    @lastModified 2026-03-10
    """
    bounds_arr = np.array(bounds, dtype=float)

    if np.any(bounds_arr <= 0):
        return False

    # Produkt der Schranken muss ≥ 1 sein
    product = np.prod(bounds_arr)
    return product >= 1.0


def minkowski_bound(discriminant: int, degree: int, num_complex_pairs: int) -> float:
    """
    @brief Berechnet die Minkowski-Schranke für die Klassenzahl eines Zahlkörpers.
    @description
        Die Minkowski-Schranke M_K gibt an, dass jede Idealklasse eines algebraischen
        Zahlkörpers K vom Grad n einen ganzen Idealteiler mit Norm ≤ M_K enthält:
            M_K = (n!/nⁿ) · (4/π)^r₂ · √|Δ_K|
        wobei:
        - n = [K:Q] der Grad des Zahlkörpers
        - r₂ = Anzahl der komplexen Einbettungspaare
        - Δ_K = Diskriminante des Zahlkörpers

        Zur Berechnung der Klassenzahl genügt es, alle Primzahlen p ≤ M_K zu prüfen.

    @param discriminant Diskriminante Δ_K des Zahlkörpers (kann negativ sein)
    @param degree Grad n = [K:Q] des Zahlkörpers
    @param num_complex_pairs Anzahl r₂ der komplexen Einbettungspaare
    @return Minkowski-Schranke M_K (obere Schranke für Idealklassenvertreter)
    @lastModified 2026-03-10
    """
    n = degree
    r2 = num_complex_pairs

    # Fakultät n!
    n_factorial = math.factorial(n)

    # Formel: M_K = (n!/nⁿ) · (4/π)^r₂ · √|Δ|
    bound = (n_factorial / (n ** n)) * ((4.0 / math.pi) ** r2) * math.sqrt(abs(discriminant))

    return bound


# ============================================================
# Klasse: QuadraticForm (Binäre quadratische Form)
# ============================================================

class QuadraticForm:
    """
    @class QuadraticForm
    @brief Repräsentiert eine binäre quadratische Form f(x,y) = ax² + bxy + cy².
    @description
        Binäre quadratische Formen sind Polynome zweiter Ordnung in zwei Variablen.
        Sie spielen eine zentrale Rolle in der algebraischen Zahlentheorie.

        Wichtige Invarianten:
        - Diskriminante: D = b² - 4ac
        - Äquivalenz: f ~ g falls f durch GL₂(ℤ)-Transformation in g überführbar
        - Reduzierte Form: |b| ≤ a ≤ c (Gauss'sche Normalform)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, a: int, b: int, c: int):
        """
        @brief Initialisiert die quadratische Form f(x,y) = ax² + bxy + cy².

        @param a Koeffizient von x²
        @param b Koeffizient von xy
        @param c Koeffizient von y²
        @lastModified 2026-03-10
        """
        self.a = int(a)
        self.b = int(b)
        self.c = int(c)

    def discriminant(self) -> int:
        """
        @brief Berechnet die Diskriminante D = b² - 4ac.
        @description
            Die Diskriminante klassifiziert die quadratische Form:
            - D < 0: positiv oder negativ definit (alle Werte gleiches Vorzeichen)
            - D = 0: parabolisch (entartet)
            - D > 0: indefinit (nimmt positive und negative Werte an)

        @return Diskriminante D = b² - 4ac (ganze Zahl)
        @lastModified 2026-03-10
        """
        return self.b * self.b - 4 * self.a * self.c

    def evaluate(self, x: int, y: int) -> int:
        """
        @brief Berechnet den Wert f(x,y) = ax² + bxy + cy².

        @param x Ganzzahliger x-Wert
        @param y Ganzzahliger y-Wert
        @return Wert der quadratischen Form an (x, y)
        @lastModified 2026-03-10
        """
        return self.a * x * x + self.b * x * y + self.c * y * y

    def is_positive_definite(self) -> bool:
        """
        @brief Prüft ob die Form positiv definit ist.
        @description
            f ist positiv definit ⟺ a > 0 und D = b² - 4ac < 0.
            Dann gilt f(x,y) > 0 für alle (x,y) ≠ (0,0).

        @return True wenn positiv definit
        @lastModified 2026-03-10
        """
        return self.a > 0 and self.discriminant() < 0

    def is_reduced(self) -> bool:
        """
        @brief Prüft ob die Form Gauss-reduziert ist.
        @description
            Eine positiv definite Form ist reduziert wenn:
                |b| ≤ a ≤ c
            mit der Zusatzbedingung: b ≥ 0 falls |b| = a oder a = c.
            Jede Äquivalenzklasse enthält genau eine reduzierte Form.

        @return True wenn Gauss-reduziert
        @lastModified 2026-03-10
        """
        a, b, c = self.a, self.b, self.c
        # Grundbedingung: |b| ≤ a ≤ c
        if not (abs(b) <= a <= c):
            return False
        # Zusatzbedingung bei Gleichheit
        if abs(b) == a or a == c:
            return b >= 0
        return True

    def reduce(self) -> 'QuadraticForm':
        """
        @brief Reduziert die Form zur Gauss'schen Normalform.
        @description
            Algorithmus zur Reduktion einer positiv definiten binären quadratischen Form:
            1. Falls |b| > a: Ersetze b durch b mod (2a) im Bereich (-a, a]
            2. Falls a > c: Tausche a und c, b → -b
            3. Wiederhole bis reduziert

            Die reduzierte Form ist die kanonische Form der Äquivalenzklasse.

        @return Reduzierte quadratische Form
        @raises ValueError wenn die Form nicht positiv definit ist
        @lastModified 2026-03-10
        """
        if not self.is_positive_definite():
            raise ValueError("Reduktion nur für positiv definite Formen definiert")

        a, b, c = self.a, self.b, self.c
        max_iter = 1000

        for _ in range(max_iter):
            # Schritt 1: Reduziere |b| ≤ a durch Translation x → x + k·y
            if abs(b) > a:
                # b → b mod 2a im Bereich (-a, a]
                b_new = b % (2 * a)
                if b_new > a:
                    b_new -= 2 * a
                # Neue Form: (a, b_new, c') mit c' = (b² - D) / (4a) - (b_new² - b²)/(4a)
                # Berechne c aus Diskriminante: c = (b² - D) / (4a)
                D = b * b - 4 * a * c  # Diskriminante bleibt gleich
                c_new = (b_new * b_new - D) // (4 * a)
                b = b_new
                c = c_new

            # Schritt 2: Falls a > c, tausche
            if a > c:
                a, c = c, a
                b = -b
                continue

            # Prüfe Abbruchbedingung
            if abs(b) <= a <= c:
                # Zusatzbedingung bei Gleichheit
                if (abs(b) == a or a == c) and b < 0:
                    b = -b
                break

        return QuadraticForm(a, b, c)

    def equivalent_forms(self, max_coeff: int = 10) -> List['QuadraticForm']:
        """
        @brief Findet äquivalente Formen durch GL₂(ℤ)-Transformationen.
        @description
            Zwei Formen f und g sind äquivalent wenn eine Matrix M ∈ GL₂(ℤ) existiert mit
            f(M·(x,y)) = g(x,y). Sucht durch Enumeration von Transformationsmatrizen.

            Transformationsmatrizen mit det=±1 werden ausprobiert:
            [[p,q],[r,s]] mit |p|,|q|,|r|,|s| ≤ 2 und |ps-qr| = 1

        @param max_coeff Maximaler Absolutwert der Koeffizienten in den Formen
        @return Liste äquivalenter Formen (eindeutige Paare (a,b,c))
        @lastModified 2026-03-10
        """
        D = self.discriminant()
        found = set()
        result = []

        # Einfache Transformationen ausprobieren: SL₂(ℤ)-Generatoren T und S
        # T: (x,y) → (x+y, y), S: (x,y) → (-y, x)
        # Erzeuge Transformationen aus Kompositionen

        search_range = range(-3, 4)
        for p in search_range:
            for q in search_range:
                for r in search_range:
                    for s in search_range:
                        # Determinante muss ±1 sein (GL₂(ℤ))
                        det = p * s - q * r
                        if abs(det) != 1:
                            continue

                        # Transformierte Form: f'(x,y) = f(px+qy, rx+sy)
                        # = a(px+qy)² + b(px+qy)(rx+sy) + c(rx+sy)²
                        a2 = self.a * p * p + self.b * p * r + self.c * r * r
                        b2 = 2 * self.a * p * q + self.b * (p * s + q * r) + 2 * self.c * r * s
                        c2 = self.a * q * q + self.b * q * s + self.c * s * s

                        if max(abs(a2), abs(b2), abs(c2)) <= max_coeff:
                            key = (a2, b2, c2)
                            if key not in found:
                                found.add(key)
                                result.append(QuadraticForm(a2, b2, c2))

        return result

    def represents(self, n: int, max_search: int = 100) -> bool:
        """
        @brief Prüft ob die Form die Zahl n darstellt.
        @description
            Eine Form f repräsentiert n wenn f(x,y) = n für ganze Zahlen x, y.
            Sucht durch Enumeration von x, y ∈ [-max_search, max_search].

        @param n Darzustellende Zahl
        @param max_search Suchbereich für x und y
        @return True wenn f(x,y) = n lösbar ist
        @lastModified 2026-03-10
        """
        search_range = range(-max_search, max_search + 1)
        for x in search_range:
            for y in search_range:
                if self.evaluate(x, y) == n:
                    return True
        return False

    def __repr__(self) -> str:
        """@brief Textdarstellung der quadratischen Form."""
        return f"QuadraticForm({self.a}x² + {self.b}xy + {self.c}y²)"

    def __eq__(self, other) -> bool:
        """@brief Vergleich zweier quadratischer Formen auf Gleichheit der Koeffizienten."""
        if not isinstance(other, QuadraticForm):
            return False
        return self.a == other.a and self.b == other.b and self.c == other.c


# ============================================================
# Klasse: GaussComposition (Gauß'sche Komposition)
# ============================================================

class GaussComposition:
    """
    @class GaussComposition
    @brief Gauß'sche Komposition binärer quadratischer Formen.
    @description
        Die Gauß'sche Komposition verleiht der Menge der Äquivalenzklassen von
        quadratischen Formen mit fester Diskriminante eine Gruppenstruktur.
        Diese Gruppe ist die Klassengruppe Cl(D) des quadratischen Zahlkörpers ℚ(√D).

        Die Klassenzahl h(D) = |Cl(D)| misst die "Komplexität" des Zahlkörpers.
        h(D) = 1 ⟺ ℤ[√D] ist faktoriell (eindeutige Primfaktorzerlegung).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def reduced_forms(discriminant: int) -> List[QuadraticForm]:
        """
        @brief Findet alle reduzierten binären quadratischen Formen mit gegebener Diskriminante.
        @description
            Eine reduzierte positiv definite Form (a,b,c) mit D = b²-4ac < 0 erfüllt:
                |b| ≤ a ≤ c und (|b| = a oder a = c → b ≥ 0)
            Suchbereich: a ≤ √(|D|/3) (da 4ac = b² - D ≥ -D, also ac ≥ -D/4, und a ≤ c → a² ≤ -D/3)

        @param discriminant Diskriminante D (muss negativ sein für positiv definite Formen)
        @return Liste aller reduzierten quadratischen Formen
        @lastModified 2026-03-10
        """
        D = discriminant
        forms = []

        if D >= 0:
            # Für indefinite Formen: vereinfachte Suche
            # Reduzierte indefinite Formen haben andere Bedingungen
            # Hier: einfache Enumeration für kleine D
            max_a = int(math.sqrt(abs(D))) + 2 if D != 0 else 1
            for a in range(-max_a, max_a + 1):
                if a == 0:
                    continue
                for b in range(-abs(a), abs(a) + 1):
                    # Berechne c aus D = b² - 4ac
                    if 4 * a != 0 and (b * b - D) % (4 * a) == 0:
                        c = (b * b - D) // (4 * a)
                        f = QuadraticForm(a, b, c)
                        if f.discriminant() == D:
                            forms.append(f)
            return forms

        # D < 0: Positiv definite Formen
        # Suchbereich für a: 1 ≤ a ≤ √(|D|/3)
        max_a = int(math.sqrt(-D / 3)) + 1

        for a in range(1, max_a + 1):
            # b läuft von -a bis a
            for b in range(-a, a + 1):
                # Berechne c aus D = b² - 4ac → c = (b² - D) / (4a)
                numerator = b * b - D
                if numerator <= 0 or numerator % (4 * a) != 0:
                    continue
                c = numerator // (4 * a)

                # Prüfe Reduktionsbedingungen: a ≤ c
                if c < a:
                    continue

                # Zusatzbedingung: bei |b|=a oder a=c muss b ≥ 0
                if (abs(b) == a or a == c) and b < 0:
                    continue

                f = QuadraticForm(a, b, c)
                if f.discriminant() == D:
                    forms.append(f)

        return forms

    @staticmethod
    def class_number(discriminant: int) -> int:
        """
        @brief Berechnet die Klassenzahl h(D) der Formen mit Diskriminante D.
        @description
            Die Klassenzahl h(D) ist die Anzahl der Äquivalenzklassen von primitiven
            binären quadratischen Formen mit Diskriminante D.
            Für D < 0: h(D) = Anzahl der reduzierten Formen.

            Bekannte Werte (Gauss):
            h(-3) = 1, h(-4) = 1, h(-7) = 1, h(-8) = 1, h(-11) = 1,
            h(-12) = 1, h(-16) = 1, h(-27) = 1, h(-28) = 1, h(-43) = 1,
            h(-67) = 1, h(-163) = 1 (Heegner-Zahlen).

        @param discriminant Diskriminante D
        @return Klassenzahl h(D) ≥ 1
        @lastModified 2026-03-10
        """
        forms = GaussComposition.reduced_forms(discriminant)
        return max(1, len(forms))

    @staticmethod
    def compose(form1: QuadraticForm, form2: QuadraticForm) -> QuadraticForm:
        """
        @brief Gauß'sche Komposition zweier quadratischer Formen (Dirichlet-Methode).
        @description
            Die Dirichlet-Komposition zweier Formen f₁ = (a₁, b₁, c₁) und f₂ = (a₂, b₂, c₂)
            mit gleicher Diskriminante D:

            1. Berechne n = gcd(a₁, a₂, (b₁+b₂)/2)
            2. Setze a = a₁·a₂/n², b = b₁ (oder b₂ nach CRT), c = (b²-D)/(4a)

            Vereinfacht (für primitive Formen mit gcd(a₁,a₂)=1):
            - a₃ = a₁·a₂
            - b₃ durch Chinesischen Restsatz aus b₁ (mod 2a₁) und b₂ (mod 2a₂)

        @param form1 Erste quadratische Form
        @param form2 Zweite quadratische Form
        @return Zusammengesetzte Form (reduziert)
        @raises ValueError wenn Diskriminanten verschieden sind
        @lastModified 2026-03-10
        """
        D1 = form1.discriminant()
        D2 = form2.discriminant()

        if D1 != D2:
            raise ValueError(f"Diskriminanten müssen gleich sein: D1={D1}, D2={D2}")

        a1, b1, c1 = form1.a, form1.b, form1.c
        a2, b2, c2 = form2.a, form2.b, form2.c
        D = D1

        # Dirichlet-Komposition: vereinfachte Version
        # n = gcd(a1, a2, (b1+b2)//2)
        if (b1 + b2) % 2 == 0:
            n = math.gcd(math.gcd(a1, a2), (b1 + b2) // 2)
        else:
            n = math.gcd(a1, a2)

        # Neues a: a₃ = a₁·a₂/n²
        a3 = (a1 * a2) // (n * n)

        # Neues b via Chinesischem Restsatz:
        # b₃ ≡ b₁ (mod 2a₁/n) und b₃ ≡ b₂ (mod 2a₂/n)
        # Einfache Näherung: Mittelwert falls gleich, sonst CRT
        if b1 == b2:
            b3 = b1
        else:
            # Chinesischer Restsatz für b
            # b₃ ≡ b₁ (mod 2·a₁) und b₃ ≡ b₂ (mod 2·a₂)
            mod1 = 2 * a1
            mod2 = 2 * a2
            g = math.gcd(mod1, mod2)

            if (b2 - b1) % g != 0:
                # Kein gemeinsames b möglich: verwende b1
                b3 = b1
            else:
                # CRT-Lösung
                lcm_val = mod1 * mod2 // g
                # Erweiterte Euler: mod1 * x ≡ 1 (mod mod2/g)
                m2g = mod2 // g
                try:
                    inv = pow(mod1 // g, -1, m2g)
                    b3 = (b1 + mod1 * ((b2 - b1) // g * inv % m2g)) % lcm_val
                    # Wähle b3 im Bereich (-a3, a3]
                    while b3 > a3:
                        b3 -= 2 * a3
                    while b3 <= -a3:
                        b3 += 2 * a3
                except (ValueError, ZeroDivisionError):
                    b3 = b1

        # Neues c aus Diskriminante: c₃ = (b₃² - D) / (4·a₃)
        numerator = b3 * b3 - D
        if a3 != 0 and numerator % (4 * a3) == 0:
            c3 = numerator // (4 * a3)
        else:
            # Fallback: Näherung
            c3 = (b3 * b3 - D + 2 * a3) // (4 * a3) if a3 != 0 else 1

        composed = QuadraticForm(a3, b3, c3)

        # Zur reduzierten Form reduzieren
        if composed.is_positive_definite():
            return composed.reduce()
        return composed


# ============================================================
# Klasse: SpherePackingTheory (Kugelpackungstheorie)
# ============================================================

class SpherePackingTheory:
    """
    @class SpherePackingTheory
    @brief Kugelpackungstheorie: Kusszahlen, Packungsdichten, bekannte optimale Packungen.
    @description
        Die Kugelpackungstheorie untersucht, wie viele gleich große Kugeln im ℝⁿ
        raumeffizient angeordnet werden können.

        Bekannte Ergebnisse:
        - dim=1: triviale Packung, Dichte = 1
        - dim=2: hexagonale Packung (Thue 1892), Dichte = π/(2√3) ≈ 0.9069
        - dim=3: FCC/HCP-Packung (Kepler-Vermutung, Hales 1998), Dichte = π/(3√2) ≈ 0.7405
        - dim=8: E₈-Gitter (optimale Packung, Viazovska 2016)
        - dim=24: Leech-Gitter (optimale Packung, Viazovska et al. 2016)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    # Bekannte Kusszahlen (Anzahl Kugeln, die eine Kugel berühren können)
    KISSING_NUMBERS = {
        1: 2,      # Trivial: 2 Nachbarn auf der Geraden
        2: 6,      # Hexagonal: 6 Nachbarn in der Ebene
        3: 12,     # FCC/HCP: 12 Nachbarn im Raum (bewiesen: Schütte & van der Waerden, 1953)
        4: 24,     # D₄-Gitter (bewiesen: Musin, 2003)
        8: 240,    # E₈-Gitter
        24: 196560  # Leech-Gitter
    }

    @staticmethod
    def kissing_number(dim: int) -> Optional[int]:
        """
        @brief Gibt die bekannte Kusszahl für die Dimension dim zurück.
        @description
            Die Kusszahl τ(n) ist die maximale Anzahl nicht-überlappender Einheitskugeln,
            die eine zentrale Einheitskugel im ℝⁿ berühren können.

            Bekannte exakte Werte: dim ∈ {1, 2, 3, 4, 8, 24}.
            Für andere Dimensionen sind nur Schranken bekannt.

        @param dim Dimension (1, 2, 3, 4, 8 oder 24 für exakte Werte)
        @return Kusszahl oder None wenn nicht bekannt
        @lastModified 2026-03-10
        """
        return SpherePackingTheory.KISSING_NUMBERS.get(dim, None)

    @staticmethod
    def packing_density_upper_bound(dim: int) -> float:
        """
        @brief Blichfeldt-Schranke für die Packungsdichte in Dimension n.
        @description
            Die Blichfeldt-Schranke (1929) gibt die triviale obere Schranke δ ≤ 1.
            Für die verbesserte Kabatiansky-Levenshtein-Schranke (1978):
                ln δ(n) ≤ -0.5990 · n + o(n)
            Diese Version gibt die triviale obere Schranke 1.0 zurück,
            da die exakten Formeln komplex sind.

        @param dim Dimension n ≥ 1
        @return Obere Schranke für die Packungsdichte (≤ 1.0)
        @lastModified 2026-03-10
        """
        if dim <= 0:
            raise ValueError(f"Dimension muss positiv sein, aber dim={dim}")

        # Triviale Blichfeldt-Schranke: δ ≤ 1
        # Verbesserung: Für dim ≥ 2 gilt die Kabatiansky-Levenshtein-Schranke
        if dim == 1:
            return 1.0

        # Kabatiansky-Levenshtein (asymptotisch):
        # δ(n) ≤ 2^{-0.5990·n}
        kl_bound = 2 ** (-0.5990 * dim)

        # Aber die triviale Schranke δ ≤ 1 ist nie verletzt
        return min(1.0, kl_bound)

    @staticmethod
    def hexagonal_packing_density() -> float:
        """
        @brief Optimale Kreispackungsdichte im ℝ² (hexagonal).
        @description
            Die hexagonale Packung ist die dichteste Kreispackung in der Ebene.
            Dichte = π / (2√3) ≈ 0.90690
            Bewiesen von Thue (1892) und Fejes Tóth (1940).

        @return Hexagonale Packungsdichte π/(2√3)
        @lastModified 2026-03-10
        """
        return math.pi / (2 * math.sqrt(3))

    @staticmethod
    def fcc_packing_density() -> float:
        """
        @brief Optimale Kugelpackungsdichte im ℝ³ (FCC/HCP, Kepler-Problem).
        @description
            Die FCC (Face-Centered Cubic) und HCP (Hexagonal Close-Packing)
            erreichen die optimale Dichte π/(3√2) ≈ 0.74048.
            Keplerscher Vermutung (1611), bewiesen von Hales (1998/2014).

        @return FCC-Packungsdichte π/(3√2)
        @lastModified 2026-03-10
        """
        return math.pi / (3 * math.sqrt(2))

    @staticmethod
    def leech_lattice_info() -> dict:
        """
        @brief Informationen über das Leech-Gitter (Dimension 24).
        @description
            Das Leech-Gitter Λ₂₄ ist das bekannteste 24-dimensionale Gitter mit:
            - Dimension: 24
            - Kusszahl: 196560 (maximal in dim=24, Viazovska et al. 2016)
            - Packungsdichte: π¹²/(12!) ≈ 0.00193 (optimal in dim=24)
            - Diskriminante: 1 (unimodulär)
            - Automorphismengruppe: 2·Co₁ (Conway-Gruppe, |Co₁| ≈ 4·10^18)
            - Verbindung zur Monster-Gruppe via Monstrous Moonshine

        @return Dictionary mit Eigenschaften des Leech-Gitters
        @lastModified 2026-03-10
        """
        # Packungsdichte des Leech-Gitters: π¹²/12!
        density = (math.pi ** 12) / math.factorial(12)

        return {
            "dimension": 24,
            "kissing_number": 196560,
            "packing_density": density,
            "determinant": 1,
            "is_unimodular": True,
            "is_even": True,
            "automorphism_group_order": 8315553613086720000,  # |2·Co₁|
            "description": "Leech-Gitter Λ₂₄: Optimale Kugelpackung in dim=24 (Viazovska et al. 2016)",
            "min_norm": 4  # Minimale Norm der kürzesten Vektoren
        }


# ============================================================
# Standalone-Funktionen
# ============================================================

def hermite_constant(n: int) -> float:
    """
    @brief Gibt die Hermite-Konstante γₙ für Dimension n zurück.
    @description
        Die Hermite-Konstante γₙ ist definiert als:
            γₙ = max_{Λ} (λ₁(Λ)² / det(Λ)^{2/n})
        d.h. der maximale Quotient aus dem quadrierten ersten Sukzessivminimum
        und der auf 2/n potenzierten Gitterdeterminante.

        Bekannte exakte Werte:
        n:  1    2    3    4    5    6    7    8
        γₙ: 1   4/3  2   4   8   64/3  64   2

        Formel: γₙ = (2/π)·Γ(n/2+1)^{2/n} (Minkowski-Schranke, asymptotisch)

    @param n Dimension (1 ≤ n ≤ 8 für exakte Werte; sonst Näherung)
    @return Hermite-Konstante γₙ
    @raises ValueError wenn n < 1
    @lastModified 2026-03-10
    """
    if n < 1:
        raise ValueError(f"Dimension n muss ≥ 1 sein, aber n={n}")

    # Exakte bekannte Werte
    exact_values = {
        1: 1.0,
        2: 4.0 / 3.0,
        3: 2.0,
        4: 4.0,
        5: 8.0,
        6: 64.0 / 3.0,
        7: 64.0,
        8: 2.0  # E₈: γ₈ = 2
    }

    if n in exact_values:
        return exact_values[n]

    # Näherung für n > 8: Minkowski-Schranke
    # γₙ ≤ (2/π)·Γ(n/2+1)^{2/n}
    gamma_val = math.gamma(n / 2 + 1)
    approx = (2.0 / math.pi) * (gamma_val ** (2.0 / n))
    return approx


def lattice_from_quadratic_form(a: int, b: int, c: int) -> 'Lattice':
    """
    @brief Konstruiert ein Gitter aus einer binären quadratischen Form.
    @description
        Jede positive definite binäre quadratische Form f(x,y) = ax² + bxy + cy²
        definiert ein Gitter im ℝ² mit Gram-Matrix:
            G = [[a, b/2], [b/2, c]]
        Die Gitterbasis B wird via Cholesky-Zerlegung G = Bᵀ·B berechnet:
            B = [[√a, 0], [b/(2√a), √(c - b²/(4a))]]

    @param a Koeffizient von x² (muss > 0 sein)
    @param b Koeffizient von xy
    @param c Koeffizient von y²
    @return Lattice-Objekt mit der entsprechenden Gitterbasis
    @raises ValueError wenn die Form nicht positiv definit ist
    @lastModified 2026-03-10
    """
    form = QuadraticForm(a, b, c)
    if not form.is_positive_definite():
        raise ValueError("Gitter aus quadratischer Form nur für positiv definite Formen definiert")

    # Gram-Matrix: G = [[a, b/2], [b/2, c]]
    G = np.array([[float(a), float(b) / 2], [float(b) / 2, float(c)]])

    # Cholesky-Zerlegung: G = Lᵀ·L (L ist untere Dreiecksmatrix)
    # B₁₁ = √a, B₂₁ = b/(2√a), B₂₂ = √(c - b²/(4a))
    b11 = math.sqrt(a)
    b21 = b / (2.0 * b11)
    b22_sq = c - b21 * b21
    if b22_sq < 0:
        raise ValueError(f"Ungültige quadratische Form: {a}x² + {b}xy + {c}y² ist nicht positiv definit")
    b22 = math.sqrt(b22_sq)

    # Basis: Zeilen sind Gittervektoren
    basis = np.array([[b11, 0.0], [b21, b22]])

    return Lattice(basis)


def successive_minima_bound(det_lattice: float, dim: int) -> float:
    """
    @brief Berechnet die Minkowski-Schranke für das Produkt der Sukzessivminima.
    @description
        Minkowskis zweiter Satz (1896):
            λ₁·λ₂·...·λₙ ≤ γₙ^{n/2} · det(Λ)^{1/n} · √n
        Vereinfacht (Minkowski-Schranke für λ₁):
            λ₁(Λ) ≤ √(γₙ) · det(Λ)^{1/n}

        Diese Funktion gibt die Schranke für das geometrische Mittel der Sukzessivminima zurück:
            (λ₁·...·λₙ)^{1/n} ≤ √(γₙ) · det(Λ)^{1/n}

    @param det_lattice Gitterdeterminante |det(Λ)| (> 0)
    @param dim Dimension n des Gitters (≥ 1)
    @return Obere Schranke für (λ₁·...·λₙ)^{1/n}
    @raises ValueError wenn det_lattice ≤ 0 oder dim < 1
    @lastModified 2026-03-10
    """
    if det_lattice <= 0:
        raise ValueError(f"Gitterdeterminante muss positiv sein, aber det={det_lattice}")
    if dim < 1:
        raise ValueError(f"Dimension muss ≥ 1 sein, aber dim={dim}")

    # Hermite-Konstante γₙ
    gamma_n = hermite_constant(dim)

    # Schranke: √(γₙ) · det(Λ)^{1/n}
    bound = math.sqrt(gamma_n) * (abs(det_lattice) ** (1.0 / dim))

    return bound
