"""
@file functional_analysis.py
@brief Funktionalanalysis: normierte Räume, Banach-Räume, Hilbert-Räume,
       lineare Operatoren, Spektraltheorie, Funktionalräume und Fixpunktsätze.
@description
    Dieses Modul implementiert die zentralen Konzepte der Funktionalanalysis:

    - Normierte Räume (NormedSpace): Vektorräume mit einer Norm ‖·‖
    - Banach-Räume (BanachSpace): vollständige normierte Räume
    - Hilbert-Räume (HilbertSpace): Banach-Räume mit Skalarprodukt
    - Lineare Operatoren (LinearOperator): beschränkte lineare Abbildungen
    - Spektraltheorie: Spektrum, Spektralradius, Spektralsatz
    - Fredholm-Alternative: Lösungsbedingung für Ax=b
    - Funktionenräume: C[a,b], Sobolev-Räume W^{k,p}, L²[a,b]
    - Klassische Sätze: Hahn-Banach, Offene-Abbildung, Gleichmäßige Beschränktheit
    - Fixpunktsätze: Banach, Schauder

    Die Funktionalanalysis bildet die theoretische Grundlage für:
    - Partielle Differentialgleichungen (PDE)
    - Quantenmechanik (Hilbert-Raum-Formalismus)
    - Numerische Analysis und Approximationstheorie

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import numpy as np
import math
from typing import Callable


# ============================================================
# Normierte Räume
# ============================================================

class NormedSpace:
    """
    @brief Normierter Raum (V, ‖·‖).
    @description
        Ein normierter Raum ist ein Vektorraum V über ℝ (oder ℂ) zusammen mit
        einer Norm ‖·‖ : V → ℝ≥0, die folgende Axiome erfüllt:

        1. Positivität:       ‖v‖ ≥ 0 für alle v ∈ V
        2. Definitheit:       ‖v‖ = 0  ⟺  v = 0
        3. Homogenität:       ‖αv‖ = |α|·‖v‖ für alle α ∈ ℝ, v ∈ V
        4. Dreiecksungleichung: ‖u+v‖ ≤ ‖u‖ + ‖v‖

        Jeder normierte Raum ist automatisch ein metrischer Raum mit d(u,v) = ‖u-v‖.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, vectors: list, norm_func: Callable):
        """
        @brief Initialisiert den normierten Raum.
        @param vectors  Liste von Vektoren (jeweils als list[float])
        @param norm_func Callable, das einen Vektor auf seinen Normwert (float) abbildet
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Vektoren als numpy-Arrays speichern für einfache Berechnung
        self.vectors = [np.array(v, dtype=float) for v in vectors]
        # Norm-Funktion als Attribut speichern
        self._norm_func = norm_func

    def norm(self, v) -> float:
        """
        @brief Berechnet die Norm ‖v‖ eines Vektors.
        @description
            Wendet die im Konstruktor übergebene Normfunktion an.
            Mathematisch: ‖v‖ ∈ ℝ≥0

        @param v Vektor (list oder np.ndarray)
        @return Normwert als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return float(self._norm_func(np.array(v, dtype=float)))

    def distance(self, u, v) -> float:
        """
        @brief Berechnet den Abstand d(u,v) = ‖u-v‖.
        @description
            Jede Norm induziert eine Metrik via d(u,v) = ‖u-v‖.
            Diese Metrik macht den Raum zu einem metrischen Raum.

        @param u Erster Vektor
        @param v Zweiter Vektor
        @return Abstand als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        u_arr = np.array(u, dtype=float)
        v_arr = np.array(v, dtype=float)
        # Metrik d(u,v) = ‖u - v‖
        return self.norm(u_arr - v_arr)

    def is_normed_space(self) -> bool:
        """
        @brief Prüft ob die gespeicherten Vektoren die Norm-Axiome erfüllen.
        @description
            Überprüft anhand der gespeicherten Vektoren:
            - Positivität ‖v‖ ≥ 0
            - Dreiecksungleichung ‖u+v‖ ≤ ‖u‖+‖v‖ (mit numerischer Toleranz)
            - Homogenität ‖αv‖ = |α|·‖v‖

        @return True wenn alle Axiome für die gespeicherten Vektoren gelten
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        tol = 1e-10  # numerische Toleranz
        for v in self.vectors:
            # Axiom 1: Positivität
            if self.norm(v) < -tol:
                return False
            # Axiom 2: Definitheit (Nullvektor hat Norm 0)
            zero = np.zeros_like(v)
            if abs(self.norm(zero)) > tol:
                return False
        # Dreiecksungleichung und Homogenität an Paaren prüfen
        for i, u in enumerate(self.vectors):
            for v in self.vectors[i:]:
                # Dreiecksungleichung: ‖u+v‖ ≤ ‖u‖+‖v‖
                if self.norm(u + v) > self.norm(u) + self.norm(v) + tol:
                    return False
            # Homogenität: ‖2v‖ = 2·‖v‖
            if abs(self.norm(2.0 * u) - 2.0 * self.norm(u)) > tol:
                return False
        return True

    def unit_ball(self, n_points: int = 100) -> list:
        """
        @brief Erzeugt Punkte auf dem Rand der Einheitskugel B(0,1).
        @description
            Die Einheitskugel (Einheitsball) ist definiert als:
                B(0,1) = {v ∈ V : ‖v‖ ≤ 1}
            Hier werden Punkte auf dem Rand ‖v‖ = 1 im ℝ² erzeugt.

        @param n_points Anzahl der Punkte auf dem Rand
        @return Liste von Punkten [x, y] mit ‖[x,y]‖ = 1
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        result = []
        # Gleichmäßig verteilte Winkel im Einheitskreis
        for i in range(n_points):
            theta = 2.0 * math.pi * i / n_points
            # Punkt auf Einheitskreis im euklidischen Sinne
            p = np.array([math.cos(theta), math.sin(theta)])
            # Auf Einheitssphäre der gegebenen Norm normieren
            n = self.norm(p)
            if n > 1e-15:
                result.append((p / n).tolist())
        return result


class BanachSpace(NormedSpace):
    """
    @brief Banach-Raum: vollständiger normierter Raum.
    @description
        Ein Banach-Raum ist ein normierter Raum, in dem jede Cauchy-Folge konvergiert.

        Cauchy-Folge: (xₙ) mit ‖xₘ-xₙ‖ → 0 für m,n → ∞
        Vollständigkeit: jede Cauchy-Folge hat einen Grenzwert in V.

        Klassische Beispiele:
        - (ℝⁿ, ‖·‖_p) für alle 1 ≤ p ≤ ∞
        - ℓᵖ = {(xₙ) : Σ|xₙ|ᵖ < ∞} mit ‖·‖_p
        - Lᵖ[a,b] = {f : ∫|f|ᵖ < ∞} mit ‖·‖_p
        - C[a,b] mit Supremumsnorm ‖·‖_∞

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def cauchy_sequence_test(self, sequence: list, tol: float = 1e-10) -> bool:
        """
        @brief Prüft ob eine Folge eine Cauchy-Folge ist.
        @description
            Eine Folge (xₙ) heißt Cauchy-Folge, wenn gilt:
                ∀ε > 0 ∃N : ∀m,n > N : ‖xₘ - xₙ‖ < ε

            Praktisch: prüft ob aufeinanderfolgende Glieder ab Mitte der Folge
            hinreichend nahe beieinander liegen.

        @param sequence Liste von Vektoren (als list oder np.ndarray)
        @param tol      Toleranzgrenze für "hinreichend klein"
        @return True wenn die Folge eine Cauchy-Folge ist
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        n = len(sequence)
        if n < 2:
            return True
        # Prüfe ob ab der Mitte der Folge alle Abstände < tol
        half = n // 2
        for i in range(half, n):
            for j in range(half, n):
                xi = np.array(sequence[i], dtype=float)
                xj = np.array(sequence[j], dtype=float)
                if self.distance(xi, xj) >= tol:
                    return False
        return True

    def is_separable(self) -> bool:
        """
        @brief Gibt an ob der Raum (basierend auf den Vektoren) separabel ist.
        @description
            Ein normierter Raum heißt separabel, wenn er eine abzählbare dichte
            Teilmenge enthält. Alle endlich-dimensionalen normierten Räume sind
            separabel. ℓᵖ (1≤p<∞) und Lᵖ[a,b] sind separabel; ℓ∞ nicht.

            Diese Implementierung gibt True zurück, da alle hier betrachteten
            endlich-dimensionalen Räume separabel sind.

        @return True (alle hier verwendeten Räume sind endlich-dimensional und damit separabel)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Endlich-dimensionale Räume sind immer separabel
        return True

    def lp_norm(self, v: list, p: float) -> float:
        """
        @brief Berechnet die ℓᵖ-Norm eines Vektors.
        @description
            Die ℓᵖ-Norm ist definiert als:
                ‖v‖_p = (Σᵢ |vᵢ|ᵖ)^{1/p}  für 1 ≤ p < ∞
                ‖v‖_∞ = max_i |vᵢ|          für p = ∞

            Wichtige Spezialfälle:
            - p=1: ‖v‖₁ = Σ|vᵢ| (Betragssumme / Mannigfaltigkeitsnorm)
            - p=2: ‖v‖₂ = √(Σvᵢ²) (euklidische Norm)
            - p=∞: ‖v‖_∞ = max|vᵢ| (Tschebyschew-Norm)

        @param v Vektor als Liste von Zahlen
        @param p Norm-Parameter (1 ≤ p ≤ ∞, float('inf') für ∞-Norm)
        @return ℓᵖ-Normwert als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        v_arr = np.array(v, dtype=float)
        if p == float('inf'):
            # ℓ∞-Norm: Maximum der Beträge
            return float(np.max(np.abs(v_arr)))
        if p <= 0:
            raise ValueError(f"p muss >= 1 sein, erhalten: {p}")
        # ℓᵖ-Norm: (Σ|vᵢ|ᵖ)^{1/p}
        return float(np.sum(np.abs(v_arr) ** p) ** (1.0 / p))


class HilbertSpace(BanachSpace):
    """
    @brief Hilbert-Raum: vollständiger normierter Raum mit Skalarprodukt.
    @description
        Ein Hilbert-Raum ist ein Banach-Raum H, auf dem ein Skalarprodukt
        ⟨·,·⟩ : H × H → ℝ (oder ℂ) definiert ist, das die Norm induziert:
            ‖v‖ = √⟨v,v⟩

        Skalarprodukt-Axiome:
        1. Linearität:      ⟨αu+βv, w⟩ = α⟨u,w⟩ + β⟨v,w⟩
        2. Symmetrie:       ⟨u,v⟩ = ⟨v,u⟩ (reell) bzw. konjugiert (komplex)
        3. Positivität:     ⟨v,v⟩ ≥ 0
        4. Definitheit:     ⟨v,v⟩ = 0 ⟺ v = 0

        Wichtige Ungleichungen:
        - Cauchy-Schwarz: |⟨u,v⟩| ≤ ‖u‖·‖v‖
        - Parallelogramm-Gesetz: ‖u+v‖² + ‖u-v‖² = 2(‖u‖²+‖v‖²)

        Klassische Beispiele:
        - ℝⁿ mit ⟨u,v⟩ = Σuᵢvᵢ (Standard-Skalarprodukt)
        - ℓ² = {(xₙ) : Σ|xₙ|² < ∞} mit ⟨x,y⟩ = Σxₙyₙ
        - L²[a,b] mit ⟨f,g⟩ = ∫ₐᵇ f(x)g(x)dx

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, vectors: list, inner_product: Callable):
        """
        @brief Initialisiert den Hilbert-Raum mit gegebenem Skalarprodukt.
        @param vectors       Liste von Vektoren
        @param inner_product Callable ⟨u,v⟩: (array, array) → float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Skalarprodukt speichern
        self._inner_product = inner_product
        # Norm aus Skalarprodukt ableiten: ‖v‖ = √⟨v,v⟩
        def ip_norm(v):
            return math.sqrt(max(0.0, inner_product(v, v)))
        # Elternklasse mit induzierter Norm initialisieren
        super().__init__(vectors, ip_norm)

    def inner_product(self, u, v) -> float:
        """
        @brief Berechnet das Skalarprodukt ⟨u,v⟩.
        @description
            Das Skalarprodukt ist die zentrale Operation des Hilbert-Raums.
            Es erlaubt Konzepte wie Orthogonalität und Winkel.

            Cauchy-Schwarz-Ungleichung: |⟨u,v⟩| ≤ ‖u‖·‖v‖

        @param u Erster Vektor
        @param v Zweiter Vektor
        @return Skalarprodukt ⟨u,v⟩ als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return float(self._inner_product(np.array(u, dtype=float),
                                         np.array(v, dtype=float)))

    def is_orthogonal(self, u, v) -> bool:
        """
        @brief Prüft ob u und v orthogonal sind (⟨u,v⟩ = 0).
        @description
            Zwei Vektoren u, v ∈ H heißen orthogonal (u ⊥ v), wenn:
                ⟨u,v⟩ = 0

            Orthogonalität verallgemeinert das geometrische Konzept
            der Rechtwinkligkeit auf beliebige Hilbert-Räume.

        @param u Erster Vektor
        @param v Zweiter Vektor
        @return True wenn ⟨u,v⟩ ≈ 0
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return abs(self.inner_product(u, v)) < 1e-10

    def gram_schmidt_orthonormalize(self, basis: list) -> list:
        """
        @brief Gram-Schmidt-Orthonormalisierung einer Basis.
        @description
            Das Gram-Schmidt-Verfahren erzeugt aus einer linear unabhängigen
            Menge {v₁,...,vₙ} eine Orthonormalbasis {e₁,...,eₙ} via:

                u₁ = v₁,                    e₁ = u₁/‖u₁‖
                uₖ = vₖ - Σⱼ<ₖ ⟨vₖ,eⱼ⟩eⱼ,  eₖ = uₖ/‖uₖ‖

            Ergebnis: ⟨eᵢ,eⱼ⟩ = δᵢⱼ (Kronecker-Delta)

        @param basis Liste von linear unabhängigen Vektoren
        @return Liste orthonormaler Vektoren (np.ndarray)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        orthonormal = []
        for v in basis:
            v_arr = np.array(v, dtype=float)
            # Orthogonalprojektion auf alle bisherigen ONB-Vektoren subtrahieren
            u = v_arr.copy()
            for e in orthonormal:
                # u = v - Σ ⟨v,e⟩·e
                u = u - self.inner_product(u, e) * np.array(e)
            # Normieren: e = u/‖u‖
            norm_u = self.norm(u)
            if norm_u > 1e-12:
                orthonormal.append((u / norm_u).tolist())
        return orthonormal

    def projection(self, v, onto: list) -> list:
        """
        @brief Orthogonalprojektion von v auf einen Unterraum.
        @description
            Sei U = span{e₁,...,eₙ} ein Unterraum mit ONB {e₁,...,eₙ}.
            Die Orthogonalprojektion von v auf U ist:
                P_U(v) = Σᵢ ⟨v,eᵢ⟩·eᵢ

            Eigenschaften:
            - P_U(v) ∈ U ist der eindeutige nächste Punkt zu v in U
            - v - P_U(v) ⊥ U (Orthogonalitätsprinzip)
            - ‖v - P_U(v)‖ = dist(v, U) (minimaler Abstand)

        @param v    Zu projizierender Vektor
        @param onto Liste von Basisvektoren des Unterraums (nicht notwendig ONB)
        @return Projektionsvektor als list
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        v_arr = np.array(v, dtype=float)
        # Erst Gram-Schmidt anwenden, falls Basis nicht orthonormal ist
        onb = self.gram_schmidt_orthonormalize(onto)
        # Projektion: P_U(v) = Σ ⟨v,eᵢ⟩·eᵢ
        proj = np.zeros_like(v_arr)
        for e in onb:
            e_arr = np.array(e, dtype=float)
            proj = proj + self.inner_product(v_arr, e_arr) * e_arr
        return proj.tolist()

    def pythagoras(self, u, v) -> bool:
        """
        @brief Prüft den Satz des Pythagoras für orthogonale Vektoren.
        @description
            Für u ⊥ v gilt der verallgemeinerte Satz des Pythagoras:
                ‖u + v‖² = ‖u‖² + ‖v‖²

            Dies ist äquivalent zu u ⊥ v, d.h. ⟨u,v⟩ = 0.
            Der Beweis folgt direkt:
                ‖u+v‖² = ⟨u+v,u+v⟩ = ‖u‖²+2⟨u,v⟩+‖v‖² = ‖u‖²+‖v‖²

        @param u Erster Vektor
        @param v Zweiter Vektor
        @return True wenn ‖u+v‖² = ‖u‖² + ‖v‖² (d.h. u ⊥ v)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        u_arr = np.array(u, dtype=float)
        v_arr = np.array(v, dtype=float)
        lhs = self.norm(u_arr + v_arr) ** 2
        rhs = self.norm(u_arr) ** 2 + self.norm(v_arr) ** 2
        return abs(lhs - rhs) < 1e-8

    def riesz_representation(self, functional: Callable, basis: list) -> list:
        """
        @brief Berechnet den Riesz-Repräsentanten eines linearen Funktionals.
        @description
            Riesz-Darstellungssatz (Riesz, 1907):
            Zu jedem stetigen linearen Funktional f ∈ H* existiert genau ein
            y ∈ H mit:
                f(x) = ⟨x, y⟩  für alle x ∈ H

            Der Repräsentant y wird über die ONB {e₁,...,eₙ} berechnet:
                y = Σᵢ f(eᵢ)·eᵢ  (da ⟨y,eᵢ⟩ = f(eᵢ) nötig)

        @param functional Lineares Funktional f: H → ℝ
        @param basis      Basis des Raums (wird orthonormalisiert)
        @return Riesz-Repräsentant y als list
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Gram-Schmidt: ONB erzeugen
        onb = self.gram_schmidt_orthonormalize(basis)
        # y = Σ f(eᵢ)·eᵢ (denn ⟨y,eᵢ⟩ = f(eᵢ))
        if not onb:
            return []
        dim = len(onb[0])
        y = np.zeros(dim, dtype=float)
        for e in onb:
            e_arr = np.array(e, dtype=float)
            # Funktionalwert am Basisvektor
            coeff = float(functional(e_arr))
            y = y + coeff * e_arr
        return y.tolist()

    def parseval_identity(self, v: list, onb: list) -> bool:
        """
        @brief Prüft die Parseval-Gleichung: ‖v‖² = Σ|⟨v,eᵢ⟩|².
        @description
            Die Parseval-Gleichung (auch Parseval-Identität) gilt für
            eine vollständige Orthonormalbasis (ONB) {eᵢ} eines Hilbert-Raums:
                ‖v‖² = Σᵢ |⟨v,eᵢ⟩|²

            Sie ist das Analogon zum Satz des Pythagoras in unendlicher Dimension.
            In endlicher Dimension gilt sie genau dann, wenn {eᵢ} eine ONB ist.

            Beziehung zur Fourier-Reihe: Die Fourier-Koeffizienten ⟨v,eᵢ⟩ sind
            optimal in dem Sinne, dass die Parseval-Gleichung minimalen Fehler liefert.

        @param v   Vektor
        @param onb Orthonormalbasis (Liste von Vektoren)
        @return True wenn ‖v‖² ≈ Σ|⟨v,eᵢ⟩|²
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        v_arr = np.array(v, dtype=float)
        norm_sq = self.norm(v_arr) ** 2
        # Parseval-Summe: Σ|⟨v,eᵢ⟩|²
        parseval_sum = sum(self.inner_product(v_arr, np.array(e, dtype=float)) ** 2
                           for e in onb)
        return abs(norm_sq - parseval_sum) < 1e-8

    def bessel_inequality(self, v: list, orthonormal_set: list) -> bool:
        """
        @brief Prüft die Bessel-Ungleichung: Σ|⟨v,eᵢ⟩|² ≤ ‖v‖².
        @description
            Die Bessel-Ungleichung gilt für jede orthonormale Menge {e₁,...,eₙ}
            (nicht notwendig vollständig) in einem Hilbert-Raum:
                Σᵢ |⟨v,eᵢ⟩|² ≤ ‖v‖²

            Gleichheit gilt genau dann, wenn {eᵢ} eine vollständige ONB ist
            (dann ist es die Parseval-Gleichung).

            Interpretation: Die Teilsumme der Fourier-Koeffizienten ist durch
            die Norm des Vektors beschränkt.

        @param v               Vektor
        @param orthonormal_set Orthonormale Menge (Liste von Vektoren)
        @return True wenn Bessel-Ungleichung erfüllt ist
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        v_arr = np.array(v, dtype=float)
        norm_sq = self.norm(v_arr) ** 2
        # Bessel-Summe: Σ|⟨v,eᵢ⟩|²
        bessel_sum = sum(self.inner_product(v_arr, np.array(e, dtype=float)) ** 2
                         for e in orthonormal_set)
        # Ungleichung: Σ|⟨v,eᵢ⟩|² ≤ ‖v‖²  (mit Toleranz für Rundungsfehler)
        return bessel_sum <= norm_sq + 1e-8


# ============================================================
# Lineare Operatoren
# ============================================================

class LinearOperator:
    """
    @brief Linearer Operator T: X → Y zwischen normierten Räumen.
    @description
        Ein Operator T: X → Y heißt linear, wenn:
            T(αu + βv) = α·T(u) + β·T(v)  für alle u,v ∈ X, α,β ∈ ℝ

        Er heißt beschränkt (stetig), wenn:
            ‖T‖ = sup{‖Tv‖_Y / ‖v‖_X : v ≠ 0} < ∞

        In endlicher Dimension sind lineare Operatoren immer beschränkt
        und können als Matrizen dargestellt werden.

        Wichtige Klassen:
        - Selbstadjungierte Operatoren: T* = T (Verallg. symmetrischer Matrizen)
        - Unitäre Operatoren: T*T = I (normerhaltend)
        - Projektionen: T² = T
        - Kompakte Operatoren: Einheitskugel → präkompakte Menge

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, matrix: np.ndarray, domain_norm: str = 'l2',
                 range_norm: str = 'l2'):
        """
        @brief Initialisiert den linearen Operator durch eine Matrix.
        @param matrix      Darstellungsmatrix des Operators (np.ndarray)
        @param domain_norm Norm auf dem Definitionsbereich ('l1','l2','linf')
        @param range_norm  Norm auf dem Wertebereich ('l1','l2','linf')
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self.matrix = np.array(matrix, dtype=float)
        self.domain_norm = domain_norm
        self.range_norm = range_norm

    def _compute_norm(self, v: np.ndarray, norm_type: str) -> float:
        """
        @brief Hilfsmethode: berechnet Norm eines Vektors.
        @param v         Vektor als np.ndarray
        @param norm_type Normtyp: 'l1', 'l2', 'linf'
        @return Normwert als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        if norm_type == 'l1':
            return float(np.sum(np.abs(v)))
        elif norm_type == 'linf':
            return float(np.max(np.abs(v)))
        else:
            # Standard: euklidische L2-Norm
            return float(np.linalg.norm(v))

    def apply(self, v: list) -> list:
        """
        @brief Wendet den Operator T auf einen Vektor an: y = T·v.
        @param v Eingangsvektor (list oder np.ndarray)
        @return Ausgabevektor T·v als list
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        v_arr = np.array(v, dtype=float)
        return (self.matrix @ v_arr).tolist()

    def is_bounded(self, test_vectors: list) -> bool:
        """
        @brief Prüft ob der Operator (auf gegebenen Testvektoren) beschränkt ist.
        @description
            Ein Operator T: X → Y heißt beschränkt, wenn:
                ‖T‖ = sup{‖Tv‖_Y : ‖v‖_X = 1} < ∞

            In endlicher Dimension sind alle linearen Operatoren beschränkt.
            Diese Methode dient der Demonstration des Konzepts.

        @param test_vectors Liste von Testvektoren
        @return True wenn ‖Tv‖/‖v‖ für alle Testvektoren endlich ist
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        for v in test_vectors:
            v_arr = np.array(v, dtype=float)
            v_norm = self._compute_norm(v_arr, self.domain_norm)
            if v_norm < 1e-15:
                continue
            tv = self.matrix @ v_arr
            tv_norm = self._compute_norm(tv, self.range_norm)
            # Verhältnis ‖Tv‖/‖v‖ muss endlich sein
            ratio = tv_norm / v_norm
            if not math.isfinite(ratio):
                return False
        return True

    def operator_norm(self) -> float:
        """
        @brief Berechnet die Operatornorm ‖T‖ = sup{‖Tv‖/‖v‖ : v≠0}.
        @description
            Für endlich-dimensionale Operatoren (Matrizen) gilt:
            - ‖T‖_2 (euklidisch) = größter Singulärwert von T
            - ‖T‖_1 = max Spaltenbetragssumme
            - ‖T‖_∞ = max Zeilenbetragssumme

            Die Operatornorm ist fundamental für die Stabilitätsanalyse.

        @return Operatornorm als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        if self.domain_norm == 'l2' and self.range_norm == 'l2':
            # ‖T‖₂ = größter Singulärwert
            singular_values = np.linalg.svd(self.matrix, compute_uv=False)
            return float(singular_values[0])
        elif self.domain_norm == 'l1' and self.range_norm == 'l1':
            # ‖T‖₁ = max Spaltenbetragssumme
            return float(np.max(np.sum(np.abs(self.matrix), axis=0)))
        elif self.domain_norm == 'linf' and self.range_norm == 'linf':
            # ‖T‖_∞ = max Zeilenbetragssumme
            return float(np.max(np.sum(np.abs(self.matrix), axis=1)))
        else:
            # Fallback: spektrale Norm
            singular_values = np.linalg.svd(self.matrix, compute_uv=False)
            return float(singular_values[0])

    def is_compact(self) -> bool:
        """
        @brief Prüft ob der Operator kompakt ist.
        @description
            Ein linearer Operator T: X → Y heißt kompakt, wenn das Bild
            jeder beschränkten Folge eine konvergente Teilfolge enthält.

            In endlicher Dimension sind ALLE beschränkten linearen Operatoren
            kompakt (Satz von Bolzano-Weierstraß).

            Der Rang des Operators bestimmt die Kompaktheit: endlicher Rang → kompakt.

        @return True (alle endlich-dim. Operatoren sind kompakt)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Endlich-dimensionale Operatoren sind immer kompakt
        return True

    def adjoint(self) -> 'LinearOperator':
        """
        @brief Berechnet den adjungierten Operator T*.
        @description
            Der adjungierte Operator T*: Y → X ist definiert durch:
                ⟨Tv, w⟩_Y = ⟨v, T*w⟩_X  für alle v ∈ X, w ∈ Y

            Für Matrizen: T* = Tᵀ (reell) bzw. T̄ᵀ (komplex, konjugiert transponiert).

            Der adjungierte Operator verallgemeinert die transponierte Matrix
            auf unendlich-dimensionale Hilbert-Räume.

        @return Adjungierter Operator T* als LinearOperator
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # T* = Tᵀ (für reelle Matrizen)
        return LinearOperator(self.matrix.T.copy(),
                               self.range_norm, self.domain_norm)

    def is_self_adjoint(self) -> bool:
        """
        @brief Prüft ob T = T* (selbstadjungiert / symmetrisch).
        @description
            T heißt selbstadjungiert, wenn T* = T, d.h.:
                ⟨Tu, v⟩ = ⟨u, Tv⟩  für alle u, v

            Wichtige Eigenschaften selbstadjungierter Operatoren:
            - Alle Eigenwerte sind reell
            - Eigenvektoren zu verschiedenen Eigenwerten sind orthogonal
            - Spektralsatz: vollständige Orthonormalbasis aus Eigenvektoren

        @return True wenn ‖T - T*‖ < Toleranz
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return np.allclose(self.matrix, self.matrix.T, atol=1e-10)

    def is_unitary(self) -> bool:
        """
        @brief Prüft ob T unitär ist (T*T = TT* = I).
        @description
            Ein Operator T heißt unitär, wenn:
                T*T = I  und  TT* = I

            Unitäre Operatoren erhalten das Skalarprodukt:
                ⟨Tu, Tv⟩ = ⟨u, v⟩  für alle u, v

            Beispiele: Rotationen, Spiegelungen, Fourier-Transformation.

        @return True wenn T*T ≈ I
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        n = self.matrix.shape[0]
        if self.matrix.shape[0] != self.matrix.shape[1]:
            return False
        # Prüfe T*T = Tᵀ·T = I
        product = self.matrix.T @ self.matrix
        return np.allclose(product, np.eye(n), atol=1e-10)

    def is_normal(self) -> bool:
        """
        @brief Prüft ob T normal ist (T*T = TT*).
        @description
            T heißt normal, wenn T*T = TT*.
            Alle selbstadjungierten und unitären Operatoren sind normal.

            Spektralsatz für normale Operatoren: jeder normale kompakte
            Operator hat eine vollständige ONB aus Eigenvektoren.

        @return True wenn ‖T*T - TT*‖ < Toleranz
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        TT_star = self.matrix @ self.matrix.T
        T_star_T = self.matrix.T @ self.matrix
        return np.allclose(TT_star, T_star_T, atol=1e-10)

    def is_projection(self) -> bool:
        """
        @brief Prüft ob T eine Projektion ist (T² = T).
        @description
            T heißt Projektion (Idempotent), wenn:
                T² = T·T = T

            Orthogonalprojektion: T² = T und T = T* (selbstadjungiert).
            Beispiel: Projektion auf einen Unterraum U.

        @return True wenn ‖T² - T‖ < Toleranz
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        T_sq = self.matrix @ self.matrix
        return np.allclose(T_sq, self.matrix, atol=1e-10)

    def kernel(self) -> list:
        """
        @brief Berechnet den Kern (Nullraum) von T.
        @description
            ker(T) = {v ∈ X : Tv = 0}

            Der Kern ist ein abgeschlossener Unterraum von X.
            dim ker(T) = Defekt von T.

            Berechnung via SVD: Vektoren zu Singulärwerten ≈ 0 spannen den Kern auf.

        @return Liste von Basisvektoren des Kerns (als list)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # SVD-basierte Kernberechnung
        _, s, Vt = np.linalg.svd(self.matrix, full_matrices=True)
        # Toleranz für numerische Null
        tol = max(self.matrix.shape) * np.finfo(float).eps * (s[0] if len(s) > 0 else 1.0)
        # Zeilen von Vt zu Singulärwerten < tol bilden den Kern
        kernel_vectors = []
        n_cols = self.matrix.shape[1]
        for i in range(n_cols):
            if i >= len(s) or s[i] < tol:
                kernel_vectors.append(Vt[i].tolist())
        return kernel_vectors

    def range_space(self) -> list:
        """
        @brief Berechnet eine Basis des Bildraums (Range) von T.
        @description
            range(T) = {Tv : v ∈ X} = Spaltenraum der Matrix

            dim range(T) = Rang(T) (Dimensionssatz: dim X = dim ker(T) + dim range(T))

            Berechnung via SVD: Spalten von U zu Singulärwerten > Toleranz.

        @return Liste von Basisvektoren des Bildraums (als list)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # SVD-basierte Bildraum-Berechnung
        U, s, _ = np.linalg.svd(self.matrix, full_matrices=True)
        tol = max(self.matrix.shape) * np.finfo(float).eps * (s[0] if len(s) > 0 else 1.0)
        # Spalten von U zu Singulärwerten > tol bilden den Bildraum
        range_vectors = []
        for i, si in enumerate(s):
            if si > tol:
                range_vectors.append(U[:, i].tolist())
        return range_vectors


class CompactOperator(LinearOperator):
    """
    @brief Kompakter linearer Operator: Einheitskugel → präkompakte Menge.
    @description
        Ein beschränkter linearer Operator K: X → Y heißt kompakt, wenn
        das Bild K(B) der Einheitskugel B = {v : ‖v‖ ≤ 1} präkompakt
        (d.h. jede Folge hat konvergente Teilfolge) in Y ist.

        Wichtige Eigenschaften:
        - In endlicher Dimension: alle beschränkten Operatoren kompakt
        - Spektrum: nur 0 als möglicher Häufungspunkt von Eigenwerten
        - Eigenwerte: abzählbar viele, mit Null als einzigem Häufungspunkt
        - Hilbert-Schmidt-Klasse: K mit Σᵢⱼ|kᵢⱼ|² < ∞

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def spectral_decomposition(self) -> dict:
        """
        @brief Berechnet die Spektralzerlegung des kompakten Operators.
        @description
            Für kompakte selbstadjungierte Operatoren T gilt der Spektralsatz:
                T = Σₙ λₙ ⟨·, eₙ⟩ eₙ

            wobei (λₙ) die (nach |λₙ| geordneten) Eigenwerte und
            {eₙ} eine zugehörige ONB aus Eigenvektoren sind.

        @return dict mit 'eigenvalues', 'eigenvectors', 'decomposition_valid'
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Symmetrisieren für selbstadjungierten Fall
        A = 0.5 * (self.matrix + self.matrix.T)
        # Eigenwerte und Eigenvektoren berechnen
        eigenvalues, eigenvectors = np.linalg.eigh(A)
        # Nach absteigendem Absolutwert sortieren
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        return {
            'eigenvalues': eigenvalues.tolist(),
            'eigenvectors': [eigenvectors[:, i].tolist()
                             for i in range(eigenvectors.shape[1])],
            'decomposition_valid': True
        }


# ============================================================
# Spektraltheorie
# ============================================================

def spectrum(operator: LinearOperator) -> dict:
    """
    @brief Berechnet das Spektrum σ(T) eines linearen Operators.
    @description
        Das Spektrum eines linearen Operators T: X → X ist:
            σ(T) = {λ ∈ ℂ : (T - λI)⁻¹ existiert nicht oder unbeschränkt}

        Zerlegung des Spektrums:
        - Punktspektrum σ_p(T):  Eigenwerte (T - λI) nicht injektiv
          d.h. ∃v≠0 : Tv = λv
        - Stetiges Spektrum σ_c(T): (T-λI)⁻¹ existiert, unbeschränkt, dicht def.
        - Residualspektrum σ_r(T): (T-λI)⁻¹ existiert, nicht dicht def.

        In endlicher Dimension: σ(T) = σ_p(T) (nur Eigenwerte).

    @param operator LinearOperator
    @return dict mit 'point_spectrum', 'eigenvalues', 'spectral_radius'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    A = operator.matrix
    # Eigenwerte berechnen (komplex, da Matrix nicht notwendig symmetrisch)
    eigenvalues = np.linalg.eigvals(A)
    # In endlicher Dimension = Punktspektrum = Eigenwerte
    point_spectrum = eigenvalues.tolist()
    return {
        'point_spectrum': point_spectrum,
        'eigenvalues': eigenvalues.tolist(),
        'spectral_radius': float(np.max(np.abs(eigenvalues))),
        'continuous_spectrum': [],   # In endl. Dim. leer
        'residual_spectrum': []      # In endl. Dim. leer
    }


def spectral_radius(operator: LinearOperator) -> float:
    """
    @brief Berechnet den Spektralradius r(T) = sup{|λ| : λ ∈ σ(T)}.
    @description
        Der Spektralradius eines Operators T ist:
            r(T) = sup{|λ| : λ ∈ σ(T)}

        Gelfand-Formel: r(T) = lim_{n→∞} ‖Tⁿ‖^{1/n}

        Für normale Operatoren (T*T = TT*) gilt: r(T) = ‖T‖.
        Der Spektralradius ist kleiner oder gleich der Operatornorm.

    @param operator LinearOperator
    @return Spektralradius als float
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    eigenvalues = np.linalg.eigvals(operator.matrix)
    return float(np.max(np.abs(eigenvalues)))


def spectral_theorem_demo(self_adjoint_matrix: np.ndarray) -> dict:
    """
    @brief Demonstration des Spektralsatzes für selbstadjungierte Operatoren.
    @description
        Spektralsatz (für kompakte selbstadjungierte Operatoren):
        Sei T: H → H ein kompakter selbstadjungierter Operator auf einem
        Hilbert-Raum H. Dann gibt es eine ONB {eₙ} aus Eigenvektoren von T und
        reelle Eigenwerte λₙ mit λₙ → 0, so dass:
            T = Σₙ λₙ ⟨·, eₙ⟩ eₙ  (Spektralzerlegung)

        Wichtige Folgerungen:
        - Alle Eigenwerte sind reell
        - Eigenvektoren zu verschiedenen Eigenwerten sind orthogonal
        - H = ⊕ₙ span{eₙ} (Hilbert-Summe)

    @param self_adjoint_matrix Symmetrische Matrix (np.ndarray)
    @return dict mit Eigenwerten, ONB, Rekonstruktionsfehler
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    A = np.array(self_adjoint_matrix, dtype=float)
    # Symmetrisieren
    A = 0.5 * (A + A.T)
    # Spektralzerlegung: A = Q·Λ·Qᵀ (Eigenzerlegung)
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # Rekonstruktion: A_rec = Σ λₙ eₙ eₙᵀ
    n = len(eigenvalues)
    A_reconstructed = sum(
        eigenvalues[i] * np.outer(eigenvectors[:, i], eigenvectors[:, i])
        for i in range(n)
    )
    reconstruction_error = float(np.linalg.norm(A - A_reconstructed))
    return {
        'eigenvalues': eigenvalues.tolist(),
        'eigenvectors': [eigenvectors[:, i].tolist() for i in range(n)],
        'reconstruction_error': reconstruction_error,
        'is_valid': reconstruction_error < 1e-10
    }


def fredholm_alternative(A: np.ndarray, b: np.ndarray) -> dict:
    """
    @brief Fredholm-Alternative: Lösungsbedingung für Ax = b.
    @description
        Fredholm-Alternative (Fredholm, 1903):
        Sei A: X → Y ein Fredholm-Operator. Dann gilt:

        Genau eine der folgenden Aussagen ist wahr:
        1. Ax = b hat eine Lösung (für alle b ⊥ ker(A*))
        2. A*y = 0 hat eine nicht-triviale Lösung

        Konkret: Ax = b ist lösbar  ⟺  b ⊥ ker(A*)
                 d.h. ⟨b, w⟩ = 0 für alle w mit A*w = 0

        Fredholm-Index: ind(A) = dim ker(A) - dim ker(A*)

    @param A Matrix (np.ndarray, m×n)
    @param b Rechte Seite (np.ndarray, m)
    @return dict mit Lösbarkeit, Kern, Index
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    # Kern von A (Nullraum)
    _, s_A, Vt_A = np.linalg.svd(A, full_matrices=True)
    tol = max(A.shape) * np.finfo(float).eps * (s_A[0] if len(s_A) > 0 else 1.0)
    n_cols = A.shape[1]
    ker_A = [Vt_A[i] for i in range(n_cols) if i >= len(s_A) or s_A[i] < tol]
    # Kern von A* = Aᵀ
    _, s_At, Vt_At = np.linalg.svd(A.T, full_matrices=True)
    tol_t = max(A.T.shape) * np.finfo(float).eps * (s_At[0] if len(s_At) > 0 else 1.0)
    n_cols_t = A.T.shape[1]
    ker_At = [Vt_At[i] for i in range(n_cols_t) if i >= len(s_At) or s_At[i] < tol_t]
    # Lösbarkeit: b ⊥ ker(A*)
    solvable = True
    for w in ker_At:
        if abs(np.dot(b, w)) > 1e-10:
            solvable = False
            break
    # Fredholm-Index: ind(A) = dim ker(A) - dim ker(A*)
    fredholm_index = len(ker_A) - len(ker_At)
    # Lösung berechnen (falls lösbar)
    solution = None
    if solvable:
        x, residuals, rank, _ = np.linalg.lstsq(A, b, rcond=None)
        solution = x.tolist()
    return {
        'solvable': solvable,
        'solution': solution,
        'ker_A': [v.tolist() for v in ker_A],
        'ker_A_star': [v.tolist() for v in ker_At],
        'fredholm_index': fredholm_index
    }


# ============================================================
# Funktionenräume
# ============================================================

class ContinuousFunctions:
    """
    @brief Raum C[a,b] der stetigen Funktionen auf [a,b].
    @description
        C[a,b] ist der Vektorraum aller stetigen Funktionen f: [a,b] → ℝ,
        ausgestattet mit der Supremumsnorm (Maximumnorm):
            ‖f‖_∞ = max{|f(x)| : x ∈ [a,b]}

        Eigenschaften:
        - Banach-Raum (vollständig bzgl. ‖·‖_∞)
        - KEIN Hilbert-Raum (Parallelogramm-Gesetz gilt nicht)
        - C[a,b] liegt dicht in L²[a,b]
        - Stone-Weierstraß: Polynome liegen dicht in C[a,b]

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, a: float = 0.0, b: float = 1.0):
        """
        @brief Initialisiert den Raum C[a,b].
        @param a Linke Intervallgrenze
        @param b Rechte Intervallgrenze
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self.a = a
        self.b = b

    def sup_norm(self, f_values: list) -> float:
        """
        @brief Berechnet die Supremumsnorm ‖f‖_∞ = max|f(x)|.
        @description
            Die Supremumsnorm (L∞-Norm) ist die natürliche Norm für C[a,b]:
                ‖f‖_∞ = sup{|f(x)| : x ∈ [a,b]}

            Da [a,b] kompakt und f stetig, wird das Supremum angenommen
            (Satz von Weierstraß): ‖f‖_∞ = max{|f(x)| : x ∈ [a,b]}.

        @param f_values Funktionswerte f(xᵢ) an diskreten Punkten
        @return Supremumsnorm als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return float(np.max(np.abs(f_values)))

    def is_dense_in_L2(self) -> bool:
        """
        @brief Bestätigt: C[a,b] liegt dicht in L²[a,b].
        @description
            Da jede stetige Funktion auf [a,b] quadratintegrierbar ist und
            da stetige Funktionen bzgl. der L²-Norm dicht in L²[a,b] liegen
            (d.h. jede L²-Funktion kann beliebig genau durch stetige Funktionen
            approximiert werden), gilt: C[a,b] ⊂ L²[a,b] dicht.

            Dies folgt aus dem Lusin'schen Satz und der Regularität.

        @return True (C[a,b] liegt stets dicht in L²[a,b])
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Theoretische Aussage: C[a,b] ist dicht in L²[a,b] bzgl. ‖·‖₂
        return True

    def stone_weierstrass_demo(self) -> dict:
        """
        @brief Demonstration des Satzes von Stone-Weierstraß.
        @description
            Stone-Weierstraß-Satz (Weierstraß 1885, Stone 1948):
            Die Polynomfunktionen P[a,b] liegen dicht in C[a,b] bzgl. ‖·‖_∞.
            D.h.: Zu jedem f ∈ C[a,b] und ε > 0 existiert ein Polynom p mit:
                ‖f - p‖_∞ < ε

            Konstruktiv: Bernstein-Polynome Bₙf konvergieren gleichmäßig gegen f.
            B_n f(x) = Σₖ₌₀ⁿ f(k/n) · C(n,k) · xᵏ · (1-x)^{n-k}

        @return dict mit Demonstrationsdaten (Polynom-Approximation von sin)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Bernstein-Approximation von f(x) = sin(2πx) auf [0,1]
        x_points = np.linspace(self.a, self.b, 100)
        f_values = np.sin(2 * math.pi * x_points)
        # Bernstein-Polynom 10. Grades
        n = 10
        # Stützstellen auf [a,b] skalieren
        f_support = np.sin(2 * math.pi * np.linspace(0, 1, n + 1))
        b_approx = np.zeros_like(x_points)
        for k in range(n + 1):
            # Bernstein-Basispolynom: B_{n,k}(t) = C(n,k) t^k (1-t)^{n-k}
            t = (x_points - self.a) / (self.b - self.a)  # auf [0,1] normieren
            binom_coeff = math.comb(n, k)
            b_approx += f_support[k] * binom_coeff * (t ** k) * ((1 - t) ** (n - k))
        error = float(np.max(np.abs(f_values - b_approx)))
        return {
            'function': 'sin(2πx)',
            'degree': n,
            'approximation_error': error,
            'dense': True,
            'statement': 'Polynome liegen dicht in C[a,b]'
        }


class SobolevSpace:
    """
    @brief Sobolev-Raum W^{k,p}(Ω).
    @description
        Der Sobolev-Raum W^{k,p}(Ω) besteht aus Funktionen f ∈ Lᵖ(Ω),
        deren schwache Ableitungen bis zur Ordnung k ebenfalls in Lᵖ(Ω) liegen:
            W^{k,p}(Ω) = {f ∈ Lᵖ(Ω) : D^α f ∈ Lᵖ(Ω) für |α| ≤ k}

        Sobolev-Norm:
            ‖f‖_{W^{k,p}} = (Σ_{|α|≤k} ‖D^α f‖_p^p)^{1/p}

        Wichtige Spezialfälle:
        - H^k = W^{k,2}: Hilbert-Raum, zentral für schwache Formulierungen von PDE
        - H¹₀(Ω): Funktionen in H¹ mit Nullrandwerten (für Dirichlet-Problem)

        Sobolev-Einbettungssatz (Sobolev 1938):
        W^{k,p}(Ω) ↪ C^m(Ω) wenn k - n/p > m  (n = Dimension von Ω)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, k: int, p: float, domain: tuple):
        """
        @brief Initialisiert den Sobolev-Raum W^{k,p}(Ω).
        @param k      Differenzierbarkeitsordnung (k ≥ 0)
        @param p      Integrierbarkeitsexponent (1 ≤ p ≤ ∞)
        @param domain Intervall (a, b) als Tupel
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self.k = k
        self.p = p
        self.a, self.b = domain

    def sobolev_norm(self, f_values: list, x_values: list, k: int) -> float:
        """
        @brief Berechnet die Sobolev-Norm ‖f‖_{W^{k,p}}.
        @description
            ‖f‖_{W^{k,p}} = (Σⱼ₌₀ᵏ ‖D^j f‖_p^p)^{1/p}

            Numerisch: Ableitungen werden durch finite Differenzen approximiert.
            Die Lᵖ-Norm wird via diskrete Quadratur berechnet.

        @param f_values Funktionswerte f(xᵢ)
        @param x_values Stützstellen xᵢ
        @param k        Maximale Ableitungsordnung
        @return Sobolev-Norm als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        f = np.array(f_values, dtype=float)
        x = np.array(x_values, dtype=float)
        dx = x[1] - x[0] if len(x) > 1 else 1.0
        total = 0.0
        current = f.copy()
        for j in range(k + 1):
            if self.p == float('inf'):
                # ‖D^j f‖_∞ = max|D^j f|
                lp_norm = float(np.max(np.abs(current)))
            else:
                # ‖D^j f‖_p^p = (Σ|D^j f(xᵢ)|^p · dx)
                lp_norm_p = float(np.sum(np.abs(current) ** self.p) * dx)
                lp_norm = lp_norm_p ** (1.0 / self.p)
            total += lp_norm ** self.p if self.p != float('inf') else lp_norm
            # Ableitung durch finite Differenz approximieren
            if j < k:
                current = np.gradient(current, dx)
        if self.p == float('inf'):
            return total
        return float(total ** (1.0 / self.p))

    def embedding_theorem(self) -> dict:
        """
        @brief Information zum Sobolev-Einbettungssatz.
        @description
            Sobolev-Einbettungssatz: Sei Ω ⊂ ℝⁿ offen und beschränkt.
                W^{k,p}(Ω) ↪ C^m(Ω̄)  wenn k - n/p > m  (m ∈ ℕ₀)
                W^{k,p}(Ω) ↪ W^{j,q}(Ω) für geeignete j, q

            Für n=1 (Intervall): W^{1,p}([a,b]) ↪ C([a,b]) wenn p > 1.
            Für H^1(Ω) = W^{1,2}(Ω): Einbettung in L^{2*} mit 2* = 2n/(n-2).

        @return dict mit Einbettungsbedingungen und Räumen
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # n=1 (eindimensionales Intervall)
        n = 1
        # W^{k,p} ↪ C^m wenn k - n/p > m
        embeddings = []
        for m in range(self.k):
            if self.k - n / self.p > m:
                embeddings.append(f'W^{{{self.k},{self.p}}} ↪ C^{m}')
        return {
            'space': f'W^{{{self.k},{self.p}}}({self.a},{self.b})',
            'dimension': n,
            'embeddings': embeddings,
            'continuous_embedding': self.k - n / self.p > 0,
            'is_hilbert': self.p == 2
        }


class L2Space:
    """
    @brief L²[a,b]-Raum: quadratintegrierbare Funktionen.
    @description
        L²[a,b] = {f messbar : ∫ₐᵇ |f(x)|² dx < ∞}

        Skalarprodukt: ⟨f, g⟩ = ∫ₐᵇ f(x)·g(x) dx
        Norm: ‖f‖₂ = √⟨f,f⟩ = √(∫ₐᵇ |f|² dx)

        L²[a,b] ist ein separabler Hilbert-Raum.
        Orthonormalbasis (Fourier): {1/√(b-a), √(2/(b-a))·cos(nπx/(b-a)), ...}

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, a: float = 0.0, b: float = 1.0):
        """
        @brief Initialisiert L²[a,b].
        @param a Linke Intervallgrenze
        @param b Rechte Intervallgrenze
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self.a = a
        self.b = b

    def inner_product(self, f: list, g: list, x: list) -> float:
        """
        @brief Berechnet das L²-Skalarprodukt ⟨f,g⟩ = ∫ₐᵇ f·g dx.
        @description
            Numerische Integration via Trapezregel:
                ⟨f,g⟩ ≈ Σᵢ f(xᵢ)·g(xᵢ)·Δxᵢ

        @param f Funktionswerte f(xᵢ)
        @param g Funktionswerte g(xᵢ)
        @param x Stützstellen xᵢ
        @return L²-Skalarprodukt als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        f_arr = np.array(f, dtype=float)
        g_arr = np.array(g, dtype=float)
        x_arr = np.array(x, dtype=float)
        # Trapezregel für ∫ f·g dx (np.trapezoid ab NumPy 2.0, np.trapz für ältere)
        trap = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        return float(trap(f_arr * g_arr, x_arr))

    def norm(self, f: list, x: list) -> float:
        """
        @brief Berechnet die L²-Norm ‖f‖₂ = √(∫ₐᵇ |f|² dx).
        @param f Funktionswerte f(xᵢ)
        @param x Stützstellen xᵢ
        @return L²-Norm als float
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return math.sqrt(max(0.0, self.inner_product(f, f, x)))

    def orthonormal_basis_fourier(self, n_terms: int, a: float, b: float) -> list:
        """
        @brief Erzeugt die Fourier-Orthonormalbasis von L²[a,b].
        @description
            Orthonormalbasis von L²[a,b]:
                e₀(x) = 1/√(b-a)
                eₙ(x) = √(2/(b-a)) · cos(nπ(x-a)/(b-a))  für n ≥ 1

            Diese bilden eine vollständige ONB von L²[a,b]:
                ⟨eₘ, eₙ⟩ = δₘₙ

        @param n_terms Anzahl der Basiselemente (inkl. Konstante)
        @param a       Linkes Intervallende
        @param b       Rechtes Intervallende
        @return Liste von ONB-Elementen als Callables
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        L = b - a
        basis_funcs = []
        # e₀ = konstant 1/√L
        basis_funcs.append(lambda x, L=L: np.ones_like(np.array(x)) / math.sqrt(L))
        # eₙ = √(2/L) · cos(nπ(x-a)/L) für n=1,...,n_terms-1
        for n in range(1, n_terms):
            def make_cosine(n, L, a):
                return lambda x: math.sqrt(2.0 / L) * np.cos(n * math.pi * (np.array(x) - a) / L)
            basis_funcs.append(make_cosine(n, L, a))
        return basis_funcs

    def projection_onto_polynomials(self, f: list, x: list, degree: int) -> list:
        """
        @brief Projiziert f auf den Unterraum der Polynome bis zum Grad `degree`.
        @description
            Die Orthogonalprojektion von f ∈ L²[a,b] auf P_d[a,b] ist das
            Polynom p* mit minimalem L²-Abstand zu f:
                p* = argmin_{deg(p)≤d} ‖f - p‖₂

            Berechnung via Gram-Schmidt auf Legendre-Polynomen oder via
            Normalengleichungen: Gram·c = rhs, wobei:
                Gram[i,j] = ⟨xⁱ, xʲ⟩, rhs[j] = ⟨f, xʲ⟩

        @param f      Funktionswerte f(xᵢ)
        @param x      Stützstellen xᵢ
        @param degree Polynomgrad
        @return Werte des Projektionspolynoms an den Stützstellen
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        f_arr = np.array(f, dtype=float)
        x_arr = np.array(x, dtype=float)
        d = degree + 1  # Anzahl der Koeffizienten
        # NumPy 2.0+: np.trapezoid statt np.trapz
        trap = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        # Gram-Matrix G[i,j] = ⟨xⁱ, xʲ⟩ = ∫ xⁱ xʲ dx
        G = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                G[i, j] = trap(x_arr ** i * x_arr ** j, x_arr)
        # Rechte Seite: rhs[j] = ⟨f, xʲ⟩ = ∫ f · xʲ dx
        rhs = np.array([trap(f_arr * x_arr ** j, x_arr) for j in range(d)])
        # Löse Normalengleichungen G·c = rhs
        try:
            c = np.linalg.lstsq(G, rhs, rcond=None)[0]
        except np.linalg.LinAlgError:
            c = np.zeros(d)
        # Polynomwerte berechnen: p(x) = Σ cⱼ xʲ
        p_values = sum(c[j] * x_arr ** j for j in range(d))
        return p_values.tolist()


# ============================================================
# Wichtige Sätze der Funktionalanalysis
# ============================================================

def hahn_banach_theorem_demo() -> dict:
    """
    @brief Demonstration des Hahn-Banach-Satzes.
    @description
        Hahn-Banach-Satz (Hahn 1927, Banach 1929):
        Sei X ein normierter Raum, U ⊆ X ein Unterraum, f: U → ℝ ein
        beschränktes lineares Funktional. Dann existiert eine Fortsetzung
        F: X → ℝ mit:
            F|_U = f  und  ‖F‖_X* = ‖f‖_U*

        Geometrische Form (Trennungssatz):
        Zwei konvexe Mengen C, D mit C ∩ D = ∅ können durch eine Hyperebene
        getrennt werden: ∃f ∈ X*, c ∈ ℝ : f(x) ≤ c ≤ f(y) für alle x∈C, y∈D.

        Folgerungen:
        - X* separiert Punkte: x≠y → ∃f : f(x)≠f(y)
        - Reichtum des Dualraums X*
        - Schwache Topologie wohldefiniert

    @return dict mit Demonstration und Schlussfolgerungen
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Demonstration: Funktional auf Unterraum fortsetzen
    # Beispiel: ℝ³, Unterraum U = {(x,0,0)}, f(x,0,0) = 2x
    # Fortsetzung auf ℝ³: F(x,y,z) = 2x (Norm gleich 2)
    u_functional_value = 2.0  # f(1,0,0) = 2
    u_norm = 2.0              # ‖f‖_U* = 2 (Norm des Funktionals auf U)
    # Fortsetzung F(x,y,z) = 2x hat gleiche Norm auf ℝ³ (bzgl. ‖·‖₂: ‖F‖=2)
    extension_norm = 2.0
    return {
        'theorem': 'Hahn-Banach',
        'subspace': 'U = {(x,0,0) : x ∈ ℝ}',
        'functional_on_U': 'f(x,0,0) = 2x',
        'extension_to_X': 'F(x,y,z) = 2x',
        'original_norm': u_norm,
        'extension_norm': extension_norm,
        'norm_preserved': abs(u_norm - extension_norm) < 1e-10,
        'consequences': [
            'Dualraum X* ist reich (trennt Punkte)',
            'Schwache Topologie ist Hausdorff',
            'Reflexivität: X ↪ X** isometrisch'
        ]
    }


def open_mapping_theorem_demo() -> dict:
    """
    @brief Demonstration des Satzes von der offenen Abbildung (Banach-Schauder).
    @description
        Satz von der offenen Abbildung (Banach 1929):
        Sei T: X → Y eine surjektive beschränkte lineare Abbildung zwischen
        Banach-Räumen X und Y. Dann ist T offen:
            T(U) offen für alle offenen U ⊆ X

        Äquivalente Formulierung: ∃c > 0 : B_Y(0,c) ⊆ T(B_X(0,1))
        (Bild der Einheitskugel enthält eine Kugel um 0)

        Korollar (Umkehrsatz): Ist T bijektiv und beschränkt, so ist T⁻¹ ebenfalls
        beschränkt (d.h. T ist ein topologischer Isomorphismus).

    @return dict mit Theoremdaten und Beispiel
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Beispiel: T: ℝ² → ℝ², T(x,y) = (2x+y, x+y) (invertierbar)
    T = np.array([[2.0, 1.0], [1.0, 1.0]])
    det_T = float(np.linalg.det(T))
    is_surjective = abs(det_T) > 1e-10
    # Inverse berechnen (falls invertierbar)
    T_inv = np.linalg.inv(T) if is_surjective else None
    inv_bounded = True  # In endlicher Dimension immer beschränkt
    return {
        'theorem': 'Offene-Abbildung (Banach-Schauder)',
        'operator': T.tolist(),
        'is_surjective': is_surjective,
        'det': det_T,
        'inverse_exists': is_surjective,
        'inverse_bounded': inv_bounded,
        'operator_norm': float(np.linalg.svd(T, compute_uv=False)[0]),
        'inverse_norm': float(np.linalg.svd(T_inv, compute_uv=False)[0]) if T_inv is not None else None,
        'consequence': 'Bijektive beschränkte Operatoren zwischen Banach-Räumen haben beschränkte Inverse'
    }


def closed_graph_theorem_demo() -> dict:
    """
    @brief Demonstration des Satzes vom abgeschlossenen Graphen.
    @description
        Satz vom abgeschlossenen Graphen (Banach):
        Sei T: X → Y linear zwischen Banach-Räumen. Dann:
            Graph(T) = {(x, Tx) : x ∈ X} ⊆ X × Y abgeschlossen
            ⟺  T ist beschränkt (stetig)

        Graph(T) abgeschlossen bedeutet:
            xₙ → x und Txₙ → y  ⟹  Tx = y

        Praktische Bedeutung: Man muss nur die Abgeschlossenheit des Graphen
        prüfen (oft einfacher), um Beschränktheit zu folgern.

    @return dict mit Theoremdaten und Beispiel
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Beispiel: T: ℝ² → ℝ², Differentiationsoperator (als Matrix)
    A = np.array([[1.0, 2.0], [3.0, 4.0]])
    op = LinearOperator(A)
    # Prüfe Graph-Abgeschlossenheit numerisch:
    # Folge xₙ → x, Axₙ → y ?= Ax
    x = np.array([1.0, 0.0])
    y_target = A @ x
    # Prüfe: Axₙ → Ax für xₙ = x + (1/n)·e₁
    errors = []
    for n in range(1, 11):
        xn = x + np.array([1.0 / n, 0.0])
        axn = A @ xn
        errors.append(float(np.linalg.norm(axn - y_target)))
    graph_closed = errors[-1] < errors[0]  # Fehler nimmt ab
    return {
        'theorem': 'Satz vom abgeschlossenen Graphen',
        'operator': A.tolist(),
        'graph_closed': graph_closed,
        'operator_bounded': True,  # In endl. Dim. immer wahr
        'convergence_errors': errors,
        'equivalence': 'Graph(T) abgeschlossen ⟺ T beschränkt'
    }


def uniform_boundedness_principle_demo() -> dict:
    """
    @brief Demonstration des Gleichmäßigen Beschränktheitsprinzips (Banach-Steinhaus).
    @description
        Satz von Banach-Steinhaus (1927):
        Sei X ein Banach-Raum, Y ein normierter Raum.
        Sei {Tₐ : α ∈ I} eine Familie beschränkter linearer Operatoren T_α: X → Y.

        Falls punktweise beschränkt: sup_α ‖Tₐx‖ < ∞ für alle x ∈ X
        Dann gleichmäßig beschränkt: sup_α ‖Tₐ‖ < ∞

        Gegenbeispiel in unvollständigem Raum: nicht möglich!
        → Vollständigkeit (Banach-Raum) ist entscheidend.

    @return dict mit Demonstration
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Familie von Operatoren Tₙ: ℝ² → ℝ², Tₙ = (1/n)·A mit A fest
    A = np.array([[1.0, 0.0], [0.0, 1.0]])
    operators = [LinearOperator((1.0 / (n + 1)) * A) for n in range(10)]
    # Punktweise Norm: ‖Tₙ x‖ für x = [1, 0]
    x_test = np.array([1.0, 0.0])
    pointwise_norms = [float(np.linalg.norm(op.apply(x_test))) for op in operators]
    # Operatornormen
    op_norms = [op.operator_norm() for op in operators]
    pointwise_bounded = max(pointwise_norms) < float('inf')
    uniformly_bounded = max(op_norms) < float('inf')
    return {
        'theorem': 'Gleichmäßiges Beschränktheitsprinzip (Banach-Steinhaus)',
        'family_size': len(operators),
        'pointwise_bounded': pointwise_bounded,
        'uniformly_bounded': uniformly_bounded,
        'max_pointwise': max(pointwise_norms),
        'max_op_norm': max(op_norms),
        'implication': 'punktweise beschränkt ⟹ gleichmäßig beschränkt (Banach-Steinhaus)'
    }


def weak_convergence_demo() -> dict:
    """
    @brief Demonstration schwacher vs. starker Konvergenz in ℓ².
    @description
        Schwache Konvergenz (xₙ ⇀ x):
            f(xₙ) → f(x) für alle f ∈ X*

        In Hilbert-Räumen (Riesz): f(xₙ) = ⟨xₙ, y⟩ → ⟨x, y⟩ für alle y.

        Klassisches Gegenbeispiel in ℓ²:
        eₙ = (0,...,0,1,0,...) (n-ter Einheitsvektor)
        - Schwach: ⟨eₙ, y⟩ = yₙ → 0 für alle y ∈ ℓ²  (eₙ ⇀ 0)
        - Stark: ‖eₙ - 0‖ = 1 ≠ 0  (keine starke Konvergenz)

        Also: schwach konvergent ≠ stark konvergent!

    @return dict mit Gegenbeispiel und Auswertung
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Einheitsvektoren eₙ in ℝ¹⁰ (Ausschnitt von ℓ²)
    dim = 10
    unit_vectors = [np.eye(dim)[n] for n in range(dim)]
    # Testvektor y = (1, 1/2, 1/3, ..., 1/10)
    y = np.array([1.0 / (n + 1) for n in range(dim)])
    # Schwache Funktionalwerte: ⟨eₙ, y⟩ = yₙ → 0
    weak_values = [float(np.dot(e, y)) for e in unit_vectors]
    # Starke Normen: ‖eₙ‖ = 1 (konstant, keine Konvergenz zu 0)
    strong_norms = [float(np.linalg.norm(e)) for e in unit_vectors]
    weak_converges = weak_values[-1] < weak_values[0]
    strong_converges = all(abs(n - 1.0) < 0.01 for n in strong_norms)  # bleibt bei 1
    return {
        'example': 'Einheitsvektoren eₙ in ℓ²',
        'weak_values': weak_values,      # ⟨eₙ, y⟩ = yₙ → 0 (schwach)
        'strong_norms': strong_norms,    # ‖eₙ‖ = 1 (stark: keine Konvergenz zu 0)
        'weakly_converges_to_zero': weak_converges,
        'strongly_converges_to_zero': not strong_converges,
        'conclusion': 'schwach konvergent ≠ stark konvergent'
    }


def riesz_representation_theorem() -> dict:
    """
    @brief Riesz-Darstellungssatz für Hilbert-Räume.
    @description
        Riesz-Darstellungssatz (Riesz 1907, Fréchet 1907):
        Sei H ein Hilbert-Raum. Dann ist die Abbildung:
            J: H → H*, y ↦ f_y  mit  f_y(x) = ⟨x, y⟩

        ein isometrischer Isomorphismus (antilinear für ℂ).

        Folgerungen:
        - H* ≅ H (Hilbert-Räume sind reflexiv: H** ≅ H)
        - ‖f_y‖_{H*} = ‖y‖_H
        - Zu jedem f ∈ H* existiert eindeutiges y ∈ H mit f(·) = ⟨·, y⟩

    @return dict mit Demonstration
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Beispiel: H = ℝ³, Funktional f(x) = 2x₁ + 3x₂ - x₃
    # Riesz-Repräsentant: y = (2, 3, -1)
    y = np.array([2.0, 3.0, -1.0])
    def functional(x):
        return float(np.dot(x, y))
    # Teste an verschiedenen Vektoren
    test_vectors = [
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([1.0, 1.0, 1.0])
    ]
    verifications = []
    for x in test_vectors:
        f_x = functional(x)
        inner_x_y = float(np.dot(x, y))  # ⟨x, y⟩
        verifications.append({
            'x': x.tolist(),
            'f(x)': f_x,
            '<x,y>': inner_x_y,
            'equal': abs(f_x - inner_x_y) < 1e-10
        })
    return {
        'theorem': 'Riesz-Darstellungssatz',
        'riesz_representative': y.tolist(),
        'functional': 'f(x) = 2x₁ + 3x₂ - x₃',
        'norm_functional': float(np.linalg.norm(y)),  # ‖f‖ = ‖y‖
        'verifications': verifications,
        'isomorphism': 'H* ≅ H via ⟨·, y⟩'
    }


def compact_operator_theory() -> dict:
    """
    @brief Theorie kompakter Operatoren.
    @description
        Theorie kompakter Operatoren (Riesz, Schauder, Hilbert-Schmidt):

        1. Spektraleigenschaften:
           - σ(K) ohne {0} besteht nur aus Eigenwerten
           - 0 ist einziger möglicher Häufungspunkt von σ(K)
           - Jeder Eigenwert ≠ 0 hat endliche Vielfachheit

        2. Hilbert-Schmidt-Operatoren:
           K ist Hilbert-Schmidt ⟺ ‖K‖_HS = (Σᵢⱼ |kᵢⱼ|²)^{1/2} < ∞
           Alle Hilbert-Schmidt-Operatoren sind kompakt.

        3. Schur-Zerlegung: Kompaktes K in Hilbert-Raum hat ONB
           aus (verallgemeinerten) Eigenvektoren.

    @return dict mit Demonstrationsdaten
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Beispiel: kompakter Operator als Hilbert-Schmidt-Operator
    K = np.array([[1.0, 0.5, 0.25],
                  [0.5, 0.25, 0.125],
                  [0.25, 0.125, 0.0625]])
    op = CompactOperator(K)
    # Hilbert-Schmidt-Norm ‖K‖_HS = √(Σᵢⱼ |kᵢⱼ|²)
    hs_norm = float(np.sqrt(np.sum(K ** 2)))
    # Spektrum
    spec = spectrum(op)
    # Spektralzerlegung
    decomp = op.spectral_decomposition()
    # Häufungspunkte der Eigenwerte
    eigenvalues = np.array(spec['eigenvalues'])
    return {
        'operator': K.tolist(),
        'hilbert_schmidt_norm': hs_norm,
        'is_compact': op.is_compact(),
        'eigenvalues': spec['eigenvalues'],
        'spectral_radius': spec['spectral_radius'],
        'spectral_decomposition_valid': decomp['decomposition_valid'],
        'accumulation_point': 0.0,
        'properties': [
            r'σ(K) \ {0} = Punktspektrum',
            '0 ist einziger Häufungspunkt',
            'Endliche Vielfachheit jedes Eigenwerts ≠ 0',
            f'Hilbert-Schmidt-Norm: {hs_norm:.4f}'
        ]
    }


# ============================================================
# Fixpunktsätze
# ============================================================

def banach_fixed_point(f: Callable, x0: float, tol: float = 1e-10,
                        max_iter: int = 1000) -> dict:
    """
    @brief Banach-Fixpunktsatz: Kontraktion hat eindeutigen Fixpunkt.
    @description
        Banachscher Fixpunktsatz (1922):
        Sei (X, d) ein vollständiger metrischer Raum und T: X → X eine
        Kontraktion mit Kontraktionskonstante q ∈ [0,1):
            d(T(x), T(y)) ≤ q·d(x,y) für alle x,y ∈ X

        Dann:
        1. T hat genau einen Fixpunkt x* ∈ X mit T(x*) = x*
        2. Die Iteration xₙ₊₁ = T(xₙ) konvergiert für jeden Startpunkt x₀
        3. A-priori-Schranke: d(xₙ, x*) ≤ qⁿ/(1-q) · d(x₁,x₀)
        4. A-posteriori-Schranke: d(xₙ, x*) ≤ q/(1-q) · d(xₙ,xₙ₋₁)

    @param f        Kontraktion T: ℝ → ℝ
    @param x0       Startwert
    @param tol      Konvergenztoleranz
    @param max_iter Maximale Iterationsanzahl
    @return dict mit Fixpunkt, Konvergenzhistorie, Kontraktionskonstante
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    history = [x0]
    x = x0
    converged = False
    for i in range(max_iter):
        x_new = f(x)
        history.append(x_new)
        # Konvergenzprüfung: |xₙ₊₁ - xₙ| < tol
        if abs(x_new - x) < tol:
            converged = True
            x = x_new
            break
        x = x_new
    # Kontraktionskonstante schätzen: q ≈ |xₙ-xₙ₋₁|/|xₙ₋₁-xₙ₋₂|
    q_estimate = None
    if len(history) >= 3:
        d1 = abs(history[-1] - history[-2])
        d2 = abs(history[-2] - history[-3])
        if d2 > 1e-15:
            q_estimate = d1 / d2
    return {
        'fixed_point': x,
        'converged': converged,
        'iterations': len(history) - 1,
        'history': history[:10],       # Nur erste 10 Iterationen
        'final_residual': abs(f(x) - x),
        'contraction_constant_estimate': q_estimate
    }


def schauder_fixed_point_demo() -> dict:
    """
    @brief Demonstration des Schauder-Fixpunktsatzes.
    @description
        Schauder-Fixpunktsatz (1930):
        Sei K ⊆ X eine nichtleere, konvexe, kompakte Teilmenge eines
        normierten Raums X. Sei T: K → K stetig. Dann hat T mindestens
        einen Fixpunkt.

        Verallgemeinert den Brouwer-Fixpunktsatz auf unendlich-dimensionale
        Banach-Räume (wobei K kompakt sein muss, was in unendl. Dim. selten ist).

        Anwendungen:
        - Existenzbeweise für ODE-Lösungen (Peano-Satz)
        - Existenz von Gleichgewichten in der Spieltheorie
        - PDE-Existenztheorie

        Beispiel: T: [0,1] → [0,1], T(x) = x² (stetig, K=[0,1] kompakt/konvex)
        Fixpunkte: x = 0 und x = 1.

    @return dict mit Demonstration und Fixpunktnachweis
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Demonstration: T(x) = x² auf K = [0, 1]
    # Fixpunkte: T(x) = x ⟺ x² = x ⟺ x(x-1) = 0 ⟺ x = 0 oder x = 1
    def T(x):
        return x ** 2
    # K = [0,1] ist konvex und kompakt
    k_convex = True
    k_compact = True
    # Fixpunkte numerisch finden
    fixed_points = []
    for x0 in np.linspace(0, 1, 20):
        result = banach_fixed_point(T, x0, tol=1e-12, max_iter=500)
        fp = round(result['fixed_point'], 8)
        if fp not in [round(p, 8) for p in fixed_points]:
            if abs(T(result['fixed_point']) - result['fixed_point']) < 1e-8:
                fixed_points.append(result['fixed_point'])
    return {
        'theorem': 'Schauder-Fixpunktsatz',
        'operator': 'T(x) = x² auf K = [0,1]',
        'K_convex': k_convex,
        'K_compact': k_compact,
        'T_continuous': True,
        'fixed_points_found': fixed_points,
        'expected_fixed_points': [0.0, 1.0],
        'applications': [
            'Peano-Existenzsatz für ODE',
            'Nash-Gleichgewichte (Spieltheorie)',
            'PDE-Existenztheorie'
        ]
    }
