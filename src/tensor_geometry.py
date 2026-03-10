"""
@file tensor_geometry.py
@brief Tensor- und Differentialgeometrie – Tensoren, Krümmung, Riemannsche Geometrie.
@description
    Implementiert grundlegende Konzepte der Tensor- und Differentialgeometrie:

    - Tensoren: Skalare, Vektoren, Matrizen, allgemeine (p,q)-Tensoren
    - Tensoroperationen: Kontraktion, äußeres Produkt, Verjüngung
    - Metriktensor: Euklidische und nicht-euklidische Metriken
    - Christoffel-Symbole: Verbindungskoeffizienten der Levi-Civita-Verbindung
    - Riemannscher Krümmungstensor
    - Ricci-Tensor und Ricci-Skalar
    - Geodäten: Kürzeste Verbindungen auf gekrümmten Flächen
    - Klassische Flächen: Sphäre, Torus, Hyperboloid, Sattelform

    Bezug zur Physik:
    - Allgemeine Relativitätstheorie: g_{μν} beschreibt Raumzeit-Krümmung
    - Einstein-Gleichungen: G_{μν} = 8πG T_{μν}
    - Geodätengleichung: d²x^μ/dτ² + Γ^μ_{αβ} dx^α/dτ dx^β/dτ = 0

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import functools
import numpy as np
from typing import Callable


# ---------------------------------------------------------------------------
# Klasse: Tensor
# ---------------------------------------------------------------------------

class Tensor:
    """
    Allgemeiner (p,q)-Tensor: p kontravariante, q kovariante Indizes.

    Ein Tensor T^{i₁...iₚ}_{j₁...j_q} ist eine multilineare Abbildung.

    Typen nach Rang:
    - (0,0)-Tensor: Skalar
    - (1,0)-Tensor: Kontravariantervektor (Pfeil)
    - (0,1)-Tensor: Kovarianter Vektor / 1-Form (Gradient)
    - (1,1)-Tensor: Lineare Abbildung / Matrix
    - (0,2)-Tensor: Bilinearform (z.B. Metriktensor)
    - (1,3)-Tensor: Riemannscher Krümmungstensor R^i_{jkl}

    Koordinatentransformation (Einstein-Summenkonvention):
        T'^{i'}_{j'} = (∂x^{i'}/∂x^i)(∂x^j/∂x^{j'}) T^i_j

    Beispiel:
        >>> import numpy as np
        >>> T = Tensor([[1, 2], [3, 4]], contravariant=1, covariant=1)
        >>> T.trace()
        5.0
    """

    def __init__(
        self,
        components: "np.ndarray | list",
        contravariant: int = 0,
        covariant: int = 0,
    ) -> None:
        """
        Initialisiert einen (p,q)-Tensor.

        @param components: Numpy-Array oder verschachtelte Liste mit Tensor-Komponenten.
        @param contravariant: Anzahl kontravarianter Indizes p (hochgestellte Indizes).
        @param covariant: Anzahl kovarianter Indizes q (tiefgestellte Indizes).
        @lastModified 2026-03-10
        """
        # Daten als float64-Array speichern
        self.data = np.array(components, dtype=float)
        self.contravariant = contravariant   # p (obere Indizes)
        self.covariant = covariant           # q (untere Indizes)
        self.rank = contravariant + covariant
        self.shape = self.data.shape

    def __add__(self, other: "Tensor") -> "Tensor":
        """
        Tensoraddition: (T + S)^{i}_{j} = T^{i}_{j} + S^{i}_{j}.

        Gleicher Rang und gleiche Form erforderlich.

        @param other: Zweiter Tensor (gleicher Typ).
        @return: Summen-Tensor.
        @raises ValueError: Bei unterschiedlichen Formen.
        @lastModified 2026-03-10
        """
        if self.shape != other.shape:
            raise ValueError(
                f"Tensoraddition erfordert gleiche Form: {self.shape} vs {other.shape}"
            )
        return Tensor(
            self.data + other.data,
            contravariant=self.contravariant,
            covariant=self.covariant,
        )

    def __sub__(self, other: "Tensor") -> "Tensor":
        """
        Tensorsubtraktion: (T - S)^{i}_{j} = T^{i}_{j} - S^{i}_{j}.

        @param other: Zweiter Tensor (gleicher Typ).
        @return: Differenz-Tensor.
        @raises ValueError: Bei unterschiedlichen Formen.
        @lastModified 2026-03-10
        """
        if self.shape != other.shape:
            raise ValueError(
                f"Tensorsubtraktion erfordert gleiche Form: {self.shape} vs {other.shape}"
            )
        return Tensor(
            self.data - other.data,
            contravariant=self.contravariant,
            covariant=self.covariant,
        )

    def __mul__(self, scalar: float) -> "Tensor":
        """
        Skalarmultiplikation: (λT)^{i}_{j} = λ · T^{i}_{j}.

        @param scalar: Reelle Zahl.
        @return: Skalierter Tensor.
        @lastModified 2026-03-10
        """
        return Tensor(
            self.data * float(scalar),
            contravariant=self.contravariant,
            covariant=self.covariant,
        )

    def __rmul__(self, scalar: float) -> "Tensor":
        """Rechtsseitige Skalarmultiplikation (λ * T)."""
        return self.__mul__(scalar)

    def outer_product(self, other: "Tensor") -> "Tensor":
        """
        Äußeres (tensorielles) Produkt T ⊗ S.

        Definition:
            (T ⊗ S)^{i₁i₂}_{j₁j₂} = T^{i₁}_{j₁} · S^{i₂}_{j₂}

        Der Rang des Ergebnis-Tensors ist (p₁+p₂, q₁+q₂).
        Geometrisch: kombiniert zwei Tensorräume zu einem größeren.

        @param other: Zweiter Tensor beliebigen Rangs.
        @return: Tensor mit erhöhtem Rang.
        @lastModified 2026-03-10

        Beispiel:
            >>> import numpy as np
            >>> u = Tensor([1.0, 2.0], contravariant=1)
            >>> v = Tensor([3.0, 4.0], contravariant=1)
            >>> T = u.outer_product(v)
            >>> T.shape
            (2, 2)
        """
        # np.outer funktioniert nur für 1D; np.tensordot mit axes=0 ist allgemeiner
        result = np.tensordot(self.data, other.data, axes=0)
        return Tensor(
            result,
            contravariant=self.contravariant + other.contravariant,
            covariant=self.covariant + other.covariant,
        )

    def contract(self, i: int, j: int) -> "Tensor":
        """
        Verjüngung (Kontraktion) über das Index-Paar (i, j).

        Setzt den i-ten oberen Index gleich dem j-ten unteren Index
        und summiert (Einstein-Summenkonvention).

        Beispiel Spur: T^i_j → T^k_k = Spur(T) (Rang (1,1) → (0,0))

        @param i: Position des kontravarianten Index (0-basiert, unter den p oberen).
        @param j: Position des kovarianten Index (0-basiert, unter den q unteren).
        @return: Kontrahierter Tensor mit Rang (p-1, q-1).
        @raises ValueError: Bei ungültigem Rang (p=0 oder q=0).
        @lastModified 2026-03-10
        """
        if self.contravariant == 0 or self.covariant == 0:
            raise ValueError(
                "Kontraktion erfordert mindestens einen kontravarianten "
                "und einen kovarianten Index."
            )
        # Der i-te obere Index liegt an Achse i
        # Der j-te untere Index liegt an Achse (contravariant + j)
        axis_upper = i
        axis_lower = self.contravariant + j

        # np.trace summiert über zwei Achsen
        result = np.trace(self.data, axis1=axis_upper, axis2=axis_lower)
        return Tensor(
            result,
            contravariant=self.contravariant - 1,
            covariant=self.covariant - 1,
        )

    def trace(self) -> float:
        """
        Spur eines (1,1)-Tensors (Matrix): Σ_i T^i_i.

        Entspricht der Summe der Diagonalelemente einer Matrix.
        Ist eine Basis-unabhängige Invariante.

        @return: Spur als Skalar.
        @raises ValueError: Falls der Tensor kein (1,1)-Tensor ist.
        @lastModified 2026-03-10
        """
        if self.contravariant != 1 or self.covariant != 1:
            raise ValueError(
                f"Spur nur für (1,1)-Tensor definiert, nicht ({self.contravariant},{self.covariant})."
            )
        return float(np.trace(self.data))

    def symmetrize(self) -> "Tensor":
        """
        Symmetrisierung eines (0,2)-Tensors: T_{(ij)} = ½(T_{ij} + T_{ji}).

        Zerlegt T in seinen symmetrischen Anteil.
        Jeder Tensor lässt sich eindeutig zerlegen: T = T_{(ij)} + T_{[ij]}.

        @return: Symmetrisierter (0,2)-Tensor.
        @raises ValueError: Falls nicht (0,2)-Tensor mit quadratischer Form.
        @lastModified 2026-03-10
        """
        if self.data.ndim != 2:
            raise ValueError("Symmetrisierung nur für 2D-Arrays (0,2)-Tensoren definiert.")
        sym = 0.5 * (self.data + self.data.T)
        return Tensor(sym, contravariant=self.contravariant, covariant=self.covariant)

    def antisymmetrize(self) -> "Tensor":
        """
        Antisymmetrisierung eines (0,2)-Tensors: T_{[ij]} = ½(T_{ij} - T_{ji}).

        Zerlegt T in seinen antisymmetrischen Anteil.
        Für eine antisymmetrische Matrix gilt: T_{ii} = 0.

        @return: Antisymmetrisierter (0,2)-Tensor.
        @raises ValueError: Falls nicht 2D-Array.
        @lastModified 2026-03-10
        """
        if self.data.ndim != 2:
            raise ValueError("Antisymmetrisierung nur für 2D-Arrays (0,2)-Tensoren definiert.")
        antisym = 0.5 * (self.data - self.data.T)
        return Tensor(antisym, contravariant=self.contravariant, covariant=self.covariant)

    def is_symmetric(self, tol: float = 1e-10) -> bool:
        """
        Prüft ob ein (0,2)-Tensor symmetrisch ist: T_{ij} = T_{ji}.

        Metriktensoren und der Ricci-Tensor sind symmetrisch.

        @param tol: Numerische Toleranz für den Vergleich.
        @return: True falls symmetrisch.
        @lastModified 2026-03-10
        """
        if self.data.ndim != 2:
            return False
        return bool(np.allclose(self.data, self.data.T, atol=tol))

    def is_antisymmetric(self, tol: float = 1e-10) -> bool:
        """
        Prüft ob ein (0,2)-Tensor antisymmetrisch ist: T_{ij} = -T_{ji}.

        Das Levi-Civita-Symbol und der Faraday-Tensor (EM-Feld) sind antisymmetrisch.

        @param tol: Numerische Toleranz.
        @return: True falls antisymmetrisch.
        @lastModified 2026-03-10
        """
        if self.data.ndim != 2:
            return False
        return bool(np.allclose(self.data, -self.data.T, atol=tol))

    def __repr__(self) -> str:
        return f"Tensor({self.contravariant},{self.covariant})\n{self.data}"


# ---------------------------------------------------------------------------
# Klasse: MetricTensor
# ---------------------------------------------------------------------------

class MetricTensor:
    """
    Metriktensor g_{ij} einer Riemannschen Mannigfaltigkeit.

    Der Metriktensor definiert die lokale Geometrie:
    - Abstände: ds² = g_{ij} dx^i dx^j
    - Winkel zwischen Vektoren: cos θ = g(u,v) / (|u| · |v|)
    - Volumenelement: dV = √|det(g)| dx¹ ∧ ... ∧ dxⁿ

    Eigenschaften:
    - Symmetrisch: g_{ij} = g_{ji}
    - Nicht-degeneriert: det(g) ≠ 0
    - Für Riemannsche Metrik: positiv definit (alle Eigenwerte > 0)
    - Für Lorentz-Metrik (ART): Signatur (-,+,+,+)

    Beispiel (euklidische Ebene):
        >>> MetricTensor([[1, 0], [0, 1]]).determinant()
        1.0
    """

    def __init__(self, components: "list[list[float]] | np.ndarray") -> None:
        """
        Initialisiert den Metriktensor.

        @param components: n×n Matrix der Metrikkomponenten g_{ij}.
        @raises ValueError: Falls nicht quadratisch oder singulär.
        @lastModified 2026-03-10
        """
        self.g = np.array(components, dtype=float)
        if self.g.ndim != 2 or self.g.shape[0] != self.g.shape[1]:
            raise ValueError("Metriktensor muss eine quadratische Matrix sein.")
        self.n = self.g.shape[0]
        # Inverse Metrik g^{ij} für Index-Hebung
        self.g_inv = np.linalg.inv(self.g)

    def inner_product(self, u: "list[float]", v: "list[float]") -> float:
        """
        Inneres Produkt zweier Vektoren: <u,v>_g = g_{ij} u^i v^j.

        Im euklidischen Fall (g = I): entspricht dem Standard-Skalarprodukt.
        Auf einer Sphäre hängt das Produkt vom Punkt ab (da g punktabhängig ist).

        @param u: Erster Vektor als Liste [u¹, u², ..., uⁿ].
        @param v: Zweiter Vektor als Liste [v¹, v², ..., vⁿ].
        @return: Skalar g_{ij} u^i v^j.
        @lastModified 2026-03-10
        """
        u_arr = np.array(u, dtype=float)
        v_arr = np.array(v, dtype=float)
        # g_{ij} u^i v^j = uᵀ · g · v
        return float(u_arr @ self.g @ v_arr)

    def norm(self, v: "list[float]") -> float:
        """
        Norm eines Vektors: |v|_g = √(g_{ij} v^i v^j).

        @param v: Vektor als Liste.
        @return: Nicht-negative reelle Zahl.
        @raises ValueError: Falls das innere Produkt negativ ist (Lorentz-Metrik).
        @lastModified 2026-03-10
        """
        ip = self.inner_product(v, v)
        if ip < 0:
            # Lorentz-Metrik: Norm kann imaginär sein
            return float(math.sqrt(abs(ip)))
        return float(math.sqrt(ip))

    def raise_index(self, tensor: Tensor, index: int) -> Tensor:
        """
        Indexhebung mit der inversen Metrik: T^i = g^{ij} T_j.

        Wandelt einen kovarianten Index in einen kontravarianten um.
        Geometrisch: wandelt eine 1-Form in einen Vektor um (musical isomorphism ♯).

        @param tensor: Tensor mit mindestens einem kovarianten Index.
        @param index: Position des kovarianten Index (0-basiert).
        @return: Tensor mit einem zusätzlichen kontravarianten Index.
        @lastModified 2026-03-10
        """
        # Kontraktion: g^{ij} T_j über Achse 'index'
        result = np.tensordot(self.g_inv, tensor.data, axes=([1], [index]))
        # Index wieder an die richtige Stelle bringen
        result = np.moveaxis(result, 0, index)
        return Tensor(
            result,
            contravariant=tensor.contravariant + 1,
            covariant=tensor.covariant - 1,
        )

    def lower_index(self, tensor: Tensor, index: int) -> Tensor:
        """
        Indexsenkung mit der Metrik: T_i = g_{ij} T^j.

        Wandelt einen kontravarianten Index in einen kovarianten um.
        Geometrisch: wandelt einen Vektor in eine 1-Form um (musical isomorphism ♭).

        @param tensor: Tensor mit mindestens einem kontravarianten Index.
        @param index: Position des kontravarianten Index (0-basiert).
        @return: Tensor mit einem zusätzlichen kovarianten Index.
        @lastModified 2026-03-10
        """
        result = np.tensordot(self.g, tensor.data, axes=([1], [index]))
        result = np.moveaxis(result, 0, index)
        return Tensor(
            result,
            contravariant=tensor.contravariant - 1,
            covariant=tensor.covariant + 1,
        )

    def determinant(self) -> float:
        """
        Determinante des Metriktensors: det(g).

        Wichtig für:
        - Volumenelement: dV = √|det(g)| d^n x
        - Nicht-Degeneriertheit: det(g) ≠ 0
        - Vorzeichen gibt Signatur an (+: Riemannsch, -: Lorentzsich bei ungeradem n)

        @return: det(g) als reelle Zahl.
        @lastModified 2026-03-10
        """
        return float(np.linalg.det(self.g))

    def volume_element(self) -> float:
        """
        Volumenelement: √|det(g)|.

        Dieser Faktor tritt bei Integration auf Mannigfaltigkeiten auf:
        ∫_M f dV = ∫ f(x) √|det(g(x))| dx¹...dxⁿ

        @return: √|det(g)| als nicht-negative reelle Zahl.
        @lastModified 2026-03-10
        """
        return float(math.sqrt(abs(self.determinant())))

    def is_flat(self, point: "list[float]", tol: float = 1e-10) -> bool:
        """
        Prüft ob die Metrik (lokal) flach ist.

        Für konstante Metriken ist der Riemann-Tensor immer null → flach.
        Diese Methode prüft nur ob die Metrik konstant ist (hinreichend für Flachheit).

        @param point: Koordinatenpunkt (für zukünftige Erweiterung mit metrischer Funktion).
        @param tol: Toleranz.
        @return: True falls Metrik konstant (kein punkt-abhängiger Term).
        @lastModified 2026-03-10
        """
        # Für ein MetricTensor-Objekt mit fester Matrix: immer flach
        # (punktabhängige Metriken werden über metric_func-Funktionen dargestellt)
        return True


# ---------------------------------------------------------------------------
# Hilfsfunktion: Numerische Ableitung des Metriktensors
# ---------------------------------------------------------------------------

def _metric_derivative(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    axis: int,
    h: float = 1e-5,
) -> np.ndarray:
    """
    Numerische partielle Ableitung der Metrik nach der axis-ten Koordinate.

    Verwendet den zentralen Differenzenquotienten:
        ∂g/∂x^k ≈ (g(x+h·eₖ) - g(x-h·eₖ)) / (2h)

    @param metric_func: Funktion point → g_{ij}(point) als numpy-Matrix.
    @param point: Koordinatenpunkt.
    @param axis: Index der Koordinate (0-basiert).
    @param h: Schrittweite.
    @return: n×n Matrix der partiellen Ableitungen ∂g_{ij}/∂x^axis.
    @lastModified 2026-03-10
    """
    p = list(point)
    # Vorwärtspunkt x + h·eₖ
    p_plus = list(p)
    p_plus[axis] += h
    # Rückwärtspunkt x - h·eₖ
    p_minus = list(p)
    p_minus[axis] -= h
    # Zentraler Differenzenquotient (O(h²)-Genauigkeit)
    return (metric_func(p_plus) - metric_func(p_minus)) / (2.0 * h)


# ---------------------------------------------------------------------------
# Christoffel-Symbole
# ---------------------------------------------------------------------------

def christoffel_symbols(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-5,
) -> np.ndarray:
    """
    Berechnet Christoffel-Symbole Γ^k_{ij} numerisch.

    Die Christoffel-Symbole (Levi-Civita-Verbindung) beschreiben,
    wie Koordinatenbasisvektoren sich von Punkt zu Punkt ändern.
    Sie sind KEIN Tensor (transformieren sich anders als Tensoren).

    Formel:
        Γ^k_{ij} = ½ g^{kl} (∂_i g_{jl} + ∂_j g_{il} - ∂_l g_{ij})

    Symmetrie: Γ^k_{ij} = Γ^k_{ji} (für torsionsfreie Verbindungen)

    Im flachen Raum mit kartesischen Koordinaten: Γ^k_{ij} = 0 überall.
    Auf einer Sphäre: Γ^θ_{φφ} = -sin(θ)cos(θ), Γ^φ_{θφ} = cos(θ)/sin(θ).

    @param metric_func: Funktion point → g_{ij}(point) als numpy-Matrix.
    @param point: Koordinatenpunkt [x¹, x², ..., xⁿ].
    @param h: Schrittweite für numerische Ableitung.
    @return: Array Γ[k,i,j] der Form (n,n,n).
    @lastModified 2026-03-10
    """
    # Metrik und ihre Inverse am gegebenen Punkt
    g = metric_func(point)
    n = g.shape[0]
    g_inv = np.linalg.inv(g)

    # Alle partiellen Ableitungen ∂_k g_{ij} vorberechnen
    # dg[k, i, j] = ∂g_{ij}/∂x^k
    dg = np.zeros((n, n, n))
    for k in range(n):
        dg[k] = _metric_derivative(metric_func, point, k, h)

    # Christoffel-Symbole berechnen: Γ^k_{ij}
    Gamma = np.zeros((n, n, n))
    for k in range(n):
        for i in range(n):
            for j in range(n):
                # Summe über l: ½ g^{kl} (∂_i g_{jl} + ∂_j g_{il} - ∂_l g_{ij})
                for l in range(n):
                    Gamma[k, i, j] += 0.5 * g_inv[k, l] * (
                        dg[i, j, l]   # ∂_i g_{jl}
                        + dg[j, i, l]  # ∂_j g_{il}
                        - dg[l, i, j]  # ∂_l g_{ij}
                    )
    return Gamma


# ---------------------------------------------------------------------------
# Riemannscher Krümmungstensor
# ---------------------------------------------------------------------------

def riemann_tensor(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-4,
) -> np.ndarray:
    """
    Riemannscher Krümmungstensor R^l_{kij} numerisch.

    Der Riemann-Tensor misst die Nicht-Kommutativität des kovarianten Ableitens.
    Anschaulich: Wenn man einen Vektor parallel entlang einer geschlossenen Kurve
    transportiert, dreht er sich um einen Winkel proportional zu R.

    Formel:
        R^l_{kij} = ∂_i Γ^l_{jk} - ∂_j Γ^l_{ik}
                    + Γ^l_{im} Γ^m_{jk} - Γ^l_{jm} Γ^m_{ik}

    Symmetrien des Riemann-Tensors:
    - R_{lkij} = -R_{lkji}  (antisymmetrisch in i,j)
    - R_{lkij} = -R_{klij}  (antisymmetrisch in k,l)
    - R_{lkij} = R_{ijlk}   (Paarvertauschung)
    - Bianchi-Identität: R^l_{k[ij;m]} = 0

    Im flachen Raum: R^l_{kij} = 0 überall.

    @param metric_func: Funktion point → g_{ij}(point).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite für numerische Ableitungen der Christoffel-Symbole.
    @return: Array R[l,k,i,j] der Form (n,n,n,n).
    @lastModified 2026-03-10
    """
    g = metric_func(point)
    n = g.shape[0]

    # Manueller Cache für Christoffel-Symbole: gleiche Punkte werden nicht neu berechnet.
    # Schlüssel: Tupel aus Punkt-Koordinaten (hashbar), Wert: numpy-Array Gamma[k,i,j].
    # Spart bei der numerischen Ableitung (2*n zusätzliche Aufrufe pro Riemann-Tensor)
    # erheblich Rechenzeit, wenn metric_func teuer ist.
    _gamma_cache: dict[tuple, np.ndarray] = {}

    def Gamma_at(pt: "list[float]") -> np.ndarray:
        """Christoffel-Symbole am Punkt pt, gecacht nach Koordinaten-Tupel."""
        key = tuple(round(c, 15) for c in pt)   # Tupel ist hashbar, auf 15 Stellen runden
        if key not in _gamma_cache:
            _gamma_cache[key] = christoffel_symbols(metric_func, list(pt), h=h * 0.1)
        return _gamma_cache[key]

    # Christoffel-Symbole und ihre Ableitungen
    Gamma = Gamma_at(list(point))

    # Numerische Ableitungen ∂_i Γ^l_{jk} via zentralem Differenzenquotienten
    # dGamma[i, l, j, k] = ∂_i Γ^l_{jk}
    dGamma = np.zeros((n, n, n, n))
    for i in range(n):
        p_plus = list(point)
        p_plus[i] += h
        p_minus = list(point)
        p_minus[i] -= h
        G_plus = Gamma_at(p_plus)
        G_minus = Gamma_at(p_minus)
        # Ableitung ∂_i von Gamma[l,j,k]
        dGamma[i] = (G_plus - G_minus) / (2.0 * h)

    # Riemann-Tensor: R^l_{kij}
    R = np.zeros((n, n, n, n))
    for l in range(n):
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    # Ableitungsterme
                    term1 = dGamma[i, l, j, k]   # ∂_i Γ^l_{jk}
                    term2 = dGamma[j, l, i, k]   # ∂_j Γ^l_{ik}
                    # Quadratische Terme: Γ^l_{im} Γ^m_{jk} - Γ^l_{jm} Γ^m_{ik}
                    quad1 = sum(Gamma[l, i, m] * Gamma[m, j, k] for m in range(n))
                    quad2 = sum(Gamma[l, j, m] * Gamma[m, i, k] for m in range(n))
                    R[l, k, i, j] = term1 - term2 + quad1 - quad2
    return R


# ---------------------------------------------------------------------------
# Ricci-Tensor
# ---------------------------------------------------------------------------

def ricci_tensor(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-4,
) -> np.ndarray:
    """
    Ricci-Tensor R_{ij} = R^k_{ikj} (Kontraktion des Riemann-Tensors).

    Der Ricci-Tensor ist eine "Komprimierung" des Riemann-Tensors und
    misst wie das Volumen einer kleinen geodätischen Kugel von dem
    flachen Fall abweicht.

    Ist symmetrisch: R_{ij} = R_{ji} (folgt aus Bianchi-Identität).

    Zentral in Einsteins Feldgleichungen:
        G_{ij} = R_{ij} - ½ g_{ij} R = 8πG/c⁴ · T_{ij}

    @param metric_func: Funktion point → g_{ij}(point).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite.
    @return: Array R[i,j] der Form (n,n), symmetrische Matrix.
    @lastModified 2026-03-10
    """
    g = metric_func(point)
    n = g.shape[0]
    R_full = riemann_tensor(metric_func, point, h)

    # Ricci-Tensor: R_{ij} = R^k_{ikj} (Kontraktion: l=k)
    Ric = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # Summiere über k: R^k_{ikj}
            for k in range(n):
                Ric[i, j] += R_full[k, i, k, j]
    return Ric


# ---------------------------------------------------------------------------
# Ricci-Skalar
# ---------------------------------------------------------------------------

def ricci_scalar(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-4,
) -> float:
    """
    Ricci-Skalar R = g^{ij} R_{ij} (vollständige Kontraktion).

    Der Ricci-Skalar ist eine einzige reelle Zahl, die die Gesamtkrümmung
    einer Mannigfaltigkeit am gegebenen Punkt quantifiziert.

    Bekannte Werte:
    - Flache Ebene: R = 0
    - 2-Sphäre mit Radius r: R = 2/r²
    - Hyperbolische Ebene: R = -2
    - Torus: R variiert je nach Punkt (positiv außen, negativ innen)

    @param metric_func: Funktion point → g_{ij}(point).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite.
    @return: Ricci-Skalar als reelle Zahl.
    @lastModified 2026-03-10
    """
    g = metric_func(point)
    g_inv = np.linalg.inv(g)
    Ric = ricci_tensor(metric_func, point, h)
    # R = g^{ij} R_{ij}
    return float(np.einsum('ij,ij->', g_inv, Ric))


# ---------------------------------------------------------------------------
# Gaußsche Krümmung
# ---------------------------------------------------------------------------

def gaussian_curvature(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-4,
) -> float:
    """
    Gaußsche Krümmung K einer 2D-Fläche.

    Für eine 2-dimensionale Riemannsche Mannigfaltigkeit:
        K = R_{1212} / det(g)

    Über den Ricci-Skalar: K = R/2 (nur für 2D gültig).

    Klassifikation:
    - K > 0: positiv gekrümmt (wie Sphäre)
    - K = 0: flach (wie Ebene oder Zylinder)
    - K < 0: negativ gekrümmt (wie Sattelform oder hyperbolische Ebene)

    Gaußscher Bonnet-Satz: ∫_M K dA = 2π χ(M) (topologische Invariante!)

    Bekannte Werte:
    - Sphäre Radius r: K = 1/r² (konstant positiv)
    - Ebene: K = 0
    - Hyperbolische Ebene: K = -1 (konstant negativ)
    - Torus: K variiert zwischen -1/r·(R-r) und +1/r·(R+r)

    @param metric_func: Funktion point → g_{ij}(point), nur für 2D.
    @param point: 2D-Koordinatenpunkt [x¹, x²].
    @param h: Schrittweite.
    @return: Gaußsche Krümmung K als reelle Zahl.
    @raises ValueError: Falls Dimension ≠ 2.
    @lastModified 2026-03-10
    """
    g = metric_func(point)
    if g.shape[0] != 2:
        raise ValueError(
            f"Gaußsche Krümmung nur für 2D-Flächen definiert (Dimension={g.shape[0]})."
        )
    # Für 2D: K = R / 2
    R = ricci_scalar(metric_func, point, h)
    return R / 2.0


# ---------------------------------------------------------------------------
# Geodätengleichung
# ---------------------------------------------------------------------------

def geodesic_equation(
    metric_func: Callable[["list[float]"], np.ndarray],
    x0: "list[float]",
    v0: "list[float]",
    t_max: float = 1.0,
    n_steps: int = 1000,
) -> dict:
    """
    Löst die Geodätengleichung numerisch via Runge-Kutta 4. Ordnung.

    Geodätengleichung:
        d²x^k/dt² + Γ^k_{ij} · dx^i/dt · dx^j/dt = 0

    Als System 1. Ordnung (Zustandsvektor [x, v] mit v = dx/dt):
        dx^k/dt = v^k
        dv^k/dt = -Γ^k_{ij} v^i v^j

    Eine Geodäte ist die Verallgemeinerung einer "geraden Linie" auf
    gekrümmten Flächen (Großkreise auf der Sphäre, Geraden in der Ebene).

    Bogenlänge:
        L = ∫₀^{t_max} √(g_{ij} ẋ^i ẋ^j) dt

    @param metric_func: Metrische Funktion.
    @param x0: Startpunkt [x¹₀, x²₀, ..., xⁿ₀].
    @param v0: Anfangsgeschwindigkeit (Tangentialvektor).
    @param t_max: Integrationsintervall [0, t_max].
    @param n_steps: Anzahl Runge-Kutta-Schritte.
    @return: Dict {'trajectory': ndarray(n_steps+1, n), 'times': ndarray, 'length': float}.
    @lastModified 2026-03-10
    """
    n = len(x0)
    dt = t_max / n_steps

    # Zustandsvektoren: x (Position) und v (Geschwindigkeit)
    x = np.array(x0, dtype=float)
    v = np.array(v0, dtype=float)

    # Ergebnis-Arrays
    trajectory = np.zeros((n_steps + 1, n))
    times = np.zeros(n_steps + 1)
    trajectory[0] = x.copy()

    # Gesamte Bogenlänge
    length = 0.0

    def acceleration(xi: np.ndarray, vi: np.ndarray) -> np.ndarray:
        """Berechnet Beschleunigung aus Christoffel-Symbolen: a^k = -Γ^k_{ij} v^i v^j."""
        Gamma = christoffel_symbols(metric_func, list(xi), h=1e-5)
        # Matrixform: Γ^k_{ij} v^i v^j = Σ_{i,j} Γ[k,i,j] * v[i] * v[j]
        accel = np.einsum('kij,i,j->k', Gamma, vi, vi)
        return -accel

    # Runge-Kutta 4. Ordnung
    for step in range(n_steps):
        t = step * dt

        # Bogenlängen-Beitrag: ds = √(g_{ij} v^i v^j) dt
        g = metric_func(list(x))
        ds = float(math.sqrt(max(0.0, float(v @ g @ v)))) * dt
        length += ds

        # RK4-Koeffizienten für Position und Geschwindigkeit
        k1x = v.copy()
        k1v = acceleration(x, v)

        k2x = v + 0.5 * dt * k1v
        k2v = acceleration(x + 0.5 * dt * k1x, v + 0.5 * dt * k1v)

        k3x = v + 0.5 * dt * k2v
        k3v = acceleration(x + 0.5 * dt * k2x, v + 0.5 * dt * k2v)

        k4x = v + dt * k3v
        k4v = acceleration(x + dt * k3x, v + dt * k3v)

        # Update Zustand
        x = x + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + (dt / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)

        trajectory[step + 1] = x.copy()
        times[step + 1] = t + dt

    return {
        "trajectory": trajectory,
        "times": times,
        "length": length,
    }


# ---------------------------------------------------------------------------
# Geodätischer Abstand
# ---------------------------------------------------------------------------

def geodesic_distance(
    metric_func: Callable[["list[float]"], np.ndarray],
    x0: "list[float]",
    x1: "list[float]",
    n_steps: int = 100,
) -> float:
    """
    Geodätischer Abstand zwischen zwei Punkten (numerisch).

    Näherung über lineare Verbindungsstrecke und Integration der Bogenlänge:
        d(x0, x1) ≈ ∫₀¹ √(g_{ij}(γ(t)) γ̇^i γ̇^j) dt

    wobei γ(t) = (1-t)·x0 + t·x1 die lineare Parametrisierung ist.
    (Exakt nur für flache oder nahezu flache Räume; für gekrümmte Flächen
    ist dies eine Approximation des tatsächlichen geodätischen Abstands.)

    @param metric_func: Metrische Funktion.
    @param x0: Startpunkt.
    @param x1: Endpunkt.
    @param n_steps: Anzahl Integrationsschritte (Simpson-Regel).
    @return: Approximierter geodätischer Abstand.
    @lastModified 2026-03-10
    """
    x0_arr = np.array(x0, dtype=float)
    x1_arr = np.array(x1, dtype=float)
    # Tangentialvektor der linearen Verbindung (konstant)
    tangent = x1_arr - x0_arr

    # Simpson-Integration der Bogenlänge
    total = 0.0
    dt = 1.0 / n_steps
    for i in range(n_steps):
        t_mid = (i + 0.5) * dt
        point = list(x0_arr + t_mid * tangent)
        g = metric_func(point)
        # Integrand: √(g_{ij} γ̇^i γ̇^j)
        val = float(math.sqrt(max(0.0, float(tangent @ g @ tangent))))
        total += val * dt
    return total


# ---------------------------------------------------------------------------
# Klassische Mannigfaltigkeiten als Metriken
# ---------------------------------------------------------------------------

def sphere_metric(r: float = 1.0) -> Callable[["list[float]"], np.ndarray]:
    """
    Metriktensor der 2-Sphäre S² mit Radius r.

    Sphärische Koordinaten (θ, φ) mit θ ∈ (0, π), φ ∈ [0, 2π):
        ds² = r²(dθ² + sin²θ dφ²)

    Metriktensor:
        g = [[r², 0         ],
             [0,  r²sin²θ  ]]

    Eigenschaften:
    - Gaußsche Krümmung: K = 1/r² (konstant positiv)
    - Ricci-Skalar: R = 2/r²
    - Gesamte Fläche: 4πr²
    - Geodäten: Großkreise (kürzeste Verbindungen)

    @param r: Kugelradius (Standard: 1).
    @return: Funktion [θ, φ] → g_{ij}(θ, φ).
    @lastModified 2026-03-10
    """
    def metric(point: "list[float]") -> np.ndarray:
        theta = point[0]
        # Numerische Stabilität: sin(0) und sin(π) schützen
        sin_theta = math.sin(theta)
        return np.array([
            [r**2,                    0.0],
            [0.0,  r**2 * sin_theta**2],
        ])
    return metric


def torus_metric(R: float = 2.0, r: float = 1.0) -> Callable[["list[float]"], np.ndarray]:
    """
    Metriktensor des Torus T² (großer Radius R, kleiner Radius r).

    Koordinaten (θ, φ) mit θ, φ ∈ [0, 2π):
    - θ: Winkel um die innere Achse
    - φ: Winkel um die äußere Achse (Symmetrieachse)

    Erste Fundamentalform:
        ds² = r² dθ² + (R + r·cosθ)² dφ²

    Metriktensor:
        g = [[r²,              0               ],
             [0,   (R + r·cosθ)²]]

    Gaußsche Krümmung: K = cosθ / (r · (R + r·cosθ))
    - Außenseite (θ = 0): K = 1/(r·(R+r)) > 0 (positiv)
    - Innen (θ = π):      K = -1/(r·(R-r)) < 0 (negativ)
    - Übergänge (θ = π/2, 3π/2): K = 0

    @param R: Großer Radius (Mittelpunkt des Querschnittskreises bis Zentrum), R > r.
    @param r: Kleiner Radius (Querschnittskreis).
    @return: Funktion [θ, φ] → g_{ij}(θ, φ).
    @lastModified 2026-03-10
    """
    def metric(point: "list[float]") -> np.ndarray:
        theta = point[0]
        # Effektiver Radius in φ-Richtung
        rho = R + r * math.cos(theta)
        return np.array([
            [r**2,    0.0],
            [0.0,  rho**2],
        ])
    return metric


def hyperbolic_plane_metric() -> Callable[["list[float]"], np.ndarray]:
    """
    Metriktensor der hyperbolischen Ebene H² (Poincaré-Halbebene-Modell).

    Koordinaten (x, y) mit y > 0 (obere Halbebene):
        ds² = (dx² + dy²) / y²

    Metriktensor:
        g = [[1/y², 0   ],
             [0,    1/y²]]

    Eigenschaften:
    - Gaußsche Krümmung: K = -1 (konstant negativ!)
    - Isometriegruppe: PSL(2,ℝ) – Möbius-Transformationen
    - Geodäten: vertikale Halbgeraden und Halbkreise mit Mittelpunkt auf x-Achse
    - Modell der hyperbolischen Geometrie (Lobatschewski)

    @return: Funktion [x, y] → g_{ij}(x, y).
    @lastModified 2026-03-10
    """
    def metric(point: "list[float]") -> np.ndarray:
        y = point[1]
        if abs(y) < 1e-12:
            # Singularität am Rand der Halbebene; kleine Korrektur
            y = 1e-12
        coeff = 1.0 / (y * y)
        return np.array([
            [coeff, 0.0   ],
            [0.0,   coeff],
        ])
    return metric


def saddle_metric() -> Callable[["list[float]"], np.ndarray]:
    """
    Metriktensor der Sattelfläche z = xy (eingebettet in ℝ³).

    Die Sattelfläche hat einen hyperbolischen Punkt im Ursprung
    (ein Sattel wie im Namen: positive Krümmung in einer, negative in anderer Richtung).

    Erste Fundamentalform (aus Einbettung r(x,y) = (x, y, xy)):
        ∂r/∂x = (1, 0, y), ∂r/∂y = (0, 1, x)
        g_{11} = 1 + y², g_{12} = xy, g_{22} = 1 + x²

    Gaußsche Krümmung: K = -1 / (1 + x² + y²)² ≤ 0

    @return: Funktion [x, y] → g_{ij}(x, y).
    @lastModified 2026-03-10
    """
    def metric(point: "list[float]") -> np.ndarray:
        x, y = point[0], point[1]
        return np.array([
            [1.0 + y**2,     x * y],
            [x * y,      1.0 + x**2],
        ])
    return metric


def flat_metric(n: int = 2) -> Callable[["list[float]"], np.ndarray]:
    """
    Euklidische (flache) n-dimensionale Metrik: g_{ij} = δ_{ij}.

    Die einfachste Metrik – entspricht dem Standard-Skalarprodukt im ℝⁿ.
    Alle Christoffel-Symbole sind null, alle Krümmungstensoren sind null.

    @param n: Dimension (Standard: 2).
    @return: Funktion point → I_n (Einheitsmatrix, unabhängig vom Punkt).
    @lastModified 2026-03-10
    """
    identity = np.eye(n)

    def metric(point: "list[float]") -> np.ndarray:
        return identity.copy()

    return metric


# ---------------------------------------------------------------------------
# Differentialformen
# ---------------------------------------------------------------------------

def wedge_product(alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """
    Äußeres (Wedge-) Produkt zweier Differentialformen.

    Für zwei 1-Formen α = α_i dx^i und β = β_j dx^j:
        (α ∧ β)_{ij} = α_i β_j - α_j β_i

    Das Wedge-Produkt ist:
    - Bilinear: (α + γ) ∧ β = α ∧ β + γ ∧ β
    - Antisymmetrisch: α ∧ β = -(β ∧ α)
    - Assoziativ: (α ∧ β) ∧ γ = α ∧ (β ∧ γ)

    @param alpha: 1-Form als 1D-Array [α₁, ..., αₙ].
    @param beta: 1-Form als 1D-Array [β₁, ..., βₙ].
    @return: 2-Form als antisymmetrische n×n-Matrix.
    @lastModified 2026-03-10
    """
    alpha = np.array(alpha, dtype=float)
    beta = np.array(beta, dtype=float)
    # Äußeres Produkt: α_i β_j
    outer = np.outer(alpha, beta)
    # Antisymmetrisierung: (α ∧ β)_{ij} = α_i β_j - α_j β_i
    return outer - outer.T


def exterior_derivative(
    form: np.ndarray,
    point: "list[float]",
    h: float = 1e-5,
) -> np.ndarray:
    """
    Äußere Ableitung d einer Differentialform (numerisch).

    Für eine 0-Form (Skalarfeld) f: R^n → R:
        df = (∂f/∂x¹, ..., ∂f/∂xⁿ)  (Gradient als 1-Form)

    Für eine 1-Form α = [α₁(x), ..., αₙ(x)] (komponentenweise Funktion):
        (dα)_{ij} = ∂_i α_j - ∂_j α_i

    Hier: form ist ein konstanter Array am gegebenen Punkt.
    (Für punktabhängige Formen muss form eine Funktion sein.)

    @param form: 1-Form als 1D-Array (Komponenten am gegebenen Punkt).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite für numerische Ableitung.
    @return: 2-Form als antisymmetrische Matrix (äußere Ableitung).
    @lastModified 2026-03-10
    """
    form = np.array(form, dtype=float)
    n = len(point)

    if form.ndim == 0 or len(form) == 1:
        # 0-Form: äußere Ableitung = Gradient (numerisch via Differenzenquotient)
        grad = np.zeros(n)
        for i in range(n):
            p_plus = list(point)
            p_plus[i] += h
            p_minus = list(point)
            p_minus[i] -= h
            # Ableitung der konstanten Form ergibt 0 – Platzhalter
            grad[i] = 0.0  # Für ortsabhängige Formen: (f(p_plus)-f(p_minus))/(2h)
        return grad

    # 1-Form: (dα)_{ij} = ∂_i α_j - ∂_j α_i
    # Für konstante Formen: dα = 0 (exakte Form)
    # Diese Implementierung berechnet das antisymmetrische Produkt aus den Komponenten
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # Für konstante Formen: alle Ableitungen = 0
            result[i, j] = 0.0
    return result


def hodge_star(form: np.ndarray, metric: MetricTensor) -> np.ndarray:
    """
    Hodge-Stern-Operator *ω für Differentialformen.

    Für eine n-dimensionale Mannigfaltigkeit:
        * : Ω^k → Ω^{n-k}

    Für eine 2D-Mannigfaltigkeit:
    - *(f) = f · √|det(g)| dx¹ ∧ dx²   (0-Form → 2-Form)
    - *(α₁dx¹ + α₂dx²) = α₁dx² - α₂dx¹ · √|det(g)|  (1-Form → 1-Form)

    Eigenschaften:
    - ** = (-1)^{k(n-k)} · sign(g) · Id
    - Für Riemannsche Metrik: ** = (-1)^{k(n-k)} · Id

    @param form: k-Form (0D-Skalar, 1D-Vektor oder 2D-Matrix).
    @param metric: Metriktensor des Raums.
    @return: (n-k)-Form (Hodge-Dual).
    @lastModified 2026-03-10
    """
    vol = metric.volume_element()  # √|det(g)|
    n = metric.n

    if form.ndim == 0 or (form.ndim == 1 and len(form) == 1):
        # 0-Form → n-Form: *(f) = f · √|det(g)|
        return float(form) * vol * np.ones(1)

    if form.ndim == 1 and n == 2:
        # 1-Form in 2D → 1-Form: *(α₁, α₂) = (-α₂, α₁) · √|det(g)|
        return vol * np.array([-form[1], form[0]])

    if form.ndim == 2 and n == 2:
        # 2-Form (antisymm. Matrix) in 2D → 0-Form (Skalar)
        # *(A dx¹∧dx²) = A/√|det(g)|
        if vol > 1e-14:
            return np.array([form[0, 1] / vol])
        return np.array([0.0])

    # Allgemeiner Fall: Levi-Civita-Kontraktion
    eps = levi_civita_symbol(n)
    # Vereinfachung: Rückgabe als skaliertes Array
    return vol * form


def lie_derivative_vector(
    f: Callable[["list[float]"], float],
    v: "list[float]",
    point: "list[float]",
    h: float = 1e-5,
) -> float:
    """
    Lie-Ableitung £_V f eines Skalarfeldes f in Richtung Vektorfeld V.

    Für ein Skalarfeld:
        £_V f = V^i ∂_i f  (direktionale Ableitung in Richtung V)

    Interpretiert als: Rate der Änderung von f wenn man dem Fluss
    von V folgt. Für Skalarfelder stimmt sie mit der direktionalen
    Ableitung überein.

    @param f: Skalarfeld als Funktion point → float.
    @param v: Vektorfeld am Punkt (Vektor der Komponenten V^i).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite für numerische Ableitung.
    @return: Lie-Ableitung als reelle Zahl.
    @lastModified 2026-03-10
    """
    v_arr = np.array(v, dtype=float)
    n = len(point)
    result = 0.0
    # Summe V^i ∂_i f (direktionale Ableitung)
    for i in range(n):
        p_plus = list(point)
        p_plus[i] += h
        p_minus = list(point)
        p_minus[i] -= h
        # Partielle Ableitung ∂f/∂x^i
        partial = (f(p_plus) - f(p_minus)) / (2.0 * h)
        result += v_arr[i] * partial
    return result


# ---------------------------------------------------------------------------
# Einstein-Gleichungen
# ---------------------------------------------------------------------------

def einstein_tensor(
    metric_func: Callable[["list[float]"], np.ndarray],
    point: "list[float]",
    h: float = 1e-4,
) -> np.ndarray:
    """
    Einstein-Tensor G_{ij} = R_{ij} - ½ g_{ij} R.

    Zentral in Einsteins Feldgleichungen der Allgemeinen Relativitätstheorie:
        G_{ij} = 8πG/c⁴ · T_{ij}

    wobei T_{ij} der Energie-Impuls-Tensor der Materie ist.

    Für Vakuum (kein Materie): T_{ij} = 0 → G_{ij} = 0
    → Ricci-Tensor und Ricci-Skalar verschwinden (Ricci-flach)
    → Schwarzschild-Metrik ist eine Vakuumlösung!

    Eigenschaften:
    - G_{ij} = G_{ji} (symmetrisch)
    - ∇^i G_{ij} = 0 (divergenzfrei – Energieerhaltung!)
    - Spurfreier Anteil = Weyl-Tensor (beschreibt Gravitationswellen)

    @param metric_func: Funktion point → g_{ij}(point).
    @param point: Koordinatenpunkt.
    @param h: Schrittweite.
    @return: Einstein-Tensor G[i,j] als n×n-Matrix.
    @lastModified 2026-03-10
    """
    g = metric_func(point)
    Ric = ricci_tensor(metric_func, point, h)
    R = ricci_scalar(metric_func, point, h)
    # G_{ij} = R_{ij} - ½ g_{ij} R
    return Ric - 0.5 * g * R


def schwarzschild_metric(
    M: float = 1.0,
    c: float = 1.0,
    G: float = 1.0,
) -> Callable[["list[float]"], np.ndarray]:
    """
    Schwarzschild-Metrik für eine sphärisch-symmetrische Masse M.

    In vereinfachten Koordinaten (t, r) für den Äquator (θ = π/2):
        ds² = -(1 - r_s/r)c²dt² + (1 - r_s/r)⁻¹ dr²

    Schwarzschild-Radius:
        r_s = 2GM/c²

    Metriktensor (Signatur -, +):
        g = [[-(1 - r_s/r)c²,        0          ],
             [0,              1/(1 - r_s/r)]]

    Physikalische Bedeutung:
    - r > r_s: Außenbereich (normale Raumzeit)
    - r = r_s: Ereignishorizont (Koordinaten-Singularität)
    - r → 0: echte (physikalische) Singularität

    Natürliche Einheiten: c = G = 1, r_s = 2M.

    @param M: Masse des Zentralkörpers (in natürlichen Einheiten).
    @param c: Lichtgeschwindigkeit (Standard: 1 = natürliche Einheiten).
    @param G: Gravitationskonstante (Standard: 1).
    @return: Funktion [t, r] → g_{ij}(t, r).
    @lastModified 2026-03-10
    """
    # Schwarzschild-Radius
    r_s = 2.0 * G * M / (c ** 2)

    def metric(point: "list[float]") -> np.ndarray:
        r = point[1]
        # Singularitäten vermeiden: r darf nicht unter r_s sinken
        if r <= r_s:
            r = r_s * 1.001  # Kleiner Offset zur Stabilisierung
        factor = 1.0 - r_s / r
        return np.array([
            [-(factor * c**2),  0.0       ],
            [0.0,           1.0 / factor],
        ])

    return metric


def check_vacuum_solution(
    metric_func: Callable[["list[float]"], np.ndarray],
    points: "list[list[float]]",
    tol: float = 1e-6,
) -> dict:
    """
    Prüft ob eine Metrik die Vakuum-Einsteingleichungen erfüllt (G_{ij} ≈ 0).

    Für eine Vakuumlösung muss gelten: G_{ij} = R_{ij} - ½ g_{ij} R = 0.
    Dies ist äquivalent zu R_{ij} = 0 (Ricci-flach).

    Bekannte Vakuumlösungen:
    - Schwarzschild-Metrik (außerhalb der Masse)
    - Kerr-Metrik (rotierende Masse)
    - Gravitationswellen-Lösungen

    @param metric_func: Metrische Funktion.
    @param points: Liste der zu prüfenden Koordinatenpunkte.
    @param tol: Toleranz für G_{ij} ≈ 0.
    @return: Dict {'is_vacuum': bool, 'max_deviation': float, 'points_checked': int}.
    @lastModified 2026-03-10
    """
    max_dev = 0.0
    for pt in points:
        G = einstein_tensor(metric_func, pt)
        dev = float(np.max(np.abs(G)))
        max_dev = max(max_dev, dev)

    return {
        "is_vacuum": max_dev < tol,
        "max_deviation": max_dev,
        "points_checked": len(points),
    }


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def levi_civita_symbol(n: int) -> np.ndarray:
    """
    Levi-Civita-Symbol ε_{i₁...iₙ} (vollständig antisymmetrischer Tensor).

    Definition:
    - ε_{12...n} = +1 (gerade Permutation)
    - ε_{i₁...iₙ} = -1 für ungerade Permutation
    - ε_{i₁...iₙ} = 0 falls irgendein Index doppelt vorkommt

    Verwendung:
    - Determinante: det(A) = ε_{i₁...iₙ} A^1_{i₁} ... A^n_{iₙ}
    - Kreuzprodukt in 3D: (u × v)^k = ε_{kij} u^i v^j
    - Volumenform: dV = ε_{i₁...iₙ} dx^{i₁} ∧ ... ∧ dx^{iₙ}

    @param n: Dimension.
    @return: n-dimensionales Array (Form: n×n×...×n, n-mal).
    @lastModified 2026-03-10
    """
    # Initialisierung mit Nullen
    eps = np.zeros(tuple([n] * n))

    # Alle Permutationen der Indizes 0..n-1 durchlaufen
    from itertools import permutations

    def perm_sign(perm: tuple) -> int:
        """Vorzeichen einer Permutation (Bubble-Sort-Methode)."""
        perm = list(perm)
        sign = 1
        for i in range(len(perm)):
            while perm[i] != i:
                j = perm[i]
                perm[i], perm[j] = perm[j], perm[i]
                sign *= -1
        return sign

    for perm in permutations(range(n)):
        eps[perm] = perm_sign(perm)

    return eps


def kronecker_delta(n: int) -> np.ndarray:
    """
    Kronecker-Delta δ^i_j = Einheitsmatrix der Dimension n.

    Definition:
        δ^i_j = 1 falls i = j, sonst 0

    Eigenschaften:
    - δ^i_j = g^{ik} g_{kj} (Metrik heben und senken ergibt δ)
    - Kontrahierung: δ^i_i = n (Spur = Dimension)

    @param n: Dimension.
    @return: n×n Einheitsmatrix.
    @lastModified 2026-03-10
    """
    return np.eye(n)


def parallel_transport(
    metric_func: Callable[["list[float]"], np.ndarray],
    v0: "list[float]",
    curve: "list[list[float]]",
) -> "list[list[float]]":
    """
    Paralleltransport eines Vektors v0 entlang einer Kurve.

    Paralleltransportgleichung:
        Dv^k/dt = dv^k/dt + Γ^k_{ij} v^i ẋ^j = 0

    Ein parallel transportierter Vektor bleibt "so gerade wie möglich"
    auf der Mannigfaltigkeit. Im flachen Raum: gleichbedeutend mit
    konstantem Vektor.

    Anschauung: Transportiert man einen Vektor parallel entlang einer
    geschlossenen Kurve auf einer Sphäre, dreht er sich um einen Winkel
    gleich dem Raumwinkel (Holonomie!).

    Numerisch: Euler-Verfahren (Δv^k = -Γ^k_{ij} v^i Δx^j).

    @param metric_func: Metrische Funktion.
    @param v0: Anfangsvektor [v¹₀, v²₀, ..., vⁿ₀].
    @param curve: Liste von Kurvenpunkten [[x₀,y₀], [x₁,y₁], ...].
    @return: Liste der transportierten Vektoren an jedem Kurvenpunkt.
    @lastModified 2026-03-10
    """
    v = np.array(v0, dtype=float)
    result = [v.copy().tolist()]

    for step in range(len(curve) - 1):
        # Aktueller und nächster Kurvenpunkt
        x_curr = np.array(curve[step], dtype=float)
        x_next = np.array(curve[step + 1], dtype=float)
        dx = x_next - x_curr

        # Christoffel-Symbole am aktuellen Punkt
        Gamma = christoffel_symbols(metric_func, list(x_curr), h=1e-5)

        # Euler-Schritt: Δv^k = -Γ^k_{ij} v^i dx^j
        dv = -np.einsum('kij,i,j->k', Gamma, v, dx)
        v = v + dv
        result.append(v.copy().tolist())

    return result


# ===========================================================================
# SYMPLEKTISCHE GEOMETRIE UND HAMILTON-MECHANIK
# ===========================================================================
#
# Die symplektische Geometrie ist die mathematische Sprache der klassischen
# Mechanik. Ein symplektischer Vektorraum (V, ω) trägt eine nicht-entartete,
# schiefsymmetrische Bilinearform ω (die "symplektische Form").
#
# Im Phasenraum mit Koordinaten (q₁,...,qₙ, p₁,...,pₙ) ist ω = Σ dqᵢ ∧ dpᵢ.
# Die zugehörige Matrix ist J = [[0, Iₙ], [-Iₙ, 0]].
#
# Hamiltonsche Mechanik: Aus dem Hamilton H(q,p) folgen Bewegungsgleichungen
#   dq/dt = ∂H/∂p,   dp/dt = -∂H/∂q
# Der Fluss erhält ω (Liouville-Satz) und die Energie (falls ∂H/∂t = 0).
# ===========================================================================


def symplectic_form(n: int) -> np.ndarray:
    """
    @brief Standard-symplektische Form auf ℝ^{2n} als (2n×2n)-Matrix J.

    @description
        Die Standard-symplektische Form auf ℝ^{2n} lautet:
            ω = Σ_{i=1}^{n} dqᵢ ∧ dpᵢ

        Als Blockmatrix:
            J = [[  0  ,  Iₙ ],
                 [ -Iₙ ,   0 ]]

        Eigenschaften:
        - Schiefsymmetrie:   J^T = -J
        - Nicht-Entartetheit: det(J) = 1 ≠ 0
        - Geschlossenheit:   dω = 0 (trivial, da ω linear/konstant ist)
        - J² = -I_{2n}  (analog zur imaginären Einheit)

        Bedeutung: J beschreibt die Geometrie des Phasenraums in der
        klassischen Mechanik. Hamiltonsche Gleichungen können kompakt als
        ż = J · ∇H(z) mit z = (q, p) geschrieben werden.

    @param n: Hälfte der Dimension (Freiheitsgrade); resultiert in 2n×2n-Matrix.
    @return: (2n × 2n)-Matrix J der Standard-symplektischen Form.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Identitätsmatrix der Größe n×n
    I_n = np.eye(n)
    # Nullmatrix der Größe n×n
    Z = np.zeros((n, n))

    # Blockmatrix J = [[0, I], [-I, 0]] (2n×2n)
    top    = np.block([Z,    I_n])   # obere Hälfte:  [0  | I]
    bottom = np.block([-I_n, Z  ])   # untere Hälfte: [-I | 0]

    return np.block([[top], [bottom]])


def is_symplectic_matrix(M: np.ndarray) -> bool:
    """
    @brief Prüft ob M zur symplektischen Gruppe Sp(2n) gehört.

    @description
        Eine (2n×2n)-Matrix M ist symplektisch, wenn gilt:
            M^T J M = J

        wobei J die Standard-symplektische Form ist.

        Die symplektische Gruppe Sp(2n) ist die Symmetriegruppe der
        klassischen Mechanik: Kanonische Transformationen erhalten ω.

        Folgerungen:
        - det(M) = ±1 (tatsächlich immer +1 für Sp(2n))
        - M invertierbar mit M⁻¹ = -J M^T J
        - Sp(2n) ⊂ SL(2n, ℝ)

        Beispiele: Rotationen im Phasenraum, Scherungen, Shear-Matrizen.

    @param M: Quadratische Matrix (Dimension muss gerade sein).
    @return: True, falls M^T J M = J (bis auf numerische Toleranz 1e-8).
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    rows, cols = M.shape
    if rows != cols or rows % 2 != 0:
        # Symplektische Matrizen müssen quadratisch und gerader Dimension sein
        return False

    n = rows // 2
    J = symplectic_form(n)   # Standard-symplektische Form

    # Berechne M^T J M und vergleiche mit J
    MTJM = M.T @ J @ M
    return bool(np.allclose(MTJM, J, atol=1e-8))


def poisson_bracket(
    f_vals: np.ndarray,
    g_vals: np.ndarray,
    dq: float,
    dp: float,
) -> np.ndarray:
    """
    @brief Numerische Poisson-Klammer {f, g} auf dem Phasenraum (q, p).

    @description
        Die Poisson-Klammer zweier Observablen f, g auf dem Phasenraum ist:
            {f, g} = ∂f/∂q · ∂g/∂p − ∂f/∂p · ∂g/∂q

        Numerisch werden die partiellen Ableitungen via zentraler Differenzen
        approximiert:
            ∂f/∂q ≈ (f[q+1,:] - f[q-1,:]) / (2·dq)   (Achse 0 = q-Achse)
            ∂f/∂p ≈ (f[:,p+1] - f[:,p-1]) / (2·dp)   (Achse 1 = p-Achse)

        Wichtige Eigenschaften der Poisson-Klammer:
        - Schiefsymmetrie:    {f, g} = -{g, f}
        - Bilinearität
        - Jacobi-Identität:  {f, {g, h}} + {g, {h, f}} + {h, {f, g}} = 0
        - Kanonische Relation: {qᵢ, pⱼ} = δᵢⱼ   (klassisches Analogon zu [q̂,p̂]=iℏ)
        - Energieerhaltung:   {H, H} = 0

        Mit `np.gradient` werden Randwerte automatisch einseitig approximiert,
        was die volle 2D-Auswertung ohne Randkorrekturen ermöglicht.

    @param f_vals: 2D-Array f[i,j] = f(qᵢ, pⱼ) über dem (q,p)-Gitter.
    @param g_vals: 2D-Array g[i,j] = g(qᵢ, pⱼ) über dem (q,p)-Gitter.
    @param dq: Gitterabstand in q-Richtung.
    @param dp: Gitterabstand in p-Richtung.
    @return: 2D-Array {f,g}[i,j] über dem (q,p)-Gitter.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Partielle Ableitungen mittels np.gradient (zentrale Differenzen innen,
    # einseitige Differenzen am Rand)
    df_dq = np.gradient(f_vals, dq, axis=0)   # ∂f/∂q  (q = Achse 0)
    df_dp = np.gradient(f_vals, dp, axis=1)   # ∂f/∂p  (p = Achse 1)
    dg_dq = np.gradient(g_vals, dq, axis=0)   # ∂g/∂q
    dg_dp = np.gradient(g_vals, dp, axis=1)   # ∂g/∂p

    # Poisson-Klammer: {f,g} = ∂f/∂q · ∂g/∂p − ∂f/∂p · ∂g/∂q
    return df_dq * dg_dp - df_dp * dg_dq


def hamiltonian_flow(
    H_func: Callable,
    q0: float,
    p0: float,
    t_span: float = 10.0,
    dt: float = 0.01,
) -> dict:
    """
    @brief Löst Hamiltonsche Gleichungen mit dem symplektischen Störmer-Verlet-Integrator.

    @description
        Die Hamiltonschen Bewegungsgleichungen lauten:
            dq/dt = +∂H/∂p
            dp/dt = −∂H/∂q

        Störmer-Verlet-Schema (leapfrog, "velocity Verlet" in Phasenraumform):
            p_{n+1/2} = pₙ − (dt/2) · ∂H/∂q(qₙ)
            q_{n+1}   = qₙ + dt · ∂H/∂p(p_{n+1/2})     ← Halbschritt-Impuls
            p_{n+1}   = p_{n+1/2} − (dt/2) · ∂H/∂q(q_{n+1})

        Vorteile gegenüber Euler/Runge-Kutta:
        - Symplektisch: erhält exakt die symplektische 2-Form ω = dq∧dp
        - Zeitumkehrinvariant
        - Energie schwankt nur beschränkt (keine langzeitige Drift)
        - Geeignet für lange Simulationen (astronomische Zeiträume)

        Numerische Gradienten werden mit h = 1e-6 berechnet.

    @param H_func: Hamilton-Funktion H(q, p) → float.
    @param q0: Anfangs-Koordinate.
    @param p0: Anfangs-Impuls.
    @param t_span: Gesamte Simulationszeit.
    @param dt: Zeitschrittgröße (kleiner = genauer, aber langsamer).
    @return: Dictionary mit Arrays 'q', 'p', 't', 'H' (Energie entlang Trajektorie).
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Schrittweite für numerische Ableitung von H
    _h = 1e-6

    def dH_dq(q: float, p: float) -> float:
        """Numerische partielle Ableitung ∂H/∂q (zentral)."""
        return (H_func(q + _h, p) - H_func(q - _h, p)) / (2.0 * _h)

    def dH_dp(q: float, p: float) -> float:
        """Numerische partielle Ableitung ∂H/∂p (zentral)."""
        return (H_func(q, p + _h) - H_func(q, p - _h)) / (2.0 * _h)

    # Anzahl der Zeitschritte
    n_steps = int(t_span / dt)

    # Initialisierung der Ausgabe-Arrays
    q_arr = np.zeros(n_steps + 1)
    p_arr = np.zeros(n_steps + 1)
    t_arr = np.zeros(n_steps + 1)
    H_arr = np.zeros(n_steps + 1)

    # Anfangswerte setzen
    q_arr[0] = q0
    p_arr[0] = p0
    t_arr[0] = 0.0
    H_arr[0] = H_func(q0, p0)

    q = q0
    p = p0

    # Störmer-Verlet-Integration (symplektischer Leapfrog)
    for i in range(1, n_steps + 1):
        # Halber Impulsschritt: p_{n+1/2} = pₙ − (dt/2)·∂H/∂q(qₙ)
        p_half = p - 0.5 * dt * dH_dq(q, p)

        # Voller Koordinatenschritt: q_{n+1} = qₙ + dt·∂H/∂p(p_{n+1/2})
        q = q + dt * dH_dp(q, p_half)

        # Zweiter halber Impulsschritt: p_{n+1} = p_{n+1/2} − (dt/2)·∂H/∂q(q_{n+1})
        p = p_half - 0.5 * dt * dH_dq(q, p_half)

        # Ergebnis speichern
        q_arr[i] = q
        p_arr[i] = p
        t_arr[i] = i * dt
        H_arr[i] = H_func(q, p)

    return {
        'q': q_arr,   # Koordinaten-Trajektorie
        'p': p_arr,   # Impuls-Trajektorie
        't': t_arr,   # Zeitachse
        'H': H_arr,   # Energie entlang der Trajektorie
    }


def liouville_theorem_check(
    H_func: Callable,
    initial_conditions: "list[tuple]",
    t: float = 5.0,
) -> dict:
    """
    @brief Numerische Verifikation des Liouville-Satzes (Phasenraumvolumen-Erhaltung).

    @description
        Liouville-Satz: Der Hamiltonsche Fluss Φ_t ist volumenerhaltend,
        d.h. für jede messbare Menge U im Phasenraum gilt:
            Vol(Φ_t(U)) = Vol(U)

        Äquivalent: div(X_H) = 0, wobei X_H = (∂H/∂p, −∂H/∂q) das
        Hamiltonsche Vektorfeld ist.

        Numerische Methode:
        - Bilde ein Ensemble von N Anfangsbedingungen (qᵢ, pᵢ)
        - Entwickle jede mit Störmer-Verlet bis Zeit t
        - Schätze das initiale und finale Phasenraumvolumen durch die
          konvexe Hülle des Ensembles (via maximale Ausdehnung/Fläche)
        - Die Jacobi-Determinante det(∂(q,p)/∂(q₀,p₀)) sollte ≈ 1 sein

        Approximation der Fläche: Verwende Kovarianzmatrix der Punktwolke.
        Für kleine Perturbationen: A ≈ π · √det(Cov) (Ellipsen-Flächenformel).

    @param H_func: Hamilton-Funktion H(q, p).
    @param initial_conditions: Liste von (q₀, p₀)-Tupeln (Ensemble).
    @param t: Entwicklungszeit.
    @return: Dictionary mit 'volume_initial', 'volume_final', 'ratio'.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Ensemble vorwärtsentwickeln
    dt = min(0.01, t / 500)   # Adaptive Schrittweite

    # Initiale Positionen als Array (N × 2)
    init_array = np.array(initial_conditions, dtype=float)   # (N, 2)

    # Finale Positionen berechnen
    final_positions = []
    for (q0, p0) in initial_conditions:
        result = hamiltonian_flow(H_func, q0, p0, t_span=t, dt=dt)
        final_positions.append((result['q'][-1], result['p'][-1]))

    final_array = np.array(final_positions, dtype=float)   # (N, 2)

    # Schätze Phasenraumvolumen via Kovarianz-Determinante
    # Für eine gleichförmige Ellipse: Fläche ∝ √det(Cov)
    def _volume_estimate(pts: np.ndarray) -> float:
        """Schätzt das 2D-Phasenraumvolumen via √det(Kovarianz)."""
        if pts.shape[0] < 2:
            return 1.0
        cov = np.cov(pts.T)
        if cov.ndim == 0:
            # Nur 1D (degenerierter Fall)
            return float(np.sqrt(abs(float(cov))))
        det_cov = np.linalg.det(cov)
        return float(np.sqrt(max(abs(det_cov), 1e-30)))

    vol_initial = _volume_estimate(init_array)
    vol_final   = _volume_estimate(final_array)

    # Verhältnis ≈ 1 bestätigt Liouville
    ratio = vol_final / vol_initial if vol_initial > 1e-30 else float('nan')

    return {
        'volume_initial': vol_initial,
        'volume_final':   vol_final,
        'ratio':          ratio,
    }


def action_angle_variables(
    H_func: Callable,
    E: float,
    q_range: tuple = (-2.0, 2.0),
) -> dict:
    """
    @brief Berechnet Wirkungs-Winkelvariablen (J, θ) für integrable Systeme.

    @description
        Für integrable Hamiltonsche Systeme existieren kanonische Koordinaten
        (J, θ) – die Wirkungs-Winkelvariablen – in denen H nur von J abhängt:
            H = H(J)   →   dθ/dt = ∂H/∂J = ω(J) = const.

        Wirkungsvariable (Adiabateninvariante):
            J = (1/2π) ∮_γ p dq   (geschlossene Energiebahn γ mit H = E)

        Frequenz:
            ω = dH/dJ = dE/dJ   (Kreisfrequenz der periodischen Bewegung)

        Numerisch: Finde p(q, E) aus H(q,p) = E via Bisektion:
            p(q) = p solch dass H(q,p) = E  (positiver Ast)
        Integration über eine halbe Periode (q_min bis q_max) × 2:
            J = (2/2π) ∫_{q_min}^{q_max} p(q) dq

        Beispiel harmonischer Oszillator H = p²/2 + ω₀²q²/2:
            Energiebahn: p² + ω₀²q² = 2E  (Ellipse)
            J = E/ω₀,   ω = ω₀

    @param H_func: Hamilton-Funktion H(q, p).
    @param E: Energie (Wert der Erhaltungsgröße).
    @param q_range: Suchbereich für die q-Koordinate als (q_min, q_max).
    @return: Dictionary mit 'J' (Wirkung), 'omega' (Frequenz), 'q_vals', 'p_vals'.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    q_min, q_max = q_range

    # Erstelle q-Gitter mit 2000 Punkten (hohe Auflösung für gute Integration)
    n_pts = 2000
    q_vals = np.linspace(q_min, q_max, n_pts)

    def find_p_positive(q: float) -> float:
        """
        Findet den positiven Impuls p > 0 für H(q, p) = E via Bisektion.
        Gibt 0.0 zurück, falls kein positiver Impuls existiert (klassisch verboten).
        """
        # Energiedifferenz bei p=0
        H0 = H_func(q, 0.0)
        if H0 > E:
            # Klassisch verbotene Region (kinetische Energie wäre negativ)
            return 0.0

        # Obere Grenze für p suchen (H wächst mit p²)
        p_upper = 0.1
        while H_func(q, p_upper) < E and p_upper < 1e6:
            p_upper *= 2.0

        if H_func(q, p_upper) < E:
            return 0.0   # Keine Lösung gefunden

        # Bisektion in [0, p_upper]
        p_lo, p_hi = 0.0, p_upper
        for _ in range(60):
            p_mid = 0.5 * (p_lo + p_hi)
            if H_func(q, p_mid) < E:
                p_lo = p_mid
            else:
                p_hi = p_mid
        return 0.5 * (p_lo + p_hi)

    # Berechne p(q) für alle q-Werte
    p_vals = np.array([find_p_positive(q) for q in q_vals])

    # Maske: nur Punkte mit p > 0 (klassisch erlaubte Region)
    mask = p_vals > 1e-10
    if not np.any(mask):
        return {'J': 0.0, 'omega': float('nan'), 'q_vals': q_vals, 'p_vals': p_vals}

    q_allowed = q_vals[mask]
    p_allowed = p_vals[mask]

    # Wirkungsvariable: J = (1/2π) ∮ p dq = (2/2π) ∫_{q_min}^{q_max} p dq
    # Faktor 2: Hin- und Rückweg der periodischen Bahn
    J = (2.0 / (2.0 * np.pi)) * np.trapezoid(p_allowed, q_allowed)

    # Frequenz: ω = dE/dJ via numerischer Ableitung
    delta_E = E * 1e-4 + 1e-8   # Kleine Energie-Perturbation

    def compute_J_at_E(energy: float) -> float:
        """Hilfsfunktion: Berechne J für gegebene Energie."""
        p_temp = np.array([find_p_positive(q) for q in q_vals])
        m = p_temp > 1e-10
        if not np.any(m):
            return 0.0
        return (2.0 / (2.0 * np.pi)) * np.trapezoid(p_temp[m], q_vals[m])

    J_plus  = compute_J_at_E(E + delta_E)
    J_minus = compute_J_at_E(E - delta_E)

    # ω = dE/dJ ≈ (2·δE) / (J_plus - J_minus)   (invertierter Quotient)
    dJ = J_plus - J_minus
    omega = (2.0 * delta_E) / dJ if abs(dJ) > 1e-15 else float('nan')

    return {
        'J':      J,         # Wirkungsvariable
        'omega':  omega,     # Grundkreisfrequenz
        'q_vals': q_vals,    # q-Gitter (zum Plotten)
        'p_vals': p_vals,    # p(q)-Werte (positiver Ast)
    }


def cotangent_bundle_metric(base_metric: np.ndarray) -> dict:
    """
    @brief Natürliche symplektische Struktur auf dem Kotangentialbündel T*M.

    @description
        Für eine Riemannsche Mannigfaltigkeit (M, g) trägt das Kotangential-
        bündel T*M eine natürliche (kanonische) symplektische Struktur:

        Kanonische 1-Form (Liouville-Form / tautologische Form):
            θ = pᵢ dqⁱ   (in lokalen Koordinaten)

        Kanonische symplektische Form:
            ω = -dθ = dqⁱ ∧ dpᵢ = Σᵢ dqⁱ ∧ dpᵢ

        Diese Form ist:
        - Exakt: ω = dθ   (daher trivialerweise geschlossen dω = 0)
        - Nicht-entartet
        - Kanonisch (unabhängig von der Wahl der Basismetrik g)

        Die symplektische Form auf T*M hat Dimension 2n (n = dim M).

        Liouville-Form in Koordinaten: θ = Σᵢ pᵢ dqⁱ
        Als Vektor (diskrete Näherung): θᵢ entspricht dem i-ten Impuls.

    @param base_metric: n×n Riemannsche Basismetrik g_{ij} auf M.
    @return: Dictionary mit 'symplectic_form' (2n×2n), 'is_exact' (bool),
             'liouville_form' (1D-Array der Länge 2n).
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n = base_metric.shape[0]   # Dimension der Basis-Mannigfaltigkeit

    # Kanonische symplektische Form auf T*M: ω = dq ∧ dp (Standard-J für n DOF)
    omega = symplectic_form(n)   # (2n × 2n)-Matrix

    # Liouville-Form θ = pᵢ dqⁱ: Als 2n-Vektor mit 1 in p-Komponenten (n bis 2n-1)
    liouville = np.zeros(2 * n)
    liouville[n:] = 1.0   # pᵢ-Koeffizienten (symbolisch = 1 für Einheitsimpulse)

    return {
        'symplectic_form': omega,
        'is_exact':        True,     # ω = dθ ist exakt (und damit geschlossen)
        'liouville_form':  liouville,
    }


def darboux_theorem_check(omega: np.ndarray, point: np.ndarray) -> dict:
    """
    @brief Numerische Überprüfung der Voraussetzungen des Darboux-Satzes.

    @description
        Darboux-Satz (Grundsatz der symplektischen Geometrie):
            Jede symplektische Mannigfaltigkeit (M, ω) ist lokal diffeomorph
            zu (ℝ^{2n}, ω_std), d.h. es existieren lokale Koordinaten (q, p),
            sodass ω = Σᵢ dqⁱ ∧ dpᵢ.

        Im Gegensatz zur Riemannschen Geometrie (wo die Krümmung ein lokales
        Invariant ist) hat die symplektische Geometrie kein lokales Invariant:
        Alle symplektischen Mannigfaltigkeiten gleicher Dimension sind lokal äquivalent!

        Numerische Überprüfung (für eine konstante symplektische Form ω auf ℝ^{2n}):
        1. Nicht-Entartetheit: det(ω) ≠ 0
        2. Geschlossenheit: dω = 0 (für konstante Matrizen trivial erfüllt)
        3. Schiefsymmetrie: ω^T = -ω

        Falls alle drei Bedingungen erfüllt sind, garantiert der Darboux-Satz,
        dass lokale Koordinaten existieren, die ω auf Normalform bringen.

    @param omega: (2n×2n)-Matrix der symplektischen Form.
    @param point: Referenzpunkt auf der Mannigfaltigkeit (für zukünftige Erweiterungen).
    @return: Dictionary mit 'is_non_degenerate', 'is_closed', 'is_antisymmetric',
             'darboux_applicable', 'det_omega'.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    rows, cols = omega.shape
    if rows != cols:
        return {
            'is_non_degenerate':  False,
            'is_closed':          False,
            'is_antisymmetric':   False,
            'darboux_applicable': False,
            'det_omega':          0.0,
        }

    # (1) Nicht-Entartetheit: det(ω) ≠ 0
    det_omega = float(np.linalg.det(omega))
    is_non_degenerate = abs(det_omega) > 1e-10

    # (2) Geschlossenheit: dω = 0
    # Für eine konstante Matrix (keine Ortsabhängigkeit) ist dω trivialerweise 0.
    is_closed = True   # Exakte 2-Form auf flachem Raum ist immer geschlossen

    # (3) Schiefsymmetrie: ω^T = -ω
    is_antisymmetric = bool(np.allclose(omega + omega.T, 0, atol=1e-10))

    # Darboux anwendbar: alle drei Bedingungen erfüllt und Dimension gerade
    darboux_applicable = (
        is_non_degenerate
        and is_closed
        and is_antisymmetric
        and rows % 2 == 0
    )

    return {
        'is_non_degenerate':  is_non_degenerate,
        'is_closed':          is_closed,
        'is_antisymmetric':   is_antisymmetric,
        'darboux_applicable': darboux_applicable,
        'det_omega':          det_omega,
    }
