"""
@file operator_algebras.py
@brief Operatoralgebren: C*-Algebren, von-Neumann-Algebren, Spektraltheorie.
@description
    Implementierung von:
    - C*-Algebren (M_n(ℂ) als kanonisches Beispiel)
    - von-Neumann-Algebren (Bikommutant-Satz, Faktoren)
    - Projektionen (Ordnungstheorie auf B(H))
    - Gelfand-Isomorphismus (kommutative C*-Algebren)
    - GNS-Konstruktion (Gelfand-Naimark-Segal)
    - Murray-von Neumann Klassifikation der Faktoren
    - K-Theorie für C*-Algebren
    - Cuntz-Algebren
    - Spurklasse-Operatoren

    Hinweis: Alle Operatoren werden durch endlichdimensionale
    Matrizen (numpy arrays) repräsentiert. Dies ermöglicht
    numerische Verifikation der algebraischen Eigenschaften.

@author Michael Fuhrmann
@date 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
from typing import List, Optional


# =============================================================================
# C*-Algebra
# =============================================================================

class CStarAlgebra:
    """
    @brief C*-Algebra, repräsentiert als Matrizenalgebra M_n(ℂ).
    @description
        Eine C*-Algebra A ist ein Banach-Algebra mit Involution * mit:
        1. (a*)* = a  (Involution ist Involution)
        2. (ab)* = b*a*  (anti-Multiplikativität)
        3. ‖a*a‖ = ‖a‖²  (C*-Identität)

        Das kanonische Beispiel ist die Algebra der n×n-Matrizen M_n(ℂ)
        mit konjugierter Transposition als Involution und Operatornorm.

    @author Michael Fuhrmann
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, elements: List[np.ndarray]):
        """
        @brief Initialisiert die C*-Algebra mit einer Liste von Erzeugern.
        @param elements Liste von quadratischen komplexen Matrizen (Erzeuger)
        @lastModified 2026-03-10
        """
        if not elements:
            raise ValueError("Mindestens ein Erzeuger erforderlich.")

        self.elements = [np.array(e, dtype=complex) for e in elements]
        # Dimensionsprüfung
        self.dim = self.elements[0].shape[0]
        for e in self.elements:
            if e.shape != (self.dim, self.dim):
                raise ValueError(
                    f"Alle Elemente müssen quadratisch mit Dimension {self.dim} sein."
                )

    def multiply(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        @brief Matrizenmultiplikation in der Algebra.
        @param a Erste Matrix
        @param b Zweite Matrix
        @return Produkt a·b
        @lastModified 2026-03-10
        """
        return a @ b

    def star(self, a: np.ndarray) -> np.ndarray:
        """
        @brief Involution (adjungierte Operation / konjugierte Transposition).
        @description
            a* = ā^T  (komplex konjugiert + transponiert).
            Für Matrizen gilt: ‖a*a‖ = ‖a‖² (C*-Identität).
        @param a Matrix in der Algebra
        @return a* (adjungierte Matrix)
        @lastModified 2026-03-10
        """
        return a.conj().T

    def norm(self, a: np.ndarray) -> float:
        """
        @brief Operatornorm ‖a‖ = σ_max(a) (größter Singulärwert).
        @description
            Für Matrizen: ‖a‖ = max { ‖av‖/‖v‖ : v ≠ 0 } = σ_max(a).
            Dies ist die Submultiplikative Norm: ‖ab‖ ≤ ‖a‖·‖b‖.
        @param a Matrix
        @return Operatornorm
        @lastModified 2026-03-10
        """
        return float(np.linalg.norm(a, ord=2))

    def is_cstar_algebra(self, sample: List[np.ndarray], tol: float = 1e-10) -> bool:
        """
        @brief Überprüft die C*-Identität ‖a*a‖ = ‖a‖² für eine Probe von Elementen.
        @param sample Liste von Testmatrizen
        @param tol    Toleranz für numerischen Vergleich
        @return True, falls C*-Identität erfüllt
        @lastModified 2026-03-10
        """
        for a in sample:
            a = np.array(a, dtype=complex)
            # C*-Identität: ‖a*a‖ = ‖a‖²
            lhs = self.norm(self.multiply(self.star(a), a))
            rhs = self.norm(a) ** 2
            if abs(lhs - rhs) > tol * max(1.0, abs(rhs)):
                return False
        return True

    def is_commutative(self, tol: float = 1e-10) -> bool:
        """
        @brief Prüft, ob alle Erzeuger kommutieren (ab = ba für alle Paare).
        @param tol Toleranz
        @return True, falls die Algebra kommutativ ist
        @lastModified 2026-03-10
        """
        for i, a in enumerate(self.elements):
            for j, b in enumerate(self.elements):
                if i >= j:
                    continue
                commutator = self.multiply(a, b) - self.multiply(b, a)
                if self.norm(commutator) > tol:
                    return False
        return True

    def spectrum(self, a: np.ndarray, tol: float = 1e-10) -> list:
        """
        @brief Berechnet das Spektrum (Eigenwerte) eines Elements.
        @description
            Spektrum: σ(a) = {λ ∈ ℂ : (a - λI) nicht invertierbar}
            = Menge der Eigenwerte (endlichdimensional).
        @param a   Matrix (Element der Algebra)
        @param tol Toleranz zum Zusammenfassen naher Eigenwerte
        @return Liste der verschiedenen Eigenwerte
        @lastModified 2026-03-10
        """
        a = np.array(a, dtype=complex)
        eigenvalues = np.linalg.eigvals(a)
        # Nahe Eigenwerte zusammenfassen (runden auf tol)
        unique_eigs = []
        for ev in eigenvalues:
            is_new = True
            for ue in unique_eigs:
                if abs(ev - ue) < tol:
                    is_new = False
                    break
            if is_new:
                unique_eigs.append(ev)
        return unique_eigs

    def spectral_radius(self, a: np.ndarray) -> float:
        """
        @brief Spektralradius r(a) = max { |λ| : λ ∈ σ(a) }.
        @description
            Für normale Elemente (a*a = aa*): r(a) = ‖a‖.
            Allgemein: r(a) = lim_{n→∞} ‖aⁿ‖^{1/n}.
        @param a Matrix
        @return Spektralradius
        @lastModified 2026-03-10
        """
        a = np.array(a, dtype=complex)
        eigs = np.linalg.eigvals(a)
        return float(np.max(np.abs(eigs)))


# =============================================================================
# von-Neumann-Algebra
# =============================================================================

class VonNeumannAlgebra(CStarAlgebra):
    """
    @brief von-Neumann-Algebra M ⊆ B(H): schwach abgeschlossene *-Unteralgebra.
    @description
        Bikommutant-Satz (von Neumann, 1929):
        M ist eine von-Neumann-Algebra genau dann, wenn M = M'' gilt,
        wobei M' = {b ∈ B(H) : ba = ab ∀a ∈ M} der Kommutant ist.

        In endlichen Dimensionen: M = Matrizenalgebra, die unter
        Adjunktion und schwachem Abschluss (= normaler Abschluss) abgeschlossen ist.

    @author Michael Fuhrmann
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def commutant(
        self,
        subset: List[np.ndarray],
        H_dim: Optional[int] = None,
        tol: float = 1e-8
    ) -> List[np.ndarray]:
        """
        @brief Berechnet den Kommutanten S' = {T ∈ B(H) : TS = ST ∀S ∈ subset}.
        @description
            Für endlichdimensionale Räume: Löse [T, S] = TS - ST = 0
            als lineares Gleichungssystem für die Einträge von T.
            Basis der Lösungen = Erzeuger des Kommutanten.

        @param subset  Liste von Matrizen (Teilmenge)
        @param H_dim   Dimension des Hilbert-Raums (Standard: dim der Elemente)
        @param tol     Toleranz für Nullraumbestimmung
        @return        Liste von Basismatrizen des Kommutanten
        @lastModified 2026-03-10
        """
        if H_dim is None:
            H_dim = self.dim

        n = H_dim
        n2 = n * n

        # [T, S] = TS - ST = 0 als lineares System für vec(T)
        # Nutze: vec(TS) = (S^T ⊗ I) vec(T), vec(ST) = (I ⊗ S) vec(T)
        # [T, S] = 0 ⟺ [(S^T ⊗ I) - (I ⊗ S)] vec(T) = 0
        A_rows = []
        for S in subset:
            S = np.array(S, dtype=complex)
            # Kommutator-Matrix für vec(T)
            C = np.kron(S.T, np.eye(n)) - np.kron(np.eye(n), S)
            A_rows.append(C)

        if not A_rows:
            # Kein Constraint: gesamtes B(H)
            basis = [
                np.eye(n, dtype=complex) * (1 if i == j else 0)
                for i in range(n) for j in range(n)
            ]
            return basis

        # Zusammenfügen und Nullraum bestimmen
        A = np.vstack(A_rows)
        _, s, Vh = np.linalg.svd(A)
        null_mask = s < tol if len(s) > 0 else np.ones(n2, dtype=bool)
        # Nullraum: Zeilen von Vh, zu denen keine große Singulärwerte gehören
        null_dim = np.sum(s < tol)
        # Vollständiger Nullraum aus den letzten Zeilen von Vh
        null_vectors = Vh[n2 - null_dim:, :] if null_dim > 0 else Vh[-1:, :]

        # Basismatrizen des Kommutanten
        basis = []
        for v in null_vectors:
            T = v.reshape(n, n)
            basis.append(T)

        return basis

    def bicommutant(
        self,
        subset: List[np.ndarray],
        tol: float = 1e-8
    ) -> List[np.ndarray]:
        """
        @brief Berechnet den Bikommutanten S'' = (S')'.
        @param subset  Liste von Matrizen
        @param tol     Toleranz
        @return        Basismatrizen des Bikommutanten
        @lastModified 2026-03-10
        """
        # Erst den Kommutanten berechnen
        commutant_basis = self.commutant(subset, tol=tol)
        # Dann den Kommutanten des Kommutanten
        return self.commutant(commutant_basis, tol=tol)

    def is_von_neumann_algebra(self, H_dim: Optional[int] = None, tol: float = 1e-6) -> bool:
        """
        @brief Prüft, ob M = M'' gilt (Bikommutant-Satz).
        @description
            Vergleicht den linearen Span der Erzeuger mit dem Bikommutanten.
            In endlichen Dimensionen: M ist von-Neumann-Algebra ⟺ M = M''.

        @param H_dim Dimension des Hilbert-Raums
        @param tol   Toleranz
        @return True, falls M = M''
        @lastModified 2026-03-10
        """
        if H_dim is None:
            H_dim = self.dim

        bicomm = self.bicommutant(self.elements, tol=tol)
        # Prüfe, ob jedes Element des Bikommutanten im Span von M liegt
        # (Grober Test: Dimension des Spans vergleichen)
        # Span von M als Vektoren
        m_vecs = np.vstack([e.flatten() for e in self.elements])
        bc_vecs = np.vstack([e.flatten() for e in bicomm]) if bicomm else np.zeros((1, H_dim**2))

        rank_m = np.linalg.matrix_rank(m_vecs, tol=tol)
        rank_bc = np.linalg.matrix_rank(bc_vecs, tol=tol)
        # Zusätzlich: Jeder Bikommutant-Basisvektor sollte im Span von M liegen
        if rank_m != rank_bc:
            return False
        return True

    def center(self, tol: float = 1e-8) -> List[np.ndarray]:
        """
        @brief Berechnet das Zentrum Z(M) = M ∩ M'.
        @description
            Das Zentrum sind alle Elemente, die mit allen anderen kommutieren:
            Z(M) = {z ∈ M : za = az ∀a ∈ M}.
        @param tol Toleranz
        @return Basismatrizen des Zentrums
        @lastModified 2026-03-10
        """
        # Zentrum = Elemente im Schnitt von M und M'
        # Vereinfachung: Kommutant der gesamten Algebra berechnen
        return self.commutant(self.elements, tol=tol)

    def is_factor(self, tol: float = 1e-6) -> bool:
        """
        @brief Prüft, ob M ein Faktor ist (Zentrum = {λI}).
        @description
            M ist ein Faktor genau dann, wenn das Zentrum trivial ist,
            d.h. Z(M) = {λI : λ ∈ ℂ}.
        @param tol Toleranz
        @return True, falls M ein Faktor ist
        @lastModified 2026-03-10
        """
        center_basis = self.center(tol=tol)
        n = self.dim
        I = np.eye(n, dtype=complex)

        # Prüfe, ob alle Zentrum-Basiselemente Vielfache der Einheitsmatrix sind
        for z in center_basis:
            # z sollte λI sein für ein λ ∈ ℂ
            # Normierter Vergleich: z / Spur(z) * n ≈ I?
            tr = np.trace(z)
            if abs(tr) < tol:
                # z ≈ 0 (kann trotzdem im Zentrum sein)
                if self.norm(z) > tol:
                    return False
            else:
                # z - (tr/n)·I sollte 0 sein
                diff = z - (tr / n) * I
                if self.norm(diff) > tol * self.norm(z):
                    return False

        return True


# =============================================================================
# Projektionen in B(H)
# =============================================================================

class Projection(CStarAlgebra):
    """
    @brief Projektionen in B(H): P = P* = P² (orthogonale Projektionen).
    @description
        Projektionen definieren eine Ordnung auf den Unterräumen von H:
        P ≤ Q ⟺ PQ = P ⟺ Bild(P) ⊆ Bild(Q).

        Partiell-isometrische Operatoren: V mit V*V = P (Projektion).

    @author Michael Fuhrmann
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def is_projection(self, P: np.ndarray, tol: float = 1e-8) -> bool:
        """
        @brief Prüft, ob P eine orthogonale Projektion ist: P = P* = P².
        @param P   Zu prüfende Matrix
        @param tol Toleranz
        @return True, falls P eine orthogonale Projektion ist
        @lastModified 2026-03-10
        """
        P = np.array(P, dtype=complex)
        # Idempotenz: P² = P
        idempotent = self.norm(self.multiply(P, P) - P) < tol
        # Selbstadjungiertheit: P* = P
        self_adjoint = self.norm(self.star(P) - P) < tol
        return idempotent and self_adjoint

    def range_space(self, P: np.ndarray, tol: float = 1e-8) -> np.ndarray:
        """
        @brief Berechnet eine Orthonormalbasis des Bildraums von P.
        @description
            Bild(P) = {Pv : v ∈ H} = Eigenraum zum Eigenwert 1.
        @param P   Projektionsmatrix
        @param tol Toleranz für die Nullraumbestimmung
        @return    Orthonormalbasis des Bildraums (Spalten der Rückgabe)
        @lastModified 2026-03-10
        """
        P = np.array(P, dtype=complex)
        # Singulärwertzerlegung: Bildraum = linkssinguläre Vektoren zu σ > tol
        U, s, _ = np.linalg.svd(P)
        # Anzahl der signifikanten Singulärwerte = Rang(P)
        rank = np.sum(s > tol)
        return U[:, :rank]

    def orthogonal_complement(self, P: np.ndarray) -> np.ndarray:
        """
        @brief Berechnet das orthogonale Komplement I - P.
        @description
            Das orthogonale Komplement von Bild(P) ist Bild(I-P) = Kern(P).
        @param P Projektionsmatrix
        @return  I - P (Projektion auf das Komplement)
        @lastModified 2026-03-10
        """
        P = np.array(P, dtype=complex)
        n = P.shape[0]
        return np.eye(n, dtype=complex) - P

    def partial_isometry(self, V: np.ndarray, tol: float = 1e-8) -> bool:
        """
        @brief Prüft, ob V eine partielle Isometrie ist: V*V ist eine Projektion.
        @description
            V ist eine partielle Isometrie ⟺ V*V ist eine Projektion (initiale Projektion).
            Äquivalent: VV*V = V.
        @param V   Zu prüfende Matrix
        @param tol Toleranz
        @return True, falls V eine partielle Isometrie ist
        @lastModified 2026-03-10
        """
        V = np.array(V, dtype=complex)
        # V*V ist eine Projektion?
        VstarV = self.multiply(self.star(V), V)
        return self.is_projection(VstarV, tol=tol)


# =============================================================================
# Gelfand-Isomorphismus
# =============================================================================

def gelfand_transform_demo() -> dict:
    """
    @brief Gelfand-Isomorphismus für kommutative C*-Algebren.
    @description
        Satz (Gelfand, 1943): Jede kommutative C*-Algebra mit Einheit ist
        isometrisch *-isomorph zu C(X), wobei X der Raum der Charaktere
        (multiplikativen Funktionale) ist.

        Beispiel: Algebra der Diagonalmatrizen D_n ⊂ M_n(ℂ).
        - D_n ist kommutativ (Diagonalmatrizen kommutieren)
        - Charaktere χ_k : D_n → ℂ, χ_k(diag(d₁,...,dₙ)) = d_k
        - Gelfand-Raum X = {χ₁,...,χₙ} ≅ {1,...,n}
        - D_n ≅ C({1,...,n}) = ℂⁿ (als C*-Algebren)

    @return dict mit Demonstrationsergebnis
    @lastModified 2026-03-10
    """
    n = 4  # Dimension

    # Diagonale C*-Algebra erzeugen
    generators = [np.diag(np.eye(n)[k]) for k in range(n)]
    algebra = CStarAlgebra(generators)

    # Prüfe Kommutativität
    is_comm = algebra.is_commutative()

    # Zufällige Diagonalmatrix
    d = np.diag(np.array([1.0 + 2j, 3.0 - 1j, -2.0 + 0j, 4.0 + 1j]))

    # Gelfand-Transformation: ĝ(χ_k) = d_k (Diagonaleintrag k)
    gelfand_values = np.diag(d)

    # C*-Norm = max Betrag der Gelfand-Werte = Operatornorm
    sup_norm = np.max(np.abs(gelfand_values))
    op_norm = algebra.norm(d)

    return {
        'description': 'Diagonalmatrizen D_n ≅ C({1,...,n})',
        'dimension': n,
        'is_commutative': is_comm,
        'gelfand_values': gelfand_values.tolist(),
        'sup_norm': float(sup_norm),
        'operator_norm': float(op_norm),
        'norms_equal': bool(abs(sup_norm - op_norm) < 1e-10),
    }


# =============================================================================
# GNS-Konstruktion
# =============================================================================

def gns_construction_demo() -> dict:
    """
    @brief GNS-Konstruktion (Gelfand-Naimark-Segal).
    @description
        Satz (GNS): Jede C*-Algebra A mit treuer positiver linearer
        Funktion φ (Zustand) trägt eine Hilbert-Raum-Darstellung π: A → B(H).

        Konstruktion:
        1. Inneres Produkt: ⟨a,b⟩ = φ(a*b)
        2. H = Vervollständigung von A bzgl. ⟨·,·⟩
        3. π(a)b = ab (Linksmultiplikation)

        Beispiel: A = M_2(ℂ), φ = normierte Spur φ(a) = Tr(a)/2.
        H = M_2(ℂ) mit Hilbert-Schmidt-Skalarprodukt ⟨A,B⟩_HS = Tr(A*B)/2.

    @return dict mit GNS-Darstellungsresultaten
    @lastModified 2026-03-10
    """
    # Für M_2(ℂ) mit normierter Spur als Zustand
    n = 2

    # Basis von M_2(ℂ): Einheitsmatrizen E_ij
    basis = []
    for i in range(n):
        for j in range(n):
            E = np.zeros((n, n), dtype=complex)
            E[i, j] = 1.0
            basis.append(E)

    # Zustand: φ(a) = Tr(a)/n (normierte Spur)
    def state(a):
        return np.trace(a) / n

    # Hilbert-Schmidt-Skalarprodukt ⟨a,b⟩ = φ(a*b) = Tr(a†b)/n
    def inner_product(a, b):
        return np.trace(a.conj().T @ b) / n

    # Gram-Matrix der Basis
    m = len(basis)
    gram = np.zeros((m, m), dtype=complex)
    for i, ei in enumerate(basis):
        for j, ej in enumerate(basis):
            gram[i, j] = inner_product(ei, ej)

    # Operatornorm im GNS-Hilbert-Raum (via Hilbert-Schmidt-Norm)
    test_op = np.array([[1, 2], [3, 4]], dtype=complex)
    hs_norm = float(np.sqrt(np.real(np.trace(test_op.conj().T @ test_op) / n)))
    op_norm = float(np.linalg.norm(test_op, ord=2))

    # Die GNS-Darstellung π: A → B(H) ist treu und isometrisch
    # Isometrie: ‖π(a)‖_B(H) = ‖a‖_A
    return {
        'description': 'GNS für M_2(ℂ) mit normierter Spur',
        'basis_size': m,
        'gram_matrix': gram.tolist(),
        'gram_matrix_rank': int(np.linalg.matrix_rank(gram)),
        'hs_norm_of_test_op': hs_norm,
        'operator_norm_of_test_op': op_norm,
        'state_of_identity': float(np.real(state(np.eye(n)))),
    }


# =============================================================================
# Faktoren-Klassifikation (Murray-von Neumann)
# =============================================================================

def factors_classification() -> dict:
    """
    @brief Murray-von Neumann Klassifikation der Faktoren.
    @description
        Faktoren sind von-Neumann-Algebren mit trivialem Zentrum.
        Klassifikation nach Art der Projektionen:

        Typ I_n:    M_n(ℂ), endlichdimensional
                    Minimalprojektion hat relative Spurwert 1/n
        Typ I_∞:   B(H) für separablen H, unendlich-dimensional
        Typ II_1:  Endliche Faktoren: ∃ normierte Spur Tr mit Tr(I)=1
        Typ II_∞:  Halbbeschränkt: Tensor II_1 ⊗ I_∞
        Typ III:   Rein unendlich: keine treue normale halbendliche Spur

        Hinweis: II_1 und III existieren nur für unendlichdimensionale H.

    @return dict mit Klassifikationsbeispielen
    @lastModified 2026-03-10
    """
    results = {}

    # --- Typ I_n: M_n(ℂ) ---
    for n in [2, 3]:
        # Erzeuger: alle Einheitsmatrizen E_ij
        generators = []
        for i in range(n):
            for j in range(n):
                E = np.zeros((n, n), dtype=complex)
                E[i, j] = 1.0
                generators.append(E)

        vna = VonNeumannAlgebra(generators)
        # Normierte Spur: Tr_n(a) = Tr(a)/n
        I_n = np.eye(n, dtype=complex)
        tr_I = np.real(np.trace(I_n)) / n  # = 1 (normiert)

        # Minimale Projektion (Rang-1-Projektion)
        P_min = np.zeros((n, n), dtype=complex)
        P_min[0, 0] = 1.0
        tr_P = np.real(np.trace(P_min)) / n  # = 1/n

        results[f'Type_I_{n}'] = {
            'description': f'M_{n}(ℂ): Faktor Typ I_{n}',
            'dimension': n,
            'is_factor': vna.is_factor(),
            'normalized_trace_of_I': float(tr_I),
            'trace_of_minimal_projection': float(tr_P),
            'dim_of_generators': len(generators),
        }

    # --- Typ I_1 (Skalare): ℂ = M_1(ℂ) ---
    results['Type_I_1'] = {
        'description': 'ℂ = M_1(ℂ): Faktor Typ I_1 (trivial)',
        'dimension': 1,
        'normalized_trace_of_I': 1.0,
    }

    # --- Typ II_1 (konzeptuell, endlichdimensionale Approximation) ---
    # Hyperendlicher Faktor R ≈ lim_{n→∞} M_2^{⊗n}(ℂ)
    # Approximation: M_4(ℂ) = M_2 ⊗ M_2 als II_1-ähnlich
    n4 = 4
    results['Type_II_1_approx'] = {
        'description': 'M_4(ℂ) ≅ M_2 ⊗ M_2: Approximation des hyperendlichen II_1-Faktors R',
        'dimension': n4,
        'note': (
            'Echter II_1-Faktor ist unendlichdimensional. '
            'Approximation: M_4(ℂ) mit normierter Spur Tr(A)/4.'
        ),
        'normalized_trace_of_I': 1.0,
    }

    # --- Typ III (konzeptuell) ---
    results['Type_III'] = {
        'description': 'Typ III Faktoren: keine halbendliche Spur, rein unendlich',
        'examples': [
            'Powers-Faktoren R_λ (0 < λ < 1): unendliche Tensorprodukte',
            'Araki-Woods Faktoren',
            'Connes-Klassifikation (1973/76): III_0, III_λ (0<λ<1), III_1',
        ],
        'note': (
            'Nicht durch endlichdimensionale Matrizen darstellbar. '
            'Erfordern Modulartheorie (Tomita-Takesaki).'
        ),
    }

    return results


# =============================================================================
# K-Theorie für C*-Algebren
# =============================================================================

def k_theory_intro() -> dict:
    """
    @brief Einführung in die K-Theorie für C*-Algebren.
    @description
        K-Theorie ordnet C*-Algebren abelsche Gruppen zu:

        K₀(A): Grothendieck-Gruppe der stabilen Äquivalenzklassen von Projektionen.
               - Für A = M_n(ℂ): K₀(M_n(ℂ)) ≅ ℤ (erzeugt von [I_n/n])
               - Für A = C(X): K₀(C(X)) ≅ K⁰(X) (topologische K-Theorie)

        K₁(A): Unitäre Äquivalenzklassen in M_∞(A)⁺.
               - Für A = M_n(ℂ): K₁(M_n(ℂ)) = 0
               - Für A = C(S¹): K₁(C(S¹)) ≅ ℤ (Windungszahl)

        Bott-Periodizität: Kₙ₊₂(A) ≅ Kₙ(A).

    @return dict mit K-Theorie-Informationen
    @lastModified 2026-03-10
    """
    # Projektionen in M_2(ℂ) klassifizieren
    n = 2

    # Alle Rang-1-Projektionen P in M_2(ℂ): P = P* = P², Rang(P)=1
    # Sie sind alle unitär äquivalent: P ~ Q falls UP U* = Q für unitäres U
    # → [P] = 1 in K₀(M_2(ℂ)) ≅ ℤ

    P1 = np.array([[1, 0], [0, 0]], dtype=complex)  # Standard-Rang-1-Projektion
    P2 = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)  # andere Rang-1-Projektion

    # Spur als K₀-Klassifikation: Tr(P) = Rang(P) = K₀-Klasse
    k0_P1 = int(round(np.real(np.trace(P1))))
    k0_P2 = int(round(np.real(np.trace(P2))))

    # Unitäre Äquivalenz: P1 ~ P2 (da gleicher Rang)
    # Finde U mit U P1 U* = P2
    # U = Drehung im 2D: U = [[cos θ, -sin θ], [sin θ, cos θ]] mit θ=π/4
    theta = np.pi / 4
    U = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ], dtype=complex)
    P1_rotated = U @ P1 @ U.conj().T
    unitary_equiv_error = float(np.linalg.norm(P1_rotated - P2))

    # Bott-Periodizität: K₀ ≅ K₂ (nur konzeptuell)
    # Für A = ℂ: K₀(ℂ) ≅ ℤ, K₁(ℂ) = 0, K₂(ℂ) ≅ ℤ,...

    return {
        'K0_description': 'K₀(A) = Grothendieck-Gruppe der Projektionen',
        'K1_description': 'K₁(A) = Unitäre Äquivalenzklassen (stabilisiert)',
        'bott_periodicity': 'Kₙ₊₂(A) ≅ Kₙ(A) (Bott-Periodizität)',
        'examples': {
            'K0_Mn_C': 'K₀(M_n(ℂ)) ≅ ℤ, erzeugt von [Rang-1-Projektion]',
            'K1_Mn_C': 'K₁(M_n(ℂ)) = 0',
            'K0_C(S1)': 'K₀(C(S¹)) ≅ ℤ',
            'K1_C(S1)': 'K₁(C(S¹)) ≅ ℤ (Windungszahl)',
        },
        'P1_K0_class': k0_P1,
        'P2_K0_class': k0_P2,
        'unitary_equivalence_error': unitary_equiv_error,
        'P1_P2_unitarily_equivalent': bool(unitary_equiv_error < 1e-10),
        'index_map': 'Exponentialfolge: 0 → K₁(I) → K₁(A) → K₁(A/I) → K₀(I) → K₀(A) → K₀(A/I) → 0',
    }


# =============================================================================
# Cuntz-Algebren
# =============================================================================

def cuntz_algebra_demo() -> dict:
    """
    @brief Demonstration der Cuntz-Algebra O_n (endlichdimensionale Approximation).
    @description
        Cuntz-Algebren O_n (Cuntz, 1977):
        Erzeugt von n Isometrien S₁,...,Sₙ mit:
        - Sᵢ*Sᵢ = I  (Sᵢ ist eine Isometrie)
        - ΣᵢSᵢSᵢ* = I  (Vollständigkeitsrelation)
        - Sᵢ*Sⱼ = δᵢⱼ I  (Isometrien mit orthogonalen Bildern)

        O_n ist einfach und rein unendlich (Typ III₁).
        O₂ ist die einfachste echte Cuntz-Algebra.

        Endlichdimensionale Approximation: Isometrien in M_n²(ℂ)
        durch Permutationsmatrizen (Toy-Modell).

    @return dict mit Cuntz-Algebra-Eigenschaften
    @lastModified 2026-03-10
    """
    # n=2: O₂ in einer endlichdimensionalen Approximation
    # Wir bauen ein Toy-Modell in M_4(ℂ):
    # S₁, S₂: Isometrien von ℂ² nach ℂ⁴, sodass S₁S₁* + S₂S₂* = I₄

    # S₁: ℂ² → ℂ⁴, x ↦ (x₁, x₂, 0, 0)
    S1 = np.zeros((4, 2), dtype=complex)
    S1[0, 0] = 1.0
    S1[1, 1] = 1.0

    # S₂: ℂ² → ℂ⁴, x ↦ (0, 0, x₁, x₂)
    S2 = np.zeros((4, 2), dtype=complex)
    S2[2, 0] = 1.0
    S2[3, 1] = 1.0

    # Prüfe Isometrie-Bedingungen
    # Sᵢ*Sᵢ = I₂
    s1_star_s1 = S1.conj().T @ S1
    s2_star_s2 = S2.conj().T @ S2
    iso_error_1 = float(np.linalg.norm(s1_star_s1 - np.eye(2)))
    iso_error_2 = float(np.linalg.norm(s2_star_s2 - np.eye(2)))

    # Vollständigkeitsrelation: S₁S₁* + S₂S₂* = I₄
    completeness = S1 @ S1.conj().T + S2 @ S2.conj().T
    completeness_error = float(np.linalg.norm(completeness - np.eye(4)))

    # Orthogonalität der Bilder: S₁*S₂ = 0
    s1_star_s2 = S1.conj().T @ S2
    orthogonality_error = float(np.linalg.norm(s1_star_s2))

    # Cuntz-Relationen erfüllt ≈ 0?
    return {
        'description': 'Toy-Modell von O₂ in M_4(ℂ): S₁, S₂: ℂ²→ℂ⁴',
        'n': 2,
        'S1_shape': list(S1.shape),
        'S2_shape': list(S2.shape),
        'S1_star_S1_error': iso_error_1,
        'S2_star_S2_error': iso_error_2,
        'completeness_error': completeness_error,
        'orthogonality_error': orthogonality_error,
        'cuntz_relations_satisfied': bool(
            iso_error_1 < 1e-10 and iso_error_2 < 1e-10
            and completeness_error < 1e-10
            and orthogonality_error < 1e-10
        ),
        'property_simple': 'O_n ist einfach (kein nichttriviales abgeschlossenes Ideal)',
        'property_purely_infinite': 'O_n ist rein unendlich (Typ III₁-Faktor)',
        'K_theory': 'K₀(O_n) ≅ ℤ/(n-1)ℤ, K₁(O_n) = 0',
    }


# =============================================================================
# Spurklasse-Operatoren
# =============================================================================

def trace_class_operators(H_dim: int = 4) -> dict:
    """
    @brief Spurklasse-Operatoren S₁(H) ⊂ B(H).
    @description
        Ein beschränkter Operator T ∈ B(H) ist in der Spurklasse S₁(H), falls:
        Tr(|T|) = Tr(√(T*T)) < ∞

        Eigenschaften:
        - S₁(H) ⊂ kompakte Operatoren ⊂ B(H)
        - Spur: Tr(T) = Σ⟨Teₙ, eₙ⟩ (absolut konvergent)
        - Zyklizität: Tr(AB) = Tr(BA)
        - Norm: ‖T‖₁ = Tr(|T|) = Σ σᵢ(T) (Summe der Singulärwerte)
        - Dualität: S₁(H)* ≅ B(H) (als Banach-Raum)

        Schatten-p-Klassen: Sₚ = {T : Σ σᵢᵖ < ∞}
        - S₁: Spurklasse, S₂: Hilbert-Schmidt-Klasse, S∞=K(H): kompakte Operatoren

    @param H_dim Dimension des Hilbert-Raums
    @return dict mit Spurklassen-Eigenschaften
    @lastModified 2026-03-10
    """
    n = H_dim

    # Zufällige positiv semidefinite Matrix als Spurklassen-Operator
    rng = np.random.default_rng(42)
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    # T = A†A ist positiv semidefinit (automatisch in S₁ für endliches dim)
    T = A.conj().T @ A / n

    # Singulärwerte
    sigma = np.linalg.svd(T, compute_uv=False)

    # Spur-Norm ‖T‖₁ = Σ σᵢ
    trace_norm = float(np.sum(sigma))

    # Spur Tr(T)
    trace_T = float(np.real(np.trace(T)))

    # Für positive T: Tr(T) = ‖T‖₁ (Spur = Spur-Norm)
    trace_equals_trace_norm = bool(abs(trace_T - trace_norm) < 1e-8)

    # Zyklizität: Tr(AB) = Tr(BA)
    B = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    B = B / np.linalg.norm(B)
    tr_AB = float(np.real(np.trace(T @ B)))
    tr_BA = float(np.real(np.trace(B @ T)))
    cyclicity_error = abs(tr_AB - tr_BA)

    # Hilbert-Schmidt-Norm ‖T‖₂ = √(Σ σᵢ²) = √Tr(T†T)
    hs_norm = float(np.sqrt(np.real(np.trace(T.conj().T @ T))))

    # Operatornorm ‖T‖ = σ_max
    op_norm = float(sigma[0])

    # Ungleichungen: ‖T‖ ≤ ‖T‖₂ ≤ ‖T‖₁
    inequalities_hold = bool(op_norm <= hs_norm + 1e-10 and hs_norm <= trace_norm + 1e-10)

    return {
        'description': f'Spurklasse-Operator in M_{n}(ℂ)',
        'singular_values': sigma.tolist(),
        'trace_norm': trace_norm,
        'trace': trace_T,
        'trace_equals_trace_norm_for_positive_T': trace_equals_trace_norm,
        'tr_AB': tr_AB,
        'tr_BA': tr_BA,
        'cyclicity_error': float(cyclicity_error),
        'cyclicity_holds': bool(cyclicity_error < 1e-8),
        'hilbert_schmidt_norm': hs_norm,
        'operator_norm': op_norm,
        'norm_inequalities_hold': inequalities_hold,
        'schatten_classes': {
            'S1_trace_class': '‖T‖₁ = Σσᵢ < ∞',
            'S2_hilbert_schmidt': '‖T‖₂ = (Σσᵢ²)^{1/2} < ∞',
            'S_inf_compact': '‖T‖ = σ_max',
        },
    }
