"""
@file homological_algebra.py
@brief Homologische Algebra: Kettenkomplexe, Homologie, Ext/Tor.
@description
    Dieses Modul implementiert grundlegende Konzepte der homologischen Algebra:

    - ChainComplex: Kettenkomplexe mit Randoperatoren ∂ₙ: Cₙ → Cₙ₋₁
    - CochainComplex: Duale Kokettenkomplexe (Kohomologie)
    - Homologie H_n = ker(∂_n) / im(∂_{n+1}) via Smith-Normalform
    - Simpliziale Komplexe → Kettenkomplex
    - Exakte Sequenzen und das Schlangen-Lemma
    - Freie Auflösungen (von ℤ/nℤ)
    - Ext und Tor Gruppen für zyklische Moduln
    - Universeller Koeffizienten-Satz

    Das fundamentale Axiom: ∂∘∂ = 0
    (Das Bild von ∂_{n+1} liegt stets im Kern von ∂_n)

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import numpy as np
import math
from typing import Optional


# ---------------------------------------------------------------------------
# Hilfsfunktionen: Smith-Normalform über ℤ
# ---------------------------------------------------------------------------

def _smith_normal_form(M: list) -> tuple:
    """
    Berechnet die Smith-Normalform einer ganzzahligen Matrix.

    Jede ganzzahlige Matrix M kann geschrieben werden als:
        M = P · D · Q
    wobei P, Q invertierbare ganzzahlige Matrizen sind und
    D eine Diagonalmatrix mit d₁ | d₂ | ... | dᵣ (Teilerfolge).

    Die Diagonaleinträge dᵢ sind die Elementarteiler.

    @param M   Ganzzahlige Matrix (verschachtelte Liste oder 2D-Array)
    @return    Tupel (D, rank, elementary_divisors) mit:
               D = Diagonalmatrix (als numpy-Array)
               rank = Rang der Matrix
               elementary_divisors = Liste der nicht-null Diagonaleinträge
    @lastModified 2026-03-10
    """
    if len(M) == 0:
        return np.array([[]]), 0, []

    # Arbeite mit einer ganzzahligen Kopie
    A = np.array(M, dtype=int)
    m, n = A.shape

    row = 0
    col = 0
    # Elementarteiler sammeln (werden am Ende ausgelesen)

    while row < m and col < n:
        # Suche das kleinste Nicht-Null-Element in der Untermatrix A[row:,col:]
        # Wähle als Pivotelement das Element mit dem kleinsten Absolutwert
        submatrix = A[row:, col:]
        nonzero_mask = submatrix != 0

        if not nonzero_mask.any():
            col += 1
            continue

        # Finde Position des kleinsten NichtNull-Elements
        nonzero_values = np.abs(submatrix[nonzero_mask])
        min_val = nonzero_values.min()
        positions = np.argwhere(nonzero_mask & (np.abs(submatrix) == min_val))
        ri, ci = positions[0]
        ri += row
        ci += col

        # Bringe Pivotelement nach (row, col) durch Zeilen-/Spaltenvertauschung
        A[[row, ri]] = A[[ri, row]]
        A[:, [col, ci]] = A[:, [ci, col]]

        # Mache Pivotelement positiv
        if A[row, col] < 0:
            A[row] = -A[row]

        # Eliminiere alle anderen Einträge in Zeile und Spalte via euklidischem Algorithmus
        changed = True
        while changed:
            changed = False

            # Eliminiere Einträge in Spalte col (außer Zeile row)
            for i in range(row + 1, m):
                if A[i, col] != 0:
                    # Euklidische Division: A[i,col] = q * A[row,col] + r
                    q = A[i, col] // A[row, col]
                    A[i] -= q * A[row]
                    if A[i, col] != 0:
                        # Noch nicht fertig: tausche und wiederhole
                        A[[row, i]] = A[[i, row]]
                        if A[row, col] < 0:
                            A[row] = -A[row]
                        changed = True

            # Eliminiere Einträge in Zeile row (außer Spalte col)
            for j in range(col + 1, n):
                if A[row, j] != 0:
                    q = A[row, j] // A[row, col]
                    A[:, j] -= q * A[:, col]
                    if A[row, j] != 0:
                        A[:, [col, j]] = A[:, [j, col]]
                        if A[row, col] < 0:
                            A[row] = -A[row]
                        changed = True

        # Prüfe ob Pivot alle anderen Einträge in seinem Block teilt
        # (für echte Smith-Normalform nötig)
        all_divisible = True
        for i in range(row + 1, m):
            for j in range(col + 1, n):
                if A[i, j] % A[row, col] != 0:
                    all_divisible = False
                    break

        if not all_divisible:
            # Addiere problematische Zeile zu Pivotzeile und starte erneut
            for i in range(row + 1, m):
                for j in range(col + 1, n):
                    if A[i, j] % A[row, col] != 0:
                        A[row] += A[i]
                        changed = True
                        break
                if changed:
                    break
            continue

        row += 1
        col += 1

    # Extrahiere Diagonaleinträge
    diag_len = min(m, n)
    diag = [A[i, i] for i in range(min(row, diag_len))]

    # Berechne Rang (Anzahl der Nicht-Null-Diagonaleinträge)
    rank = sum(1 for d in diag if d != 0)
    elementary_divisors = [abs(d) for d in diag if d != 0]

    return A, rank, elementary_divisors


def _gcd(a: int, b: int) -> int:
    """
    Größter gemeinsamer Teiler via Euklidischem Algorithmus.

    @param a   Erste ganze Zahl
    @param b   Zweite ganze Zahl
    @return    gcd(a, b) ≥ 0
    @lastModified 2026-03-10
    """
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


# ---------------------------------------------------------------------------
# Kettenkomplex
# ---------------------------------------------------------------------------

class ChainComplex:
    """
    Kettenkomplex C_*: ... → C_{n+1} →^{∂_{n+1}} C_n →^{∂_n} C_{n-1} → ...

    Das fundamentale Axiom lautet: ∂∘∂ = 0
    d.h. jedes Bild von ∂_{n+1} liegt im Kern von ∂_n.

    Implementierung:
    - groups: {n: rank} — die Gruppe Cₙ = ℤ^rank
    - boundaries: {n: matrix} — die Randmatrix ∂_n: Cₙ → Cₙ₋₁
    """

    def __init__(self, groups: dict, boundaries: dict):
        """
        Initialisiert den Kettenkomplex.

        @param groups      Dictionary {n: rank} mit ℤ^rank in Grad n
        @param boundaries  Dictionary {n: matrix} mit Randoperator ∂_n: Cₙ → Cₙ₋₁
                           Die Matrix hat Dimension rank(Cₙ₋₁) × rank(Cₙ)
        @lastModified 2026-03-10
        """
        # Speichere Gruppen und Randoperatoren
        self.groups = {int(n): int(r) for n, r in groups.items()}
        self.boundaries = {}
        for n, mat in boundaries.items():
            self.boundaries[int(n)] = np.array(mat, dtype=int)

    def is_complex(self) -> bool:
        """
        Prüft das Kettenkomplex-Axiom: ∂_{n-1} ∘ ∂_n = 0 für alle n.

        D.h. die Komposition zweier aufeinanderfolgender Randoperatoren ist die Nullmatrix.

        @return   True wenn ∂∘∂ = 0 für alle aufeinanderfolgenden Paare
        @lastModified 2026-03-10
        """
        for n in sorted(self.boundaries.keys()):
            # Prüfe ob auch ∂_{n-1} vorhanden ist
            if (n - 1) in self.boundaries:
                d_n = self.boundaries[n]       # ∂_n: C_n → C_{n-1}
                d_n_minus_1 = self.boundaries[n - 1]  # ∂_{n-1}: C_{n-1} → C_{n-2}

                # Berechne die Komposition ∂_{n-1} ∘ ∂_n
                composition = d_n_minus_1 @ d_n

                # Prüfe ob das Ergebnis die Nullmatrix ist
                if not np.all(composition == 0):
                    return False
        return True

    def homology(self, n: int) -> dict:
        """
        Berechnet die n-te Homologiegruppe H_n = ker(∂_n) / im(∂_{n+1}).

        Algorithmus:
        1. Berechne ker(∂_n) via Nullraum der Matrix ∂_n
        2. Berechne im(∂_{n+1}) via Spaltenraum der Matrix ∂_{n+1}
        3. Teile durch: H_n ≅ ℤ^β ⊕ ℤ/d₁ℤ ⊕ ... ⊕ ℤ/dₖℤ (Struktursatz)
           Hier ist β die Betti-Zahl und dᵢ die Torsionskoeffizienten.

        @param n   Grad der Homologiegruppe
        @return    Dictionary {'rank': β, 'torsion': [d₁,...,dₖ], 'free': bool}
        @lastModified 2026-03-10
        """
        # Rang der Gruppe C_n
        rank_n = self.groups.get(n, 0)

        if rank_n == 0:
            return {'rank': 0, 'torsion': [], 'free': True}

        # Randoperator ∂_n: C_n → C_{n-1}
        d_n = self.boundaries.get(n, None)

        # Randoperator ∂_{n+1}: C_{n+1} → C_n
        d_np1 = self.boundaries.get(n + 1, None)

        # Berechne Kern von ∂_n (Nullraum)
        if d_n is not None and d_n.size > 0:
            _, rank_dn, _ = _smith_normal_form(d_n.tolist())
            # dim(ker ∂_n) = rank(C_n) - rank(∂_n)
            dim_ker = rank_n - rank_dn
        else:
            # ∂_n = 0 → Kern = C_n vollständig
            dim_ker = rank_n
            rank_dn = 0

        # Berechne Bild von ∂_{n+1} (Spaltenraum in C_n)
        if d_np1 is not None and d_np1.size > 0:
            _, rank_dnp1, _ = _smith_normal_form(d_np1.tolist())
            # Rang des Bildes = Rang der Matrix
            dim_im = rank_dnp1
        else:
            dim_im = 0
            rank_dnp1 = 0

        # Freie Betti-Zahl: β_n = dim(ker ∂_n) - dim(im ∂_{n+1})
        betti = max(0, dim_ker - dim_im)

        # Torsion: aus Smith-Normalform von ∂_{n+1}
        torsion = []
        if d_np1 is not None and d_np1.size > 0:
            _, _, elem_div = _smith_normal_form(d_np1.tolist())
            # Torsionskoeffizienten sind Elementarteiler > 1
            torsion = [d for d in elem_div if d > 1]

        return {
            'rank': betti,
            'torsion': torsion,
            'free': len(torsion) == 0
        }

    def euler_characteristic(self) -> int:
        """
        Berechnet die Euler-Charakteristik des Kettenkomplexes.

        χ = Σₙ (-1)^n · rank(Cₙ) = Σₙ (-1)^n · β_n

        Beide Formeln liefern dasselbe Ergebnis (Satz von Euler).

        @return   Euler-Charakteristik (ganze Zahl)
        @lastModified 2026-03-10
        """
        chi = 0
        for n, rank in self.groups.items():
            # Alternierend summieren: + für gerades n, - für ungerades n
            chi += ((-1) ** n) * rank
        return chi

    def betti_numbers(self) -> dict:
        """
        Berechnet alle Betti-Zahlen β_n = rank(H_n).

        Die Betti-Zahlen messen den "freien" Teil der Homologie:
        β_n = Anzahl der unabhängigen n-dimensionalen Löcher.

        @return   Dictionary {n: β_n} für alle Grade n
        @lastModified 2026-03-10
        """
        result = {}
        for n in self.groups.keys():
            h = self.homology(n)
            result[n] = h['rank']
        return result


# ---------------------------------------------------------------------------
# Kokettenkomplex
# ---------------------------------------------------------------------------

class CochainComplex:
    """
    Kokettenkomplex C^*: ... → C^{n-1} →^{d^{n-1}} C^n →^{d^n} C^{n+1} → ...

    Dual zu ChainComplex via Hom(·, ℤ).
    Kobrandoperator: d^n = (∂_{n+1})ᵀ (Transponierte des Randoperators).

    Das Axiom lautet dual: d∘d = 0.
    """

    def __init__(self, groups: dict, coboundaries: dict):
        """
        Initialisiert den Kokettenkomplex.

        @param groups        Dictionary {n: rank} mit ℤ^rank in Kograd n
        @param coboundaries  Dictionary {n: matrix} mit Kobrandoperator d^n: C^n → C^{n+1}
        @lastModified 2026-03-10
        """
        self.groups = {int(n): int(r) for n, r in groups.items()}
        self.coboundaries = {}
        for n, mat in coboundaries.items():
            self.coboundaries[int(n)] = np.array(mat, dtype=int)

    def is_complex(self) -> bool:
        """
        Prüft das Kokettenkomplex-Axiom: d^{n+1} ∘ d^n = 0 für alle n.

        @return   True wenn d∘d = 0 für alle aufeinanderfolgenden Paare
        @lastModified 2026-03-10
        """
        for n in sorted(self.coboundaries.keys()):
            if (n + 1) in self.coboundaries:
                d_n = self.coboundaries[n]       # d^n: C^n → C^{n+1}
                d_np1 = self.coboundaries[n + 1]  # d^{n+1}: C^{n+1} → C^{n+2}
                composition = d_np1 @ d_n
                if not np.all(composition == 0):
                    return False
        return True

    def cohomology(self, n: int) -> dict:
        """
        Berechnet die n-te Kohomologiegruppe H^n = ker(d^n) / im(d^{n-1}).

        @param n   Kograd der Kohomologiegruppe
        @return    Dictionary {'rank': β, 'torsion': [d₁,...], 'free': bool}
        @lastModified 2026-03-10
        """
        rank_n = self.groups.get(n, 0)
        if rank_n == 0:
            return {'rank': 0, 'torsion': [], 'free': True}

        # Kobrandoperator d^n: C^n → C^{n+1}
        d_n = self.coboundaries.get(n, None)

        # Kobrandoperator d^{n-1}: C^{n-1} → C^n
        d_nm1 = self.coboundaries.get(n - 1, None)

        # Kern von d^n
        if d_n is not None and d_n.size > 0:
            _, rank_dn, _ = _smith_normal_form(d_n.tolist())
            dim_ker = rank_n - rank_dn
        else:
            dim_ker = rank_n

        # Bild von d^{n-1}
        if d_nm1 is not None and d_nm1.size > 0:
            _, rank_dnm1, _ = _smith_normal_form(d_nm1.tolist())
            dim_im = rank_dnm1
        else:
            dim_im = 0

        betti = max(0, dim_ker - dim_im)

        # Torsion aus dem Bild
        torsion = []
        if d_nm1 is not None and d_nm1.size > 0:
            _, _, elem_div = _smith_normal_form(d_nm1.tolist())
            torsion = [d for d in elem_div if d > 1]

        return {
            'rank': betti,
            'torsion': torsion,
            'free': len(torsion) == 0
        }


# ---------------------------------------------------------------------------
# Simplizialer Komplex → Kettenkomplex
# ---------------------------------------------------------------------------

def simplicial_complex_to_chain(simplices: dict) -> ChainComplex:
    """
    Erzeugt einen Kettenkomplex aus simplizialen Daten.

    simplices: {k: [[v₀,...,vₖ], ...]} — Liste der k-Simplizes (als sortierte Tuples)

    Der Randoperator ∂ₖ bildet jeden k-Simplex auf seine (k-1)-dimensionalen Seiten ab:
        ∂[v₀, v₁, ..., vₖ] = Σᵢ (-1)^i [v₀, ..., v̂ᵢ, ..., vₖ]

    Die Matrix von ∂ₖ hat:
    - Zeilen: (k-1)-Simplizes
    - Spalten: k-Simplizes
    - Einträge: ±1 (Vorzeichen der Inzidenz)

    @param simplices   Dictionary mit k → Liste der k-Simplizes
    @return            ChainComplex mit berechneten Randoperatoren
    @lastModified 2026-03-10
    """
    # Normalisiere alle Simplizes (sortierte Tupel)
    normalized = {}
    for k, slist in simplices.items():
        normalized[int(k)] = [tuple(sorted(s)) for s in slist]

    # Erzeuge Gruppen: Cₖ = ℤ^{Anzahl der k-Simplizes}
    groups = {k: len(slist) for k, slist in normalized.items()}

    # Erzeuge Randoperatoren
    boundaries = {}
    for k in sorted(normalized.keys()):
        if k == 0:
            continue  # Kein Randoperator für 0-Simplizes

        k_simplices = normalized[k]            # Liste der k-Simplizes
        km1_simplices = normalized.get(k - 1, [])  # Liste der (k-1)-Simplizes

        if len(km1_simplices) == 0:
            continue

        # Erstelle Indexierung der (k-1)-Simplizes
        km1_index = {s: i for i, s in enumerate(km1_simplices)}

        # Randmatrix: Zeilen = (k-1)-Simplizes, Spalten = k-Simplizes
        n_rows = len(km1_simplices)
        n_cols = len(k_simplices)
        boundary_matrix = np.zeros((n_rows, n_cols), dtype=int)

        for col_j, sigma in enumerate(k_simplices):
            # Berechne den Rand von sigma
            for i in range(len(sigma)):
                # Entferne den i-ten Knoten: Seite mit Vorzeichen (-1)^i
                face = tuple(sigma[:i] + sigma[i+1:])
                sign = (-1) ** i

                # Suche die Seite in den (k-1)-Simplizes
                if face in km1_index:
                    row_i = km1_index[face]
                    boundary_matrix[row_i, col_j] = sign

        boundaries[k] = boundary_matrix

    return ChainComplex(groups, boundaries)


# ---------------------------------------------------------------------------
# Exakte Sequenzen
# ---------------------------------------------------------------------------

class ExactSequence:
    """
    Exakte Sequenz von abelschen Gruppen:
        ... → Aₙ →^{fₙ} Aₙ₋₁ → ...

    Exaktheit bedeutet: im(fₙ) = ker(fₙ₋₁) an jedem Punkt.

    Gruppen sind als ℤ^rank kodiert, Abbildungen als ganzzahlige Matrizen.
    """

    def __init__(self, groups: list, maps: list):
        """
        Initialisiert die exakte Sequenz.

        @param groups   Liste der Gruppenränge [rank(A₀), rank(A₁), ...]
        @param maps     Liste der Abbildungsmatrizen [f₁, f₂, ...]
                        (maps[i]: groups[i+1] → groups[i])
        @lastModified 2026-03-10
        """
        self.groups = [int(r) for r in groups]
        self.maps = [np.array(m, dtype=int) for m in maps]

    def is_exact(self) -> bool:
        """
        Prüft die Exaktheit der Sequenz an jeder Stelle.

        Exaktheit bei Aₙ: im(fₙ₊₁) = ker(fₙ)

        @return   True wenn exakt an allen inneren Stellen
        @lastModified 2026-03-10
        """
        # Wir prüfen an jeder inneren Stelle
        for i in range(1, len(self.maps)):
            f_next = self.maps[i]    # fₙ₊₁: Aₙ₊₁ → Aₙ
            f_curr = self.maps[i-1]  # fₙ: Aₙ → Aₙ₋₁

            # Prüfe: im(f_next) ⊆ ker(f_curr)
            # d.h. f_curr ∘ f_next = 0
            composition = f_curr @ f_next
            if not np.all(composition == 0):
                return False

            # Prüfe: Rang der Komposition entspricht Exaktheit
            _, rank_next, _ = _smith_normal_form(f_next.tolist())
            _, rank_curr, _ = _smith_normal_form(f_curr.tolist())

            n = self.groups[i]
            # dim(ker f_curr) = n - rank(f_curr)
            # dim(im f_next) = rank(f_next)
            dim_ker = n - rank_curr
            dim_im = rank_next

            if dim_ker != dim_im:
                return False

        return True

    def is_short_exact(self) -> bool:
        """
        Prüft ob die Sequenz eine kurze exakte Sequenz ist: 0 → A → B → C → 0.

        Bedingungen:
        - Genau 3 Gruppen und 2 Abbildungen
        - f₁ ist injektiv (Kern = 0)
        - f₂ ist surjektiv (Bild = C)
        - Exaktheit in der Mitte: im(f₁) = ker(f₂)

        @return   True wenn kurze exakte Sequenz
        @lastModified 2026-03-10
        """
        if len(self.groups) != 3 or len(self.maps) != 2:
            return False

        f1 = self.maps[0]  # f₁: A → B
        f2 = self.maps[1]  # f₂: B → C

        # f₁ injektiv: Rang(f₁) = Rang(A)
        _, rank_f1, _ = _smith_normal_form(f1.tolist())
        if rank_f1 != self.groups[0]:
            return False

        # f₂ surjektiv: Rang(f₂) = Rang(C)
        _, rank_f2, _ = _smith_normal_form(f2.tolist())
        if rank_f2 != self.groups[2]:
            return False

        # Exaktheit in der Mitte: im(f₁) = ker(f₂)
        composition = f2 @ f1
        if not np.all(composition == 0):
            return False

        dim_ker_f2 = self.groups[1] - rank_f2
        dim_im_f1 = rank_f1
        return dim_ker_f2 == dim_im_f1

    def split(self) -> bool:
        """
        Prüft ob die kurze exakte Sequenz 0 → A → B → C → 0 zerfällt: B ≅ A ⊕ C.

        Eine kurze exakte Sequenz zerfällt genau dann, wenn:
        - B ≅ A ⊕ C (als abelsche Gruppen)
        - Äquivalent: rank(B) = rank(A) + rank(C) ohne Torsion

        Für freie abelsche Gruppen zerfällt jede kurze exakte Sequenz.

        @return   True wenn die Sequenz zerfällt
        @lastModified 2026-03-10
        """
        if len(self.groups) != 3:
            return False

        rank_A = self.groups[0]
        rank_B = self.groups[1]
        rank_C = self.groups[2]

        # Einfaches Kriterium für freie Module: B ≅ A ⊕ C ⟺ rank(B) = rank(A) + rank(C)
        return rank_B == rank_A + rank_C


# ---------------------------------------------------------------------------
# 5-Lemma und Schlangen-Lemma
# ---------------------------------------------------------------------------

def five_lemma_check(diagram: dict) -> dict:
    """
    Überprüft das 5-Lemma für ein kommutatives Diagramm mit zwei Zeilen.

    Das 5-Lemma besagt: In einem kommutativen Diagramm
        A₁ → A₂ → A₃ → A₄ → A₅
        |    |    |    |    |
        B₁ → B₂ → B₃ → B₄ → B₅

    wenn f₁, f₂, f₄, f₅ Isomorphismen sind, dann ist auch f₃ ein Isomorphismus.

    @param diagram   Dictionary mit Struktur:
                     'top_maps': [M₁₂, M₂₃, M₃₄, M₄₅] (obere Zeile)
                     'bot_maps': [N₁₂, N₂₃, N₃₄, N₄₅] (untere Zeile)
                     'vert_maps': [f₁, f₂, f₃, f₄, f₅] (vertikale Abbildungen)
    @return          Dictionary mit Ergebnis und Analyse
    @lastModified 2026-03-10
    """
    top_maps = [np.array(m, dtype=int) for m in diagram.get('top_maps', [])]
    bot_maps = [np.array(m, dtype=int) for m in diagram.get('bot_maps', [])]
    vert_maps = [np.array(m, dtype=int) for m in diagram.get('vert_maps', [])]

    if len(vert_maps) != 5:
        return {'valid': False, 'reason': '5 vertikale Abbildungen erwartet'}

    # Prüfe ob die 4 äußeren Abbildungen Isomorphismen sind
    iso_status = []
    for i, f in enumerate(vert_maps):
        if f.shape[0] != f.shape[1]:
            iso_status.append(False)
            continue
        det = int(round(np.linalg.det(f.astype(float))))
        # Isomorphismus ⟺ |det| = 1 (für ganzzahlige Matrizen)
        iso_status.append(abs(det) == 1)

    # Prüfe Kommutatvität des Diagramms
    commutes = True
    for i in range(min(len(top_maps), len(bot_maps))):
        if i < len(vert_maps) - 1:
            # f_{i+1} ∘ top_maps[i] = bot_maps[i] ∘ f_i
            lhs = vert_maps[i + 1] @ top_maps[i] if i < len(top_maps) else None
            rhs = bot_maps[i] @ vert_maps[i] if i < len(bot_maps) else None
            if lhs is not None and rhs is not None:
                if not np.allclose(lhs, rhs):
                    commutes = False
                    break

    # 5-Lemma: f₁,f₂,f₄,f₅ Iso → f₃ Iso
    outer_isos = iso_status[0] and iso_status[1] and iso_status[3] and iso_status[4]
    f3_is_iso = iso_status[2]

    return {
        'commutes': commutes,
        'iso_status': iso_status,
        'outer_isos': outer_isos,
        'f3_is_iso': f3_is_iso,
        'five_lemma_conclusion': outer_isos and commutes,
        'verified': outer_isos and f3_is_iso
    }


def snake_lemma(A: int, B: int, C: int, Ap: int, Bp: int, Cp: int,
                maps: dict) -> dict:
    """
    Das Schlangen-Lemma (Snake Lemma).

    Gegeben ein kommutatives Diagramm:
        0 → A  →^f  B  →^g  C  → 0
            |α       |β      |γ
        0 → A' →^f' B' →^g' C' → 0

    mit exakten Zeilen, erzeugt das Schlangen-Lemma eine lange exakte Sequenz:
        0 → ker(α) → ker(β) → ker(γ) →^δ coker(α) → coker(β) → coker(γ) → 0

    Für zyklische Gruppen ℤ/nℤ berechnen wir die Kerne und Kokerne explizit.

    @param A, B, C     Ränge (Ordnungen) der oberen Zeile
    @param Ap, Bp, Cp  Ränge der unteren Zeile
    @param maps        Dictionary mit 'alpha', 'beta', 'gamma', 'f', 'g', 'fp', 'gp'
    @return            Dictionary mit der langen exakten Sequenz
    @lastModified 2026-03-10
    """
    alpha = np.array(maps.get('alpha', [[1]]), dtype=int)
    beta = np.array(maps.get('beta', [[1]]), dtype=int)
    gamma = np.array(maps.get('gamma', [[1]]), dtype=int)

    # Berechne Kern und Kokern für die senkrechten Abbildungen
    def compute_kernel_rank(M, domain_rank):
        """Rang des Kerns: dim_ker = domain_rank - rank(M)"""
        if M.size == 0:
            return domain_rank
        _, r, _ = _smith_normal_form(M.tolist())
        return max(0, domain_rank - r)

    def compute_cokernel_rank(M, codomain_rank):
        """Rang des Kokerns: dim_coker = codomain_rank - rank(M)"""
        if M.size == 0:
            return codomain_rank
        _, r, _ = _smith_normal_form(M.tolist())
        return max(0, codomain_rank - r)

    ker_alpha = compute_kernel_rank(alpha, A)
    ker_beta = compute_kernel_rank(beta, B)
    ker_gamma = compute_kernel_rank(gamma, C)

    coker_alpha = compute_cokernel_rank(alpha, Ap)
    coker_beta = compute_cokernel_rank(beta, Bp)
    coker_gamma = compute_cokernel_rank(gamma, Cp)

    # Verbindungsmorphismus δ: ker(γ) → coker(α)
    # Existenz ist durch das Schlangen-Lemma garantiert

    return {
        'ker_alpha': ker_alpha,
        'ker_beta': ker_beta,
        'ker_gamma': ker_gamma,
        'coker_alpha': coker_alpha,
        'coker_beta': coker_beta,
        'coker_gamma': coker_gamma,
        'long_exact_sequence': [
            ('0', 0),
            ('ker(α)', ker_alpha),
            ('ker(β)', ker_beta),
            ('ker(γ)', ker_gamma),
            ('δ (Verbindungsmorphismus)', None),
            ('coker(α)', coker_alpha),
            ('coker(β)', coker_beta),
            ('coker(γ)', coker_gamma),
            ('0', 0)
        ]
    }


# ---------------------------------------------------------------------------
# Freie Auflösungen
# ---------------------------------------------------------------------------

def free_resolution(n: int, relations: list) -> ChainComplex:
    """
    Erzeugt die minimale freie Auflösung von ℤ/nℤ.

    Die freie Auflösung von ℤ/nℤ ist:
        0 → ℤ →^n ℤ → ℤ/nℤ → 0

    Die Randmatrix der Multiplikation mit n ist einfach die 1×1-Matrix [n].

    @param n           Modulus (n ≥ 1)
    @param relations   Relationsmatrix (für allgemeinere Moduln)
    @return            ChainComplex der freien Auflösung
    @lastModified 2026-03-10
    """
    if n <= 0:
        raise ValueError(f"Modulus n muss positiv sein, aber {n} erhalten.")

    if n == 1:
        # ℤ/1ℤ = 0: triviale Auflösung
        return ChainComplex(
            groups={0: 0},
            boundaries={}
        )

    # Standardauflösung von ℤ/nℤ:
    # C₁ = ℤ, C₀ = ℤ
    # ∂₁: C₁ → C₀ ist die Multiplikation mit n: 1 ↦ n
    groups = {0: 1, 1: 1}

    # Randmatrix ∂₁ = [n] (1×1-Matrix)
    boundary_matrix = [[n]]

    boundaries = {1: boundary_matrix}

    return ChainComplex(groups, boundaries)


def projective_dimension(module_rank: int, relations: list) -> int:
    """
    Berechnet die projektive Dimension eines Moduls.

    Die projektive Dimension pd(M) ist die Länge der kürzesten projektiven Auflösung.

    Für ℤ/nℤ (n > 1): pd(ℤ/nℤ) = 1 (Auflösung der Länge 1)
    Für ℤ selbst: pd(ℤ) = 0 (bereits projektiv)
    Für 0: pd(0) = -∞ (konventionell -1)

    @param module_rank   Rang des Moduls (0 = torsion, 1 = zyklisch, ...)
    @param relations     Relationsmatrix (leere Liste = freier Modul)
    @return              Projektive Dimension (0, 1, oder 2 für endlich erzeugte ℤ-Moduln)
    @lastModified 2026-03-10
    """
    if module_rank == 0 and len(relations) == 0:
        # Nullmodul: Dimension -1 (Konvention)
        return -1

    if len(relations) == 0:
        # Freier Modul: projektiv, Dimension 0
        return 0

    # Prüfe ob der Modul selbst bereits projektiv (frei) ist
    R = np.array(relations, dtype=int)
    _, rank_R, elem_div = _smith_normal_form(R.tolist())

    # Torsionskoeffizienten prüfen
    has_torsion = any(d > 1 for d in elem_div)

    if not has_torsion and rank_R == module_rank:
        return 0  # Freier Modul

    if has_torsion:
        return 1  # ℤ/nℤ hat projektive Dimension 1

    return 1  # Allgemeiner endlich erzeugter ℤ-Modul


# ---------------------------------------------------------------------------
# Ext und Tor
# ---------------------------------------------------------------------------

def ext_group(n: int, m: int, degree: int) -> dict:
    """
    Berechnet Ext^k(ℤ/nℤ, ℤ/mℤ) für Grad k = degree.

    Die Ext-Gruppen messen die Verlängerungen von Moduln:

    Ergebnisse via freier Auflösung 0 → ℤ →^n ℤ → ℤ/nℤ → 0:
    - k=0: Ext⁰ = Hom(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
    - k=1: Ext¹(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
    - k≥2: Ext^k = 0 (da pd(ℤ/nℤ) = 1)

    @param n      Modulus des ersten Arguments
    @param m      Modulus des zweiten Arguments
    @param degree Grad k der Ext-Gruppe
    @return       Dictionary {'group': 'ℤ/dℤ', 'order': d, 'trivial': bool}
    @lastModified 2026-03-10
    """
    g = _gcd(n, m)

    if degree == 0:
        # Hom(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
        # Ein Homomorphismus f: ℤ/nℤ → ℤ/mℤ ist durch f(1) bestimmt,
        # wobei n·f(1) = 0 in ℤ/mℤ, d.h. m | n·f(1)
        return {
            'group': f'ℤ/{g}ℤ',
            'order': g,
            'trivial': g == 1,
            'description': f'Hom(ℤ/{n}ℤ, ℤ/{m}ℤ) ≅ ℤ/{g}ℤ'
        }
    elif degree == 1:
        # Ext¹(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
        # Berechnet aus der freien Auflösung: Anwenden von Hom(·, ℤ/mℤ) liefert
        # 0 → Hom(ℤ,ℤ/mℤ) →^{n·} Hom(ℤ,ℤ/mℤ) → ...
        # Der Kokern ist ℤ/mℤ / (n · ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
        return {
            'group': f'ℤ/{g}ℤ',
            'order': g,
            'trivial': g == 1,
            'description': f'Ext¹(ℤ/{n}ℤ, ℤ/{m}ℤ) ≅ ℤ/{g}ℤ'
        }
    else:
        # k ≥ 2: Ext^k = 0, da pd(ℤ/nℤ) = 1
        return {
            'group': '0',
            'order': 1,
            'trivial': True,
            'description': f'Ext^{degree}(ℤ/{n}ℤ, ℤ/{m}ℤ) = 0 (pd = 1)'
        }


def tor_group(n: int, m: int, degree: int) -> dict:
    """
    Berechnet Tor_k(ℤ/nℤ, ℤ/mℤ) für Grad k = degree.

    Die Tor-Gruppen messen die Torsion im Tensorprodukt:

    Ergebnisse via freier Auflösung von ℤ/nℤ:
    - k=0: Tor₀ = ℤ/nℤ ⊗ ℤ/mℤ ≅ ℤ/gcd(n,m)ℤ
    - k=1: Tor₁(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
    - k≥2: Tor_k = 0 (da pd(ℤ/nℤ) = 1)

    @param n      Modulus des ersten Arguments
    @param m      Modulus des zweiten Arguments
    @param degree Grad k der Tor-Gruppe
    @return       Dictionary {'group': 'ℤ/dℤ', 'order': d, 'trivial': bool}
    @lastModified 2026-03-10
    """
    g = _gcd(n, m)

    if degree == 0:
        # ℤ/nℤ ⊗_ℤ ℤ/mℤ ≅ ℤ/gcd(n,m)ℤ
        return {
            'group': f'ℤ/{g}ℤ',
            'order': g,
            'trivial': g == 1,
            'description': f'ℤ/{n}ℤ ⊗ ℤ/{m}ℤ ≅ ℤ/{g}ℤ'
        }
    elif degree == 1:
        # Tor₁(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ
        # Aus der freien Auflösung: Tensorieren mit ℤ/mℤ liefert
        # ℤ/mℤ →^{n·} ℤ/mℤ → ...
        # Der Kern von (n·): ℤ/mℤ → ℤ/mℤ hat Ordnung gcd(n,m)
        return {
            'group': f'ℤ/{g}ℤ',
            'order': g,
            'trivial': g == 1,
            'description': f'Tor₁(ℤ/{n}ℤ, ℤ/{m}ℤ) ≅ ℤ/{g}ℤ'
        }
    else:
        # k ≥ 2: Tor_k = 0
        return {
            'group': '0',
            'order': 1,
            'trivial': True,
            'description': f'Tor_{degree}(ℤ/{n}ℤ, ℤ/{m}ℤ) = 0 (pd = 1)'
        }


# ---------------------------------------------------------------------------
# Universeller Koeffizienten-Satz
# ---------------------------------------------------------------------------

def universal_coefficient_theorem(betti: dict, torsion: dict) -> dict:
    """
    Wendet den Universellen Koeffizienten-Satz (UCT) an.

    Der UCT besagt für die Kohomologie mit ganzzahligen Koeffizienten:

        H^n(X; ℤ) ≅ Hom(H_n(X), ℤ) ⊕ Ext¹(H_{n-1}(X), ℤ)

    Dabei gilt:
    - Hom(ℤ^β ⊕ Torsion, ℤ) ≅ ℤ^β (freier Teil)
    - Ext¹(ℤ/dℤ, ℤ) ≅ ℤ/dℤ (Torsion aus vorherigem Grad)

    Kurzform: Die Kohomologie "erbt" die Torsion einen Grad höher.

    @param betti    Dictionary {n: β_n} mit Betti-Zahlen
    @param torsion  Dictionary {n: [d₁,...,dₖ]} mit Torsionskoeffizienten
    @return         Dictionary mit Kohomologiegruppen pro Grad
    @lastModified 2026-03-10
    """
    result = {}

    # Bestimme alle vorkommenden Grade
    all_degrees = set(betti.keys()) | set(torsion.keys())
    if all_degrees:
        max_degree = max(all_degrees) + 1
    else:
        max_degree = 1

    for n in range(max_degree + 1):
        # Freier Teil: Hom(H_n, ℤ) ≅ ℤ^{β_n}
        free_rank = betti.get(n, 0)

        # Torsion aus Ext¹(H_{n-1}, ℤ):
        # Für H_{n-1} = ℤ/dℤ gilt Ext¹(ℤ/dℤ, ℤ) ≅ ℤ/dℤ
        prev_torsion = torsion.get(n - 1, [])

        # Kohomologiegruppe: H^n ≅ ℤ^{free_rank} ⊕ (⊕ ℤ/dᵢℤ)
        if free_rank == 0 and len(prev_torsion) == 0:
            group_str = '0'
        elif free_rank > 0 and len(prev_torsion) == 0:
            group_str = f'ℤ^{free_rank}' if free_rank > 1 else 'ℤ'
        elif free_rank == 0 and len(prev_torsion) > 0:
            tors_str = ' ⊕ '.join([f'ℤ/{d}ℤ' for d in prev_torsion])
            group_str = tors_str
        else:
            free_str = f'ℤ^{free_rank}' if free_rank > 1 else 'ℤ'
            tors_str = ' ⊕ '.join([f'ℤ/{d}ℤ' for d in prev_torsion])
            group_str = f'{free_str} ⊕ {tors_str}'

        result[n] = {
            'free_rank': free_rank,
            'torsion': prev_torsion,
            'group': group_str
        }

    return result
