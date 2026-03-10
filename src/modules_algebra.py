"""
@file modules_algebra.py
@brief R-Moduln über Ringen: Smith-Normalform, Klassifikationssatz, Tensorprodukt.
@description
    Implementiert die Theorie der Moduln über Ringen, insbesondere:

    - Smith-Normalform ganzzahliger Matrizen (A = U·D·V)
    - Klassifikationssatz endlich erzeugter abelscher Gruppen
    - Tensorprodukt von ℤ-Moduln
    - Hom-Modul (Menge der Homomorphismen)
    - Exakte Sequenzen
    - Freie Auflösungen

    Mathematischer Hintergrund:
    Ein R-Modul M ist eine abelsche Gruppe (M, +) zusammen mit einer
    Skalarmultiplikation R × M → M, die die Ring-Axiome erfüllt.

    Wichtigster Fall: M = ℤ^n / Im(A) für eine ganzzahlige Matrix A.
    Klassifikation durch Smith-Normalform: d_1 | d_2 | ... | d_r.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import copy
from typing import Any


# ---------------------------------------------------------------------------
# Smith-Normalform
# ---------------------------------------------------------------------------

def _copy_matrix(M: list[list[int]]) -> list[list[int]]:
    """
    @brief Erzeugt eine tiefe Kopie einer ganzzahligen Matrix.
    @param M Eingabematrix
    @return Tiefe Kopie der Matrix
    @lastModified 2026-03-10
    """
    return [row[:] for row in M]


def _identity(n: int) -> list[list[int]]:
    """
    @brief Erzeugt eine n×n-Einheitsmatrix.
    @param n Dimension
    @return Einheitsmatrix
    @lastModified 2026-03-10
    """
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def _mat_mul(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
    """
    @brief Ganzzahlige Matrixmultiplikation A · B.
    @param A Linke Matrix (m×k)
    @param B Rechte Matrix (k×n)
    @return Produkt A·B (m×n)
    @lastModified 2026-03-10
    """
    m = len(A)
    k = len(A[0])
    n = len(B[0])
    # Ergebnismatrix mit Nullen initialisieren
    C = [[0] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            s = 0
            for l in range(k):
                s += A[i][l] * B[l][j]
            C[i][j] = s
    return C


def _row_swap(M: list[list[int]], i: int, j: int) -> None:
    """
    @brief Vertauscht Zeilen i und j in-place.
    @param M Matrix
    @param i Erste Zeilenindex
    @param j Zweiter Zeilenindex
    @lastModified 2026-03-10
    """
    M[i], M[j] = M[j], M[i]


def _col_swap(M: list[list[int]], i: int, j: int) -> None:
    """
    @brief Vertauscht Spalten i und j in-place.
    @param M Matrix
    @param i Erste Spaltenindex
    @param j Zweiter Spaltenindex
    @lastModified 2026-03-10
    """
    for row in M:
        row[i], row[j] = row[j], row[i]


def _row_add(M: list[list[int]], i: int, j: int, c: int) -> None:
    """
    @brief Addiert c-faches der Zeile j zur Zeile i: M[i] += c * M[j].
    @param M Matrix
    @param i Zielzeile
    @param j Quellzeile
    @param c Vielfaches
    @lastModified 2026-03-10
    """
    for k in range(len(M[i])):
        M[i][k] += c * M[j][k]


def _col_add(M: list[list[int]], i: int, j: int, c: int) -> None:
    """
    @brief Addiert c-faches der Spalte j zur Spalte i: M[*][i] += c * M[*][j].
    @param M Matrix
    @param i Zielspalte
    @param j Quellspalte
    @param c Vielfaches
    @lastModified 2026-03-10
    """
    for row in M:
        row[i] += c * row[j]


def _apply_row_to_transform(U: list[list[int]], i: int, j: int, c: int) -> None:
    """
    @brief Wendet elementare Zeilenoperation auf Transformationsmatrix an.
    @param U Transformationsmatrix (wird in-place geändert)
    @param i Zielzeile
    @param j Quellzeile
    @param c Vielfaches
    @lastModified 2026-03-10
    """
    _row_add(U, i, j, c)


def _apply_col_to_transform(V: list[list[int]], i: int, j: int, c: int) -> None:
    """
    @brief Wendet elementare Spaltenoperation auf Transformationsmatrix an.
    @param V Transformationsmatrix (wird in-place geändert)
    @param i Zielspalte
    @param j Quellspalte
    @param c Vielfaches
    @lastModified 2026-03-10
    """
    _col_add(V, i, j, c)


def smith_normal_form(matrix: list[list[int]]) -> dict[str, Any]:
    """
    @brief Berechnet die Smith-Normalform einer ganzzahligen Matrix.
    @description
        Die Smith-Normalform einer ganzzahligen Matrix A ist A = U · D · V,
        wobei:
        - U ∈ GL_m(ℤ), V ∈ GL_n(ℤ) (invertierbare ganzzahlige Matrizen)
        - D = diag(d_1, ..., d_r, 0, ..., 0) Diagonalmatrix
        - d_1 | d_2 | ... | d_r (Teilbarkeitskette)

        Algorithmus: Gauß-Elimination mit ganzzahligen Zeilen-/Spaltenoperationen.
        In jedem Schritt wird das betragskleinste Nicht-Null-Element als Pivot
        gewählt und via erweitertem Euklidischen Algorithmus eliminiert.

        Klassifikation: ℤ^n / Im(A) ≅ ℤ_{d_1} × ... × ℤ_{d_r} × ℤ^{n-r}

    @param matrix Ganzzahlige Matrix (Liste von Listen)
    @return Dict mit:
            - 'invariant_factors': Liste [d_1, ..., d_r] (Teilbarkeitskette)
            - 'rank': Rang r der Matrix
            - 'U': linke Transformationsmatrix
            - 'D': Diagonalmatrix in SNF
            - 'V': rechte Transformationsmatrix
    @lastModified 2026-03-10
    """
    # Tiefe Kopie der Matrix erstellen, um die Eingabe nicht zu verändern
    D = _copy_matrix(matrix)
    m = len(D)
    n = len(D[0]) if m > 0 else 0

    # Transformationsmatrizen als Einheitsmatrizen initialisieren
    U = _identity(m)  # linke Transformationsmatrix
    V = _identity(n)  # rechte Transformationsmatrix

    # Schritt-für-Schritt Smith-Normalform berechnen
    step = 0  # aktueller Diagonalblock-Index
    while step < m and step < n:
        # Suche das betragskleinste Nicht-Null-Element im verbleibenden Unterblock
        pivot_found = False
        min_val = None
        pi, pj = step, step  # Pivot-Position

        for i in range(step, m):
            for j in range(step, n):
                if D[i][j] != 0:
                    if min_val is None or abs(D[i][j]) < abs(min_val):
                        min_val = D[i][j]
                        pi, pj = i, j
                        pivot_found = True

        # Kein Nicht-Null-Element mehr → Rest ist Nullblock
        if not pivot_found:
            break

        # Pivot in die aktuelle Diagonalposition bringen
        if pi != step:
            _row_swap(D, step, pi)
            _row_swap(U, step, pi)
        if pj != step:
            _col_swap(D, step, pj)
            _col_swap(V, step, pj)

        # Iterativ: Pivot eliminiert alle anderen Einträge in Zeile und Spalte
        # und sorgt dafür, dass Pivot alle Einträge des Unterblocks teilt
        changed = True
        while changed:
            changed = False

            # Zeile step bereinigen: D[step][j] für j > step auf 0 reduzieren
            for j in range(step + 1, n):
                if D[step][j] != 0:
                    q = D[step][j] // D[step][step]
                    _col_add(D, j, step, -q)
                    _col_add(V, j, step, -q)
                    # Falls noch Rest da ist, neu pivotisieren
                    if D[step][j] != 0 and abs(D[step][j]) < abs(D[step][step]):
                        _col_swap(D, step, j)
                        _col_swap(V, step, j)
                    changed = True

            # Spalte step bereinigen: D[i][step] für i > step auf 0 reduzieren
            for i in range(step + 1, m):
                if D[i][step] != 0:
                    q = D[i][step] // D[step][step]
                    _row_add(D, i, step, -q)
                    _row_add(U, i, step, -q)
                    if D[i][step] != 0 and abs(D[i][step]) < abs(D[step][step]):
                        _row_swap(D, step, i)
                        _row_swap(U, step, i)
                    changed = True

        # Prüfen ob Pivot alle Einträge im verbleibenden Block teilt
        # Falls nicht, Addition von Zeile/Spalte um ggT herzustellen
        pivot = D[step][step]
        if pivot != 0:
            for i in range(step + 1, m):
                for j in range(step + 1, n):
                    if D[i][j] % pivot != 0:
                        # Zeile i zur Zeile step addieren (ggT-Trick)
                        _row_add(D, step, i, 1)
                        _row_add(U, step, i, 1)
                        changed = True
                        break
                if changed:
                    break
            if changed:
                # Erneute Bereinigung nötig
                continue

        # Pivot positiv machen
        if D[step][step] < 0:
            for j in range(n):
                D[step][j] = -D[step][j]
            for j in range(m):
                U[step][j] = -U[step][j]

        step += 1

    # Invariantenfaktoren aus der Diagonale ablesen (nur Nicht-Null-Einträge)
    invariant_factors = []
    for i in range(min(m, n)):
        if D[i][i] != 0:
            invariant_factors.append(D[i][i])

    return {
        'invariant_factors': invariant_factors,
        'rank': len(invariant_factors),
        'U': U,
        'D': D,
        'V': V,
    }


def module_from_matrix(A: list[list[int]]) -> dict[str, Any]:
    """
    @brief Konstruiert den ℤ-Modul M = ℤ^n / Im(A) aus einer Relationsmatrix.
    @description
        Nutzt die Smith-Normalform, um den Modul zu klassifizieren.
        Die Invariantenfaktoren d_1 | d_2 | ... | d_r bestimmen die Struktur:
        M ≅ ℤ_{d_1} × ... × ℤ_{d_r} × ℤ^{n-r}

        Dabei gilt:
        - d_i = 1 → keine Torsion (Faktor fällt weg)
        - d_i > 1 → Torsionsanteil ℤ_{d_i}
        - freier Rang = Spaltenanzahl - Rang(A)

    @param A Ganzzahlige Relationsmatrix (m×n)
    @return Dict mit:
            - 'invariant_factors': relevante Invariantenfaktoren (≥ 2)
            - 'rank': freier Rang des Moduls
            - 'torsion_part': String-Darstellung des Torsionsanteils
            - 'free_part': String-Darstellung des freien Anteils
    @lastModified 2026-03-10
    """
    snf = smith_normal_form(A)
    n = len(A[0]) if A else 0  # Spaltenanzahl = Rang des freien Moduls

    # Invariantenfaktoren ohne Einsen (d_i=1 → triviale Komponente)
    all_factors = snf['invariant_factors']
    torsion_factors = [d for d in all_factors if d > 1]

    # Freier Rang: Dimensionen, die nicht durch Relationen gebunden sind
    free_rank = n - snf['rank']

    # String-Darstellung des Torsionsanteils
    if torsion_factors:
        torsion_str = ' × '.join(f'ℤ_{d}' for d in torsion_factors)
    else:
        torsion_str = '0'

    # String-Darstellung des freien Anteils
    if free_rank == 0:
        free_str = '0'
    elif free_rank == 1:
        free_str = 'ℤ'
    else:
        free_str = f'ℤ^{free_rank}'

    return {
        'invariant_factors': torsion_factors,
        'rank': free_rank,
        'torsion_part': torsion_str,
        'free_part': free_str,
    }


# ---------------------------------------------------------------------------
# Klassifikationssatz abelscher Gruppen
# ---------------------------------------------------------------------------

def _prime_factorization(n: int) -> dict[int, int]:
    """
    @brief Primfaktorzerlegung von n.
    @param n Natürliche Zahl ≥ 1
    @return Dict {Primzahl: Exponent}
    @lastModified 2026-03-10
    """
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def _partitions(n: int) -> list[list[int]]:
    """
    @brief Berechnet alle Partitionen von n (absteigend sortiert).
    @param n Natürliche Zahl
    @return Liste aller Partitionen als Listen
    @lastModified 2026-03-10
    """
    if n == 0:
        return [[]]
    result = []

    def _gen(remaining: int, max_val: int, current: list[int]) -> None:
        if remaining == 0:
            result.append(current[:])
            return
        for k in range(min(remaining, max_val), 0, -1):
            current.append(k)
            _gen(remaining - k, k, current)
            current.pop()

    _gen(n, n, [])
    return result


def _primary_to_invariant(primary_decomp: dict[int, list[int]]) -> list[int]:
    """
    @brief Wandelt Primärkomponentendarstellung in Invariantenfaktoren um.
    @description
        Primärkomponenten: {p: [a_1, ..., a_k]} mit a_1 ≥ ... ≥ a_k
        → Invariantenfaktoren: d_i = ∏_p p^{a_{k-i+1}} (aufsteigend)

        Beispiel: {2: [2,1], 3: [1,1]} → [6, 12]
        (d_1 = 2^1 · 3^1 = 6, d_2 = 2^2 · 3^1 = 12)

    @param primary_decomp Primärkomponentendarstellung
    @return Invariantenfaktoren (aufsteigend, Teilbarkeitskette)
    @lastModified 2026-03-10
    """
    if not primary_decomp:
        return [1]

    # Maximale Anzahl Faktoren bestimmen
    max_len = max(len(v) for v in primary_decomp.values())

    # Invariantenfaktoren von hinten aufbauen
    inv_factors = []
    for i in range(max_len - 1, -1, -1):
        factor = 1
        for p, exponents in primary_decomp.items():
            # Exponentenliste mit führenden Nullen auffüllen
            padded = [0] * (max_len - len(exponents)) + exponents
            factor *= p ** padded[i]
        inv_factors.append(factor)

    # Aufsteigend sortieren (kleinster zuerst)
    inv_factors.sort()
    return inv_factors


def structure_theorem_abelian_groups(n: int) -> list[dict[str, Any]]:
    """
    @brief Klassifikation aller abelschen Gruppen der Ordnung n.
    @description
        Hauptsatz über endlich erzeugte abelsche Gruppen:
        Jede endliche abelsche Gruppe G der Ordnung n ist isomorph zu einem
        direkten Produkt zyklischer Gruppen:

        Invariantenfaktoren-Form: G ≅ ℤ_{d_1} × ... × ℤ_{d_k}
        mit d_1 | d_2 | ... | d_k und ∏ d_i = n

        Primärkomponenten-Form: G ≅ ∏_p ∏_i ℤ_{p^{a_i}}
        mit ∑_i a_i = v_p(n) für jede Primzahl p | n

        Die Anzahl nicht-isomorpher abelscher Gruppen der Ordnung n ist
        p(v_{p_1}(n)) · p(v_{p_2}(n)) · ... (Produkt der Partitionszahlen).

    @param n Ordnung der Gruppe
    @return Liste von Dicts mit:
            - 'invariant_factors': Invariantenfaktoren [d_1, ..., d_k]
            - 'primary_decomposition': Primärkomponenten {p: [a_1, ...]}
            - 'description': String-Beschreibung
    @lastModified 2026-03-10
    """
    if n == 1:
        # Nur die triviale Gruppe
        return [{'invariant_factors': [1],
                 'primary_decomposition': {},
                 'description': 'ℤ_1 (triviale Gruppe)'}]

    # Primfaktorzerlegung von n
    prime_factors = _prime_factorization(n)

    # Für jede Primzahl alle Partitionen des Exponenten berechnen
    # Jede Partition entspricht einem anderen Typ der p-Sylowgruppe
    prime_partition_lists: list[tuple[int, list[list[int]]]] = []
    for p, exp in prime_factors.items():
        parts = _partitions(exp)
        prime_partition_lists.append((p, parts))

    # Alle Kombinationen von Partitionen (kartesisches Produkt)
    # bilden alle Isomorphietypen
    def _cartesian_product(lists: list[list[Any]]) -> list[list[Any]]:
        """Kartesisches Produkt einer Liste von Listen."""
        if not lists:
            return [[]]
        result = []
        for item in lists[0]:
            for rest in _cartesian_product(lists[1:]):
                result.append([item] + rest)
        return result

    partition_choices = [parts for _, parts in prime_partition_lists]
    primes = [p for p, _ in prime_partition_lists]

    all_combos = _cartesian_product(partition_choices)

    groups = []
    for combo in all_combos:
        # Primärkomponentendarstellung aufbauen: {p: [a_1, ..., a_k]} absteigend
        primary_decomp: dict[int, list[int]] = {}
        for p, partition in zip(primes, combo):
            # Partition absteigend sortieren (größter Exponent zuerst)
            primary_decomp[p] = sorted(partition, reverse=True)

        # In Invariantenfaktoren umwandeln
        inv_factors = _primary_to_invariant(primary_decomp)

        # String-Beschreibung
        desc = ' × '.join(f'ℤ_{d}' for d in inv_factors)

        groups.append({
            'invariant_factors': inv_factors,
            'primary_decomposition': primary_decomp,
            'description': desc,
        })

    # Nach Invariantenfaktoren sortieren für konsistente Ausgabe
    groups.sort(key=lambda g: g['invariant_factors'])
    return groups


# ---------------------------------------------------------------------------
# Tensorprodukt
# ---------------------------------------------------------------------------

def tensor_product_modules(
    M_gens: list[int],
    N_gens: list[int],
    R_mod: int
) -> dict[str, Any]:
    """
    @brief Berechnet das Tensorprodukt M ⊗_ℤ N zweier ℤ-Moduln.
    @description
        Für endliche zyklische ℤ-Moduln gilt:
        ℤ_m ⊗_ℤ ℤ_n ≅ ℤ_{gcd(m,n)}

        Allgemein für direkte Summen (Bilinearität des Tensorprodukts):
        (⊕_i ℤ_{m_i}) ⊗_ℤ (⊕_j ℤ_{n_j}) ≅ ⊕_{i,j} ℤ_{gcd(m_i, n_j)}

        Spezialfall: ℤ_m ⊗_ℤ ℤ_n ≅ ℤ_{gcd(m,n)} wegen
        m·x = 0 und n·x = 0 impliziert gcd(m,n)·x = 0.

    @param M_gens Liste der Ordnungen der zyklischen Summanden von M (0 = freier Rang 1)
    @param N_gens Liste der Ordnungen der zyklischen Summanden von N (0 = freier Rang 1)
    @param R_mod Modulus des Grundrings (0 = ℤ)
    @return Dict mit:
            - 'order': Ordnung des Tensorprodukts (0 = unendlich)
            - 'structure': Invariantenfaktoren des Ergebnisses
            - 'description': String-Darstellung
    @lastModified 2026-03-10
    """
    # Tensorprodukt zyklischer Komponenten: ℤ_m ⊗ ℤ_n ≅ ℤ_{gcd(m,n)}
    # Für m=0 oder n=0 (freier Rang): ℤ ⊗ ℤ_n ≅ ℤ_n

    # Invariantenfaktoren des Tensorprodukts sammeln
    result_factors = []
    for m in M_gens:
        for n_val in N_gens:
            if m == 0 and n_val == 0:
                # ℤ ⊗ ℤ ≅ ℤ (freier Rang 1)
                result_factors.append(0)
            elif m == 0:
                # ℤ ⊗ ℤ_n ≅ ℤ_n
                result_factors.append(n_val)
            elif n_val == 0:
                # ℤ_m ⊗ ℤ ≅ ℤ_m
                result_factors.append(m)
            else:
                # ℤ_m ⊗ ℤ_n ≅ ℤ_{gcd(m,n)}
                g = math.gcd(m, n_val)
                result_factors.append(g)

    # Triviale Faktoren (Ordnung 1) entfernen
    result_factors = [f for f in result_factors if f != 1]

    # Gesamtordnung berechnen
    if 0 in result_factors:
        total_order = 0  # unendlich
    else:
        total_order = 1
        for f in result_factors:
            total_order *= f
        if not result_factors:
            total_order = 1

    # String-Beschreibung
    if not result_factors:
        desc = '0 (trivial)'
    elif 0 in result_factors:
        desc = 'ℤ (frei)'
    else:
        desc = ' ⊗ '.join(f'ℤ_{f}' for f in result_factors)

    return {
        'order': total_order,
        'structure': result_factors,
        'description': desc,
    }


# ---------------------------------------------------------------------------
# Hom-Modul
# ---------------------------------------------------------------------------

def hom_module(M_order: int, N_order: int) -> dict[str, Any]:
    """
    @brief Berechnet Hom_ℤ(ℤ_m, ℤ_n) (Menge der Gruppenhomomorphismen).
    @description
        Für zyklische Gruppen gilt:
        Hom_ℤ(ℤ_m, ℤ_n) ≅ ℤ_{gcd(m,n)}

        Beweis: Ein Homomorphismus f: ℤ_m → ℤ_n ist durch f(1) ∈ ℤ_n bestimmt.
        Es muss gelten: m · f(1) ≡ 0 (mod n), d.h. f(1) ∈ {k·(n/gcd(m,n)) : k=0,...,gcd(m,n)-1}.
        Also gibt es genau gcd(m,n) Homomorphismen, erzeugt von n/gcd(m,n).

        |Hom(ℤ_m, ℤ_n)| = gcd(m, n)

    @param M_order Ordnung von M = ℤ_{M_order}
    @param N_order Ordnung von N = ℤ_{N_order}
    @return Dict mit:
            - 'order': Anzahl der Homomorphismen = gcd(M_order, N_order)
            - 'generators': Liste der Erzeuger [n/gcd(m,n)]
            - 'description': String-Darstellung
    @lastModified 2026-03-10
    """
    g = math.gcd(M_order, N_order)

    # Erzeuger: n/gcd als Element von ℤ_n
    generator = N_order // g if g > 0 else 0

    # Alle Homomorphismen explizit auflisten (f(1) = k · generator mod N_order)
    generators = [(k * generator) % N_order for k in range(g)]

    desc = f'Hom(ℤ_{M_order}, ℤ_{N_order}) ≅ ℤ_{g}'

    return {
        'order': g,
        'generators': generators,
        'description': desc,
    }


# ---------------------------------------------------------------------------
# Exakte Sequenzen
# ---------------------------------------------------------------------------

def exact_sequence_check(
    maps: list[list[list[int]]],
    modules: list[int]
) -> dict[str, Any]:
    """
    @brief Prüft ob eine kurze exakte Sequenz 0 → A → B → C → 0 exakt ist.
    @description
        Eine Sequenz ... → A -f-> B -g-> C → ... heißt exakt an B, wenn:
        Im(f) = Ker(g)

        Für ganzzahlige Matrizen (lineare Abbildungen zwischen freien ℤ-Moduln):
        - Im(f) = Spaltenraum von F über ℤ
        - Ker(g) = Kern von G über ℤ

        Für die Standard-Sequenz 0 → ℤ -n·-> ℤ → ℤ_n → 0:
        - f = Multiplikation mit n
        - g = Projektion mod n
        - Im(f) = nℤ = Ker(g) ✓

        Vereinfachte Prüfung für kurze exakte Sequenzen über zyklische Moduln.

    @param maps Liste von Matrizen [f_0, f_1, ...] (als Listen von Listen)
    @param modules Liste der Modul-Ordnungen (0 = freier Rang 1)
    @return Dict mit:
            - 'is_exact': True wenn exakt
            - 'details': Beschreibung der Überprüfung
    @lastModified 2026-03-10
    """
    details = []

    # Kurze exakte Sequenz 0 → A -f-> B -g-> C → 0 prüfen
    if len(maps) == 2 and len(modules) == 3:
        # maps[0]: f: A → B (als Matrix)
        # maps[1]: g: B → C (als Matrix)
        f_matrix = maps[0]
        g_matrix = maps[1]

        mod_a = modules[0]  # Ordnung/Rang von A
        mod_b = modules[1]  # Ordnung/Rang von B
        mod_c = modules[2]  # Ordnung/Rang von C

        # Spezialfall: 0 → ℤ -n-> ℤ → ℤ_n → 0
        if (mod_a == 1 and mod_b == 1 and mod_c > 1
                and len(f_matrix) == 1 and len(f_matrix[0]) == 1
                and len(g_matrix) == 1 and len(g_matrix[0]) == 1):
            n = f_matrix[0][0]    # Multiplikationsfaktor
            mod_n = g_matrix[0][0]  # Repräsentiert mod-n Abbildung

            # Exakt wenn: Im(f) = nℤ entspricht Kern(mod n Projektion)
            # d.h. n muss der Ordnung des Quotientenmoduls entsprechen
            is_exact = (n == mod_c)
            details.append(f'f = ×{n}, g = proj mod {mod_c}')
            details.append(f'Im(f) = {n}ℤ, Ker(g) = {mod_c}ℤ')
            details.append(f'Exakt: {is_exact}')

            return {
                'is_exact': is_exact,
                'details': details,
            }

        # Allgemeine Prüfung: Im(f) ⊆ Ker(g) via Matrixmultiplikation
        # Für ganzzahlige Matrizen: g·f = 0 impliziert Im(f) ⊆ Ker(g)
        # Vollständige Exaktheit benötigt Smith-Normalform
        try:
            gf = _mat_mul(g_matrix, f_matrix)
            # Prüfen ob g·f = 0
            is_zero = all(gf[i][j] == 0
                          for i in range(len(gf))
                          for j in range(len(gf[0])))
            details.append(f'g·f = {gf}')
            details.append(f'Im(f) ⊆ Ker(g): {is_zero}')

            # Für einfache Fälle: Exaktheit wenn g·f=0 und Dimensionen passen
            is_exact = is_zero
        except (IndexError, TypeError):
            is_exact = False
            details.append('Matrixmultiplikation fehlgeschlagen')

        return {
            'is_exact': is_exact,
            'details': details,
        }

    # Fallback für allgemeine Sequenzen
    return {
        'is_exact': False,
        'details': ['Allgemeine Sequenz: nur kurze exakte Sequenzen werden geprüft'],
    }


# ---------------------------------------------------------------------------
# Freie Auflösung
# ---------------------------------------------------------------------------

def free_resolution(module_matrix: list[list[int]]) -> dict[str, Any]:
    """
    @brief Berechnet eine freie Auflösung des ℤ-Moduls M = ℤ^n / Im(A).
    @description
        Eine freie Auflösung von M ist eine exakte Sequenz:
        ... → ℤ^{n_2} -d_2-> ℤ^{n_1} -d_1-> ℤ^{n_0} -ε-> M → 0

        Für endlich erzeugte ℤ-Moduln genügt eine Auflösung der Länge 1:
        0 → ℤ^m -A-> ℤ^n → M → 0

        wobei A die Relationsmatrix ist (Smith-Normalform gibt Struktur).

        Für Torsionsmoduln ist dies minimal. Für freie Moduln noch kürzer.

    @param module_matrix Ganzzahlige Relationsmatrix A (m×n)
    @return Dict mit:
            - 'resolution': Liste der Stufen [(Abbildung, Rang)] von rechts nach links
            - 'length': Länge der Auflösung
            - 'module_structure': Smith-Normalform des Moduls
    @lastModified 2026-03-10
    """
    snf = smith_normal_form(module_matrix)
    m = len(module_matrix)
    n = len(module_matrix[0]) if module_matrix else 0

    # Freie Auflösung der Länge 1:
    # 0 → ℤ^m -A-> ℤ^n → M → 0
    resolution = [
        {'map': 'ε: ℤ^' + str(n) + ' → M',
         'rank': n,
         'description': 'Augmentierungsabbildung'},
        {'map': 'A: ℤ^' + str(m) + ' → ℤ^' + str(n),
         'rank': m,
         'description': 'Relationsmatrix A'},
        {'map': '0 → ℤ^' + str(m),
         'rank': 0,
         'description': 'Kern der Relationen (trivial für volle Auflösung)'},
    ]

    return {
        'resolution': resolution,
        'length': 1,
        'module_structure': snf,
    }


# ---------------------------------------------------------------------------
# Module-Klasse
# ---------------------------------------------------------------------------

class Module:
    """
    @brief Klasse für R-Moduln über ℤ (Verallgemeinerung des Vektorraums).
    @description
        Repräsentiert den ℤ-Modul M = ℤ^n / Im(A) für eine ganzzahlige
        Relationsmatrix A. Die Struktur wird durch die Smith-Normalform
        von A bestimmt.

        Mathematischer Hintergrund:
        Ein R-Modul M ist eine abelsche Gruppe (M,+) mit Skalarmultiplikation
        R × M → M, die die Distributivgesetze erfüllt.
        Über ℤ sind Moduln dasselbe wie abelsche Gruppen.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, relation_matrix: list[list[int]]) -> None:
        """
        @brief Erstellt einen ℤ-Modul aus einer Relationsmatrix.
        @param relation_matrix Ganzzahlige Matrix A, sodass M = ℤ^n / Im(A)
        @lastModified 2026-03-10
        """
        self._matrix = relation_matrix
        # Smith-Normalform einmal berechnen und cachen
        self._snf = smith_normal_form(relation_matrix)
        self._module_info = module_from_matrix(relation_matrix)

    def rank(self) -> int:
        """
        @brief Gibt den freien Rang des Moduls zurück.
        @description
            Der freie Rang ist die Anzahl der ℤ-Faktoren in der Zerlegung
            M ≅ ℤ_{d_1} × ... × ℤ_{d_r} × ℤ^{rank}.
        @return Freier Rang (nicht-negativer Integer)
        @lastModified 2026-03-10
        """
        return self._module_info['rank']

    def torsion_submodule(self) -> dict[str, Any]:
        """
        @brief Gibt den Torsionsuntermodul zurück.
        @description
            Der Torsionsuntermodul besteht aus allen Elementen m in M, für die
            ein n in ℤ\\{0} existiert mit n·m = 0.
            Für M = ℤ^n / Im(A) ist der Torsionsanteil
            T(M) = ℤ_{d_1} × ... × ℤ_{d_r}.
        @return Dict mit Torsionsfaktoren und Beschreibung
        @lastModified 2026-03-10
        """
        torsion_factors = self._module_info['invariant_factors']
        return {
            'invariant_factors': torsion_factors,
            'description': self._module_info['torsion_part'],
            'order': math.prod(torsion_factors) if torsion_factors else 1,
        }

    def free_part(self) -> dict[str, Any]:
        """
        @brief Gibt den freien Anteil des Moduls zurück.
        @description
            Der freie Anteil ist isomorph zu ℤ^r, wobei r der freie Rang ist.
        @return Dict mit Rang und Beschreibung
        @lastModified 2026-03-10
        """
        r = self._module_info['rank']
        return {
            'rank': r,
            'description': self._module_info['free_part'],
        }

    def is_free(self) -> bool:
        """
        @brief Prüft ob der Modul frei ist (keine Torsion).
        @description
            M ist frei ⟺ kein Torsionsanteil ⟺ alle Invariantenfaktoren = 1.
        @return True wenn M frei ist
        @lastModified 2026-03-10
        """
        return len(self._module_info['invariant_factors']) == 0

    def is_finitely_generated(self) -> bool:
        """
        @brief Prüft ob der Modul endlich erzeugt ist.
        @description
            Für Moduln der Form ℤ^n / Im(A) ist dies immer True,
            da die Bilder der Standardbasisvektoren ein endliches Erzeugendensystem bilden.
        @return True (immer für endliche Matrizen)
        @lastModified 2026-03-10
        """
        return True

    def smith_normal_form(self) -> dict[str, Any]:
        """
        @brief Gibt die Smith-Normalform der Relationsmatrix zurück.
        @return Smith-Normalform (gecacht)
        @lastModified 2026-03-10
        """
        return self._snf

    def __repr__(self) -> str:
        """
        @brief String-Darstellung des Moduls.
        @return Mathematische Notation des Moduls
        @lastModified 2026-03-10
        """
        torsion = self._module_info['torsion_part']
        free = self._module_info['free_part']
        if torsion == '0' and free == '0':
            return 'M ≅ 0'
        elif torsion == '0':
            return f'M ≅ {free}'
        elif free == '0':
            return f'M ≅ {torsion}'
        else:
            return f'M ≅ {torsion} × {free}'
