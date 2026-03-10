"""
@file algebraic_topology.py
@brief Algebraische Topologie – simpliziale und singuläre Homologie, Kohomologie,
       Homotopiegruppen, Faserbündel, Spektralsequenzen, K-Theorie und CW-Komplexe.
@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0

Dieses Modul implementiert die grundlegenden Konzepte der algebraischen Topologie
als Python-Klassen und Funktionen. Es dient sowohl als Lernmaterial als auch als
Berechnungswerkzeug.

Mathematische Grundbegriffe:
  - Simplizialkomplex: Kombination von Simplizes (Punkte, Kanten, Dreiecke, ...)
  - Homologiegruppe H_k: Maß für k-dimensionale "Löcher" in einem Raum
  - Betti-Zahlen β_k: Ränge der freien Teile der Homologiegruppen
  - Euler-Charakteristik: χ = Σ (-1)^k · β_k

Formeln (KaTeX):
  χ(X) = \\sum_{k=0}^{\\dim X} (-1)^k \\cdot \\beta_k
  H_k(X) = \\ker(\\partial_k) / \\operatorname{im}(\\partial_{k+1})
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple, Optional, Set
from itertools import combinations

import numpy as np

# ============================================================
# Hilfsfunktionen: Smith-Normalform und Homologie
# ============================================================

def compute_smith_normal_form(
    A: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    @brief Berechnet die Smith-Normalform einer ganzzahligen Matrix.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Smith-Normalform einer ganzzahligen Matrix A ist eine Diagonalmatrix D,
    sodass D = U · A · V für invertierbare ganzzahlige Matrizen U und V gilt,
    wobei d_1 | d_2 | ... | d_r (Teilbarkeitsbedingung).

    Mathematische Bedeutung:
        Jede ganzzahlige Matrix A über ℤ hat eine Smith-Normalform D = UAV,
        wobei D = diag(d_1, ..., d_r, 0, ..., 0) mit d_i | d_{i+1}.
        Dies ermöglicht die Berechnung der Torsionsuntergruppen in der Homologie:
        H_k ≅ ℤ^β ⊕ ℤ/d_1 ⊕ ... ⊕ ℤ/d_t

    Formel (KaTeX):
        D = U \\cdot A \\cdot V, \\quad U, V \\in GL_n(\\mathbb{Z})

    @param A Ganzzahlige Matrix als np.ndarray (dtype int).
    @return Tupel (D, U, V): D Smith-Normalform, U und V Transformationsmatrizen.

    Beispiel:
    >>> A = np.array([[2, 4], [1, 3]], dtype=int)
    >>> D, U, V = compute_smith_normal_form(A)
    """
    # Sicherheitshalber Integer-Typ erzwingen
    M = np.array(A, dtype=int).copy()
    m, n = M.shape

    # Initialisierung der Transformationsmatrizen als Einheitsmatrizen
    U = np.eye(m, dtype=int)  # Linkstransformation (Zeilenoperationen)
    V = np.eye(n, dtype=int)  # Rechtstransformation (Spaltenoperationen)

    pivot_row = 0  # Aktuelle Pivot-Zeile
    pivot_col = 0  # Aktuelle Pivot-Spalte

    while pivot_row < m and pivot_col < n:
        # Suche nach dem kleinsten Nicht-Null-Element im verbleibenden Unterblock
        sub = M[pivot_row:, pivot_col:]

        # Finde Position des betragsmäßig kleinsten Nicht-Null-Eintrags
        nonzero_mask = sub != 0
        if not nonzero_mask.any():
            break  # Kein Nicht-Null-Element mehr → fertig

        # Betragsmäßig kleinsten Eintrag als Pivot wählen
        abs_sub = np.where(nonzero_mask, np.abs(sub), np.iinfo(int).max)
        min_idx = np.unravel_index(abs_sub.argmin(), sub.shape)
        ri, ci = min_idx[0] + pivot_row, min_idx[1] + pivot_col

        # Pivot in die Diagonalposition bringen (Zeilen- und Spaltentausch)
        if ri != pivot_row:
            M[[pivot_row, ri]] = M[[ri, pivot_row]]
            U[[pivot_row, ri]] = U[[ri, pivot_row]]
        if ci != pivot_col:
            M[:, [pivot_col, ci]] = M[:, [ci, pivot_col]]
            V[:, [pivot_col, ci]] = V[:, [ci, pivot_col]]

        # Iteration: Eliminiere alle anderen Einträge in Zeile und Spalte
        progress = True
        while progress:
            progress = False

            # Spaltenelimination: Eliminiere M[pivot_row, j] für j > pivot_col
            for j in range(pivot_col + 1, n):
                if M[pivot_row, j] != 0:
                    if M[pivot_row, pivot_col] == 0:
                        # Pivot ist 0, tausche Spalten
                        M[:, [pivot_col, j]] = M[:, [j, pivot_col]]
                        V[:, [pivot_col, j]] = V[:, [j, pivot_col]]
                        progress = True
                        break
                    q = M[pivot_row, j] // M[pivot_row, pivot_col]
                    M[:, j] -= q * M[:, pivot_col]
                    V[:, j] -= q * V[:, pivot_col]
                    if M[pivot_row, j] != 0:
                        # Noch nicht eliminiert: erneuter Tausch nötig
                        M[:, [pivot_col, j]] = M[:, [j, pivot_col]]
                        V[:, [pivot_col, j]] = V[:, [j, pivot_col]]
                        progress = True

            # Zeilenelimination: Eliminiere M[i, pivot_col] für i > pivot_row
            for i in range(pivot_row + 1, m):
                if M[i, pivot_col] != 0:
                    if M[pivot_row, pivot_col] == 0:
                        M[[pivot_row, i]] = M[[i, pivot_row]]
                        U[[pivot_row, i]] = U[[i, pivot_row]]
                        progress = True
                        break
                    q = M[i, pivot_col] // M[pivot_row, pivot_col]
                    M[i] -= q * M[pivot_row]
                    U[i] -= q * U[pivot_row]
                    if M[i, pivot_col] != 0:
                        M[[pivot_row, i]] = M[[i, pivot_row]]
                        U[[pivot_row, i]] = U[[i, pivot_row]]
                        progress = True

        # Vorzeichen korrigieren: Diagonaleinträge sollen positiv sein
        if M[pivot_row, pivot_col] < 0:
            M[pivot_row] = -M[pivot_row]
            U[pivot_row] = -U[pivot_row]

        pivot_row += 1
        pivot_col += 1

    return M, U, V


def homology_from_boundary_matrices(
    boundaries: Dict[int, np.ndarray],
) -> Dict[int, Dict]:
    """
    @brief Berechnet Homologiegruppen aus gegebenen Randoperatormatrizen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Aus den Randoperatoren ∂_k werden die Homologiegruppen
        H_k = ker(∂_k) / im(∂_{k+1})
    berechnet. Dabei wird die Smith-Normalform verwendet, um
    Torsionsanteile zu ermitteln.

    Formel (KaTeX):
        H_k = \\ker(\\partial_k) / \\operatorname{im}(\\partial_{k+1})
        \\beta_k = \\operatorname{rank}(\\ker \\partial_k) - \\operatorname{rank}(\\operatorname{im} \\partial_{k+1})

    @param boundaries Dict mit Randmatrizen: boundaries[k] = Matrix für ∂_k.
    @return Dict mit Homologie-Informationen pro Grad k.
    """
    if not boundaries:
        return {}

    # Dimensionen aller Kettengruppen bestimmen
    all_dims: Dict[int, int] = {}
    for k, mat in boundaries.items():
        # ∂_k: C_k → C_{k-1}, also Spalten = dim(C_k)
        all_dims[k] = mat.shape[1]
        all_dims[k - 1] = mat.shape[0]

    # Sortierte Grade
    max_k = max(all_dims.keys()) if all_dims else 0
    min_k = min(all_dims.keys()) if all_dims else 0

    result: Dict[int, Dict] = {}

    for k in range(min_k, max_k + 1):
        dim_ck = all_dims.get(k, 0)

        # Rang des Bildes von ∂_{k+1}: im(∂_{k+1}) ⊂ C_k
        rank_image_next = 0
        torsion_coefficients = []
        if (k + 1) in boundaries:
            D, _, _ = compute_smith_normal_form(boundaries[k + 1])
            # Diagonaleinträge zählen → Rang und Torsion
            diag = [D[i, i] for i in range(min(D.shape)) if D[i, i] != 0]
            rank_image_next = len(diag)
            # Torsionskoeffizienten: d_i > 1 bedeutet Torsion
            torsion_coefficients = [d for d in diag if d > 1]

        # Rang des Kerns von ∂_k
        rank_boundary_k = 0
        if k in boundaries:
            D_k, _, _ = compute_smith_normal_form(boundaries[k])
            diag_k = [D_k[i, i] for i in range(min(D_k.shape)) if D_k[i, i] != 0]
            rank_boundary_k = len(diag_k)

        # Rang des Kerns: nullity(∂_k) = dim(C_k) - rank(∂_k)
        rank_kernel = dim_ck - rank_boundary_k

        # Betti-Zahl: β_k = rank(ker ∂_k) - rank(im ∂_{k+1})
        betti = max(0, rank_kernel - rank_image_next)

        # Torsion beschreiben
        torsion_str = ""
        if torsion_coefficients:
            parts = [f"ℤ/{d}" for d in torsion_coefficients]
            torsion_str = " ⊕ ".join(parts)

        result[k] = {
            "betti": betti,
            "rank": betti,
            "torsion": torsion_str,
            "dim_chain_group": dim_ck,
            "rank_kernel": rank_kernel,
            "rank_image": rank_image_next,
            "description": f"H_{k} ≅ ℤ^{betti}" + (
                f" ⊕ {torsion_str}" if torsion_str else ""
            ),
        }

    return result


def classify_surface(genus: int, orientable: bool) -> Dict:
    """
    @brief Klassifiziert kompakte Flächen nach Geschlecht und Orientierbarkeit.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Klassifikationssatz kompakter Flächen:
      - Orientierbar: Σ_g (zusammenhängende Summe von g Tori)
        χ = 2 - 2g, π_1 = ⟨a_1,b_1,...,a_g,b_g | [a_1,b_1]···[a_g,b_g]⟩
      - Nicht-orientierbar: N_k (k Kreuzkappenröhren)
        χ = 2 - k, π_1 = ⟨c_1,...,c_k | c_1²···c_k²⟩

    Formel (KaTeX):
        \\chi(\\Sigma_g) = 2 - 2g \\quad \\text{(orientierbar)}
        \\chi(N_k) = 2 - k \\quad \\text{(nicht-orientierbar)}

    @param genus Geschlecht der Fläche (g ≥ 0 orientierbar, k ≥ 1 nicht-orientierbar).
    @param orientable True für orientierbare, False für nicht-orientierbare Flächen.
    @return Dict mit Flächeninformationen.
    """
    if orientable:
        # Orientierbare Fläche mit Geschlecht g
        euler_char = 2 - 2 * genus
        if genus == 0:
            name = "S² (2-Sphäre)"
            homology = {0: "ℤ", 1: "0", 2: "ℤ"}
            pi1 = "trivial (1)"
        elif genus == 1:
            name = "T² (Torus)"
            homology = {0: "ℤ", 1: "ℤ²", 2: "ℤ"}
            pi1 = "ℤ × ℤ"
        else:
            name = f"Σ_{genus} (Geschlecht-{genus}-Fläche)"
            homology = {
                0: "ℤ",
                1: f"ℤ^{2*genus}",
                2: "ℤ",
            }
            gens = ", ".join(
                [f"a_{i}, b_{i}" for i in range(1, genus + 1)]
            )
            rel = "·".join([f"[a_{i},b_{i}]" for i in range(1, genus + 1)])
            pi1 = f"⟨{gens} | {rel}⟩"

        return {
            "name": name,
            "orientable": True,
            "genus": genus,
            "euler_characteristic": euler_char,
            "homology": homology,
            "pi1": pi1,
            "classification": "orientierbare kompakte Fläche",
        }
    else:
        # Nicht-orientierbare Fläche (k Kreuzkappenröhren)
        k = genus  # Hier ist genus die Anzahl der Kreuzkappenröhren
        euler_char = 2 - k
        if k == 1:
            name = "RP² (reell-projektive Ebene)"
            homology = {0: "ℤ", 1: "ℤ/2", 2: "0"}
            pi1 = "ℤ/2"
        elif k == 2:
            name = "Klein-Flasche"
            homology = {0: "ℤ", 1: "ℤ ⊕ ℤ/2", 2: "0"}
            pi1 = "⟨a, b | abab⁻¹⟩"
        else:
            name = f"N_{k} ({k} Kreuzkappenröhren)"
            homology = {
                0: "ℤ",
                1: f"ℤ^{k-1} ⊕ ℤ/2",
                2: "0",
            }
            gens = ", ".join([f"c_{i}" for i in range(1, k + 1)])
            rel = "·".join([f"c_{i}²" for i in range(1, k + 1)])
            pi1 = f"⟨{gens} | {rel}⟩"

        return {
            "name": name,
            "orientable": False,
            "genus": k,
            "euler_characteristic": euler_char,
            "homology": homology,
            "pi1": pi1,
            "classification": "nicht-orientierbare kompakte Fläche",
        }


def van_kampen_free_product(g1: str, g2: str) -> str:
    """
    @brief Berechnet das freie Produkt zweier Gruppen (Seifert-van-Kampen).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Seifert-van-Kampen-Satz besagt: Wenn X = A ∪ B mit A ∩ B wegzusammenhängend,
    dann π_1(X) = π_1(A) *_{π_1(A∩B)} π_1(B) (amalgamiertes Produkt).
    Für A ∩ B einfach zusammenhängend gilt: π_1(X) = π_1(A) * π_1(B).

    Formel (KaTeX):
        \\pi_1(A \\cup B) \\cong \\pi_1(A) *_{\\pi_1(A \\cap B)} \\pi_1(B)

    @param g1 Beschreibung der ersten Gruppe (z.B. "ℤ").
    @param g2 Beschreibung der zweiten Gruppe (z.B. "ℤ").
    @return String-Beschreibung des freien Produkts.
    """
    # Spezialfall: triviale Gruppen eliminieren
    if g1 in ("1", "trivial", "{1}"):
        return g2
    if g2 in ("1", "trivial", "{1}"):
        return g1
    return f"({g1}) * ({g2})"


def lyndon_hochschild_serre_demo() -> str:
    """
    @brief Demonstriert die Lyndon-Hochschild-Serre-Spektralsequenz.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Für eine Gruppenextension 1 → N → G → Q → 1 gibt es eine Spektralsequenz
    E²_{p,q} = H_p(Q; H_q(N; ℤ)) ⟹ H_{p+q}(G; ℤ).

    Formel (KaTeX):
        E^2_{p,q} = H_p(Q; H_q(N; \\mathbb{Z})) \\Rightarrow H_{p+q}(G; \\mathbb{Z})

    @return Beschreibungsstring der Spektralsequenz.
    """
    return (
        "Lyndon-Hochschild-Serre-Spektralsequenz:\n"
        "Für 1 → N → G → Q → 1 konvergiert\n"
        "E²_{p,q} = H_p(Q; H_q(N; ℤ)) ⟹ H_{p+q}(G; ℤ).\n"
        "Beispiel: G = ℤ, N = nℤ, Q = ℤ/n:\n"
        "  E²_{0,0} = ℤ/n, E²_{1,0} = 0, ...\n"
        "  Dies rekonstruiert H_*(ℤ) = (ℤ, 0, 0, ...)."
    )


def classifying_space_demo(group: str) -> str:
    """
    @brief Beschreibt den klassifizierenden Raum BG einer Gruppe G.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der klassifizierende Raum BG ist ein topologischer Raum mit
    π_1(BG) = G und π_k(BG) = 0 für k > 1 (Eilenberg-MacLane-Raum K(G,1)).

    Formel (KaTeX):
        \\pi_1(BG) = G, \\quad \\pi_k(BG) = 0 \\text{ für } k > 1

    @param group Name der Gruppe (z.B. "ℤ", "ℤ/2", "ℤ²").
    @return Beschreibungsstring des klassifizierenden Raums.
    """
    # Bekannte klassifizierende Räume
    known: Dict[str, str] = {
        "ℤ": "BZ = S¹ (Kreis), H_k(BZ) = ℤ für k=0,1, sonst 0",
        "ℤ/2": "B(ℤ/2) = RP^∞, H_k(B(ℤ/2)) = ℤ/2 für k ungerade > 0",
        "ℤ²": "B(ℤ²) = T² (Torus), H_0=ℤ, H_1=ℤ², H_2=ℤ",
        "ℤ/n": "B(ℤ/n) = L^∞(n) (Linsenraum), H_{2k-1}=ℤ/n für k≥1",
        "1": "B(1) = * (Punkt), H_0=ℤ, H_k=0 für k>0",
    }
    if group in known:
        return known[group]
    return (
        f"BG für G={group}: K(G,1)-Eilenberg-MacLane-Raum.\n"
        f"π_1(BG)={group}, π_k(BG)=0 für k>1.\n"
        f"H_*(BG) hängt von der Gruppenstruktur von G ab."
    )


# ============================================================
# Klasse: SimplicialComplex
# ============================================================

class SimplicialComplex:
    """
    @brief Repräsentiert einen Simplizialkomplex.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein Simplizialkomplex K ist eine Menge von Simplizes (Ecken, Kanten,
    Dreiecken, Tetraedern, ...), abgeschlossen unter Seiten-Operationen.

    Formeln (KaTeX):
        \\chi(K) = \\sum_{k=0}^{\\dim K} (-1)^k \\cdot |K_k|
        \\partial_p [v_0, ..., v_p] = \\sum_{i=0}^p (-1)^i [v_0, ..., \\hat{v}_i, ..., v_p]
    """

    def __init__(self, simplices: List[tuple]) -> None:
        """
        @brief Initialisiert den Simplizialkomplex mit einer Liste von Simplizes.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param simplices Liste von Tupeln, z.B. [(0,), (1,), (0,1), (0,1,2)].
                         Tupel der Länge k+1 repräsentieren k-Simplizes.
        """
        # Alle gegebenen Simplizes und ihre Seiten sammeln
        self._simplices: Set[tuple] = set()
        for s in simplices:
            self._add_with_faces(tuple(sorted(s)))

    def _add_with_faces(self, simplex: tuple) -> None:
        """
        @brief Fügt einen Simplex und alle seine Seiten hinzu (Abschlusseigenschaft).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param simplex Tupel der Ecken des Simplex.
        """
        self._simplices.add(simplex)
        # Alle echten Teilmengen als Seiten hinzufügen
        for k in range(1, len(simplex)):
            for face in combinations(simplex, k):
                self._simplices.add(tuple(sorted(face)))

    def faces(self) -> List[tuple]:
        """
        @brief Gibt alle Simplizes (Seiten) des Komplexes zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Sortierte Liste aller Simplizes.
        """
        return sorted(self._simplices, key=lambda s: (len(s), s))

    def _simplices_of_dim(self, k: int) -> List[tuple]:
        """
        @brief Gibt alle k-Simplizes zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param k Dimension (0=Punkte, 1=Kanten, 2=Dreiecke, ...).
        @return Liste der k-Simplizes in lexikographischer Ordnung.
        """
        return sorted(
            [s for s in self._simplices if len(s) == k + 1],
            key=lambda s: s,
        )

    def euler_characteristic(self) -> int:
        """
        @brief Berechnet die Euler-Charakteristik χ = V - E + F - ...
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            \\chi = \\sum_{k=0}^{\\dim K} (-1)^k \\cdot f_k

        wobei f_k die Anzahl der k-Simplizes ist.

        @return Euler-Charakteristik als ganze Zahl.
        """
        chi = 0
        k = 0
        while True:
            sims = self._simplices_of_dim(k)
            if not sims:
                break
            chi += ((-1) ** k) * len(sims)
            k += 1
        return chi

    def boundary_matrix(self, p: int) -> np.ndarray:
        """
        @brief Berechnet den Randoperator ∂_p als Matrix.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Der Randoperator ∂_p: C_p → C_{p-1} ist definiert durch:
            ∂_p [v_0, ..., v_p] = Σ_{i=0}^p (-1)^i [v_0, ..., v̂_i, ..., v_p]

        Formel (KaTeX):
            (\\partial_p)_{\\sigma, \\tau} = \\begin{cases}
              (-1)^i & \\text{wenn } \\tau = [v_0,...,\\hat{v}_i,...,v_p] \\text{ Seite von } \\sigma \\\\
              0 & \\text{sonst}
            \\end{cases}

        @param p Dimension der Simplizes (∂_p: C_p → C_{p-1}).
        @return Randmatrix als np.ndarray (dtype int), Form: |C_{p-1}| × |C_p|.
        """
        # p-Simplizes (Spalten) und (p-1)-Simplizes (Zeilen)
        p_simplices = self._simplices_of_dim(p)
        pm1_simplices = self._simplices_of_dim(p - 1)

        if not p_simplices or not pm1_simplices:
            return np.zeros((max(1, len(pm1_simplices)), max(1, len(p_simplices))), dtype=int)

        # Index-Dictionaries für schnellen Zugriff
        row_idx = {s: i for i, s in enumerate(pm1_simplices)}

        # Randmatrix aufbauen
        mat = np.zeros((len(pm1_simplices), len(p_simplices)), dtype=int)

        for col, sigma in enumerate(p_simplices):
            # Alle Seiten von sigma berechnen
            for i, v in enumerate(sigma):
                # Seite: sigma ohne die i-te Ecke
                face = tuple(sorted(sigma[:i] + sigma[i + 1:]))
                sign = (-1) ** i
                if face in row_idx:
                    mat[row_idx[face], col] = sign

        return mat

    def betti_numbers(self) -> Dict[int, int]:
        """
        @brief Berechnet die Betti-Zahlen β_k = dim H_k.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Betti-Zahlen messen die Anzahl der k-dimensionalen "Löcher":
          β_0 = Anzahl Zusammenhangskomponenten
          β_1 = Anzahl unabhängiger Schleifen
          β_2 = Anzahl 2D-Hohlräume
          ...

        Formel (KaTeX):
            \\beta_k = \\operatorname{rank}(H_k) = \\operatorname{rank}(\\ker \\partial_k)
                       - \\operatorname{rank}(\\operatorname{im} \\partial_{k+1})

        @return Dict mit β_k-Werten.
        """
        max_dim = 0
        while self._simplices_of_dim(max_dim):
            max_dim += 1
        max_dim -= 1

        if max_dim < 0:
            return {}

        # Randmatrizen berechnen
        boundaries = {}
        for p in range(1, max_dim + 1):
            boundaries[p] = self.boundary_matrix(p)

        # Betti-Zahlen über Smith-Normalform
        result = {}
        for k in range(max_dim + 1):
            dim_ck = len(self._simplices_of_dim(k))

            # Rang von ∂_k (Bild von ∂_k in C_{k-1})
            rank_dk = 0
            if k in boundaries and boundaries[k].size > 0:
                D, _, _ = compute_smith_normal_form(boundaries[k])
                rank_dk = int(np.sum(np.diag(D[:min(D.shape)]) != 0))

            # Rang von ∂_{k+1} (Bild von ∂_{k+1} in C_k)
            rank_dk1 = 0
            if (k + 1) in boundaries and boundaries[k + 1].size > 0:
                D1, _, _ = compute_smith_normal_form(boundaries[k + 1])
                rank_dk1 = int(np.sum(np.diag(D1[:min(D1.shape)]) != 0))

            # β_k = dim(ker ∂_k) - dim(im ∂_{k+1})
            ker_dim = dim_ck - rank_dk
            betti = max(0, ker_dim - rank_dk1)
            result[k] = betti

        return result


# ============================================================
# Klasse: SimplicialHomology
# ============================================================

class SimplicialHomology:
    """
    @brief Berechnet die simpliziale Homologie eines Simplizialkomplexes.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die simpliziale Homologie nutzt die Kettengruppen C_k und Randoperatoren
    ∂_k, um die Homologiegruppen H_k = ker(∂_k) / im(∂_{k+1}) zu berechnen.
    """

    def __init__(self, complex: SimplicialComplex) -> None:
        """
        @brief Initialisiert SimplicialHomology mit einem Simplizialkomplex.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param complex Ein SimplicialComplex-Objekt.
        """
        self._complex = complex

    def chain_groups(self) -> Dict[int, int]:
        """
        @brief Gibt die Dimensionen der Kettengruppen C_k zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        C_k ist der freie abelsche Gruppe, erzeugt von den k-Simplizes.
        dim(C_k) = Anzahl der k-Simplizes.

        @return Dict: k → dim(C_k).
        """
        result = {}
        k = 0
        while True:
            sims = self._complex._simplices_of_dim(k)
            if not sims:
                break
            result[k] = len(sims)
            k += 1
        return result

    def boundary_operators(self) -> Dict[int, np.ndarray]:
        """
        @brief Gibt alle Randoperatoren als Matrizen zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Dict: k → ∂_k (Matrix).
        """
        dims = self.chain_groups()
        if not dims:
            return {}
        max_k = max(dims.keys())
        return {k: self._complex.boundary_matrix(k) for k in range(1, max_k + 1)}

    def homology_groups(self) -> Dict[int, Dict]:
        """
        @brief Berechnet alle Homologiegruppen H_k mit Rang und Torsion.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Verwendet Smith-Normalform für die vollständige Strukturbestimmung
        inklusive Torsionsanteile.

        @return Dict: k → {"rank": β_k, "torsion": str, "description": str}.
        """
        dims = self.chain_groups()
        if not dims:
            return {}
        max_k = max(dims.keys())

        # Alle Randmatrizen berechnen
        boundaries = {}
        for k in range(1, max_k + 1):
            mat = self._complex.boundary_matrix(k)
            if mat.size > 0:
                boundaries[k] = mat

        result = {}
        for k in range(max_k + 1):
            dim_ck = dims.get(k, 0)

            # Rang des Kerns von ∂_k
            rank_dk = 0
            if k in boundaries:
                D, _, _ = compute_smith_normal_form(boundaries[k])
                diag_vals = [D[i, i] for i in range(min(D.shape))]
                rank_dk = sum(1 for d in diag_vals if d != 0)

            # Rang und Torsion des Bildes von ∂_{k+1}
            rank_dk1 = 0
            torsion = []
            if (k + 1) in boundaries:
                D1, _, _ = compute_smith_normal_form(boundaries[k + 1])
                diag_vals1 = [D1[i, i] for i in range(min(D1.shape))]
                nonzero = [d for d in diag_vals1 if d != 0]
                rank_dk1 = len(nonzero)
                torsion = [abs(d) for d in nonzero if abs(d) > 1]

            ker_dim = dim_ck - rank_dk
            betti = max(0, ker_dim - rank_dk1)

            torsion_str = " ⊕ ".join(f"ℤ/{d}" for d in torsion) if torsion else ""
            desc = f"H_{k} ≅ ℤ^{betti}"
            if torsion_str:
                desc += f" ⊕ {torsion_str}"

            result[k] = {
                "rank": betti,
                "betti": betti,
                "torsion": torsion_str,
                "description": desc,
            }

        return result

    def betti_numbers(self) -> Dict[int, int]:
        """
        @brief Gibt die Betti-Zahlen als Dict zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Dict: k → β_k.
        """
        return self._complex.betti_numbers()


# ============================================================
# Klasse: SingularHomology (abstrakte Darstellung)
# ============================================================

class SingularHomology:
    """
    @brief Abstrakte Darstellung der singulären Homologie bekannter Räume.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die singuläre Homologie ist der natürliche Begriff der Homologie für
    allgemeine topologische Räume. Für bekannte Räume sind die Gruppen
    explizit tabellarisch dargestellt.
    """

    def __init__(self, space_name: str = "") -> None:
        """
        @brief Initialisiert SingularHomology.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param space_name Name des topologischen Raums (optional).
        """
        self._space_name = space_name

    def homology_of_sphere(self, n: int) -> Dict[int, str]:
        """
        @brief Gibt die singulären Homologiegruppen der n-Sphäre S^n zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H_k(S^n; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0 \\text{ oder } k = n \\\\
              0 & \\text{sonst}
            \\end{cases}

        @param n Dimension der Sphäre (n ≥ 0).
        @return Dict: k → Gruppen-Beschreibung.
        """
        if n < 0:
            raise ValueError(f"Sphärendimension muss ≥ 0 sein, erhalten: {n}")

        result = {}
        # S^0 = zwei Punkte: H_0 = ℤ², H_k = 0 für k > 0
        if n == 0:
            result[0] = "ℤ ⊕ ℤ"
            return result

        for k in range(n + 1):
            if k == 0:
                result[k] = "ℤ"
            elif k == n:
                result[k] = "ℤ"
            else:
                result[k] = "0"
        return result

    def homology_of_torus(self) -> Dict[int, str]:
        """
        @brief Gibt die singulären Homologiegruppen des 2-Torus T² zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H_k(T^2; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0, 2 \\\\
              \\mathbb{Z}^2 & k = 1 \\\\
              0 & k > 2
            \\end{cases}

        @return Dict: k → Gruppen-Beschreibung.
        """
        return {
            0: "ℤ",
            1: "ℤ ⊕ ℤ",
            2: "ℤ",
        }

    def homology_of_rp(self, n: int) -> Dict[int, str]:
        """
        @brief Gibt die singulären Homologiegruppen des reell-projektiven Raums RP^n zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H_k(\\mathbb{RP}^n; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0 \\\\
              \\mathbb{Z}/2 & k \\text{ ungerade}, 0 < k < n \\\\
              \\mathbb{Z} & k = n, n \\text{ ungerade} \\\\
              0 & \\text{sonst}
            \\end{cases}

        @param n Dimension des projektiven Raums (n ≥ 1).
        @return Dict: k → Gruppen-Beschreibung.
        """
        if n < 1:
            raise ValueError(f"RP^n erfordert n ≥ 1, erhalten: {n}")

        result: Dict[int, str] = {}
        for k in range(n + 1):
            if k == 0:
                result[k] = "ℤ"
            elif k == n and n % 2 == 1:
                # RP^n (n ungerade) ist orientierbar: H_n = ℤ
                result[k] = "ℤ"
            elif k % 2 == 1 and 0 < k < n:
                # Ungerade Dimension: ℤ/2 Torsion
                result[k] = "ℤ/2"
            else:
                result[k] = "0"
        return result

    def homology_of_cp(self, n: int) -> Dict[int, str]:
        """
        @brief Gibt die singulären Homologiegruppen des komplex-projektiven Raums CP^n zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H_k(\\mathbb{CP}^n; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0, 2, 4, ..., 2n \\\\
              0 & \\text{sonst}
            \\end{cases}

        @param n Komplex-Dimension (n ≥ 1).
        @return Dict: k → Gruppen-Beschreibung.
        """
        if n < 1:
            raise ValueError(f"CP^n erfordert n ≥ 1, erhalten: {n}")

        result: Dict[int, str] = {}
        for k in range(2 * n + 1):
            if k % 2 == 0:
                result[k] = "ℤ"
            else:
                result[k] = "0"
        return result

    def mayer_vietoris_demo(self) -> str:
        """
        @brief Demonstriert die Mayer-Vietoris-Sequenz am Beispiel S^2 = D² ∪ D².
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Mayer-Vietoris-Sequenz ist eine lange exakte Sequenz:
        ... → H_n(A∩B) → H_n(A) ⊕ H_n(B) → H_n(X) → H_{n-1}(A∩B) → ...

        Formel (KaTeX):
            \\cdots \\to H_n(A \\cap B) \\to H_n(A) \\oplus H_n(B)
            \\to H_n(X) \\xrightarrow{\\partial} H_{n-1}(A \\cap B) \\to \\cdots

        @return Beschreibungsstring der Sequenz.
        """
        return (
            "Mayer-Vietoris für S² = D²_+ ∪ D²_-, A∩B ≅ S¹:\n"
            "  ... → H_2(S¹) → H_2(D²_+)⊕H_2(D²_-) → H_2(S²) → H_1(S¹) → H_1(D²_+)⊕H_1(D²_-) → ...\n"
            "  ... → 0 → 0⊕0 → H_2(S²) → ℤ → 0⊕0 → H_1(S²) → ...\n"
            "Exaktheit erzwingt: H_2(S²) ≅ ℤ, H_1(S²) = 0.\n"
            "Allgemein: H_n(S^k) = ℤ für n=0,k; = 0 sonst."
        )


# ============================================================
# Klasse: CohomologyRing
# ============================================================

class CohomologyRing:
    """
    @brief Beschreibt den Kohomologiering eines topologischen Raums.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Kohomologiering H*(X; ℤ) trägt eine Ring-Struktur durch das Cup-Produkt:
        ∪: H^p(X) × H^q(X) → H^{p+q}(X)
    """

    def __init__(self, space_name: str = "") -> None:
        """
        @brief Initialisiert CohomologyRing.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param space_name Name des Raums.
        """
        self._space_name = space_name

    def cohomology_of_sphere(self, n: int) -> Dict[int, str]:
        """
        @brief Gibt die Kohomologiegruppen der n-Sphäre zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H^k(S^n; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0 \\text{ oder } k = n \\\\
              0 & \\text{sonst}
            \\end{cases}

        @param n Dimension der Sphäre.
        @return Dict: k → Kohomologie-Beschreibung.
        """
        if n < 1:
            raise ValueError(f"Sphärendimension muss ≥ 1 sein, erhalten: {n}")

        result: Dict[int, str] = {}
        for k in range(n + 1):
            if k == 0 or k == n:
                result[k] = "ℤ"
            else:
                result[k] = "0"
        return result

    def cohomology_of_torus(self) -> Dict[int, str]:
        """
        @brief Gibt die Kohomologiegruppen des 2-Torus T² zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Formel (KaTeX):
            H^k(T^2; \\mathbb{Z}) = \\begin{cases}
              \\mathbb{Z} & k = 0, 2 \\\\
              \\mathbb{Z}^2 & k = 1 \\\\
              0 & k > 2
            \\end{cases}

        @return Dict: k → Kohomologie-Beschreibung.
        """
        return {
            0: "ℤ",
            1: "ℤ ⊕ ℤ",
            2: "ℤ",
        }

    def cup_product_demo(self, space: str) -> str:
        """
        @brief Demonstriert das Cup-Produkt für einen gegebenen Raum.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Das Cup-Produkt macht H*(X) zu einem graduierten Ring.
        Für S^n gilt: das Cup-Produkt ist trivial (da H^k = 0 für 0 < k < n).
        Für T² = S¹ × S¹: α ∪ β = [T²] (Fundamentalklasse).

        Formel (KaTeX):
            \\cup: H^p(X) \\times H^q(X) \\to H^{p+q}(X)

        @param space Raumname ("sphere", "torus", "CP2", ...).
        @return Beschreibungsstring des Cup-Produkts.
        """
        if space.lower() in ("sphere", "s^n", "sn"):
            return (
                "Cup-Produkt auf S^n:\n"
                "H*(S^n) = ℤ[x]/(x²), deg(x) = n.\n"
                "x ∪ x = 0 für n > 0 (kein nicht-triviales Produkt)."
            )
        elif space.lower() in ("torus", "t²", "t2"):
            return (
                "Cup-Produkt auf T² = S¹ × S¹:\n"
                "H*(T²) = ℤ[α,β]/(α²,β²), deg(α)=deg(β)=1.\n"
                "α ∪ β = [T²] ∈ H²(T²) ≅ ℤ (Fundamentalklasse).\n"
                "β ∪ α = -α ∪ β (Antikommutativität für ungerade Grade)."
            )
        elif space.lower() in ("cp2", "cp^2"):
            return (
                "Cup-Produkt auf CP²:\n"
                "H*(CP²) = ℤ[x]/(x³), deg(x) = 2.\n"
                "x ∪ x = x² ∈ H⁴(CP²) ≅ ℤ (nicht-trivial!).\n"
                "Dies unterscheidet CP² von S² ∨ S⁴."
            )
        return f"Cup-Produkt für '{space}' nicht explizit implementiert."

    def poincare_duality_demo(self, space: str, dim: int) -> str:
        """
        @brief Demonstriert die Poincaré-Dualität für eine orientierbare Mannigfaltigkeit.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Poincaré-Dualität: Für eine kompakte, orientierbare n-Mannigfaltigkeit M:
            H^k(M; ℤ) ≅ H_{n-k}(M; ℤ) (grob gesprochen)

        Formel (KaTeX):
            H^k(M; \\mathbb{Z}) \\cong H_{n-k}(M; \\mathbb{Z})

        @param space Name der Mannigfaltigkeit.
        @param dim Dimension der Mannigfaltigkeit.
        @return Beschreibungsstring der Poincaré-Dualität.
        """
        return (
            f"Poincaré-Dualität für {space} (dim={dim}):\n"
            f"H^k({space}) ≅ H_{{dim-k}}({space}) für k=0,...,{dim}.\n"
            f"Beispiel: H^0 ≅ H_{dim} (Fundamentalklasse), H^{dim} ≅ H_0 ≅ ℤ.\n"
            f"Voraussetzung: {space} ist kompakt und orientierbar."
        )


# ============================================================
# Klasse: HomotopyGroups
# ============================================================

class HomotopyGroups:
    """
    @brief Beschreibt Homotopiegruppen bekannter topologischer Räume.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Homotopiegruppen π_k(X) klassifizieren k-dimensionale Schleifen:
      π_1(X) = Fundamentalgruppe (Weg-Homotopie-Klassen)
      π_k(X) = Homotopie-Klassen stetiger Abbildungen S^k → X
    """

    # Bekannte Homotopiegruppen von Sphären (Tabelle)
    # π_k(S^n): Zeile = k, Spalte = n
    _SPHERE_HOMOTOPY: Dict[Tuple[int, int], str] = {
        (0, 1): "0", (1, 1): "ℤ", (2, 1): "0", (3, 1): "0",
        (0, 2): "0", (1, 2): "0", (2, 2): "ℤ", (3, 2): "ℤ",
        (4, 2): "ℤ/2", (5, 2): "ℤ/2", (6, 2): "ℤ/12",
        (0, 3): "0", (1, 3): "0", (2, 3): "0", (3, 3): "ℤ",
        (4, 3): "ℤ/2", (5, 3): "ℤ/2", (6, 3): "ℤ/12",
        (0, 4): "0", (1, 4): "0", (2, 4): "0", (3, 4): "0",
        (4, 4): "ℤ", (5, 4): "ℤ/2", (6, 4): "ℤ/2",
    }

    def __init__(self) -> None:
        """
        @brief Initialisiert HomotopyGroups.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        pass

    def fundamental_group(self, space: str) -> str:
        """
        @brief Gibt die Fundamentalgruppe π_1 eines bekannten Raums zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bekannte Fundamentalgruppen:
          π_1(S¹) = ℤ  (Umlaufzahl)
          π_1(T²) = ℤ × ℤ  (zwei unabhängige Schleifen)
          π_1(S^n) = 1 für n ≥ 2  (einfach zusammenhängend)
          π_1(RP^n) = ℤ/2 für n ≥ 2

        @param space Raumname (z.B. "S1", "T2", "S2", "RP2", "RP3").
        @return String mit der Fundamentalgruppe.
        """
        known: Dict[str, str] = {
            "S1": "ℤ",
            "S^1": "ℤ",
            "circle": "ℤ",
            "T2": "ℤ × ℤ",
            "T^2": "ℤ × ℤ",
            "torus": "ℤ × ℤ",
            "S2": "1 (trivial)",
            "S^2": "1 (trivial)",
            "S3": "1 (trivial)",
            "S^3": "1 (trivial)",
            "S4": "1 (trivial)",
            "S^4": "1 (trivial)",
            "RP2": "ℤ/2",
            "RP^2": "ℤ/2",
            "RP3": "ℤ/2",
            "RP^3": "ℤ/2",
            "Klein": "⟨a,b | abab⁻¹⟩",
            "Moebius": "ℤ",
            "point": "1 (trivial)",
            "disk": "1 (trivial)",
        }
        return known.get(space, f"π_1({space}) nicht in Tabelle")

    def higher_homotopy_sphere(self, n: int, k: int) -> str:
        """
        @brief Gibt bekannte höhere Homotopiegruppen π_k(S^n) zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Allgemeine Regeln:
          - π_k(S^n) = 0 für k < n (Hurewicz: trivialer für k < n)
          - π_n(S^n) = ℤ für alle n ≥ 1
          - Für k > n: komplizierte Torsionsgruppen (Freudenthal-Suspension)

        Formel (KaTeX):
            \\pi_k(S^n) = 0 \\text{ für } k < n \\quad (\\text{Hurewicz-Theorem})
            \\pi_n(S^n) = \\mathbb{Z} \\quad \\text{(Abbildungsgrad)}

        @param n Dimension der Sphäre.
        @param k Dimension der Homotopiegruppe.
        @return String mit der Gruppe.
        """
        if k < 0 or n < 1:
            return "0"
        if k == 0:
            return "0" if n >= 1 else "ℤ/2"
        if k < n:
            # Hurewicz: alle π_k = 0 für k < n
            return "0"
        if k == n:
            # π_n(S^n) = ℤ immer
            return "ℤ"

        # Tabelle für bekannte Werte
        if (k, n) in self._SPHERE_HOMOTOPY:
            return self._SPHERE_HOMOTOPY[(k, n)]

        # Stabile Homotopiegruppen (n ≥ k+2): Freudenthal-Periodizität
        if n >= k + 2:
            # Stabile Homotopiegruppe π_{k-n}^s
            diff = k - n
            stable: Dict[int, str] = {
                1: "ℤ/2",
                2: "ℤ/2",
                3: "ℤ/24",
                4: "0",
                5: "0",
                6: "ℤ/2",
                7: "ℤ/240",
                8: "ℤ/2 ⊕ ℤ/2",
            }
            if diff in stable:
                return stable[diff]

        return f"π_{k}(S^{n}) unbekannt (komplex)"

    def seifert_van_kampen_demo(self) -> str:
        """
        @brief Demonstriert den Seifert-van-Kampen-Satz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Seifert-van-Kampen: Wenn X = A ∪ B und A, B, A∩B wegzusammenhängend:
            π_1(X) ≅ π_1(A) *_{π_1(A∩B)} π_1(B)

        Formel (KaTeX):
            \\pi_1(X) \\cong \\pi_1(A) *_{\\pi_1(A \\cap B)} \\pi_1(B)

        @return Beschreibungsstring.
        """
        return (
            "Seifert-van-Kampen-Satz:\n"
            "Für X = A ∪ B (A,B,A∩B wegzusammenhängend):\n"
            "  π_1(X) = π_1(A) *_{π_1(A∩B)} π_1(B)\n"
            "\nBeispiel: Torus T² = (S¹×D²) ∪ (D²×S¹), A∩B ≅ S¹×S¹:\n"
            "  π_1(A) = ℤ, π_1(B) = ℤ, π_1(A∩B) = ℤ×ℤ\n"
            "  π_1(T²) = ℤ×ℤ (aus der Amalgamierung)\n"
            "\nBeispiel: Figur-8 (Wedge S¹∨S¹):\n"
            "  π_1(S¹∨S¹) = ℤ * ℤ (freies Produkt)"
        )

    def covering_space_demo(self) -> str:
        """
        @brief Demonstriert die Überlagerungsraumtheorie.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Überlagerungsräume: p: X̃ → X mit p_*π_1(X̃) ⊂ π_1(X).
        Klassifikation: Isomorphieklassen von Überlagerungen ↔
        Untergruppen von π_1(X).

        @return Beschreibungsstring.
        """
        return (
            "Überlagerungsraumtheorie:\n"
            "• p: ℝ → S¹, t ↦ e^{2πit}: universelle Überlagerung von S¹\n"
            "  π_1(ℝ) = 1, Decktransformationsgruppe = ℤ = π_1(S¹)\n"
            "• p: S^n → RP^n (n≥1): zweifache Überlagerung\n"
            "  Decktransformationen: ℤ/2 = π_1(RP^n) für n≥2\n"
            "• Klassifikationssatz: Überlagerungen von X ↔ Untergruppen von π_1(X)\n"
            "  (für lokal einfach-zusammenhängende X)"
        )

    def long_exact_sequence_demo(self) -> str:
        """
        @brief Demonstriert die lange exakte Homotopiesequenz eines Faserraums.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Für ein Faserbündel F → E → B gilt die lange exakte Sequenz:
        ... → π_n(F) → π_n(E) → π_n(B) → π_{n-1}(F) → ... → π_0(B)

        Formel (KaTeX):
            \\cdots \\to \\pi_n(F) \\to \\pi_n(E) \\to \\pi_n(B)
            \\xrightarrow{\\partial} \\pi_{n-1}(F) \\to \\cdots

        @return Beschreibungsstring.
        """
        return (
            "Lange exakte Homotopiesequenz für F → E → B:\n"
            "  ... → π_n(F) → π_n(E) → π_n(B) → π_{n-1}(F) → ...\n"
            "\nBeispiel: Hopf-Faserung S¹ → S³ → S²:\n"
            "  ... → π_n(S¹) → π_n(S³) → π_n(S²) → π_{n-1}(S¹) → ...\n"
            "  Für n=2: π_2(S¹)=0 → π_2(S³)=0 → π_2(S²)=ℤ → π_1(S¹)=ℤ → π_1(S³)=0\n"
            "  Exaktheit: π_2(S²) ≅ ℤ (Abbildungsgrad), π_3(S²) ≅ ℤ (Hopf-Invariante)"
        )


# ============================================================
# Klasse: FiberBundle
# ============================================================

class FiberBundle:
    """
    @brief Beschreibt Faserbündel und ihre topologischen Invarianten.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein Faserbündel F → E → B besteht aus:
      - Basisraum B
      - Totalraum E
      - Faser F
      - Projektion p: E → B mit lokalem Produktcharakter
    """

    def __init__(self, base: str, fiber: str, total: str) -> None:
        """
        @brief Initialisiert FiberBundle.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param base Name des Basisraums (z.B. "S^2").
        @param fiber Name der Faser (z.B. "S^1").
        @param total Name des Totalraums (z.B. "S^3").
        """
        self.base = base
        self.fiber = fiber
        self.total = total

    def hopf_fibration_demo(self) -> Dict:
        """
        @brief Demonstriert die Hopf-Faserung S¹ → S³ → S².
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Hopf-Faserung ist das historisch erste nicht-triviale Faserbündel:
          h: S³ → S² mit Faser S¹
          [z₁:z₂] ↦ (|z₁|²-|z₂|², 2·Re(z₁z̄₂), 2·Im(z₁z̄₂))

        Formel (KaTeX):
            h: S^3 \\subset \\mathbb{C}^2 \\to \\mathbb{CP}^1 \\cong S^2
            (z_1, z_2) \\mapsto [z_1 : z_2]

        @return Dict mit Informationen zur Hopf-Faserung.
        """
        return {
            "name": "Hopf-Faserung",
            "fiber": "S¹",
            "total_space": "S³",
            "base_space": "S²",
            "map": "h(z₁,z₂) = [z₁:z₂] ∈ CP¹ ≅ S²",
            "hopf_invariant": 1,
            "pi3_s2": "ℤ (erzeugt von der Hopf-Faserung)",
            "characteristic_class": "c₁(γ) = 1 ∈ H²(S²;ℤ)",
            "description": (
                "Die Hopf-Faserung zeigt: π₃(S²) = ℤ ≠ 0.\n"
                "Jeder Kreis S¹ in S³ projiziert auf einen Punkt in S².\n"
                "Je zwei verschiedene Fasern verlinken sich einmal (Hopf-Link)."
            ),
        }

    def vector_bundle_demo(self) -> str:
        """
        @brief Demonstriert Vektorbündel und ihre Eigenschaften.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Beschreibungsstring.
        """
        return (
            "Vektorbündel:\n"
            "• Tangentialbündel TM → M: Faser = ℝⁿ, dim(M) = n\n"
            "  TM trivial ⟺ M parallelisierbar (z.B. S¹, S³, S⁷)\n"
            "• Tautologisches Bündel γ¹ → RP^n: Moebius-Band für n=1\n"
            "• Kanonisches Linienbündel L → CP^n: Hopf-Bündel für n=1\n"
            "• Stiefel-Whitney-Klassen w_i ∈ H^i(M;ℤ/2): Obstruktionen\n"
            f"  w_1(M) = 0 ⟺ M orientierbar\n"
            f"Aktuelles Bündel: {self.fiber} → {self.total} → {self.base}"
        )

    def characteristic_classes_demo(self) -> str:
        """
        @brief Demonstriert charakteristische Klassen von Vektorbündeln.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Charakteristische Klassen sind kohomologische Invarianten von Bündeln:
          - Chern-Klassen c_k ∈ H^{2k}(B; ℤ) für komplexe Bündel
          - Pontryagin-Klassen p_k ∈ H^{4k}(B; ℤ) für reelle Bündel
          - Stiefel-Whitney-Klassen w_k ∈ H^k(B; ℤ/2)

        @return Beschreibungsstring.
        """
        return (
            "Charakteristische Klassen:\n"
            "Chern-Klassen (komplexe Bündel E → B):\n"
            "  c(E) = 1 + c₁(E) + c₂(E) + ... ∈ H*(B; ℤ)\n"
            "  c₁(TCP^n) = (n+1)·h, h = c₁(L*) Hyperebenenklasse\n"
            "\nPontryagin-Klassen (reelle Bündel E → B):\n"
            "  p_k(E) = (-1)^k c_{2k}(E ⊗ ℂ) ∈ H^{4k}(B; ℤ)\n"
            "  Hirzebruch-Signaturformel: σ(M) = L(p₁,...,pₙ)[M]\n"
            "\nEuler-Klasse (orientiertes Bündel Rang r):\n"
            "  e(E) ∈ H^r(B; ℤ), e(TM)[M] = χ(M) (Gauß-Bonnet)"
        )

    def euler_class_demo(self) -> str:
        """
        @brief Demonstriert die Euler-Klasse.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Euler-Klasse e(E) ∈ H^r(B; ℤ) eines orientierten Bündels vom Rang r:
          Auswertung auf der Fundamentalklasse: e(TM)[M] = χ(M)

        Formel (KaTeX):
            e(TM)[M] = \\chi(M) = \\sum_k (-1)^k \\beta_k

        @return Beschreibungsstring.
        """
        return (
            "Euler-Klasse e(E) ∈ H^r(B; ℤ) für Rang-r-Bündel E → B:\n"
            "• e(TM)[M] = χ(M) = Euler-Charakteristik (Gauß-Bonnet-Chern)\n"
            "• e(TM) = 0 ⟺ M hat ein nirgendwo verschwindendes Vektorfeld\n"
            "• Sphären S^n: e(TS^n) = χ(S^n) = 1+(-1)^n\n"
            "  S² (χ=2): 'Haarsatz' – kein nirgendwo-verschwindenes Vektorfeld\n"
            "  S³ (χ=0): parallelisierbar"
        )


# ============================================================
# Klasse: SpectralSequence
# ============================================================

class SpectralSequence:
    """
    @brief Beschreibt Spektralsequenzen in der algebraischen Topologie.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Spektralsequenzen sind algebraische Werkzeuge zur Berechnung von
    Homologie-/Kohomologiegruppen durch iterative Approximationen.
    """

    def __init__(self, name: str) -> None:
        """
        @brief Initialisiert SpectralSequence.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param name Name der Spektralsequenz (z.B. "Serre", "Leray").
        """
        self._name = name

    def serre_spectral_sequence_demo(self) -> str:
        """
        @brief Demonstriert die Serre-Spektralsequenz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Für ein Faserbündel F → E → B mit π_1(B)-Wirkung auf H_*(F):
          E²_{p,q} = H_p(B; H_q(F; ℤ)) ⟹ H_{p+q}(E; ℤ)

        Formel (KaTeX):
            E^2_{p,q} = H_p(B; H_q(F; \\mathbb{Z})) \\Rightarrow H_{p+q}(E; \\mathbb{Z})

        @return Beschreibungsstring.
        """
        return (
            "Serre-Spektralsequenz:\n"
            "Für F → E → B (B einfach-zusammenhängend):\n"
            "  E²_{p,q} = H_p(B; ℤ) ⊗ H_q(F; ℤ) ⟹ H_{p+q}(E; ℤ)\n"
            "\nAnwendung: Hopf-Faserung S¹ → S³ → S²:\n"
            "  E²_{p,q} = H_p(S²; ℤ) ⊗ H_q(S¹; ℤ)\n"
            "  Nicht-null: E²_{0,0}=ℤ, E²_{2,0}=ℤ, E²_{0,1}=ℤ, E²_{2,1}=ℤ\n"
            "  d₂: E²_{2,1} → E²_{0,0} ist Isomorphismus (d₂[S¹×S²]=d₂[S³]=ℤ)\n"
            "  E∞ → H_*(S³): H_0=H_3=ℤ, H_1=H_2=0 ✓"
        )

    def leray_hirsch_theorem_demo(self) -> str:
        """
        @brief Demonstriert den Leray-Hirsch-Satz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Leray-Hirsch: Wenn H*(F; ℝ) von globalen Klassen erzeugt wird,
        dann H*(E; ℝ) ≅ H*(B; ℝ) ⊗ H*(F; ℝ) (als H*(B)-Modul).

        @return Beschreibungsstring.
        """
        return (
            "Leray-Hirsch-Satz:\n"
            "Für F → E → B: Falls ι*: H*(E;ℝ) → H*(F;ℝ) surjektiv,\n"
            "dann H*(E;ℝ) ≅ H*(B;ℝ) ⊗_ℝ H*(F;ℝ) (als Vektorräume).\n"
            "\nAnwendung: Flaggen-Mannigfaltigkeit, Projektivisierung von Bündeln.\n"
            "Projektivisierung P(E) → B eines Vektorbündels E:\n"
            "  H*(P(E); ℤ) ≅ H*(B; ℤ)[ξ] / (ξ^r + c₁·ξ^{r-1} + ... + c_r)\n"
            "  ξ = c₁(O(1)), r = Rang(E)."
        )

    def e2_page_description(self) -> str:
        """
        @brief Beschreibt die E₂-Seite der aktuellen Spektralsequenz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Beschreibungsstring der E₂-Seite.
        """
        return (
            f"Spektralsequenz '{self._name}' – E₂-Seite:\n"
            "E²_{p,q} ist als bigraduiertes Modul aufgebaut.\n"
            "Differentiale d_r: E^r_{p,q} → E^r_{p-r, q+r-1} (Homologie)\n"
            "E^{r+1}_{p,q} = ker(d_r) / im(d_r) (nächste Seite)\n"
            "E^∞_{p,q} = Grad-Stücke des Limits H_{p+q}(E; ℤ)\n"
            "Die Sequenz konvergiert, wenn sie ab Seite r_0 stationär ist."
        )


# ============================================================
# Klasse: KTheory
# ============================================================

class KTheory:
    """
    @brief Beschreibt K-Theorie und Bott-Periodizität.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die K-Theorie K(X) klassifiziert Vektorbündel über X:
      K(X) = Grothendieck-Gruppe der Vektorbündel über X
      K̃(X) = reduzierte K-Theorie (Kern von K(X) → K(Punkt))
    """

    def __init__(self) -> None:
        """
        @brief Initialisiert KTheory.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        pass

    def k_group_of_sphere(self, n: int) -> Dict:
        """
        @brief Gibt die reduzierte K-Theorie-Gruppen K̃(S^n) zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bott-Periodizität für K-Theorie:
          K̃(S^{2k}) = ℤ (erzeugt vom Hopf-Bündel oder Differenzbündel)
          K̃(S^{2k+1}) = 0

        Formel (KaTeX):
            \\tilde{K}(S^n) = \\begin{cases}
              \\mathbb{Z} & n \\text{ gerade} \\\\
              0 & n \\text{ ungerade}
            \\end{cases}

        @param n Dimension der Sphäre (n ≥ 0).
        @return Dict mit K-Theorie-Informationen.
        """
        if n < 0:
            raise ValueError(f"Sphärendimension muss ≥ 0 sein, erhalten: {n}")

        if n % 2 == 0:
            k_reduced = "ℤ"
            generator = "Hopf-Bündel H - 1 (n=2), allg. Bott-Element" if n == 2 else "Bott-Klasse β^{n/2}"
        else:
            k_reduced = "0"
            generator = "trivial"

        return {
            "space": f"S^{n}",
            "k_tilde": k_reduced,
            "k_unreduced": f"ℤ ⊕ {k_reduced}" if k_reduced == "ℤ" else "ℤ",
            "generator": generator,
            "bott_period": f"Bott-Periodizität: K̃(S^n) ≅ K̃(S^{{n+2}})",
        }

    def bott_periodicity_demo(self) -> str:
        """
        @brief Demonstriert die Bott-Periodizität.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bott-Periodizität (komplexe K-Theorie):
          K(X × S²) ≅ K(X) × K(X)  →  Periode 2

        Formel (KaTeX):
            \\tilde{K}(\\Sigma^2 X) \\cong \\tilde{K}(X) \\quad \\text{(Periode 2)}
            \\pi_n(BU \\times \\mathbb{Z}) \\cong \\begin{cases}
              \\mathbb{Z} & n \\text{ gerade} \\\\
              0 & n \\text{ ungerade}
            \\end{cases}

        @return Beschreibungsstring.
        """
        return (
            "Bott-Periodizitätssatz:\n"
            "Komplexe K-Theorie: K̃(Σ²X) ≅ K̃(X) (Periode 2)\n"
            "Reelle K-Theorie: K̃_ℝ(Σ⁸X) ≅ K̃_ℝ(X) (Periode 8)\n"
            "\nFolgerungen:\n"
            "  K̃(S^{2k}) = ℤ für k ≥ 0\n"
            "  K̃(S^{2k+1}) = 0 für k ≥ 0\n"
            "\nAnwendungen:\n"
            "  • Divisionsalgebren: ℝ,ℂ,ℍ,𝕆 (dim 1,2,4,8) – Hopf-Invariante Eins\n"
            "  • Parallelisierbarkeit: S^n parallelisierbar ⟺ n ∈ {1,3,7}\n"
            "  • Index-Theorie: Atiyah-Singer-Index-Satz"
        )

    def chern_character_demo(self) -> str:
        """
        @brief Demonstriert den Chern-Charakter ch: K(X) → H*(X; ℚ).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Der Chern-Charakter verbindet K-Theorie mit rationaler Kohomologie:
          ch(E) = r + c₁(E) + (c₁²-2c₂)/2 + ...

        Formel (KaTeX):
            \\text{ch}(E) = \\sum_{k=0}^{\\infty} \\frac{c_1(E)^k}{k!}
                          = r + c_1 + \\frac{c_1^2 - 2c_2}{2} + \\cdots

        @return Beschreibungsstring.
        """
        return (
            "Chern-Charakter ch: K(X) → H^{even}(X; ℚ):\n"
            "  ch(E) = r + c₁(E) + (c₁(E)²-2c₂(E))/2! + ...\n"
            "  ch(E⊕F) = ch(E) + ch(F)  (additiv)\n"
            "  ch(E⊗F) = ch(E)·ch(F)   (multiplikativ)\n"
            "\nCh ist ein Ringhomomorphismus:\n"
            "  ch: K(X)⊗ℚ → H^{even}(X; ℚ) (Isomorphismus für X endl. CW)\n"
            "\nTodd-Klasse td(TX): ch(D(E))·td(TX) = ∫_X ch(E)·td(TX) (HRR)"
        )

    def atiyah_singer_index_theorem_demo(self) -> str:
        """
        @brief Demonstriert den Atiyah-Singer-Indexsatz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Atiyah-Singer: Für einen elliptischen Differentialoperator D auf M:
          index(D) = ∫_M ch(σ(D)) · td(TM ⊗ ℂ)

        Formel (KaTeX):
            \\operatorname{ind}(D) = \\int_M \\operatorname{ch}(\\sigma(D))
                                     \\cdot \\operatorname{td}(TM_{\\mathbb{C}})

        @return Beschreibungsstring.
        """
        return (
            "Atiyah-Singer-Indexsatz:\n"
            "Für elliptischen Operator D: Γ(E) → Γ(F) auf kompakter M:\n"
            "  ind(D) = dim ker(D) - dim coker(D) (analytischer Index)\n"
            "  = ∫_M ch(σ(D)) · td(TM_ℂ) (topologischer Index)\n"
            "\nSpezialfälle:\n"
            "  • Gauß-Bonnet: ind(d+d*) = χ(M) = ∫_M e(TM)\n"
            "  • Hirzebruch-Signatur: ind(D_L²) = σ(M) = ∫_M L(p₁,...)\n"
            "  • Dirac-Operator: ind(D) = Â(M) = ∫_M Â(p₁,...) (Spin-Mannigfaltigkeit)"
        )


# ============================================================
# Klasse: CWComplex
# ============================================================

class CWComplex:
    """
    @brief Beschreibt CW-Komplexe und ihre zelluläre Homologie.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein CW-Komplex ist ein topologischer Raum, aufgebaut durch sukzessives
    Anheften von Zellen: X = X^(-1) ⊂ X^0 ⊂ X^1 ⊂ ... ⊂ X^n = X

    Formel (KaTeX):
        \\chi(X) = \\sum_k (-1)^k \\cdot \\#\\{k\\text{-Zellen}\\}
    """

    def __init__(self, cells: Dict[int, int]) -> None:
        """
        @brief Initialisiert CWComplex mit Zellenanzahlen pro Dimension.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param cells Dict: k → Anzahl der k-Zellen (z.B. {0: 1, 2: 1} für S²).
        """
        self._cells = cells

    def euler_characteristic(self) -> int:
        """
        @brief Berechnet die Euler-Charakteristik χ = Σ (-1)^k · #k-Zellen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Euler-Charakteristik.
        """
        return sum((-1) ** k * count for k, count in self._cells.items())

    def attaching_map_example(self) -> str:
        """
        @brief Gibt ein Beispiel einer Anheftungsabbildung zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @return Beschreibungsstring.
        """
        # Beispiel basierend auf den vorhandenen Zellen ableiten
        max_dim = max(self._cells.keys()) if self._cells else 0
        n_cells = sum(self._cells.values())

        if self._cells.get(0, 0) == 1 and self._cells.get(2, 0) == 1:
            return (
                "S² als CW-Komplex:\n"
                "  X^0 = {*} (ein 0-Zelle e⁰)\n"
                "  X^1 = X^0 (keine 1-Zellen)\n"
                "  X^2 = X^0 ∪_φ e² (φ: ∂D² → {*} konstant)\n"
                "  Anheftungsabbildung: φ: S¹ → {*} (konstante Abbildung)"
            )
        elif self._cells.get(0, 0) == 1 and self._cells.get(1, 0) == 2 and self._cells.get(2, 0) == 1:
            return (
                "T² (Torus) als CW-Komplex:\n"
                "  e⁰: ein Punkt, e¹_α, e¹_β: zwei 1-Zellen, e²: eine 2-Zelle\n"
                "  Anheftungsabbildung φ: ∂D² → T¹ = S¹∨S¹\n"
                "  φ repräsentiert das Wort aba⁻¹b⁻¹ (Kommutator)"
            )
        return (
            f"CW-Komplex mit {n_cells} Zellen (max. Dim. {max_dim}):\n"
            f"Zellen: {self._cells}\n"
            "Anheftungsabbildungen bestimmen die Topologie des Raums."
        )

    def cellular_homology_demo(self) -> Dict[int, int]:
        """
        @brief Berechnet die zelluläre Homologie (vereinfacht: Betti-Zahlen).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die zelluläre Homologie nutzt Anheftungsgrade, um Randoperatoren
        zu definieren. Ohne explizite Anheftungsabbildungen geben wir
        die Zellenanzahlen als obere Schranken zurück.

        Formel (KaTeX):
            H_k^{CW}(X) \\cong H_k(X) \\quad \\text{(stimmt mit sing. Homologie überein)}

        @return Dict: k → Betti-Zahl (untere Schranke aus Zellstruktur).
        """
        # Vereinfachte Berechnung: ohne Anheftungsgrade nehmen wir Zellanzahl
        result = {}
        for k, count in sorted(self._cells.items()):
            result[k] = count
        return result
