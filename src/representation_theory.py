"""
@file representation_theory.py
@brief Darstellungstheorie: Gruppen-Darstellungen, Charaktere, Schur-Lemma,
       Maschke-Satz, Burnside-Lemma, irreduzible Darstellungen.
@description
    Dieses Modul implementiert die grundlegenden Strukturen der Darstellungstheorie:

    Klassen:
    - Representation  – Darstellung ρ: G → GL(n, ℂ) einer endlichen Gruppe

    Freie Funktionen:
    - trivial_representation()       – Triviale Darstellung ρ(g) = [1]
    - regular_representation()       – Reguläre Darstellung via Cayley-Tabelle
    - character_table()              – Charaktertafel (Irreduzible × Konjugationsklassen)
    - schur_lemma()                  – Schur-Lemma: Intertwiner-Analyse
    - maschke_theorem_check()        – Maschke-Satz über vollständige Reduzibilität
    - decompose_representation()     – Zerlegung in irreduzible Teildarstellungen
    - burnside_lemma()               – Anzahl der Orbits via Fixpunktzählung
    - z2_representations()           – Alle irreduziblen Darst. von ℤ/2ℤ
    - s3_representations()           – Alle irreduziblen Darst. von S₃

    Mathematische Grundlagen:
    - Darstellung: Gruppenhomomorphismus ρ: G → GL(V)
    - Charakter: χ(g) = Tr(ρ(g))
    - Schur-Orthogonalität: ⟨χᵢ, χⱼ⟩ = δᵢⱼ
    - Burnside: |X/G| = (1/|G|) · Σ_{g∈G} |Fix(g)|
    - Maschke: Jede Darst. über ℂ ist vollständig reduzibel

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import numpy as np
from typing import Any
from itertools import product as iproduct


# =============================================================================
# KLASSE: Representation
# =============================================================================

class Representation:
    """
    Darstellung einer endlichen Gruppe G in GL(V).

    Eine Darstellung ist ein Gruppenhomomorphismus
        ρ: G → GL(n, ℂ),
    d.h. ρ(g·h) = ρ(g)·ρ(h) für alle g, h ∈ G.

    @param group_elements  Liste aller Gruppenelemente (hashbar)
    @param matrices        Dictionary {element: numpy-Matrix (n×n, komplex)}
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        group_elements: list,
        matrices: dict,
        multiplication_table: dict | None = None
    ):
        """
        Initialisiert die Darstellung.

        @param group_elements      Geordnete Liste der Gruppenelemente
        @param matrices            {g: np.ndarray} – komplexe Matrizen
        @param multiplication_table  Optional: {(g,h): g*h} für Homomorphismus-Test
        @lastModified 2026-03-10
        """
        # Gruppenelemente speichern
        self.group_elements = list(group_elements)
        # Matrizen als komplexe Arrays speichern
        self.matrices = {g: np.array(m, dtype=complex) for g, m in matrices.items()}
        # Multiplikationstabelle für Homomorphismus-Prüfung
        self.multiplication_table = multiplication_table
        # Dimension aus erster Matrix ableiten
        first = next(iter(self.matrices.values()))
        self._dim = first.shape[0]

    def dimension(self) -> int:
        """
        Gibt die Dimension n des Darstellungsraums ℂⁿ zurück.

        @return  Dimension n der Matrizen (n×n)
        @lastModified 2026-03-10
        """
        return self._dim

    def is_homomorphism(self, tol: float = 1e-9) -> bool:
        """
        Prüft, ob ρ ein Gruppenhomomorphismus ist: ρ(g·h) = ρ(g)·ρ(h).

        Benötigt eine Multiplikationstabelle. Ohne diese wird True zurückgegeben
        (keine Prüfung möglich).

        @param tol  Toleranz für numerischen Vergleich
        @return     True wenn Homomorphismus-Eigenschaft erfüllt
        @lastModified 2026-03-10
        """
        if self.multiplication_table is None:
            # Ohne Multiplikationstabelle kann nicht geprüft werden
            return True

        for g in self.group_elements:
            for h in self.group_elements:
                # Produkt g*h in der Gruppe bestimmen
                gh = self.multiplication_table.get((g, h))
                if gh is None:
                    continue
                # ρ(g*h) berechnen
                rho_gh = self.matrices[gh]
                # ρ(g) * ρ(h) berechnen
                rho_g_rho_h = self.matrices[g] @ self.matrices[h]
                # Vergleich mit Toleranz
                if not np.allclose(rho_gh, rho_g_rho_h, atol=tol):
                    return False
        return True

    def character(self) -> dict:
        """
        Berechnet den Charakter χ(g) = Tr(ρ(g)) für alle g ∈ G.

        Der Charakter ist eine Klassenfunktion und enthält vollständige
        Information über die Darstellung bis auf Äquivalenz.

        @return  Dictionary {g: χ(g)} mit komplexen Spurwerten
        @lastModified 2026-03-10
        """
        return {g: complex(np.trace(self.matrices[g])) for g in self.group_elements}

    def is_irreducible(self, tol: float = 1e-6) -> bool:
        """
        Prüft Irreduzibilität via Schur-Orthogonalitätsrelation.

        Eine Darstellung ist irreduzibel genau dann wenn:
            ⟨χ, χ⟩ = (1/|G|) · Σ_{g∈G} |χ(g)|² = 1

        @param tol  Toleranz für numerischen Vergleich mit 1
        @return     True wenn irreduzibel
        @lastModified 2026-03-10
        """
        chi = self.character()
        n = len(self.group_elements)
        # Inneres Produkt des Charakters mit sich selbst berechnen
        inner = sum(abs(chi[g]) ** 2 for g in self.group_elements) / n
        return abs(inner - 1.0) < tol

    def direct_sum(self, other: 'Representation') -> 'Representation':
        """
        Bildet die direkte Summe zweier Darstellungen: ρ₁ ⊕ ρ₂.

        Die Matrizen werden block-diagonal angeordnet:
            (ρ₁ ⊕ ρ₂)(g) = diag(ρ₁(g), ρ₂(g))

        @param other  Zweite Darstellung (gleiche Gruppe)
        @return       Direkte Summe als neue Darstellung
        @lastModified 2026-03-10
        """
        new_matrices = {}
        for g in self.group_elements:
            m1 = self.matrices[g]
            m2 = other.matrices[g]
            n1, n2 = m1.shape[0], m2.shape[0]
            # Block-diagonale Matrix aufbauen
            block = np.zeros((n1 + n2, n1 + n2), dtype=complex)
            block[:n1, :n1] = m1
            block[n1:, n1:] = m2
            new_matrices[g] = block
        return Representation(
            self.group_elements,
            new_matrices,
            self.multiplication_table
        )

    def tensor_product(self, other: 'Representation') -> 'Representation':
        """
        Bildet das Tensorprodukt zweier Darstellungen: ρ₁ ⊗ ρ₂.

        Verwendet das Kronecker-Produkt der Matrizen:
            (ρ₁ ⊗ ρ₂)(g) = ρ₁(g) ⊗ ρ₂(g)

        @param other  Zweite Darstellung (gleiche Gruppe)
        @return       Tensorprodukt als neue Darstellung
        @lastModified 2026-03-10
        """
        new_matrices = {}
        for g in self.group_elements:
            # Kronecker-Produkt der Matrizen
            new_matrices[g] = np.kron(self.matrices[g], other.matrices[g])
        return Representation(
            self.group_elements,
            new_matrices,
            self.multiplication_table
        )

    def is_equivalent(self, other: 'Representation', tol: float = 1e-6) -> bool:
        """
        Prüft Äquivalenz zweier Darstellungen via Charaktervergleich.

        Zwei Darstellungen sind äquivalent ⟺ ihre Charaktere stimmen überein:
            χ_ρ(g) = χ_σ(g) für alle g ∈ G

        (Vollständige Äquivalenz via Charaktertheorie, da Charaktere Darstellungen
        bis auf Äquivalenz klassifizieren.)

        @param other  Andere Darstellung derselben Gruppe
        @param tol    Toleranz für numerischen Vergleich
        @return       True wenn äquivalent
        @lastModified 2026-03-10
        """
        if self.dimension() != other.dimension():
            return False
        chi1 = self.character()
        chi2 = other.character()
        for g in self.group_elements:
            if abs(chi1.get(g, 0) - chi2.get(g, 0)) > tol:
                return False
        return True


# =============================================================================
# FREIE FUNKTIONEN
# =============================================================================

def trivial_representation(group_elements: list) -> 'Representation':
    """
    Erstellt die triviale Darstellung: ρ(g) = [1] für alle g ∈ G.

    Die triviale Darstellung ist stets irreduzibel und eindimensional.

    @param group_elements  Liste der Gruppenelemente
    @return                Triviale 1D-Darstellung
    @lastModified 2026-03-10
    """
    # Jedes Element wird auf die 1×1-Einheitsmatrix abgebildet
    matrices = {g: np.array([[1.0]], dtype=complex) for g in group_elements}
    return Representation(group_elements, matrices)


def regular_representation(
    group_elements: list,
    multiplication_table: list[list]
) -> 'Representation':
    """
    Erstellt die reguläre Darstellung einer Gruppe.

    Die reguläre Darstellung wirkt durch Links-Multiplikation auf ℂ[G]:
        ρ(g) · e_h = e_{g·h}

    Die Matrix [ρ(g)]_{h,k} = 1 wenn g·k = h, sonst 0.

    @param group_elements     Liste der Gruppenelemente
    @param multiplication_table  2D-Liste: mul_table[i][j] = Index von g_i * g_j
    @return                   Reguläre |G|-dimensionale Darstellung
    @lastModified 2026-03-10
    """
    n = len(group_elements)
    matrices = {}

    for i, g in enumerate(group_elements):
        # Permutationsmatrix für Links-Multiplikation mit g aufbauen
        mat = np.zeros((n, n), dtype=complex)
        for k in range(n):
            # g · e_k = e_{g*k}
            gk_idx = multiplication_table[i][k]
            mat[gk_idx, k] = 1.0
        matrices[g] = mat

    # Multiplikationstabelle als Dict aufbauen
    mult_dict = {}
    for i, g in enumerate(group_elements):
        for j, h in enumerate(group_elements):
            mult_dict[(g, h)] = group_elements[multiplication_table[i][j]]

    return Representation(group_elements, matrices, mult_dict)


def character_table(
    group_elements: list,
    conjugacy_classes: list[list],
    representations: list['Representation']
) -> list[list]:
    """
    Berechnet die Charaktertafel einer Gruppe.

    Die Charaktertafel hat:
    - Zeilen: irreduzible Darstellungen
    - Spalten: Konjugationsklassen
    - Einträge: χᵢ(Cⱼ) = Charakter der i-ten Darst. auf einem Element der j-ten Klasse

    @param group_elements      Alle Gruppenelemente
    @param conjugacy_classes   Liste von Listen; jede enthält Elemente einer Klasse
    @param representations     Liste der (irreduziblen) Darstellungen
    @return                    2D-Liste [Darst.-Index][Klassen-Index] = χᵢ(Cⱼ)
    @lastModified 2026-03-10
    """
    table = []
    for rep in representations:
        chi = rep.character()
        # Für jede Konjugationsklasse den Charakter des repräsentativen Elements nehmen
        row = []
        for cls in conjugacy_classes:
            rep_elem = cls[0]  # Repräsentant der Klasse
            row.append(chi.get(rep_elem, 0.0))
        table.append(row)
    return table


def schur_lemma(
    rho: 'Representation',
    sigma: 'Representation',
    T: np.ndarray,
    tol: float = 1e-8
) -> dict:
    """
    Wendet das Schur-Lemma auf eine Intertwiner-Matrix T an.

    Schur-Lemma: Wenn T·ρ(g) = σ(g)·T für alle g ∈ G gilt, dann:
    - Falls ρ und σ inäquivalent sind: T = 0
    - Falls ρ = σ irreduzibel und T nicht null: T = λI

    @param rho    Erste irreduzible Darstellung
    @param sigma  Zweite irreduzible Darstellung
    @param T      Kandidat-Intertwiner-Matrix
    @param tol    Toleranz für Null-Prüfung
    @return       Dict mit 'is_intertwiner', 'type' ('zero'|'scalar'|'general'), 'scalar'
    @lastModified 2026-03-10
    """
    T = np.array(T, dtype=complex)

    # Prüfen ob T tatsächlich ein Intertwiner ist: T·ρ(g) = σ(g)·T
    is_intertwiner = True
    for g in rho.group_elements:
        lhs = T @ rho.matrices[g]
        rhs = sigma.matrices[g] @ T
        if not np.allclose(lhs, rhs, atol=tol):
            is_intertwiner = False
            break

    if not is_intertwiner:
        return {
            'is_intertwiner': False,
            'type': 'none',
            'scalar': None,
            'message': 'T ist kein Intertwiner (T·ρ(g) ≠ σ(g)·T)'
        }

    # Prüfen ob T die Nullmatrix ist
    if np.allclose(T, 0, atol=tol):
        return {
            'is_intertwiner': True,
            'type': 'zero',
            'scalar': None,
            'message': 'T = 0 (ρ und σ inäquivalent, oder T ist trivial null)'
        }

    # Prüfen ob T ein Vielfaches der Identität ist
    n = T.shape[0]
    if T.shape[0] == T.shape[1]:
        # Skalar λ = T[0,0] wenn T = λI
        lam = T[0, 0]
        if np.allclose(T, lam * np.eye(n), atol=tol):
            return {
                'is_intertwiner': True,
                'type': 'scalar',
                'scalar': complex(lam),
                'message': f'T = λI mit λ = {lam:.6f} (ρ = σ irreduzibel)'
            }

    return {
        'is_intertwiner': True,
        'type': 'general',
        'scalar': None,
        'message': 'T ist Intertwiner aber kein Skalar-Vielfaches der Identität'
    }


def maschke_theorem_check(group_order: int, char_field: int = 0) -> dict:
    """
    Überprüft die Voraussetzungen des Maschke-Satzes.

    Maschke-Satz: Sei G eine endliche Gruppe, k ein Körper.
    Falls char(k) = 0 oder char(k) ∤ |G|, dann ist jede k[G]-Darstellung
    vollständig reduzibel (direkte Summe irreduzibler Darstellungen).

    @param group_order   Ordnung der Gruppe |G|
    @param char_field    Charakteristik des Körpers (0 für ℂ, ℝ, ℚ)
    @return              Dict mit 'applicable', 'reason', 'conclusion'
    @lastModified 2026-03-10
    """
    if char_field == 0:
        # Über ℂ, ℝ oder ℚ gilt der Satz immer
        return {
            'applicable': True,
            'char_field': 0,
            'group_order': group_order,
            'reason': 'char(k) = 0 (z.B. ℂ, ℝ, ℚ): teilt nicht |G|',
            'conclusion': 'Jede Darstellung ist vollständig reduzibel (direkte Summe irreduzibler Darstellungen).'
        }

    # Prüfen ob char_field die Gruppenordnung teilt
    if group_order % char_field == 0:
        return {
            'applicable': False,
            'char_field': char_field,
            'group_order': group_order,
            'reason': f'char(k) = {char_field} teilt |G| = {group_order}',
            'conclusion': 'Maschke-Satz gilt NICHT. Es können unzerlegbare, nicht-irreduzible Darstellungen existieren.'
        }
    else:
        return {
            'applicable': True,
            'char_field': char_field,
            'group_order': group_order,
            'reason': f'char(k) = {char_field} teilt nicht |G| = {group_order}',
            'conclusion': 'Maschke-Satz gilt. Jede Darstellung ist vollständig reduzibel.'
        }


def decompose_representation(
    rep: 'Representation',
    irreps: list['Representation'],
    tol: float = 1e-6
) -> dict:
    """
    Zerlegt eine Darstellung in irreduzible Teildarstellungen.

    Multiplizität der i-ten irreduziblen Darstellung:
        mᵢ = (1/|G|) · Σ_{g∈G} χ(g) · χᵢ(g)*

    Dies nutzt die Schur-Orthogonalitätsrelationen der Charaktere.

    @param rep     Zu zerlegende Darstellung
    @param irreps  Liste der irreduziblen Darstellungen
    @param tol     Toleranz für Rundung auf ganzzahlige Multiplizitäten
    @return        Dict {'multiplicities': [m₀, m₁, ...], 'decomposition': str}
    @lastModified 2026-03-10
    """
    n = len(rep.group_elements)
    chi = rep.character()
    multiplicities = []

    for irrep in irreps:
        chi_i = irrep.character()
        # Inneres Produkt ⟨χ, χᵢ⟩ = (1/|G|) Σ_g χ(g)·χᵢ(g)*
        inner = sum(
            chi[g] * np.conj(chi_i[g]) for g in rep.group_elements
        ) / n
        # Multiplizität sollte ganzzahlig sein
        m = round(inner.real)
        multiplicities.append(m)

    # Zerlegungsformel als String aufbauen
    parts = []
    for i, m in enumerate(multiplicities):
        if m > 0:
            dim = irreps[i].dimension()
            parts.append(f"{m}·ρ_{i+1}(dim={dim})")
    decomp_str = " ⊕ ".join(parts) if parts else "0"

    return {
        'multiplicities': multiplicities,
        'decomposition': f"ρ ≅ {decomp_str}",
        'total_dimension': sum(m * irreps[i].dimension() for i, m in enumerate(multiplicities))
    }


def burnside_lemma(
    group_elements: list,
    action: dict
) -> int:
    """
    Berechnet die Anzahl der Orbits mit dem Burnside-Lemma (auch Cauchy-Frobenius-Lemma).

    |X/G| = (1/|G|) · Σ_{g∈G} |Fix(g)|

    wobei Fix(g) = {x ∈ X : g·x = x} die Fixpunktmenge von g ist.

    @param group_elements  Liste der Gruppenelemente
    @param action          {g: permutation_as_list} – Permutation von X als Index-Liste
                           action[g][x] = g·x (Index des Bildes)
    @return                Anzahl der Orbits (ganzzahlig)
    @lastModified 2026-03-10
    """
    total_fixed = 0

    for g in group_elements:
        perm = action[g]
        # Fixpunkte zählen: Elemente x mit perm[x] = x
        fixed = sum(1 for x in range(len(perm)) if perm[x] == x)
        total_fixed += fixed

    # Burnside-Formel: Mittelwert der Fixpunkte
    n = len(group_elements)
    # Das Ergebnis ist immer ganzzahlig
    return total_fixed // n


# =============================================================================
# BEISPIEL-DARSTELLUNGEN
# =============================================================================

def z2_representations() -> list['Representation']:
    """
    Erstellt alle irreduziblen Darstellungen von ℤ/2ℤ = {0, 1}.

    ℤ/2ℤ hat genau 2 irreduzible Darstellungen (Grad 1 = Anzahl Konjugationsklassen):
    - Triviale Darst.: ρ(0) = [1], ρ(1) = [1]
    - Signum-Darst.:  ρ(0) = [1], ρ(1) = [-1]

    Multiplikation: 0+0=0, 0+1=1, 1+0=1, 1+1=0

    @return  Liste [triviale_Darst., signum_Darst.]
    @lastModified 2026-03-10
    """
    elements = [0, 1]
    # Multiplikationstabelle für ℤ/2ℤ: Addition mod 2
    mult_table = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0}

    # Triviale Darstellung
    trivial_mats = {
        0: np.array([[1.0]], dtype=complex),
        1: np.array([[1.0]], dtype=complex)
    }
    trivial = Representation(elements, trivial_mats, mult_table)

    # Signum-Darstellung (alternierende Darstellung)
    sign_mats = {
        0: np.array([[1.0]], dtype=complex),
        1: np.array([[-1.0]], dtype=complex)
    }
    sign = Representation(elements, sign_mats, mult_table)

    return [trivial, sign]


def s3_representations() -> list['Representation']:
    """
    Erstellt alle irreduziblen Darstellungen von S₃ (Ordnung 6).

    S₃ = {e, (12), (13), (23), (123), (132)} hat 3 Konjugationsklassen:
    - {e}
    - {(12), (13), (23)}  (Transpositionen)
    - {(123), (132)}      (3-Zyklen)

    Damit gibt es genau 3 irreduzible Darstellungen:
    1. Triviale Darst. (1D): ρ(σ) = [1]
    2. Signum-Darst.  (1D): ρ(σ) = [sgn(σ)]
    3. Standard-Darst.(2D): Orthogonale Wirkung auf das 2D-Komplement

    @return  Liste [trivial, signum, standard]
    @lastModified 2026-03-10
    """
    # Gruppenelemente als Tupel-Permutationen von {0,1,2}
    # e=(0,1,2), (12)=(1,0,2), (13)=(2,1,0), (23)=(0,2,1), (123)=(1,2,0), (132)=(2,0,1)
    e   = 'e'
    t12 = '(12)'
    t13 = '(13)'
    t23 = '(23)'
    c123 = '(123)'
    c132 = '(132)'
    elements = [e, t12, t13, t23, c123, c132]

    # Multiplikationstabelle für S₃ (Links-Komposition σ∘τ)
    mult_table = {
        (e, e): e,       (e, t12): t12,   (e, t13): t13,
        (e, t23): t23,   (e, c123): c123, (e, c132): c132,
        (t12, e): t12,   (t12, t12): e,   (t12, t13): c132,
        (t12, t23): c123,(t12, c123): t23,(t12, c132): t13,
        (t13, e): t13,   (t13, t12): c123,(t13, t13): e,
        (t13, t23): c132,(t13, c123): t12,(t13, c132): t23,
        (t23, e): t23,   (t23, t12): c132,(t23, t13): c123,
        (t23, t23): e,   (t23, c123): t13,(t23, c132): t12,
        (c123, e): c123, (c123, t12): t13,(c123, t13): t23,
        (c123, t23): t12,(c123, c123): c132,(c123, c132): e,
        (c132, e): c132, (c132, t12): t23,(c132, t13): t12,
        (c132, t23): t13,(c132, c123): e, (c132, c132): c123,
    }

    # -------------------------------------------------------------------------
    # 1) Triviale Darstellung
    # -------------------------------------------------------------------------
    trivial_mats = {g: np.array([[1.0 + 0j]]) for g in elements}
    trivial = Representation(elements, trivial_mats, mult_table)

    # -------------------------------------------------------------------------
    # 2) Signum-Darstellung
    # -------------------------------------------------------------------------
    # Vorzeichen: gerade Permutationen → +1, ungerade → -1
    signs = {e: 1, t12: -1, t13: -1, t23: -1, c123: 1, c132: 1}
    sign_mats = {g: np.array([[float(signs[g]) + 0j]]) for g in elements}
    signum = Representation(elements, sign_mats, mult_table)

    # -------------------------------------------------------------------------
    # 3) Standard-Darstellung (2D)
    # Wirkung auf ℝ²: Spiegelungen und Rotationen des gleichseitigen Dreiecks
    # -------------------------------------------------------------------------
    # Rotation um 120° und 240°
    cos120 = np.cos(2 * np.pi / 3)
    sin120 = np.sin(2 * np.pi / 3)
    rot120 = np.array([
        [cos120, -sin120],
        [sin120,  cos120]
    ], dtype=complex)
    rot240 = rot120 @ rot120  # Umgekehrte Rotation

    # Spiegelungen: (12) spiegelt an x-Achse
    refl_12 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    # (13) = refl um 60° Achse = rot120 @ refl_12 @ rot120^{-1}
    refl_13 = rot120 @ refl_12 @ np.linalg.inv(rot120)
    # (23) = refl um 120° Achse
    refl_23 = rot240 @ refl_12 @ np.linalg.inv(rot240)

    std_mats = {
        e:    np.eye(2, dtype=complex),
        t12:  refl_12,
        t13:  refl_13,
        t23:  refl_23,
        # Hinweis: (123) entspricht Rotation um 240° und (132) um 120°,
        # da Links-Komposition die Richtung der Matrizenmultiplikation umkehrt.
        c123: rot240,
        c132: rot120
    }
    standard = Representation(elements, std_mats, mult_table)

    return [trivial, signum, standard]
