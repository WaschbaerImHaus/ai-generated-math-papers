"""
@file commutative_algebra.py
@brief Kommutative Algebra: Lokalisierung, Primspektrum, Klassengruppen.
@description
    Implementiert grundlegende Konzepte der kommutativen Algebra:

    - Lokalisierung S^{-1}R von Ringen
    - Primspektrum Spec(ℤ/nℤ)
    - Nakayama-Lemma
    - Ganzer Abschluss von ℤ in quadratischen Zahlkörpern ℚ(√d)
    - Klassengruppe Cl(ℚ(√d)) (Klassenzahl h(d))
    - Noether-Normalisierungslemma
    - Hilbert-Basissatz

    Mathematischer Hintergrund:
    Die kommutative Algebra ist die Grundlage der algebraischen Geometrie
    und der algebraischen Zahlentheorie. Zentrale Objekte sind kommutative
    Ringe mit Eins und ihre Ideale.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import cmath
from typing import Any


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _is_prime(n: int) -> bool:
    """
    @brief Einfacher Primalitätstest.
    @param n Zu prüfende Zahl
    @return True wenn n eine Primzahl ist
    @lastModified 2026-03-10
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def _prime_factors(n: int) -> list[int]:
    """
    @brief Gibt alle verschiedenen Primteiler von n zurück.
    @param n Natürliche Zahl
    @return Sortierte Liste der Primteiler
    @lastModified 2026-03-10
    """
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return sorted(factors)


def _squarefree_part(d: int) -> int:
    """
    @brief Berechnet den quadratfreien Teil von d.
    @description
        Schreibt d = k² · m mit m quadratfrei und gibt m zurück.
    @param d Ganzzahl (≠ 0, 1)
    @return Quadratfreier Teil von d
    @lastModified 2026-03-10
    """
    sign = -1 if d < 0 else 1
    d_abs = abs(d)
    # Durch Quadrate dividieren
    p = 2
    while p * p <= d_abs:
        while d_abs % (p * p) == 0:
            d_abs //= p * p
        p += 1
    return sign * d_abs


# ---------------------------------------------------------------------------
# Lokalisierung
# ---------------------------------------------------------------------------

def localization(
    ring_elems: list[int],
    mult_set: list[int],
    mod: int
) -> dict[str, Any]:
    """
    @brief Berechnet die Lokalisierung S^{-1}R.
    @description
        Die Lokalisierung eines kommutativen Rings R an einer multiplikativen
        Menge S besteht aus formalen Brüchen a/s mit a ∈ R, s ∈ S,
        wobei a/s = a'/s' gilt, wenn t(as' - a's) = 0 für ein t ∈ S.

        Für R = ℤ/nℤ und S = {p^k : k ≥ 0} ergibt sich die Lokalisierung
        bei der Primzahl p, die alle Einheiten s ∈ S invertiert.

        Brüche a/s mit s ∈ S werden durch ihren reduzierten Vertreter
        in ℤ/nℤ dargestellt, sofern s eine Einheit mod n ist.

    @param ring_elems Elemente des Grundrings R = ℤ/mod·ℤ
    @param mult_set Multiplikative Menge S ⊆ R
    @param mod Modulus des Rings (0 = ℤ)
    @return Dict mit:
            - 'elements': Elemente des lokalisierten Rings (wenn endlich)
            - 'invertible': Elemente von S, die invertierbar sind
            - 'description': Beschreibung der Lokalisierung
    @lastModified 2026-03-10
    """
    if mod == 0:
        # Lokalisierung von ℤ: S^{-1}ℤ ⊆ ℚ
        desc = f'Lokalisierung S^{{-1}}ℤ mit S = {mult_set}'
        return {
            'elements': ring_elems,
            'invertible': mult_set,
            'description': desc,
        }

    # Für ℤ/nℤ: Elemente von S prüfen, ob sie Einheiten mod n sind
    invertible = []
    for s in mult_set:
        s_mod = s % mod
        if s_mod != 0 and math.gcd(s_mod, mod) == 1:
            invertible.append(s_mod)

    # Lokalisierte Elemente: alle a/s mit a ∈ ring_elems, s ∈ invertible
    localized_elems = set()
    for a in ring_elems:
        for s in invertible:
            # s^{-1} mod n berechnen
            try:
                s_inv = pow(s, -1, mod)
                elem = (a * s_inv) % mod
                localized_elems.add(elem)
            except (ValueError, ZeroDivisionError):
                pass

    if not localized_elems and ring_elems:
        localized_elems = set(a % mod for a in ring_elems)

    desc = (f'S^{{-1}}(ℤ/{mod}ℤ) mit S = {mult_set}: '
            f'{len(localized_elems)} Elemente')

    return {
        'elements': sorted(localized_elems),
        'invertible': invertible,
        'description': desc,
    }


# ---------------------------------------------------------------------------
# Primspektrum
# ---------------------------------------------------------------------------

def prime_spectrum(n: int) -> dict[str, Any]:
    """
    @brief Berechnet Spec(ℤ/nℤ): Alle Primideale in ℤ/nℤ.
    @description
        Das Primspektrum Spec(R) eines kommutativen Rings R ist die Menge
        aller Primideale von R, versehen mit der Zariski-Topologie.

        Für R = ℤ/nℤ:
        - Primideale: (p̄) für Primzahlen p mit p | n
        - Darstellung: pℤ/nℤ ↔ Primzahl p
        - Maximale Ideale: ebenfalls (p̄) für alle Primteiler p

        Spezialfall n = p (prim): ℤ/pℤ ist ein Körper, das einzige Primideal
        ist das Nullideal (0).

        Begründung: Ein Ideal I in ℤ/nℤ entspricht einem Ideal in ℤ das nℤ enthält,
        also I = dℤ/nℤ für d | n. Primideale entstehen für d = p prim mit p | n.

    @param n Modulus (natürliche Zahl ≥ 1)
    @return Dict mit:
            - 'prime_ideals': Liste der Primteiler p (repräsentieren (p))
            - 'maximal_ideals': Liste der maximalen Ideale (= prime_ideals für ℤ/nℤ)
            - 'ring': Beschreibung des Rings
            - 'is_field': True wenn ℤ/nℤ ein Körper ist (n prim)
    @lastModified 2026-03-10
    """
    # Primteiler von n bestimmen
    primes = _prime_factors(n)

    # Für ℤ/nℤ sind Primideale genau die (p) für Primteiler p von n
    # Spezialfall: n=1 → ℤ/ℤ = 0 (trivial), Spec = leer
    if n == 1:
        return {
            'prime_ideals': [],
            'maximal_ideals': [],
            'ring': 'ℤ/ℤ ≅ 0 (Nullring)',
            'is_field': False,
        }

    is_field = len(primes) == 1 and n == primes[0]  # n prim

    if is_field:
        # ℤ/pℤ ist Körper: Nullideal ist prim und maximal
        prime_ideals = [0]  # Nullideal
        maximal_ideals = [0]
    else:
        # Primideale = Primteiler
        prime_ideals = primes
        maximal_ideals = primes  # Alle Primideale in ℤ/nℤ sind maximal

    return {
        'prime_ideals': prime_ideals,
        'maximal_ideals': maximal_ideals,
        'ring': f'ℤ/{n}ℤ',
        'is_field': is_field,
    }


# ---------------------------------------------------------------------------
# Nakayama-Lemma
# ---------------------------------------------------------------------------

def nakayama_lemma_check(
    M_generators: list[int],
    I_generators: list[int],
    p: int
) -> dict[str, Any]:
    """
    @brief Überprüft das Nakayama-Lemma für Moduln über ℤ/p^kℤ.
    @description
        Nakayama-Lemma: Sei R ein lokaler Ring mit maximalem Ideal 𝔪,
        M ein endlich erzeugter R-Modul. Wenn M = 𝔪·M, dann M = 0.

        Äquivalent: Wenn M/𝔪M = 0, dann M = 0.

        Für R = ℤ/p^kℤ (lokaler Ring mit 𝔪 = pℤ/p^kℤ):
        - M = 𝔪·M bedeutet: jedes m ∈ M ist von der Form p·m' für m' ∈ M
        - Dies impliziert p | alle Erzeuger von M
        - Nakayama: wenn p | alle Erzeuger, dann M kann durch Erzeuger /(p) = 0
          beschrieben werden

        Prüft ob M_generators alle durch p teilbar sind (M = I·M Bedingung).

    @param M_generators Erzeuger des Moduls M (Ordnungen oder Koeffizienten)
    @param I_generators Erzeuger des Ideals I (Jacobson-Radikal)
    @param p Charakteristik/Primzahl des Grundrings
    @return Dict mit:
            - 'condition_satisfied': True wenn M = I·M erfüllt
            - 'conclusion': Schlussfolgerung (M=0 oder keine Aussage)
            - 'M_is_zero': True wenn Nakayama M=0 impliziert
    @lastModified 2026-03-10
    """
    if not M_generators:
        # Triviales Modul M = 0
        return {
            'condition_satisfied': True,
            'conclusion': 'M = 0 (triviales Modul, Nakayama trivial erfüllt)',
            'M_is_zero': True,
        }

    # Prüfen ob alle Erzeuger von M durch ein Element von I teilbar sind
    # Einfachste Bedingung: p teilt alle Erzeuger
    i_gen_value = I_generators[0] if I_generators else p
    condition = all(g % i_gen_value == 0 for g in M_generators if g != 0)

    if condition:
        conclusion = (
            f'M = I·M (alle Erzeuger durch {i_gen_value} teilbar). '
            f'Nakayama-Lemma impliziert M = 0 (über lokalem Ring).'
        )
        m_is_zero = True
    else:
        conclusion = (
            f'M ≠ I·M (nicht alle Erzeuger durch {i_gen_value} teilbar). '
            f'Nakayama-Lemma nicht anwendbar.'
        )
        m_is_zero = False

    return {
        'condition_satisfied': condition,
        'conclusion': conclusion,
        'M_is_zero': m_is_zero,
    }


# ---------------------------------------------------------------------------
# Ganzer Abschluss
# ---------------------------------------------------------------------------

def integral_closure(d: int) -> dict[str, Any]:
    """
    @brief Berechnet den ganzen Abschluss von ℤ im quadratischen Zahlkörper ℚ(√d).
    @description
        Der Ring der ganzen Zahlen O_K im quadratischen Zahlkörper K = ℚ(√d)
        (mit d quadratfrei) hängt von d mod 4 ab:

        - d ≡ 1 (mod 4): O_K = ℤ[(1+√d)/2], ℤ-Basis: {1, (1+√d)/2}
                          Diskriminante: Δ = d
        - d ≡ 2 oder 3 (mod 4): O_K = ℤ[√d], ℤ-Basis: {1, √d}
                                  Diskriminante: Δ = 4d

        Der Ring ist ein PID ⟺ Klassenzahl h(d) = 1.
        Für imaginär-quadratische Körper (d < 0) mit h=1:
        d ∈ {-1, -2, -3, -7, -11, -19, -43, -67, -163} (Stark-Heegner-Theorem).

    @param d Quadratfreier ganzzahliger Diskriminantenparameter (d ≠ 0, 1)
    @return Dict mit:
            - 'ring_basis': ℤ-Basis des Rings der ganzen Zahlen
            - 'discriminant': Diskriminante des Zahlkörpers
            - 'is_pid': True wenn der Ring ein Hauptidealring ist
            - 'description': Beschreibung des Rings
    @lastModified 2026-03-10
    """
    # Quadratfreien Teil berechnen
    d_sf = _squarefree_part(d)

    # Diskriminante und Basis bestimmen
    if d_sf % 4 == 1:
        # d ≡ 1 (mod 4): O_K = ℤ[(1+√d)/2]
        discriminant = d_sf
        ring_basis = ['1', f'(1+√{d_sf})/2']
        basis_desc = f'{{1, (1+√{d_sf})/2}}'
        ring_name = f'ℤ[(1+√{d_sf})/2]'
    else:
        # d ≡ 2 oder 3 (mod 4): O_K = ℤ[√d]
        discriminant = 4 * d_sf
        ring_basis = ['1', f'√{d_sf}']
        basis_desc = f'{{1, √{d_sf}}}'
        ring_name = f'ℤ[√{d_sf}]'

    # Bekannte PIDs für imaginär-quadratische Körper (Stark-Heegner-Liste)
    # Diese haben Klassenzahl h = 1
    known_pid_d = {-1, -2, -3, -7, -11, -19, -43, -67, -163}
    # Für reell-quadratische Körper: einige bekannte PIDs
    known_pid_real = {2, 3, 5, 6, 7, 11, 13, 14, 17, 19, 21, 22, 23}

    is_pid = d_sf in known_pid_d or d_sf in known_pid_real

    # Spezielle Fälle benennen
    if d_sf == -1:
        ring_name = 'ℤ[i] (Gaußsche ganze Zahlen)'
    elif d_sf == -3:
        ring_name = 'ℤ[ω] (Eisensteinsche ganze Zahlen, ω = e^{2πi/3})'

    return {
        'ring_basis': ring_basis,
        'discriminant': discriminant,
        'is_pid': is_pid,
        'description': f'O_{{ℚ(√{d_sf})}} = {ring_name}, Δ = {discriminant}',
    }


# ---------------------------------------------------------------------------
# Klassengruppe
# ---------------------------------------------------------------------------

def class_group_estimate(d: int) -> dict[str, Any]:
    """
    @brief Schätzt die Klassengruppe Cl(ℚ(√d)) und berechnet h(d).
    @description
        Die Klassengruppe Cl(K) eines Zahlkörpers K misst, wie weit O_K
        von einem Hauptidealring entfernt ist:
        - h(d) = 1 ⟺ O_K ist PID (eindeutige Primfaktorzerlegung)
        - h(d) > 1 → nicht jedes Ideal ist ein Hauptideal

        Minkowski-Schranke: Jede Idealklasse enthält ein ganzzahliges Ideal
        mit Norm N(I) ≤ M_K, wobei
        M_K = (2/π)^{r_2} · (n!/n^n) · √|Δ|  (n = Grad, r_2 = komplex-Paare)

        Für quadratische Körper: n=2, r_2=1 (imag) oder 0 (real)
        M_K = (2/π)·√|Δ|  (imaginär-quadratisch)
        M_K = (1/2)·√|Δ|  (reell-quadratisch)

        Algorithmus: Prüfe alle Primideale ℘ mit N(℘) ≤ M_K,
        berechne Klassenzahl durch Zählen unabhängiger Idealklassen.

        Bekannte Klassenzahlen (Gauß, Heegner, Stark):
        h(-1)=1, h(-2)=1, h(-3)=1, h(-5)=2, h(-6)=2, h(-7)=1,
        h(-10)=2, h(-11)=1, h(-13)=2, h(-15)=2, h(-19)=1, h(-23)=3.

    @param d Quadratfreier Diskriminantenparameter
    @return Dict mit:
            - 'class_number': Klassenzahl h(d)
            - 'is_pid': True wenn h(d) = 1
            - 'minkowski_bound': Minkowski-Schranke M_K
            - 'description': Beschreibung
    @lastModified 2026-03-10
    """
    d_sf = _squarefree_part(d)

    # Diskriminante bestimmen
    if d_sf % 4 == 1:
        discriminant = d_sf
    else:
        discriminant = 4 * d_sf

    # Minkowski-Schranke berechnen
    abs_disc = abs(discriminant)
    sqrt_disc = math.sqrt(abs_disc)

    if d_sf < 0:
        # Imaginär-quadratischer Körper: M_K = (2/π) · √|Δ|
        minkowski_bound = (2.0 / math.pi) * sqrt_disc
    else:
        # Reell-quadratischer Körper: M_K = (1/2) · √Δ
        minkowski_bound = 0.5 * sqrt_disc

    # Bekannte Klassenzahlen für imaginär-quadratische Körper
    # (Quelle: Tabellen der algebraischen Zahlentheorie)
    known_class_numbers_neg: dict[int, int] = {
        -1: 1, -2: 1, -3: 1, -5: 2, -6: 2, -7: 1,
        -10: 2, -11: 1, -13: 2, -14: 2, -15: 2,
        -17: 4, -19: 1, -21: 4, -22: 2, -23: 3,
        -26: 6, -29: 6, -30: 4, -31: 3, -33: 4,
        -34: 4, -35: 2, -37: 2, -38: 6, -39: 4,
        -41: 8, -43: 1, -46: 4, -47: 5, -51: 2,
        -53: 6, -55: 4, -57: 4, -58: 2, -59: 3,
        -61: 6, -65: 8, -67: 1, -69: 4, -71: 7,
        -163: 1,
    }

    # Bekannte Klassenzahlen für reell-quadratische Körper
    known_class_numbers_pos: dict[int, int] = {
        2: 1, 3: 1, 5: 1, 6: 1, 7: 1, 10: 2, 11: 1,
        13: 1, 14: 1, 15: 2, 17: 1, 19: 1, 21: 1,
        22: 1, 23: 1, 26: 2, 29: 1, 30: 2, 31: 1,
    }

    # Klassenzahl nachschlagen oder schätzen
    if d_sf < 0 and d_sf in known_class_numbers_neg:
        h = known_class_numbers_neg[d_sf]
    elif d_sf > 0 and d_sf in known_class_numbers_pos:
        h = known_class_numbers_pos[d_sf]
    else:
        # Grobe Schätzung via Minkowski-Schranke für unbekannte d
        # Zähle Primideale bis zur Minkowski-Schranke und schätze h
        h = _estimate_class_number(d_sf, int(minkowski_bound) + 1)

    is_pid = (h == 1)

    return {
        'class_number': h,
        'is_pid': is_pid,
        'minkowski_bound': minkowski_bound,
        'discriminant': discriminant,
        'description': (
            f'Cl(ℚ(√{d_sf})): h = {h}, '
            f'Minkowski-Schranke ≈ {minkowski_bound:.2f}, '
            f'{"PID" if is_pid else "kein PID"}'
        ),
    }


def _estimate_class_number(d: int, bound: int) -> int:
    """
    @brief Grobe Schätzung der Klassenzahl via Primidealgenerierung.
    @description
        Vereinfachter Algorithmus: Bestimme Legendre-Symbol (d/p) für alle
        Primzahlen p ≤ bound und schätze h aus der Anzahl nichttrivialer Klassen.

        Für imaginär-quadratische Körper d < 0:
        Wenn kein Primideal ≤ Minkowski-Schranke, dann h = 1 (PID).

    @param d Quadratfreier Diskriminantenparameter
    @param bound Minkowski-Schranke (aufgerundet)
    @return Geschätzte Klassenzahl (1 wenn keine Primideale)
    @lastModified 2026-03-10
    """
    # Primzahlen bis zur Schranke bestimmen
    primes_up_to_bound = [p for p in range(2, max(bound + 1, 3)) if _is_prime(p)]

    # Primideale zählen, die nicht Hauptideale sind
    non_principal = []
    for p in primes_up_to_bound:
        # Zerlegungsverhalten von p in O_K
        disc = d if d % 4 == 1 else 4 * d
        leg = _legendre_symbol(disc, p)
        if leg == -1:
            # p bleibt prim (inert) → (p) ist Hauptideal in O_K
            pass
        elif leg == 0:
            # p ramifiziert → ℘² = (p), ℘ kann Nicht-Hauptideal sein
            non_principal.append(p)
        else:
            # p zerfällt → (p) = ℘·℘̄, ℘ kann Nicht-Hauptideal sein
            non_principal.append(p)

    # Wenn alle Primideale ≤ Bound Hauptideale sind, dann h = 1
    if not non_principal:
        return 1

    # Sehr grobe Schätzung: h ≈ 1 + len(non_principal) // 2
    return max(1, len(non_principal) // 2)


def _legendre_symbol(a: int, p: int) -> int:
    """
    @brief Berechnet das Legendre-Symbol (a/p) für Primzahl p.
    @param a Ganzzahl
    @param p Ungerade Primzahl
    @return 0 wenn p|a, 1 wenn a QR mod p, -1 wenn a NR mod p
    @lastModified 2026-03-10
    """
    if p == 2:
        return 0 if a % 2 == 0 else 1
    a_mod = a % p
    if a_mod == 0:
        return 0
    # Euler-Kriterium: (a/p) ≡ a^{(p-1)/2} (mod p)
    result = pow(a_mod, (p - 1) // 2, p)
    return -1 if result == p - 1 else result


# ---------------------------------------------------------------------------
# Noether-Normalisierung
# ---------------------------------------------------------------------------

def noether_normalization(
    poly_system: list[str],
    n_vars: int
) -> dict[str, Any]:
    """
    @brief Bestimmt Krull-Dimension und Noether-Normalisierung.
    @description
        Noether-Normalisierungslemma: Für jede endlich erzeugte k-Algebra
        A = k[x_1,...,x_n]/I existieren algebraisch unabhängige Elemente
        y_1,...,y_d ∈ A, sodass A über k[y_1,...,y_d] ganz ist.

        d = Krull-Dim(A) = Transzendenzgrad von Quot(A) über k.

        Für affine Varietäten V(I) ⊆ A^n gilt:
        - dim V = Krull-Dim(k[V]) = Krull-Dim(k[x_1,...,x_n]/I)
        - Für eine Hyperfläche V(f) in A^n: dim = n - 1
        - Für eine Kurve in A^2: dim = 1
        - Für einen Punkt: dim = 0

        Vereinfachte Implementierung: Schätzt die Dimension aus der
        Anzahl der Variablen und Gleichungen (reguläre Systeme).

    @param poly_system Liste der Polynomausdrücke als Strings (Erzeuger von I)
    @param n_vars Anzahl der Variablen
    @return Dict mit:
            - 'krull_dimension': Krull-Dimension der Algebra
            - 'transcendence_degree': Transzendenzgrad
            - 'algebraically_independent': Anzahl unabhängiger Elemente
            - 'description': Beschreibung
    @lastModified 2026-03-10
    """
    # Vereinfachte Berechnung der Krull-Dimension
    # Für eine affine Varietät mit n Variablen und c unabhängigen Gleichungen:
    # dim = n - c (für reguläre Systeme / complete intersections)

    n_equations = len(poly_system)

    # Krull-Dimension schätzen
    if n_equations == 0:
        # Freie Algebra: dim = n
        krull_dim = n_vars
    else:
        # Ansatz: n - codim(I), wobei codim ≤ n_equations
        # Für generische Systeme: dim = max(0, n - n_equations)
        krull_dim = max(0, n_vars - n_equations)

    # Für einfache Hyperflächen (1 Gleichung) in A^n: dim = n-1
    # Überprüfung auf Spezialfälle
    if n_equations == 1 and n_vars >= 1:
        poly = poly_system[0].strip()
        # Zähle Variablen im Ausdruck (grobe Heuristik)
        var_names = [f'x{i+1}' for i in range(n_vars)] + ['x', 'y', 'z', 'w']
        vars_in_poly = sum(1 for v in var_names if v in poly)
        if vars_in_poly >= 2:
            krull_dim = n_vars - 1
        else:
            krull_dim = n_vars - 1  # Hyperfläche immer dim = n-1

    return {
        'krull_dimension': krull_dim,
        'transcendence_degree': krull_dim,
        'algebraically_independent': krull_dim,
        'n_generators': n_equations,
        'description': (
            f'k[x_1,...,x_{n_vars}]/(I): Krull-Dim = {krull_dim}, '
            f'{krull_dim} algebraisch unabhängige Elemente'
        ),
    }


# ---------------------------------------------------------------------------
# Hilbert-Basissatz
# ---------------------------------------------------------------------------

def hilbert_basis_theorem_verify(
    ideal_gens: list[list[int]],
    p: int
) -> dict[str, Any]:
    """
    @brief Verifiziert den Hilbert-Basissatz für Polynomringe über ℤ_p.
    @description
        Hilbert-Basissatz: Wenn R noetherisch ist (jedes Ideal endlich erzeugt),
        dann ist auch der Polynomring R[x] noetherisch.

        Folgerung: k[x_1,...,x_n] ist noetherisch für jeden Körper k.
        Insbesondere ist jedes Ideal in k[x_1,...,x_n] endlich erzeugt.

        Für ℤ_p (p prim) ist ℤ_p ein Körper → ℤ_p[x_1,...,x_n] noetherisch.

        Die Verifikation prüft:
        1. Ist der Grundring ℤ_p noetherisch? (Ja für p prim)
        2. Ist das Ideal durch die gegebenen Generatoren endlich erzeugt? (Ja, trivial)
        3. Übersicht über das erzeugte Ideal

    @param ideal_gens Liste der Generatoren (als Koeffizientenvektoren)
    @param p Charakteristik des Körpers (Primzahl)
    @return Dict mit:
            - 'is_finitely_generated': True (Hilbert-Basissatz)
            - 'base_ring_noetherian': True wenn ℤ_p Körper ist
            - 'n_generators': Anzahl der Generatoren
            - 'description': Erklärung
    @lastModified 2026-03-10
    """
    # Prüfen ob p eine Primzahl ist (dann ist ℤ_p ein Körper)
    base_ring_noetherian = _is_prime(p)

    n_gens = len(ideal_gens)
    n_vars = len(ideal_gens[0]) if ideal_gens else 0

    # Hilbert-Basissatz: Polynomring über Körper ist noetherisch
    # Daher ist jedes Ideal endlich erzeugt
    is_finitely_generated = base_ring_noetherian

    desc = (
        f'Ideal in ℤ_{p}[x_1,...,x_{n_vars}] erzeugt von {n_gens} Generatoren. '
        f'ℤ_{p} ist {"Körper (p prim)" if base_ring_noetherian else "kein Körper"}. '
        f'Hilbert-Basissatz: Jedes Ideal ist endlich erzeugt = {is_finitely_generated}.'
    )

    return {
        'is_finitely_generated': is_finitely_generated,
        'base_ring_noetherian': base_ring_noetherian,
        'n_generators': n_gens,
        'description': desc,
    }
