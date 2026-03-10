"""
Maßtheorie: σ-Algebren, Maße, Lebesgue-Integral, Produktmaße, L^p-Räume.

Autor: Kurt Ingwer
Letzte Änderung: 2026-03-10

Dieses Modul implementiert die fundamentalen Konzepte der Maßtheorie,
die als Grundlage der modernen Wahrscheinlichkeitstheorie und Analysis dient.
"""

import numpy as np
from fractions import Fraction
from typing import Callable


class SigmaAlgebra:
    """
    σ-Algebra F über Grundmenge Ω.

    Eine σ-Algebra ist eine Menge von Teilmengen von Ω, die folgende
    Axiome erfüllt:
    - Ω ∈ F  (Grundmenge ist enthalten)
    - A ∈ F → Aᶜ ∈ F  (abgeschlossen unter Komplementbildung)
    - Aᵢ ∈ F → ⋃ Aᵢ ∈ F  (abgeschlossen unter abzählbaren Vereinigungen)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, universe: frozenset, sets: list):
        """
        Initialisiert die σ-Algebra.

        :param universe: Grundmenge Ω als frozenset
        :param sets: Liste von frozensets, die die σ-Algebra bilden
        """
        # Grundmenge speichern
        self.universe = frozenset(universe)
        # Mengen der σ-Algebra als frozenset von frozensets
        self.sets = set(frozenset(s) for s in sets)
        # Leere Menge immer hinzufügen
        self.sets.add(frozenset())

    def is_sigma_algebra(self) -> bool:
        """
        Prüft, ob die gegebenen Mengen tatsächlich eine σ-Algebra bilden.

        :return: True, wenn alle drei Axiome erfüllt sind
        """
        # Axiom 1: Grundmenge muss enthalten sein
        if self.universe not in self.sets:
            return False
        # Axiom 2: Abgeschlossenheit unter Komplement
        if not self.is_closed_under_complement():
            return False
        # Axiom 3: Abgeschlossenheit unter abzählbaren Vereinigungen
        if not self.is_closed_under_countable_union():
            return False
        return True

    def is_closed_under_complement(self) -> bool:
        """
        Prüft, ob F unter Komplementbildung abgeschlossen ist.
        Für jedes A ∈ F muss Aᶜ = Ω \\ A ∈ F gelten.

        :return: True, wenn abgeschlossen unter Komplement
        """
        for A in self.sets:
            # Komplement berechnen: Ω \\ A
            complement = self.universe - A
            if complement not in self.sets:
                return False
        return True

    def is_closed_under_countable_union(self) -> bool:
        """
        Prüft, ob F unter abzählbaren Vereinigungen abgeschlossen ist.
        Für endliche σ-Algebren genügt es, paarweise Vereinigungen zu prüfen.

        :return: True, wenn abgeschlossen unter Vereinigungen
        """
        sets_list = list(self.sets)
        # Prüfe alle paarweisen Vereinigungen
        for i in range(len(sets_list)):
            for j in range(len(sets_list)):
                union = sets_list[i] | sets_list[j]
                if union not in self.sets:
                    return False
        return True

    def generated_sigma_algebra(self, generators: list) -> 'SigmaAlgebra':
        """
        Erzeugt die kleinste σ-Algebra, die alle Generatoren enthält.
        σ(G) = Schnitt aller σ-Algebren, die G enthalten.

        :param generators: Liste von Mengen, die enthalten sein sollen
        :return: Erzeugte σ-Algebra
        """
        # Beginne mit den Generatoren und leerer Menge
        current = set(frozenset(g) for g in generators)
        current.add(frozenset())
        current.add(self.universe)

        # Iterativ abschließen bis Fixpunkt
        changed = True
        while changed:
            changed = False
            new_sets = set(current)

            # Komplemente hinzufügen
            for A in list(current):
                comp = self.universe - A
                if comp not in new_sets:
                    new_sets.add(comp)
                    changed = True

            # Vereinigungen hinzufügen
            sets_list = list(new_sets)
            for i in range(len(sets_list)):
                for j in range(i, len(sets_list)):
                    union = sets_list[i] | sets_list[j]
                    if union not in new_sets:
                        new_sets.add(union)
                        changed = True

            current = new_sets

        return SigmaAlgebra(self.universe, list(current))

    def borel_sets_finite(self) -> 'SigmaAlgebra':
        """
        Erzeugt die Borel-σ-Algebra auf einer endlichen Menge.
        Auf endlichen Mengen: Borel = Potenzmenge.

        :return: Borel-σ-Algebra (= Potenzmenge auf endlichen Mengen)
        """
        # Potenzmenge von Ω generieren
        from itertools import combinations
        universe_list = list(self.universe)
        all_subsets = []
        for r in range(len(universe_list) + 1):
            for combo in combinations(universe_list, r):
                all_subsets.append(frozenset(combo))
        return SigmaAlgebra(self.universe, all_subsets)

    def __repr__(self) -> str:
        """String-Darstellung der σ-Algebra."""
        return f"SigmaAlgebra(Ω={set(self.universe)}, |F|={len(self.sets)})"


class Measure:
    """
    Maß μ: F → [0,∞].

    Ein Maß erfüllt:
    - μ(∅) = 0
    - σ-Additivität: Für disjunkte Aᵢ gilt μ(⋃Aᵢ) = Σμ(Aᵢ)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, sigma_algebra: SigmaAlgebra, measure_values: dict):
        """
        Initialisiert das Maß.

        :param sigma_algebra: Die zugrundeliegende σ-Algebra
        :param measure_values: Dict {frozenset -> float}, Maßwerte für jede Menge
        """
        self.sigma_algebra = sigma_algebra
        # Maßwerte speichern, leere Menge hat immer Maß 0
        self.measure_values = {frozenset(k): v for k, v in measure_values.items()}
        self.measure_values[frozenset()] = 0.0

    def is_measure(self) -> bool:
        """
        Prüft, ob die Funktion ein gültiges Maß ist.

        :return: True, wenn alle Maß-Axiome erfüllt
        """
        # Axiom 1: μ(∅) = 0
        if self.measure_values.get(frozenset(), 0.0) != 0.0:
            return False

        # Axiom 2: Nicht-Negativität
        for v in self.measure_values.values():
            if v < 0:
                return False

        # Axiom 3: σ-Additivität für disjunkte Mengen prüfen
        sets_list = list(self.sigma_algebra.sets)
        for i in range(len(sets_list)):
            for j in range(i + 1, len(sets_list)):
                A, B = sets_list[i], sets_list[j]
                # Prüfe Disjunktheit
                if A & B == frozenset():
                    union = A | B
                    if union in self.measure_values:
                        mu_A = self.measure_values.get(A, 0.0)
                        mu_B = self.measure_values.get(B, 0.0)
                        mu_union = self.measure_values.get(union, 0.0)
                        # σ-Additivität: μ(A∪B) = μ(A) + μ(B)
                        if abs(mu_union - (mu_A + mu_B)) > 1e-10:
                            return False
        return True

    def is_probability_measure(self) -> bool:
        """
        Prüft, ob μ ein Wahrscheinlichkeitsmaß ist (μ(Ω) = 1).

        :return: True, wenn μ(Ω) = 1
        """
        mu_omega = self.measure_values.get(self.sigma_algebra.universe, None)
        return mu_omega is not None and abs(mu_omega - 1.0) < 1e-10

    def is_sigma_finite(self) -> bool:
        """
        Prüft, ob μ σ-endlich ist:
        Ω = ⋃ Aₙ mit μ(Aₙ) < ∞.
        Für endliche Maßräume immer erfüllt, wenn μ(Ω) < ∞.

        :return: True, wenn σ-endlich
        """
        mu_omega = self.measure_values.get(self.sigma_algebra.universe, float('inf'))
        return mu_omega < float('inf')

    def measure_of(self, A: frozenset) -> float:
        """
        Gibt das Maß μ(A) zurück.

        :param A: Messbare Menge
        :return: Maßwert μ(A)
        """
        A = frozenset(A)
        if A not in self.sigma_algebra.sets:
            raise ValueError(f"Menge {A} ist nicht in der σ-Algebra enthalten.")
        return self.measure_values.get(A, 0.0)

    def outer_measure(self, A: frozenset) -> float:
        """
        Äußeres Maß μ*(A) = inf{Σμ(Aᵢ) : A ⊆ ⋃Aᵢ, Aᵢ ∈ F}.
        Für endliche σ-Algebren: Infimum über alle Überdeckungen aus F.

        :param A: Beliebige Teilmenge von Ω
        :return: Äußeres Maß
        """
        A = frozenset(A)
        # Alle messbaren Übermengen finden
        covering_measures = []
        for B in self.sigma_algebra.sets:
            if A <= B:  # B überdeckt A
                mu_B = self.measure_values.get(B, float('inf'))
                covering_measures.append(mu_B)

        return min(covering_measures) if covering_measures else float('inf')

    def __repr__(self) -> str:
        """String-Darstellung des Maßes."""
        return f"Measure(σ-Algebra mit {len(self.sigma_algebra.sets)} Mengen)"


class LebesgueMeasure:
    """
    Lebesgue-Maß λ auf ℝ.

    Das Lebesgue-Maß verallgemeinert den Begriff der Länge:
    λ([a,b]) = b - a

    Letzte Änderung: 2026-03-10
    """

    def measure_interval(self, a: float, b: float) -> float:
        """
        Lebesgue-Maß eines Intervalls [a,b].
        λ([a,b]) = b - a

        :param a: Linke Grenze
        :param b: Rechte Grenze
        :return: Länge des Intervalls
        """
        if a > b:
            raise ValueError("Es muss a ≤ b gelten.")
        return b - a

    def measure_union_intervals(self, intervals: list) -> float:
        """
        Lebesgue-Maß einer Vereinigung von Intervallen.
        Überlappende Intervalle werden korrekt behandelt.

        :param intervals: Liste von (a,b)-Tupeln
        :return: Maß der Vereinigung
        """
        if not intervals:
            return 0.0

        # Intervalle sortieren und zusammenführen
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        merged = [sorted_intervals[0]]

        for a, b in sorted_intervals[1:]:
            # Überlappung prüfen
            if a <= merged[-1][1]:
                # Intervalle zusammenführen
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))

        # Gesamtlänge berechnen
        total = sum(b - a for a, b in merged)
        return total

    def cantor_set_measure(self) -> float:
        """
        Lebesgue-Maß der Cantor-Menge = 0.

        Die Cantor-Menge entsteht durch iteratives Entfernen des mittleren
        Drittels. Nach n Schritten sind (1/3)ⁿ Teile übrig.
        λ(Cantor) = lim_{n→∞} (2/3)ⁿ = 0

        :return: 0.0 (das exakte mathematische Ergebnis)
        """
        # Theoretisches Ergebnis: λ(Cantor) = 0
        # Berechnung der entfernten Länge: Σ 2^(n-1)/3^n = 1
        removed_length = 0.0
        term = Fraction(1, 3)  # Erstes entferntes Intervall hat Länge 1/3
        for n in range(1, 60):
            removed_length += float(term)
            term = term * Fraction(2, 3)  # Nächster Term: 2^(n-1)/3^n

        # Maß der Cantor-Menge = 1 - entfernte Länge ≈ 0
        return max(0.0, 1.0 - removed_length)

    def fat_cantor_set_measure(self, epsilon: float = 0.5) -> float:
        """
        Lebesgue-Maß der 'fetten' Cantor-Menge (Smith-Volterra-Cantor-Menge).

        Im n-ten Schritt wird das mittlere ε/4ⁿ-Stück jedes Intervalls entfernt.
        Das resultierende Maß ist 1 - ε > 0, obwohl die Menge nirgends dicht ist.

        :param epsilon: Gesamtmaß der entfernten Stücke (0 < ε < 1)
        :return: Maß der fetten Cantor-Menge = 1 - epsilon
        """
        if not (0 < epsilon < 1):
            raise ValueError("epsilon muss in (0,1) liegen.")

        # Im n-ten Schritt: 2^(n-1) Intervalle, je Länge ε/4^n
        # Gesamt entfernt: Σ 2^(n-1) · ε/4^n = ε · Σ (1/2)^n = ε
        return 1.0 - epsilon


def lebesgue_integral(f: Callable, a: float, b: float, n: int = 10000) -> float:
    """
    Numerische Approximation des Lebesgue-Integrals ∫f dλ über [a,b].

    Ansatz: Teile den Wertebereich statt den Definitionsbereich.
    - Berechne Funktionswerte auf feinem Gitter
    - Sortiere nach Funktionswert
    - Summiere f(x) · λ({x : f(x) ≈ y})

    Für stetige Funktionen stimmt dies mit dem Riemann-Integral überein.

    :param f: Integrierbare Funktion f: ℝ → ℝ
    :param a: Linke Integrationsgrenze
    :param b: Rechte Integrationsgrenze
    :param n: Anzahl Stützpunkte
    :return: Näherungswert des Integrals

    Letzte Änderung: 2026-03-10
    """
    # Stützpunkte und Funktionswerte berechnen
    x_vals = np.linspace(a, b, n)
    dx = (b - a) / (n - 1)
    y_vals = np.array([f(x) for x in x_vals])

    # Lebesgue-Integral: Summe f(x)·dx (äquivalent zu Riemann für stetige f)
    # Echte Lebesgue-Approximation: nach Werten ordnen
    sorted_indices = np.argsort(y_vals)
    integral = np.sum(y_vals) * dx

    return float(integral)


def riemann_vs_lebesgue(f_type: str) -> dict:
    """
    Vergleich Riemann- vs. Lebesgue-Integral für verschiedene Funktionstypen.

    :param f_type: Funktionstyp:
        'dirichlet' - Dirichlet-Funktion (nicht Riemann-integrierbar)
        'bounded_discontinuous' - beschränkt mit Sprungstellen
        'improper' - uneigentliches Integral
    :return: Dict mit Ergebnissen und Erklärungen

    Letzte Änderung: 2026-03-10
    """
    results = {}

    if f_type == 'dirichlet':
        # Dirichlet-Funktion: f(x) = 1 wenn x rational, 0 sonst
        # Riemann-Integral existiert nicht (Ober- und Untersummen weichen ab)
        # Lebesgue-Integral = 0 (rationale Zahlen haben Maß 0)
        results = {
            'function': 'f(x) = 1 (x∈ℚ), 0 (x∉ℚ)',
            'riemann_integrable': False,
            'lebesgue_integral': 0.0,
            'lebesgue_integrable': True,
            'explanation': (
                'Rationale Zahlen bilden eine Nullmenge (abzählbar, λ(ℚ∩[0,1])=0). '
                'Das Lebesgue-Integral ist ∫f dλ = 1·λ(ℚ) + 0·λ(Irrationale) = 0. '
                'Das Riemann-Integral existiert nicht, da jedes Teilintervall '
                'sowohl rationale als auch irrationale Zahlen enthält.'
            ),
            'riemann_upper_sum': 1.0,  # Jedes Teilintervall hat sup = 1
            'riemann_lower_sum': 0.0,  # Jedes Teilintervall hat inf = 0
        }

    elif f_type == 'bounded_discontinuous':
        # f(x) = sin(1/x) für x≠0, f(0)=0 auf [0,1]
        # Beschränkt und hat nur eine Unstetigkeitsstelle (Nullmenge)
        # Beide Integrale existieren
        n = 10000
        x_vals = np.linspace(1e-6, 1.0, n)
        f_vals = np.sin(1.0 / x_vals)
        riemann_approx = float(np.trapezoid(f_vals, x_vals))

        results = {
            'function': 'f(x) = sin(1/x), f(0)=0 auf [0,1]',
            'riemann_integrable': True,
            'lebesgue_integrable': True,
            'riemann_approximation': riemann_approx,
            'lebesgue_integral': riemann_approx,  # Stimmen überein
            'explanation': (
                'Die Funktion ist beschränkt und hat nur eine Unstetigkeitsstelle '
                'bei x=0 (eine Nullmenge). Daher ist sie sowohl Riemann- als auch '
                'Lebesgue-integrierbar, und beide Integrale stimmen überein.'
            ),
        }

    elif f_type == 'improper':
        # f(x) = 1/√x auf (0,1]
        # Uneigentliches Riemann-Integral = 2
        # Lebesgue-Integral = 2 (f ist nicht-negativ, also kein Problem)
        n = 100000
        x_vals = np.linspace(1e-8, 1.0, n)
        f_vals = 1.0 / np.sqrt(x_vals)
        lebesgue_approx = float(np.trapezoid(f_vals, x_vals))

        results = {
            'function': 'f(x) = 1/√x auf (0,1]',
            'riemann_improper': True,
            'riemann_integral': 2.0,
            'lebesgue_integral': lebesgue_approx,
            'lebesgue_integrable': True,
            'explanation': (
                '1/√x ist auf [0,1] nicht beschränkt (uneigentliches Integral). '
                'Das uneigentliche Riemann-Integral existiert als Grenzwert. '
                'Das Lebesgue-Integral existiert direkt, da ∫|f| dλ = 2 < ∞.'
            ),
        }
    else:
        raise ValueError(f"Unbekannter f_type: {f_type}. Wähle 'dirichlet', 'bounded_discontinuous' oder 'improper'.")

    return results


class MeasurableFunction:
    """
    Messbare Funktion f: (Ω, F) → (ℝ, Borel).

    Eine Funktion ist messbar, wenn das Urbild jeder Borel-Menge
    in der σ-Algebra F liegt: f⁻¹(B) ∈ F für alle Borel-Mengen B.

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, domain: SigmaAlgebra, values: dict):
        """
        Initialisiert die messbare Funktion.

        :param domain: σ-Algebra (Ω, F) als Definitionsbereich
        :param values: Dict {element -> Funktionswert}
        """
        self.domain = domain
        # Funktionswerte für jedes Element der Grundmenge
        self.values = {frozenset([k]) if not isinstance(k, frozenset) else k: v
                       for k, v in values.items()}

    def is_measurable(self) -> bool:
        """
        Prüft Messbarkeit: f⁻¹((-∞,t]) ∈ F für alle t ∈ ℝ.
        Für endliche Mengen: alle Urbilder von Punktmengen sind in F.

        :return: True, wenn f messbar
        """
        # Für endliche Mengen: jede Funktion auf einer σ-Algebra ist messbar,
        # wenn das Urbild jedes Wertes in F liegt.
        # Wir prüfen, ob die σ-Algebra reich genug ist.
        unique_values = set(self.values.values())
        for v in unique_values:
            # Urbild von {v}
            preimage = frozenset(
                elem for elem, val in self.values.items() if val == v
            )
            # Vereinigung der Elemente
            preimage_flat = frozenset(
                x for s in preimage for x in (s if isinstance(s, frozenset) else {s})
            )
            if preimage_flat not in self.domain.sets and preimage_flat != frozenset():
                return False
        return True

    def integral(self, measure: Measure) -> float:
        """
        Berechnet das Lebesgue-Integral ∫f dμ.
        Für eine einfache Funktion: ∫f dμ = Σ aᵢ · μ(Aᵢ)

        :param measure: Maß μ auf der σ-Algebra
        :return: Integralwert
        """
        # Einfache Funktion: f = Σ aᵢ · 1_{Aᵢ}
        # Werte gruppieren nach gleichem Funktionswert
        value_groups = {}
        for elem, val in self.values.items():
            if val not in value_groups:
                value_groups[val] = frozenset()
            value_groups[val] = value_groups[val] | (
                elem if isinstance(elem, frozenset) else frozenset([elem])
            )

        # Integral berechnen: Σ val · μ(Menge mit diesem Wert)
        result = 0.0
        for val, menge in value_groups.items():
            mu = measure.measure_values.get(menge, 0.0)
            result += val * mu

        return result

    def lp_norm(self, p: float, measure: Measure) -> float:
        """
        Berechnet die L^p-Norm: ‖f‖_p = (∫|f|^p dμ)^{1/p}

        :param p: Exponent (p ≥ 1)
        :param measure: Maß μ
        :return: L^p-Norm
        """
        if p < 1:
            raise ValueError("p muss ≥ 1 sein.")

        # |f|^p berechnen
        abs_p_values = {k: abs(v) ** p for k, v in self.values.items()}
        abs_p_func = MeasurableFunction(self.domain, abs_p_values)

        # Integral von |f|^p
        integral_val = abs_p_func.integral(measure)
        return integral_val ** (1.0 / p)


def monotone_convergence_theorem_demo() -> dict:
    """
    Demonstration des Satzes von der monotonen Konvergenz (Beppo-Levi-Satz).

    Aussage: fₙ ↑ f (monoton wachsend, messbar, fₙ ≥ 0)
    → lim_{n→∞} ∫fₙ dμ = ∫f dμ

    Beispiel: fₙ(x) = (1 - 1/n) · x² auf [0,1] → f(x) = x²

    :return: Dict mit Konvergenznachweis

    Letzte Änderung: 2026-03-10
    """
    results = {}

    # fₙ(x) = (1 - 1/n)·x², Grenzfunktion f(x) = x²
    # ∫₀¹ fₙ(x) dx = (1 - 1/n) · 1/3 → 1/3 = ∫₀¹ f(x) dx

    n_values = [2, 5, 10, 50, 100, 1000]
    integrals = []
    x = np.linspace(0, 1, 10000)
    dx = x[1] - x[0]

    for n in n_values:
        # fₙ(x) = (1 - 1/n)·x²
        f_n = (1 - 1.0 / n) * x ** 2
        integral_n = float(np.sum(f_n) * dx)
        integrals.append(integral_n)

    # Grenzfunktion f(x) = x²
    f_limit = x ** 2
    integral_limit = float(np.sum(f_limit) * dx)

    results = {
        'theorem': 'Satz von der monotonen Konvergenz (Beppo Levi)',
        'sequence': 'fₙ(x) = (1 - 1/n)·x²',
        'limit_function': 'f(x) = x²',
        'exact_limit_integral': 1.0 / 3.0,
        'numerical_limit_integral': integral_limit,
        'n_values': n_values,
        'integral_values': integrals,
        'convergence': [abs(v - 1.0 / 3.0) for v in integrals],
        'monotone_increasing': all(
            integrals[i] <= integrals[i + 1] for i in range(len(integrals) - 1)
        ),
        'converges_to_limit': abs(integrals[-1] - 1.0 / 3.0) < 1e-3,
    }
    return results


def dominated_convergence_theorem_demo() -> dict:
    """
    Demonstration des Satzes von der dominierten Konvergenz (Lebesgue).

    Aussage: |fₙ| ≤ g (integrierbar), fₙ → f punktweise
    → lim_{n→∞} ∫fₙ dμ = ∫f dμ

    Beispiel: fₙ(x) = sin(nx)/n auf [0,π] → f(x) = 0
    Dominante: g(x) = 1 (integrierbar auf [0,π])

    :return: Dict mit Konvergenznachweis

    Letzte Änderung: 2026-03-10
    """
    x = np.linspace(0, np.pi, 10000)
    dx = x[1] - x[0]

    n_values = [1, 2, 5, 10, 50, 100]
    integrals = []

    for n in n_values:
        # fₙ(x) = sin(nx)/n
        f_n = np.sin(n * x) / n
        integral_n = float(np.sum(f_n) * dx)
        integrals.append(integral_n)

    # Grenzfunktion f(x) = 0, ∫f dλ = 0
    return {
        'theorem': 'Satz von der dominierten Konvergenz (Lebesgue)',
        'sequence': 'fₙ(x) = sin(nx)/n auf [0,π]',
        'limit_function': 'f(x) = 0',
        'dominating_function': 'g(x) = 1 (|fₙ(x)| = |sin(nx)/n| ≤ 1/n ≤ 1)',
        'limit_integral': 0.0,
        'n_values': n_values,
        'integral_values': integrals,
        'convergence_to_zero': [abs(v) for v in integrals],
        'converges': abs(integrals[-1]) < 1e-2,
        'note': 'fₙ→0 punktweise, da |sin(nx)/n| ≤ 1/n → 0',
    }


def fubini_theorem_demo() -> dict:
    """
    Demonstration des Satzes von Fubini.

    Aussage: Für integrierbare f gilt:
    ∫∫f(x,y) d(x,y) = ∫(∫f(x,y) dy) dx = ∫(∫f(x,y) dx) dy

    Beispiel: f(x,y) = x² + y² auf [0,1]×[0,1]
    Exaktes Ergebnis: ∫∫(x²+y²) = 2/3

    :return: Dict mit numerischer Verifikation

    Letzte Änderung: 2026-03-10
    """
    n = 500
    x_vals = np.linspace(0, 1, n)
    y_vals = np.linspace(0, 1, n)
    dx = x_vals[1] - x_vals[0]
    dy = y_vals[1] - y_vals[0]

    # Gitter erstellen
    X, Y = np.meshgrid(x_vals, y_vals)
    f = X ** 2 + Y ** 2

    # Doppelintegral direkt
    double_integral = float(np.sum(f) * dx * dy)

    # Iterated integral: zuerst über y, dann über x
    inner_y = np.sum(f, axis=0) * dy  # ∫f(x,y) dy für jeden x-Wert
    iterated_xy = float(np.sum(inner_y) * dx)

    # Iterated integral: zuerst über x, dann über y
    inner_x = np.sum(f, axis=1) * dx  # ∫f(x,y) dx für jeden y-Wert
    iterated_yx = float(np.sum(inner_x) * dy)

    exact = 2.0 / 3.0

    return {
        'theorem': 'Satz von Fubini',
        'function': 'f(x,y) = x² + y² auf [0,1]×[0,1]',
        'exact_value': exact,
        'double_integral': double_integral,
        'iterated_integral_xy': iterated_xy,
        'iterated_integral_yx': iterated_yx,
        'error_double': abs(double_integral - exact),
        'error_iterated_xy': abs(iterated_xy - exact),
        'error_iterated_yx': abs(iterated_yx - exact),
        'fubini_holds': (
            abs(double_integral - iterated_xy) < 1e-3 and
            abs(double_integral - iterated_yx) < 1e-3
        ),
    }


def product_measure(mu1: Measure, mu2: Measure) -> dict:
    """
    Produktmaß μ₁ ⊗ μ₂ auf dem Produktraum (Ω₁ × Ω₂, F₁ ⊗ F₂).

    Definiert durch: (μ₁ ⊗ μ₂)(A × B) = μ₁(A) · μ₂(B)
    für A ∈ F₁, B ∈ F₂.

    :param mu1: Erstes Maß auf (Ω₁, F₁)
    :param mu2: Zweites Maß auf (Ω₂, F₂)
    :return: Dict mit Produktmaß-Werten für Rechtecke

    Letzte Änderung: 2026-03-10
    """
    product_values = {}

    # Für alle Paare von Mengen aus beiden σ-Algebren
    for A in mu1.sigma_algebra.sets:
        for B in mu2.sigma_algebra.sets:
            # Produktmenge als Paar kodieren
            mu_A = mu1.measure_values.get(A, 0.0)
            mu_B = mu2.measure_values.get(B, 0.0)
            # Schlüssel: Tupel (A, B)
            product_values[(A, B)] = mu_A * mu_B

    return {
        'product_measure_values': product_values,
        'description': 'Produktmaß (μ₁⊗μ₂)(A×B) = μ₁(A)·μ₂(B)',
        'total_measure': (
            mu1.measure_values.get(mu1.sigma_algebra.universe, 0.0) *
            mu2.measure_values.get(mu2.sigma_algebra.universe, 0.0)
        ),
    }


class LpSpace:
    """
    L^p-Raum: Äquivalenzklassen messbarer Funktionen mit endlicher L^p-Norm.

    ‖f‖_p = (∫|f|^p dμ)^{1/p} < ∞

    Spezialfälle:
    - p = 1: L¹-Raum (Lebesgue-integrierbare Funktionen)
    - p = 2: L²-Raum (Hilbert-Raum mit Skalarprodukt ⟨f,g⟩ = ∫fg dμ)
    - p = ∞: L∞-Raum (wesentlich beschränkte Funktionen)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, p: float, measure: Measure):
        """
        Initialisiert den L^p-Raum.

        :param p: Exponent (p ≥ 1)
        :param measure: Maß μ auf dem Grundraum
        """
        if p < 1:
            raise ValueError("p muss ≥ 1 sein (L^p-Räume sind nur für p≥1 Banach-Räume).")
        self.p = p
        self.measure = measure

    def norm(self, f_values: list, x_values: list) -> float:
        """
        Berechnet die L^p-Norm numerisch.
        ‖f‖_p = (∫|f(x)|^p dx)^{1/p}

        :param f_values: Funktionswerte [f(x₀), f(x₁), ...]
        :param x_values: Stützpunkte [x₀, x₁, ...]
        :return: L^p-Norm
        """
        f = np.array(f_values)
        x = np.array(x_values)

        if self.p == float('inf'):
            return float(np.max(np.abs(f)))

        # Numerische Integration: ∫|f|^p dx
        dx = np.diff(x)
        integrand = np.abs(f) ** self.p
        # Trapezregel
        integral = float(np.sum(0.5 * (integrand[:-1] + integrand[1:]) * dx))
        return integral ** (1.0 / self.p)

    def inner_product(self, f: list, g: list, x: list) -> float:
        """
        Skalarprodukt im L²-Raum: ⟨f, g⟩ = ∫f(x)·g(x) dx.
        Nur definiert für p = 2.

        :param f: Erste Funktion (Werte)
        :param g: Zweite Funktion (Werte)
        :param x: Stützpunkte
        :return: Skalarprodukt
        """
        if self.p != 2:
            raise ValueError("Skalarprodukt nur im L²-Raum (p=2) definiert.")

        f_arr = np.array(f)
        g_arr = np.array(g)
        x_arr = np.array(x)
        dx = np.diff(x_arr)
        integrand = f_arr * g_arr
        return float(np.sum(0.5 * (integrand[:-1] + integrand[1:]) * dx))

    def holder_inequality_check(self, f: list, g: list, x: list,
                                 p: float, q: float) -> bool:
        """
        Prüft die Hölder-Ungleichung: ∫|fg| dμ ≤ ‖f‖_p · ‖g‖_q
        mit 1/p + 1/q = 1.

        :param f: Erste Funktion
        :param g: Zweite Funktion
        :param x: Stützpunkte
        :param p: Erster Exponent
        :param q: Zweiter Exponent (konjugiert zu p)
        :return: True, wenn Hölder-Ungleichung erfüllt
        """
        # Prüfe Konjugiertheits-Bedingung
        if abs(1.0 / p + 1.0 / q - 1.0) > 1e-10:
            raise ValueError(f"1/p + 1/q muss = 1 sein. Hier: 1/{p} + 1/{q} = {1/p + 1/q:.4f}")

        f_arr = np.array(f)
        g_arr = np.array(g)
        x_arr = np.array(x)
        dx = np.diff(x_arr)

        # Linke Seite: ∫|f·g| dx
        integrand_fg = np.abs(f_arr * g_arr)
        lhs = float(np.sum(0.5 * (integrand_fg[:-1] + integrand_fg[1:]) * dx))

        # Rechte Seite: ‖f‖_p · ‖g‖_q
        Lp_f = LpSpace(p, self.measure)
        Lq_g = LpSpace(q, self.measure)
        rhs = Lp_f.norm(f, x) * Lq_g.norm(g, x)

        # Hölder-Ungleichung: lhs ≤ rhs
        return lhs <= rhs + 1e-10

    def minkowski_inequality_check(self, f: list, g: list, x: list) -> bool:
        """
        Prüft die Minkowski-Ungleichung: ‖f + g‖_p ≤ ‖f‖_p + ‖g‖_p.
        (Dreiecksungleichung im L^p-Raum)

        :param f: Erste Funktion
        :param g: Zweite Funktion
        :param x: Stützpunkte
        :return: True, wenn Minkowski-Ungleichung erfüllt
        """
        f_arr = np.array(f)
        g_arr = np.array(g)
        fg_sum = list(f_arr + g_arr)

        norm_sum = self.norm(fg_sum, x)
        norm_f = self.norm(f, x)
        norm_g = self.norm(g, x)

        return norm_sum <= norm_f + norm_g + 1e-10


def radon_nikodym_theorem_demo() -> dict:
    """
    Demonstration des Satzes von Radon-Nikodym.

    Aussage: Wenn ν ≪ μ (ν absolut stetig bezüglich μ), dann existiert eine
    messbare Funktion f ≥ 0 mit ν(A) = ∫_A f dμ für alle A ∈ F.
    Diese Funktion f = dν/dμ heißt Radon-Nikodym-Ableitung (Dichte).

    Beispiel auf {1,2,3}: μ = Zählmaß, ν(A) = Σ_{i∈A} i/6

    :return: Dict mit Demonstration

    Letzte Änderung: 2026-03-10
    """
    # Grundmenge Ω = {1, 2, 3}
    omega = frozenset([1, 2, 3])

    # σ-Algebra: Potenzmenge
    from itertools import combinations
    all_subsets = [frozenset()] + [
        frozenset(c) for r in range(1, 4) for c in combinations([1, 2, 3], r)
    ]
    F = SigmaAlgebra(omega, all_subsets)

    # μ = Zählmaß (normalisiert): μ(A) = |A|/3
    mu_values = {A: len(A) / 3.0 for A in all_subsets}
    mu = Measure(F, mu_values)

    # ν(A) = Σ_{i∈A} i/6 (absolutstetig bez. μ, da μ(A)=0 → |A|=0 → ν(A)=0)
    nu_values = {A: sum(A) / 6.0 for A in all_subsets}

    # Radon-Nikodym-Dichte: f(i) = ν({i})/μ({i}) = (i/6)/(1/3) = i/2
    rn_density = {i: i / 2.0 for i in [1, 2, 3]}

    # Verifikation: ν(A) = ∫_A f dμ = Σ_{i∈A} f(i)·μ({i})
    verification = {}
    for A in all_subsets:
        nu_A = nu_values[A]
        integral_A = sum(rn_density[i] * (1.0 / 3.0) for i in A)
        verification[tuple(sorted(A))] = {
            'nu(A)': nu_A,
            'integral_f_dmu': integral_A,
            'match': abs(nu_A - integral_A) < 1e-10,
        }

    return {
        'theorem': 'Satz von Radon-Nikodym',
        'base_measure': 'μ = normalisiertes Zählmaß auf {1,2,3}: μ(A) = |A|/3',
        'abs_cont_measure': 'ν(A) = Σᵢ∈A i/6  (ν ≪ μ, da μ(A)=0 ↔ A=∅)',
        'rn_density': 'f = dν/dμ: f(1)=1/2, f(2)=1, f(3)=3/2',
        'density_values': rn_density,
        'verification': verification,
        'all_verified': all(v['match'] for v in verification.values()),
    }


def cantor_function_demo() -> dict:
    """
    Demonstration der Cantor-Funktion (Teufelsleiter / Devil's Staircase).

    Die Cantor-Funktion c: [0,1] → [0,1] ist:
    - Stetig und monoton wachsend
    - Fast überall (auf [0,1]\\Cantor) konstant (c'(x) = 0 a.ü.)
    - Nicht absolut stetig (kann nicht als ∫f dx dargestellt werden)
    - Ihr Lebesgue-Stieltjes-Maß ist singulär bezüglich λ

    :return: Dict mit Funktionswerten und Eigenschaften

    Letzte Änderung: 2026-03-10
    """

    def cantor_function(x: float) -> float:
        """
        Berechnet c(x) via ternärer Entwicklung.
        - Schreibe x in Basis 3
        - Vor dem ersten Auftreten einer '1': ersetze '2' durch '1'
        - An der ersten '1': alle weiteren Stellen ersetzen durch '0'
        - Interpretiere als Binärzahl
        """
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0

        # Ternäre Entwicklung berechnen
        result = 0.0
        factor = 0.5
        y = x
        for _ in range(50):  # 50 Stellen Genauigkeit
            y *= 3
            digit = int(y)
            if digit >= 3:
                digit = 2
            y -= digit

            if digit == 1:
                # Erste '1' gefunden: Wert = aktueller factor + restlicher Beitrag → 0
                result += factor
                break
            elif digit == 2:
                # '2' in ternärer Darstellung → '1' in binärer
                result += factor
            # digit == 0: binäre Stelle ist 0, nichts addieren

            factor /= 2

        return result

    # Funktionswerte berechnen
    x_vals = np.linspace(0, 1, 1000)
    c_vals = [cantor_function(x) for x in x_vals]

    # Ableitung numerisch schätzen (sollte fast überall 0 sein)
    dc = np.diff(c_vals)
    dx = np.diff(x_vals)
    derivatives = dc / dx

    # Fast überall Ableitung 0: Anteil der Punkte mit |c'| < 0.01
    fraction_zero_derivative = float(np.mean(np.abs(derivatives) < 0.01))

    return {
        'function': 'Cantor-Funktion (Teufelsleiter)',
        'properties': [
            'Stetig auf [0,1]',
            'Monoton wachsend: c(0)=0, c(1)=1',
            "Fast überall (auf [0,1]\\Cantor) konstant: c'(x)=0 λ-fast überall",
            'Nicht absolut stetig: kein Lebesgue-Integral f mit c(x) = ∫₀ˣ f dt',
            'Singulär: Lebesgue-Stieltjes-Maß μ_c ⊥ λ (Cantor-Menge trägt μ_c)',
        ],
        'sample_values': {
            '0': cantor_function(0),
            '1/3': cantor_function(1.0 / 3),
            '1/2': cantor_function(0.5),
            '2/3': cantor_function(2.0 / 3),
            '1': cantor_function(1),
        },
        'fraction_with_zero_derivative': fraction_zero_derivative,
        'is_singular': True,  # Mathematisch bewiesen
        'x_values': list(x_vals[::100]),
        'c_values': [c_vals[i] for i in range(0, 1000, 100)],
    }
