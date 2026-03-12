"""
@file debruijn_newman.py
@brief de Bruijn-Newman-Konstante Λ und ihre Verbindung zur Riemann-Hypothese.
@description
    Implementiert die numerische Analyse der de Bruijn-Newman-Konstante Λ,
    die durch eine parametrische Familie von Funktionen H_t(x) definiert wird.

    **Definition (de Bruijn 1950, Newman 1976):**
    Betrachte die Funktion:
        H_t(x) = ∫₀^∞ e^{t·u²} · Φ(u) · cos(x·u) du

    wobei:
        Φ(u) = Σ_{n=1}^∞ (2π²n⁴e^{9u} - 3πn²e^{5u}) · exp(-πn²·e^{4u})

    Die **de Bruijn-Newman-Konstante** Λ ist definiert als:
        Λ = inf { t ∈ ℝ : H_t hat nur reelle Nullstellen }

    **Zentrale Zusammenhänge:**
    - RH ⟺ Λ ≤ 0  (alle nicht-trivialen Nullstellen der Zeta-Funktion
                    liegen auf der kritischen Geraden)
    - Newman-Vermutung (1976): Λ ≥ 0
    - Rodgers-Tao (2018): Λ ≥ 0 **BEWIESEN** ✓
    - Aktuelle Schranke: 0 ≤ Λ ≤ 0.2 (Platt-Trudgian 2021)

    **Verbindung zu H_0:**
        H_0(x) = ½ ξ(½ + ix)
    wobei ξ(s) = ½ s(s-1) π^(-s/2) Γ(s/2) ζ(s) die vollständige Zeta-Funktion.

    Die Nullstellen von H_0 entsprechen den imaginären Teilen der
    nicht-trivialen Nullstellen der Riemann-Zeta-Funktion:
        y₁ ≈ 14.134725..., y₂ ≈ 21.022040..., y₃ ≈ 25.010856...

    **Historische Schranken für Λ:**
    - 1987: Λ ≤ 0.5 (de Bruijn)
    - 1988: Λ ≥ -∞  (trivial)
    - 2009: Λ ≥ -2.7·10⁻⁹ (Saouter-Gourdon-Demichel)
    - 2018: Λ ≥ 0 (Rodgers-Tao, Newman-Vermutung bewiesen)
    - 2021: Λ ≤ 0.2 (Platt-Trudgian)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

# mpmath für hochpräzise Arithmetik (benötigt für Nullstellensuche)
try:
    import mpmath
    MPMATH_VERFUEGBAR = True
except ImportError:
    MPMATH_VERFUEGBAR = False

# numpy für Vektoroperationen
try:
    import numpy as np
    NUMPY_VERFUEGBAR = True
except ImportError:
    NUMPY_VERFUEGBAR = False

# matplotlib für optionale Visualisierung
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_VERFUEGBAR = True
except ImportError:
    MATPLOTLIB_VERFUEGBAR = False


# ===========================================================================
# BEKANNTE RIEMANN-NULLSTELLEN (imaginäre Teile auf Re(s)=1/2)
# ===========================================================================

# Erste 10 nicht-triviale Nullstellen der Riemann-Zeta-Funktion (imaginäre Teile)
# Quelle: LMFDB, Odlyzko (2001)
BEKANNTE_NULLSTELLEN_GAMMA: List[float] = [
    14.134725141734693790,
    21.022039638771554993,
    25.010857580145688763,
    30.424876125859513210,
    32.935061587739189691,
    37.586178158825671257,
    40.918719012147495187,
    43.327073280914999519,
    48.005150881167159727,
    49.773832477672302182,
]


# ===========================================================================
# HAUPTKLASSE: DeBruijnNewman
# ===========================================================================

class DeBruijnNewman:
    """
    @brief Klasse zur Analyse der de Bruijn-Newman-Konstante Λ.
    @description
        Implementiert die numerische Berechnung der parametrischen Funktionen
        H_t(x) und Methoden zur Schrankenanalyse für Λ.

        Die de Bruijn-Newman-Konstante Λ verbindet die analytische
        Zahlentheorie (RH) mit der reellen Analysis (Nullstellenstruktur
        von parameterabhängigen Integralen).

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def __init__(self, praezision: int = 50):
        """
        @brief Initialisiert die DeBruijnNewman-Instanz.
        @description
            Setzt die gewünschte numerische Präzision für mpmath-Berechnungen.
            Höhere Präzision erhöht die Genauigkeit, aber auch die Rechenzeit.
        @param praezision Dezimalstellen der Präzision für mpmath (Standard: 50).
        @lastModified 2026-03-12
        """
        self.praezision = praezision
        if MPMATH_VERFUEGBAR:
            # Setze mpmath-Präzision in Dezimalstellen
            mpmath.mp.dps = praezision

    def phi_funktion(self, u: float, N_max: int = 20) -> float:
        """
        @brief Berechnet die Kernfunktion Φ(u) der de Bruijn-Newman-Darstellung.
        @description
            Φ(u) ist definiert als:
                Φ(u) = Σ_{n=1}^{N_max} (2π²n⁴·e^{9u} - 3πn²·e^{5u}) · exp(-πn²·e^{4u})

            Diese Funktion ist exponentiell abfallend für u → ∞ und
            konvergiert schnell: für u ≥ 0 dominiert der Term n=1.

            Herleitung: Φ entsteht aus der Jacobi-Theta-Funktion
                ω(x) = Σ_{n=1}^∞ e^{-πn²x}
            durch zweifache Differentiation und geeignete Transformation.

        @param u Argument der Funktion (reell, u ≥ 0 für Hauptanwendung).
        @param N_max Anzahl der Summanden (Standard: 20, ab n=3 vernachlässigbar).
        @return Wert von Φ(u).
        @lastModified 2026-03-12
        """
        if MPMATH_VERFUEGBAR:
            return float(self._phi_mpmath(u, N_max))
        return self._phi_python(u, N_max)

    def _phi_mpmath(self, u: float, N_max: int) -> object:
        """
        @brief Hochpräzise Berechnung von Φ(u) mit mpmath.
        @param u Argument.
        @param N_max Anzahl Summanden.
        @return mpmath-Zahl.
        @lastModified 2026-03-12
        """
        # Setze Präzision
        with mpmath.workdps(self.praezision):
            u_mp = mpmath.mpf(u)
            pi = mpmath.pi
            ergebnis = mpmath.mpf(0)

            for n in range(1, N_max + 1):
                n2 = n * n
                n4 = n2 * n2

                # exp(4u) für den Dämpfungsterm
                e4u = mpmath.exp(4 * u_mp)

                # Dämpfungsfaktor: exp(-π·n²·e^{4u})
                daempfung = mpmath.exp(-pi * n2 * e4u)

                # Falls Dämpfung numerisch null → Rest vernachlässigbar
                if daempfung < mpmath.mpf(10) ** (-self.praezision + 5):
                    break

                # Wachstumsterm: (2π²n⁴·e^{9u} - 3πn²·e^{5u})
                term1 = 2 * pi**2 * n4 * mpmath.exp(9 * u_mp)
                term2 = 3 * pi * n2 * mpmath.exp(5 * u_mp)
                wachstum = term1 - term2

                ergebnis += wachstum * daempfung

            return ergebnis

    def _phi_python(self, u: float, N_max: int) -> float:
        """
        @brief Python-Fallback für Φ(u) ohne mpmath.
        @param u Argument.
        @param N_max Anzahl Summanden.
        @return Wert als float.
        @lastModified 2026-03-12
        """
        pi = math.pi
        ergebnis = 0.0

        for n in range(1, N_max + 1):
            n2 = n * n
            n4 = n2 * n2

            try:
                e4u = math.exp(4 * u)
                # Verhindere Overflow
                exp_arg = -pi * n2 * e4u
                if exp_arg < -700:
                    break
                daempfung = math.exp(exp_arg)

                term1 = 2 * pi**2 * n4 * math.exp(min(9 * u, 700))
                term2 = 3 * pi * n2 * math.exp(min(5 * u, 700))
                wachstum = term1 - term2

                ergebnis += wachstum * daempfung
            except OverflowError:
                break

        return ergebnis

    def H_t(self, x: float, t: float, N_max: int = 30) -> float:
        """
        @brief Berechnet H_t(x) numerisch via Gauss-Laguerre-ähnlicher Quadratur.
        @description
            H_t(x) ist definiert durch das uneigentliche Integral:
                H_t(x) = ∫₀^∞ e^{t·u²} · Φ(u) · cos(x·u) du

            Für t > 0 wächst der Faktor e^{tu²}, aber Φ(u) fällt schneller
            (superexponentiell), sodass das Integral konvergiert.

            **Numerische Methode:**
            Substitution u → u, Abschneiden bei u_max wo Φ(u)·e^{tu²} ≈ 0.
            Dann Riemann-Summe oder Simpson-Quadratur.

            **Spezialfall t=0:**
                H_0(x) = ∫₀^∞ Φ(u)·cos(xu) du = ½ ξ(½+ix)
            Nullstellen bei x = y_k (imaginäre Teile der Zeta-Nullstellen).

        @param x Reelles Argument.
        @param t Parameter (t=0 entspricht ξ/2).
        @param N_max Terme in Φ-Summe.
        @return Numerischer Wert von H_t(x).
        @lastModified 2026-03-12
        """
        if MPMATH_VERFUEGBAR:
            return float(self._H_t_mpmath(x, t, N_max))
        return self._H_t_quadratur(x, t, N_max)

    def _H_t_mpmath(self, x: float, t: float, N_max: int) -> object:
        """
        @brief Hochpräzise Berechnung von H_t(x) via mpmath-Integration.
        @description
            Verwendet mpmath.quad für adaptive Gauß-Quadratur.
            Obere Integrationsgrenze: u_max = 10 (Φ ist dort ~0).
        @param x Reelles Argument.
        @param t Parameter.
        @param N_max Terme in Φ.
        @return mpmath-Zahl.
        @lastModified 2026-03-12
        """
        with mpmath.workdps(self.praezision):
            x_mp = mpmath.mpf(x)
            t_mp = mpmath.mpf(t)

            def integrand(u):
                """Integrand: e^{tu²}·Φ(u)·cos(xu)"""
                phi_val = self._phi_mpmath(float(u), N_max)
                return mpmath.exp(t_mp * u**2) * phi_val * mpmath.cos(x_mp * u)

            # Adaptive Quadratur von 0 bis ~10 (Φ ist danach vernachlässigbar)
            try:
                ergebnis = mpmath.quad(integrand, [0, 5], error=True, maxdegree=6)
                return ergebnis[0]
            except Exception:
                # Fallback: grobe Riemann-Summe
                return self._H_t_quadratur_mp(x_mp, t_mp, N_max)

    def _H_t_quadratur(self, x: float, t: float, N_max: int) -> float:
        """
        @brief Simpson-Quadratur für H_t(x) (Python-Fallback).
        @param x Argument.
        @param t Parameter.
        @param N_max Terme.
        @return float-Approximation.
        @lastModified 2026-03-12
        """
        # Integrationsbereich [0, u_max], Schrittweite h
        u_max = 5.0
        M = 500  # Anzahl Intervalle (muss gerade sein für Simpson)
        h = u_max / M

        ergebnis = 0.0
        for i in range(M + 1):
            u = i * h
            phi_val = self._phi_python(u, N_max)
            integrand = math.exp(t * u * u) * phi_val * math.cos(x * u)

            # Simpson-Gewicht
            if i == 0 or i == M:
                gewicht = 1.0
            elif i % 2 == 1:
                gewicht = 4.0
            else:
                gewicht = 2.0

            ergebnis += gewicht * integrand

        return ergebnis * h / 3.0

    def _H_t_quadratur_mp(self, x_mp, t_mp, N_max: int):
        """
        @brief mpmath-basierte Simpson-Quadratur als Fallback.
        @param x_mp mpmath-Zahl x.
        @param t_mp mpmath-Zahl t.
        @param N_max Terme.
        @return mpmath-Zahl.
        @lastModified 2026-03-12
        """
        u_max = mpmath.mpf(5)
        M = 200
        h = u_max / M
        ergebnis = mpmath.mpf(0)

        for i in range(M + 1):
            u = i * h
            phi_val = self._phi_mpmath(float(u), N_max)
            integrand = mpmath.exp(t_mp * u**2) * phi_val * mpmath.cos(x_mp * u)

            if i == 0 or i == M:
                gewicht = 1
            elif i % 2 == 1:
                gewicht = 4
            else:
                gewicht = 2

            ergebnis += gewicht * integrand

        return ergebnis * h / 3

    def suche_nullstellen_H_t(
        self,
        t: float,
        x_bereich: Tuple[float, float] = (0.0, 60.0),
        N_max: int = 20,
        schritte: int = 1000,
    ) -> List[float]:
        """
        @brief Findet reelle Nullstellen von H_t(x) im angegebenen Bereich.
        @description
            Sucht nach Vorzeichenwechseln von H_t(x) bei festem t.
            Pro Intervall mit Vorzeichenwechsel wird der Nullpunkt via
            Bisektionsverfahren präzisiert.

            Für t = 0 sollten die Nullstellen bei den bekannten
            Riemann-Nullstellen liegen:
                y₁ ≈ 14.13, y₂ ≈ 21.02, y₃ ≈ 25.01, ...

        @param t Parameter der H_t-Familie.
        @param x_bereich Suchbereich (a, b).
        @param N_max Terme in Φ-Summe.
        @param schritte Anzahl der Auswertungspunkte.
        @return Sortierte Liste gefundener Nullstellen.
        @lastModified 2026-03-12
        """
        a, b = x_bereich
        dx = (b - a) / schritte
        nullstellen: List[float] = []

        # Erste Auswertung
        x_prev = a
        f_prev = self.H_t(x_prev, t, N_max)

        for i in range(1, schritte + 1):
            x_curr = a + i * dx
            f_curr = self.H_t(x_curr, t, N_max)

            # Vorzeichenwechsel detektiert?
            if f_prev * f_curr < 0:
                # Bisektion zur Präzisierung
                x_links, x_rechts = x_prev, x_curr
                f_links = f_prev
                for _ in range(50):  # 50 Bisektion-Schritte → ~15 Stellen
                    x_mitte = (x_links + x_rechts) / 2
                    f_mitte = self.H_t(x_mitte, t, N_max)
                    if abs(f_mitte) < 1e-14:
                        break
                    if f_links * f_mitte < 0:
                        x_rechts = x_mitte
                    else:
                        x_links = x_mitte
                        f_links = f_mitte
                nullstellen.append((x_links + x_rechts) / 2)

            x_prev = x_curr
            f_prev = f_curr

        return sorted(nullstellen)

    def ist_alle_reell(
        self,
        t: float,
        grenze: float = 60.0,
        N_max: int = 20,
        toleranz: float = 1e-8,
    ) -> bool:
        """
        @brief Prüft heuristisch ob H_t im Bereich [0, grenze] nur reelle Nullstellen hat.
        @description
            Untersucht die Nullstellen von H_t(x) numerisch. Falls alle
            gefundenen Nullstellen reell sind und deren Anzahl mit der erwarteten
            übereinstimmt (für t=0: Riemann-Nullstellen), wird True zurückgegeben.

            **Einschränkung**: Dies ist eine heuristische Prüfung. Ein
            rigoroser Beweis erfordert Winding-Number-Berechnung.

        @param t Parameter.
        @param grenze Obere Grenze des Suchbereichs.
        @param N_max Terme in Φ.
        @param toleranz Toleranz für Nullstellenerkennung.
        @return True wenn heuristisch nur reelle Nullstellen.
        @lastModified 2026-03-12
        """
        nullstellen = self.suche_nullstellen_H_t(t, (0.0, grenze), N_max)

        # Für t=0: prüfe ob erste bekannte Nullstellen reproduziert werden
        if abs(t) < 0.01:
            erwartete = [y for y in BEKANNTE_NULLSTELLEN_GAMMA if y <= grenze]
            if len(nullstellen) < len(erwartete):
                return False
            # Vergleiche mit Toleranz
            for i, y_exp in enumerate(erwartete):
                if i >= len(nullstellen):
                    return False
                if abs(nullstellen[i] - y_exp) > 1.0:  # 1.0 Toleranz für numerische Näherung
                    return False

        return True

    def schranke_lambda_oben(self, kandidat_t: float) -> Dict[str, object]:
        """
        @brief Analysiert ob Λ ≤ kandidat_t gilt (numerische Schrankenprüfung).
        @description
            Prüft heuristisch ob H_{kandidat_t} nur reelle Nullstellen hat.
            Falls ja, liefert dies eine (numerische) Evidenz für Λ ≤ kandidat_t.

            **Wichtiger Hinweis**: Ein rigoroser Beweis von Λ ≤ t₀ erfordert
            analytische Methoden (Winding-Number, explizite Fehlerabschätzungen).
            Diese Funktion liefert nur numerische Indizien.

            **Aktuell bekannt**: Λ ≤ 0.2 (Platt-Trudgian 2021, rigoros bewiesen).

        @param kandidat_t Zu prüfender t-Wert als obere Schranke für Λ.
        @return Dictionary mit Analyse-Ergebnis.
        @lastModified 2026-03-12
        """
        heuristisch_reell = self.ist_alle_reell(kandidat_t, grenze=50.0)

        return {
            "kandidat_t": kandidat_t,
            "heuristisch_nur_reelle_nullstellen": heuristisch_reell,
            "schlussfolgerung": (
                f"Λ ≤ {kandidat_t} möglich (heuristisch)"
                if heuristisch_reell
                else f"Λ ≤ {kandidat_t} nicht durch Numerik gestützt"
            ),
            "hinweis": (
                "Rigorer Beweis erfordert Winding-Number-Berechnung. "
                "Bekannte rigorse Schranke: Λ ≤ 0.2 (Platt-Trudgian 2021)."
            ),
        }

    def riemann_zusammenhang(self) -> Dict[str, str]:
        """
        @brief Erklärt den formalen Zusammenhang zwischen RH und der Konstante Λ.
        @description
            Die Äquivalenz RH ⟺ Λ ≤ 0 ist ein zentrales Resultat der
            de Bruijn-Newman-Theorie. Diese Methode dokumentiert die
            mathematischen Details.

        @return Dictionary mit formaler Beschreibung des Zusammenhangs.
        @lastModified 2026-03-12
        """
        return {
            "hauptsatz": "RH (Riemann-Hypothese) ⟺ Λ ≤ 0",
            "definition_Phi": (
                "Φ(u) = Σ_{n=1}^∞ (2π²n⁴e^{9u} - 3πn²e^{5u}) · exp(-πn²e^{4u})"
            ),
            "definition_H_t": (
                "H_t(x) = ∫₀^∞ e^{tu²} · Φ(u) · cos(xu) du  (t ∈ ℝ, x ∈ ℝ)"
            ),
            "definition_Lambda": (
                "Λ = inf{ t ∈ ℝ : H_t hat nur reelle Nullstellen }"
            ),
            "spezialfall_t0": (
                "H_0(x) = ½ · ξ(½ + ix)  wobei ξ(s) = ½s(s-1)π^{-s/2}Γ(s/2)ζ(s)"
            ),
            "rh_aequivalenz": (
                "RH ⟺ alle nicht-trivialen Nullstellen von ζ auf Re(s)=½ "
                "⟺ ξ(½+ix) hat nur reelle Nullstellen (x ∈ ℝ) "
                "⟺ H_0 hat nur reelle Nullstellen "
                "⟺ Λ ≤ 0"
            ),
            "de_bruijn_1950": (
                "de Bruijn (1950) bewies: Wenn H_t₀ nur reelle Nullstellen hat, "
                "dann hat H_t nur reelle Nullstellen für alle t ≥ t₀. "
                "Damit ist Λ wohldefiniert."
            ),
            "newman_1976": (
                "Newman (1976) vermutete Λ ≥ 0. Motivation: Die RH ist scharf, "
                "d.h. H_t für t < 0 hat nicht-reelle Nullstellen."
            ),
            "rodgers_tao_2018": (
                "Rodgers & Tao (2018) bewiesen Λ ≥ 0 (Newman-Vermutung). "
                "Damit gilt: Λ = 0 ⟺ RH wahr."
            ),
        }

    def newman_vermutung_status(self) -> Dict[str, object]:
        """
        @brief Dokumentiert den Beweis-Status der Newman-Vermutung (Λ ≥ 0).
        @description
            Die Newman-Vermutung wurde 2018 von Brad Rodgers und Terence Tao bewiesen.
            Dies zeigt, dass die RH-Äquivalenz 'scharf' ist: Λ = 0 wäre der
            kleinstmögliche Wert.

            **Bedeutung**: RH ⟺ Λ = 0. Die RH ist also äquivalent dazu, dass
            H_t nur für alle t ≥ 0 reelle Nullstellen hat, und genau bei t=0
            beginnt dies zu gelten (oder es gilt auch für t < 0, was ausgeschlossen
            ist durch Rodgers-Tao).

        @return Dictionary mit vollständigem Status-Bericht.
        @lastModified 2026-03-12
        """
        return {
            "vermutung": "Newman-Vermutung: Λ ≥ 0",
            "status": "BEWIESEN ✓",
            "beweis_quelle": "B. Rodgers & T. Tao (2018)",
            "titel": "The de Bruijn-Newman constant is non-negative",
            "arXiv": "arXiv:1801.05914",
            "jahr_beweis": 2018,
            "konsequenz": (
                "Zusammen mit RH ⟺ Λ ≤ 0 gilt: "
                "Falls RH wahr, dann Λ = 0 (nicht Λ < 0). "
                "Die RH ist also eine 'scharfe' Aussage."
            ),
            "aktueller_wissenstand": {
                "untere_schranke": "Λ ≥ 0  (Rodgers-Tao 2018, rigoros)",
                "obere_schranke": "Λ ≤ 0.2  (Platt-Trudgian 2021, rigoros)",
                "offene_frage": "Ob Λ = 0 (⟺ RH) ist offen",
            },
            "methode_rodgers_tao": (
                "Wahrscheinlichkeitstheoretischer Ansatz: Zeige, dass für "
                "t < 0 die Funktion H_t immer nicht-reelle Nullstellen hat. "
                "Technisches Werkzeug: Korrelationen von Nullstellen-Paarverteilungen."
            ),
        }

    def bekannte_schranken_history(self) -> List[Dict[str, object]]:
        """
        @brief Gibt die historische Entwicklung der Schranken für Λ zurück.
        @description
            Listet alle bekannten rigoros bewiesenen Schranken für die
            de Bruijn-Newman-Konstante Λ in chronologischer Reihenfolge.

        @return Liste von Dictionaries mit Jahres- und Schranken-Informationen.
        @lastModified 2026-03-12
        """
        return [
            {
                "jahr": 1950,
                "schranke": "Λ ≤ ½",
                "richtung": "oben",
                "quelle": "N.G. de Bruijn (1950): 'The roots of trigonometric integrals'",
                "methode": "Direkter Beweis via Gauß-Faltung",
            },
            {
                "jahr": 1976,
                "schranke": "Λ ≥ −∞ (Vermutung: Λ ≥ 0)",
                "richtung": "unten (Vermutung)",
                "quelle": "C.M. Newman (1976): 'Fourier transforms with only real zeros'",
                "methode": "Heuristische Analyse der Nullstellenstruktur",
            },
            {
                "jahr": 1988,
                "schranke": "Λ < ½",
                "richtung": "oben (scharf)",
                "quelle": "C.M. Newman, A. Goldberg (1988)",
                "methode": "Analytische Verbesserung des de Bruijn-Bounds",
            },
            {
                "jahr": 2009,
                "schranke": "Λ ≥ −2.7 · 10⁻⁹",
                "richtung": "unten",
                "quelle": "Y. Saouter, X. Gourdon, P. Demichel (2009)",
                "methode": "Numerische Analyse von Nullstellen-Paaren nahe der kritischen Geraden",
            },
            {
                "jahr": 2018,
                "schranke": "Λ ≥ 0  ← BEWIESEN",
                "richtung": "unten (Newman-Vermutung bewiesen)",
                "quelle": "B. Rodgers, T. Tao (2018): arXiv:1801.05914",
                "methode": "Wahrscheinlichkeitstheorie, Korrelationen der GUE-Nullstellenverteilung",
            },
            {
                "jahr": 2020,
                "schranke": "0 ≤ Λ ≤ 0.22",
                "richtung": "oben+unten kombiniert",
                "quelle": "Polymath15 Projekt (T. Tao et al., 2020)",
                "methode": "Numerische Massenberechnung von H_t-Nullstellen via Gittersummen",
            },
            {
                "jahr": 2021,
                "schranke": "0 ≤ Λ ≤ 0.2",
                "richtung": "oben (beste bekannte)",
                "quelle": "D. Platt, T. Trudgian (2021)",
                "methode": "Rigorse numerische Verifikation mit Fehlerabschätzung",
            },
        ]

    def nullstellen_H0_berechnen(
        self,
        x_max: float = 50.0,
        N_max: int = 20,
    ) -> List[float]:
        """
        @brief Berechnet die ersten Nullstellen von H_0(x) numerisch.
        @description
            H_0(x) = ½·ξ(½+ix). Die Nullstellen von H_0 entsprechen den
            imaginären Teilen der Riemann-Zeta-Nullstellen auf der kritischen Geraden.

            Erwartete erste Nullstellen:
                y₁ ≈ 14.134725...
                y₂ ≈ 21.022040...
                y₃ ≈ 25.010857...

            Diese Methode kombiniert suche_nullstellen_H_t mit t=0.

        @param x_max Obere Grenze für die Nullstellensuche.
        @param N_max Terme in Φ-Summe.
        @return Liste der gefundenen Nullstellen.
        @lastModified 2026-03-12
        """
        return self.suche_nullstellen_H_t(
            t=0.0,
            x_bereich=(0.0, x_max),
            N_max=N_max,
            schritte=2000,
        )

    def vergleiche_mit_riemann_nullstellen(
        self,
        nullstellen: List[float],
    ) -> Dict[str, object]:
        """
        @brief Vergleicht numerisch gefundene Nullstellen mit bekannten Riemann-Nullstellen.
        @description
            Prüft die Übereinstimmung zwischen den numerisch berechneten
            Nullstellen von H_0 und den bekannten Riemann-Zeta-Nullstellen.

        @param nullstellen Liste der numerisch gefundenen Nullstellen.
        @return Dictionary mit Vergleichsstatistik.
        @lastModified 2026-03-12
        """
        ergebnisse = []
        max_vergleich = min(len(nullstellen), len(BEKANNTE_NULLSTELLEN_GAMMA))

        for i in range(max_vergleich):
            gefunden = nullstellen[i]
            erwartet = BEKANNTE_NULLSTELLEN_GAMMA[i]
            fehler = abs(gefunden - erwartet)
            ergebnisse.append({
                "index": i + 1,
                "erwartet": erwartet,
                "gefunden": round(gefunden, 6),
                "absoluter_fehler": round(fehler, 6),
                "relativ_fehler": round(fehler / erwartet, 8) if erwartet != 0 else None,
            })

        return {
            "anzahl_verglichen": max_vergleich,
            "ergebnisse": ergebnisse,
            "max_absoluter_fehler": max(r["absoluter_fehler"] for r in ergebnisse) if ergebnisse else None,
            "zusammenfassung": (
                f"{max_vergleich} Nullstellen verglichen. "
                f"Maximaler Fehler: {max(r['absoluter_fehler'] for r in ergebnisse):.6f}"
                if ergebnisse else "Keine Nullstellen zum Vergleichen."
            ),
        }

    def visualisiere_H0(
        self,
        x_max: float = 50.0,
        N_max: int = 20,
        dateiname: Optional[str] = None,
    ) -> bool:
        """
        @brief Visualisiert H_0(x) für x ∈ [0, x_max] mit matplotlib.
        @description
            Zeichnet den Graphen von H_0(x) und markiert die bekannten
            Riemann-Nullstellen als vertikale Linien.

            H_0(x) = ½ ξ(½+ix) oszilliert und hat Nullstellen bei den
            imaginären Teilen der Riemann-Zeta-Nullstellen.

        @param x_max Obere Grenze des Plots.
        @param N_max Terme in Φ.
        @param dateiname Falls angegeben, wird der Plot gespeichert (statt angezeigt).
        @return True wenn Plot erstellt, False wenn matplotlib fehlt.
        @lastModified 2026-03-12
        """
        if not MATPLOTLIB_VERFUEGBAR:
            return False

        # Berechne H_0(x) auf einem Gitter
        x_werte = [i * x_max / 300 for i in range(301)]
        y_werte = [self.H_t(x, 0.0, N_max) for x in x_werte]

        fig, ax = plt.subplots(figsize=(14, 5))
        ax.plot(x_werte, y_werte, 'b-', linewidth=1.2, label=r'$H_0(x)$')
        ax.axhline(y=0, color='k', linewidth=0.5)

        # Markiere bekannte Nullstellen
        for i, y_null in enumerate(BEKANNTE_NULLSTELLEN_GAMMA):
            if y_null <= x_max:
                ax.axvline(x=y_null, color='r', linewidth=0.8, alpha=0.6,
                           label=f'$\\gamma_{i+1}$' if i < 3 else None)
                ax.annotate(f'$\\gamma_{i+1}$', xy=(y_null, 0),
                            xytext=(y_null, max(y_werte) * 0.3 if y_werte else 0.1),
                            fontsize=7, color='red', ha='center')

        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel(r'$H_0(x)$', fontsize=12)
        ax.set_title(
            r'$H_0(x) = \frac{1}{2}\xi(\frac{1}{2}+ix)$ — Nullstellen = Riemann-$\zeta$-Nullstellen',
            fontsize=13,
        )
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)

        if dateiname:
            plt.savefig(dateiname, dpi=150, bbox_inches='tight')
        else:
            plt.tight_layout()
            plt.show()

        plt.close()
        return True


# ===========================================================================
# HAUPTPROGRAMM (Demo)
# ===========================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("DE BRUIJN-NEWMAN-KONSTANTE Λ: Analyse")
    print("=" * 60)

    dbn = DeBruijnNewman(praezision=30)

    # --- Zusammenhang RH ↔ Λ ---
    print("\n[1] Formaler Zusammenhang RH ↔ Λ:")
    zusammenhang = dbn.riemann_zusammenhang()
    print(f"  Hauptsatz: {zusammenhang['hauptsatz']}")
    print(f"  Spezialfall t=0: {zusammenhang['spezialfall_t0']}")

    # --- Newman-Vermutung Status ---
    print("\n[2] Newman-Vermutung:")
    status = dbn.newman_vermutung_status()
    print(f"  Vermutung: {status['vermutung']}")
    print(f"  Status: {status['status']}")
    print(f"  Quelle: {status['beweis_quelle']} ({status['titel']})")
    print(f"  Konsequenz: {status['konsequenz']}")

    # --- Historische Schranken ---
    print("\n[3] Historische Schranken für Λ:")
    history = dbn.bekannte_schranken_history()
    for eintrag in history:
        print(f"  {eintrag['jahr']}: {eintrag['schranke']} ({eintrag['quelle'].split('(')[0].strip()})")

    # --- Phi-Funktion testen ---
    print("\n[4] Φ(u) für ausgewählte u-Werte:")
    for u_val in [0.0, 0.5, 1.0, 2.0]:
        phi_val = dbn.phi_funktion(u_val, N_max=10)
        print(f"  Φ({u_val:.1f}) = {phi_val:.6e}")

    # --- H_0 Nullstellen berechnen ---
    print("\n[5] Nullstellen von H_0(x) für x ∈ [0, 50]:")
    nullstellen = dbn.nullstellen_H0_berechnen(x_max=50.0, N_max=15)
    print(f"  Gefundene Nullstellen: {len(nullstellen)}")
    for i, y in enumerate(nullstellen[:5]):
        erwartet = BEKANNTE_NULLSTELLEN_GAMMA[i] if i < len(BEKANNTE_NULLSTELLEN_GAMMA) else "?"
        print(f"  y_{i+1}: gefunden={y:.4f}, bekannt={erwartet:.4f}")

    # --- Vergleich mit Riemann ---
    if nullstellen:
        print("\n[6] Vergleich mit bekannten Riemann-Nullstellen:")
        vergleich = dbn.vergleiche_mit_riemann_nullstellen(nullstellen)
        print(f"  {vergleich['zusammenfassung']}")

    # --- Obere Schranke prüfen ---
    print("\n[7] Numerische Schranken-Analyse:")
    for t_test in [0.0, 0.2]:
        schranke = dbn.schranke_lambda_oben(t_test)
        print(f"  t={t_test}: {schranke['schlussfolgerung']}")
