"""
@file riemann_siegel_ext.py
@brief Erweiterte Riemann-Siegel-Formel: Z(t)-Funktion, Gram-Punkte und Fehlerterme.
@description
    Dieses Modul implementiert die erweiterte Riemann-Siegel-Formel zur
    effizienten Berechnung der Riemann-Zeta-Funktion auf der kritischen Geraden.

    **Riemann-Siegel-Formel** (Riemann ~1852, Siegel 1932 veröffentlicht):
        Z(t) = e^{iθ(t)} · ζ(1/2 + it)

    Dabei ist Z(t) für reelle t stets reell (Hardy-Z-Funktion).

    **Siegelsche θ-Funktion**:
        θ(t) = Im(ln Γ(1/4 + it/2)) − t/2 · ln π

    **Hauptformel** (N = ⌊√(t/2π)⌋):
        Z(t) = 2 · Σ_{n=1}^{N} n^{−1/2} · cos(θ(t) − t·ln n) + R(t)

    **Fehlerterm**:
        R(t) = O(t^{−1/4})  (explizite Schranken nach Gabcke 1979)

    **Gram-Punkte** gₙ (nach J.-P. Gram, 1903):
        θ(gₙ) = n·π  ⟺  Z(gₙ) alterniert im Vorzeichen mit (−1)ⁿ
        (Gram's Law: Z(gₙ) > 0 für gerade n, Z(gₙ) < 0 für ungerade n)
        → Gilt NICHT immer (Gram-Ausnahmen / "Gram failures")

    **Gram-Blöcke**:
        Intervall [gₙ, gₙ₊₁] wird "Gram-Block der Breite k" wenn k aufeinanderfolgende
        Gram-Intervalle zu einem Block zusammengefasst werden müssen (Rosser-Regel).

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import cmath
from typing import Dict, List, Optional, Tuple

import mpmath
import numpy as np

# mpmath-Präzision für hochgenaue Berechnungen
mpmath.mp.dps = 50


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _theta_mpmath(t: float) -> float:
    """
    Berechnet die Siegelsche θ-Funktion mit mpmath für hohe Genauigkeit.

    θ(t) = Im(ln Γ(1/4 + it/2)) − t/2 · ln π

    Die θ-Funktion wird benötigt, um Z(t) als reelle Funktion darzustellen.

    @param t: Reeller Parameter (t > 0 für sinnvolle Ergebnisse)
    @return: θ(t) als float
    @lastModified: 2026-03-12
    """
    # Im(log Γ(1/4 + it/2)) − t/2 · ln π
    z = mpmath.mpc("0.25", 0) + mpmath.mpc(0, t / 2.0)
    lg = mpmath.loggamma(z)
    return float(lg.imag) - (t / 2.0) * math.log(math.pi)


def _theta_approx(t: float) -> float:
    """
    Schnelle Näherungsformel für θ(t) via Stirling-Entwicklung.

    Stirling-Asymptotik für große t:
        θ(t) ≈ t/2·ln(t/2π) − t/2 − π/8 + 1/(48t) + 7/(5760t³) + ...

    Fehler: O(t^{−5}) — ausreichend genau für t > 10.

    @param t: Reeller Parameter t > 0
    @return: Näherung von θ(t)
    @lastModified: 2026-03-12
    """
    t2pi = t / (2.0 * math.pi)
    result = (t / 2.0) * math.log(t2pi) - t / 2.0 - math.pi / 8.0
    # Asymptotische Korrekturen
    if t > 1.0:
        result += 1.0 / (48.0 * t)
        if t > 10.0:
            result += 7.0 / (5760.0 * t ** 3)
    return result


# ===========================================================================
# HAUPTKLASSE: RiemannSiegelExtended
# ===========================================================================

class RiemannSiegelExtended:
    """
    Erweiterte Implementierung der Riemann-Siegel-Formel.

    Berechnet die Hardy-Z-Funktion Z(t) = e^{iθ(t)} · ζ(1/2 + it),
    die für reelles t stets reell ist (|Z(t)| = |ζ(1/2+it)|).

    Nullstellen von Z(t) entsprechen Nullstellen von ζ(s) auf der
    kritischen Geraden Re(s) = 1/2.

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self, use_high_precision: bool = False):
        """
        Initialisiert das RiemannSiegelExtended-Objekt.

        @param use_high_precision: Wenn True, wird mpmath für alle Berechnungen
                                    verwendet (langsamer, aber genauer).
        @lastModified: 2026-03-12
        """
        self.use_high_precision = use_high_precision
        # Cache für berechnete Gram-Punkte
        self._gram_cache: Dict[int, float] = {}

    def theta(self, t: float) -> float:
        """
        Berechnet θ(t) = Im(ln Γ(1/4 + it/2)) − t/2 · ln π.

        Die θ-Funktion gibt den Phasenwinkel von ζ(1/2+it) an.
        Sie ist monoton wachsend für t > 6.29... (Minimum).

        @param t: Reeller Parameter
        @return: θ(t)
        @lastModified: 2026-03-12
        """
        if self.use_high_precision or abs(t) < 20.0:
            return _theta_mpmath(t)
        return _theta_approx(t)

    def Z(self, t: float, terms: Optional[int] = None) -> float:
        """
        Berechnet die Hardy-Z-Funktion Z(t) via Riemann-Siegel-Formel.

        Hauptformel (N = ⌊√(t/2π)⌋):
            Z(t) = 2 · Σ_{n=1}^{N} n^{−1/2} · cos(θ(t) − t·ln n) + R(t)

        Die Hauptsumme konvergiert für alle t > 0.
        Der Fehlerterm R(t) = O(t^{−1/4}) wird durch den C₀-Term approximiert.

        **Hinweis**: Die Riemann-Hypothese besagt, dass alle nicht-trivialen
        Nullstellen von ζ(s) auf Re(s) = 1/2 liegen, d.h. alle reellen
        Nullstellen von Z(t) (t > 0).

        @param t: Reeller Parameter t > 0
        @param terms: Anzahl der Summanden (Standard: ⌊√(t/2π)⌋)
        @return: Z(t) als reelle Zahl
        @lastModified: 2026-03-12
        """
        if t <= 0:
            raise ValueError(f"Z(t) ist nur für t > 0 definiert, erhalten: t={t}")

        # Anzahl der Hauptsummanden: N = ⌊√(t/2π)⌋
        N = terms if terms is not None else int(math.sqrt(t / (2.0 * math.pi)))
        N = max(N, 1)

        # Siegelsche θ-Funktion
        th = self.theta(t)

        # Hauptsumme der Riemann-Siegel-Formel
        # Σ_{n=1}^{N} n^{−1/2} · cos(θ(t) − t·ln n)
        main_sum = 0.0
        for n in range(1, N + 1):
            main_sum += math.cos(th - t * math.log(n)) / math.sqrt(n)
        main_sum *= 2.0

        # Fehlerterm: C₀-Korrektur nach Riemann-Siegel
        # p = {√(t/2π)} (gebrochener Anteil von √(t/2π))
        frac_sqrt = math.sqrt(t / (2.0 * math.pi)) - N
        # C₀-Term: cos(2π(p² − p − 1/16)) / cos(2πp) · t^{−1/4}
        # (Gabcke 1979, vereinfachte Form)
        remainder = self._remainder_term(t, N, frac_sqrt)

        return main_sum + remainder

    def _remainder_term(self, t: float, N: int, p: float) -> float:
        """
        Berechnet den Fehlerterm R(t) der Riemann-Siegel-Formel.

        Näherungsformel (Hauptterm des Fehlerterms, Gabcke 1979):
            R(t) ≈ (−1)^{N−1} · t^{−1/4} · C₀(p)

        mit C₀(p) = cos(2π(p² − p − 1/16)) / cos(2πp)

        Absolute Schranke: |R(t)| ≤ 0.053 · t^{−3/4} + 0.14 · t^{−5/4}

        @param t: Zeitparameter
        @param N: Anzahl der Hauptsummanden
        @param p: Gebrochener Anteil von √(t/2π)
        @return: Näherung des Fehlerterms
        @lastModified: 2026-03-12
        """
        # Vorzeichen alterniert mit N
        sign = (-1) ** (N - 1)
        # t^{-1/4} Faktor
        t_factor = t ** (-0.25)
        # C₀(p)-Term
        denom = math.cos(2.0 * math.pi * p)
        if abs(denom) < 1e-10:
            return 0.0  # Singularität vermeiden
        c0 = math.cos(2.0 * math.pi * (p * p - p - 1.0 / 16.0)) / denom
        return sign * t_factor * c0

    def error_bound(self, t: float) -> float:
        """
        Obere Schranke für den absoluten Fehler des Hauptterms.

        Nach Gabcke (1979):
            |R(t)| ≤ 0.053 · t^{−3/4} + 0.14 · t^{−5/4}

        Für t > 100 ist der Fehler kleiner als 0.01.

        @param t: Zeitparameter t > 0
        @return: Obere Fehlerschranke
        @lastModified: 2026-03-12
        """
        return 0.053 * t ** (-0.75) + 0.14 * t ** (-1.25)

    def find_gram_point(self, n: int, tol: float = 1e-10) -> float:
        """
        Berechnet den n-ten Gram-Punkt gₙ: θ(gₙ) = n·π.

        Gram-Punkte sind definiert als die Lösungen von θ(t) = nπ.
        Für n = 0: g₀ ≈ 17.845...
        Für n = -1: g₋₁ ≈ 9.666... (oft als erster Gram-Punkt gelistet)

        Methode: Newton-Verfahren mit Startpunkt aus Invertierung der
        Stirling-Näherung θ(t) ≈ t/2·ln(t/2π) − t/2 − π/8.

        @param n: Gram-Index (ganze Zahl, auch negativ)
        @param tol: Toleranz für Newton-Verfahren
        @return: Gram-Punkt gₙ
        @lastModified: 2026-03-12
        """
        # Cache-Prüfung
        if n in self._gram_cache:
            return self._gram_cache[n]

        # Startwert: Umkehrung der Stirling-Asymptotik
        # θ(t) ≈ n·π → t ≈ 2π·exp(W(n·π + π/8 + 1) + 1)
        # Einfache Näherung für n ≥ 0:
        target = n * math.pi
        # Schätzung: t ≈ 2πe · exp(W(...)) via Iteration
        # Startpunkt: t₀ so dass θ(t₀) ≈ target
        if n >= 0:
            # Grobe Näherung aus inverser Stirling
            t = max(2.0 * math.pi * math.e * (abs(n) + 1) ** 0.5, 10.0)
            # Bessere Schätzung
            for _ in range(5):
                t = 2.0 * math.pi * math.exp(
                    (target + math.pi / 8.0 + t / 2.0) / (t / 2.0)
                ) if t > 0.1 else 10.0
        else:
            t = max(2.0 * math.pi * math.e, 6.3)

        # Sicherstellen dass t > 0
        t = max(t, 6.0)

        # Newton-Verfahren: f(t) = θ(t) − n·π = 0
        # f'(t) ≈ d/dt θ(t) ≈ (1/2)·ln(t/2π) + kleine Korrekturen
        for _ in range(50):
            th = _theta_mpmath(t)
            f = th - target
            # Ableitung θ'(t) ≈ ln(t/2π)/2 für große t
            dt = max(t / (4.0 * math.pi), 0.1)
            # Numerische Ableitung
            dth = (_theta_mpmath(t + 1e-5) - _theta_mpmath(t - 1e-5)) / 2e-5
            if abs(dth) < 1e-15:
                break
            t_new = t - f / dth
            if t_new <= 0:
                t_new = t / 2.0
            if abs(t_new - t) < tol:
                t = t_new
                break
            t = t_new

        self._gram_cache[n] = t
        return t

    def gram_sign(self, n: int) -> int:
        """
        Gibt das erwartete Vorzeichen von Z(gₙ) gemäß Gram's Law zurück.

        Gram's Law: Z(gₙ) sollte das Vorzeichen (−1)ⁿ haben.
        - n gerade: Z(gₙ) > 0 erwartet → +1
        - n ungerade: Z(gₙ) < 0 erwartet → -1

        Dies gilt statistisch in ~73% aller Fälle (Titchmarsh).

        @param n: Gram-Index
        @return: Erwartetes Vorzeichen: +1 oder -1
        @lastModified: 2026-03-12
        """
        return 1 if n % 2 == 0 else -1

    def check_gram_law(self, n: int) -> Tuple[bool, float, float]:
        """
        Prüft ob das Gram-Gesetz für den n-ten Gram-Punkt gilt.

        Berechnet Z(gₙ) und vergleicht das Vorzeichen mit (−1)ⁿ.

        @param n: Gram-Index
        @return: Tupel (gram_law_holds, gram_point, Z_value)
        @lastModified: 2026-03-12
        """
        gn = self.find_gram_point(n)
        z_val = self.Z(gn)
        expected_sign = self.gram_sign(n)
        holds = (z_val > 0 and expected_sign > 0) or (z_val < 0 and expected_sign < 0)
        return holds, gn, z_val

    def gram_statistics(self, n_max: int = 50) -> Dict:
        """
        Berechnet Gram-Statistiken für Gram-Punkte g₀ bis g_{n_max}.

        Untersucht: Wie oft gilt Gram's Law? Wie viele Ausnahmen gibt es?
        Historisch gilt das Gesetz in ~73% der Fälle (Titchmarsh 1935).

        @param n_max: Maximaler Gram-Index
        @return: Dictionary mit Statistiken
        @lastModified: 2026-03-12
        """
        total = 0
        successes = 0
        failures = []

        for n in range(0, n_max + 1):
            try:
                holds, gn, z_val = self.check_gram_law(n)
                total += 1
                if holds:
                    successes += 1
                else:
                    failures.append((n, gn, z_val))
            except Exception:
                pass

        return {
            "total": total,
            "successes": successes,
            "failures": len(failures),
            "failure_rate": (total - successes) / total if total > 0 else 0.0,
            "failure_list": failures[:10],  # Erste 10 Ausnahmen
        }

    def find_zeros_in_interval(self, t_min: float, t_max: float,
                                num_points: int = 1000) -> List[float]:
        """
        Sucht Nullstellen von Z(t) im Intervall [t_min, t_max].

        Methode: Vorzeichenwechsel-Suche (Bisektionsverfahren nach
        Vorzeichenwechsel-Detektion im Gitter).

        @param t_min: Untere Grenze (t_min > 0)
        @param t_max: Obere Grenze
        @param num_points: Anzahl der Abtastpunkte
        @return: Liste der gefundenen Nullstellen (gerundete Werte)
        @lastModified: 2026-03-12
        """
        ts = np.linspace(t_min, t_max, num_points)
        zs = [self.Z(float(t)) for t in ts]
        zeros = []

        for i in range(len(zs) - 1):
            # Vorzeichenwechsel detektiert
            if zs[i] * zs[i + 1] < 0:
                # Bisektionsverfahren zur Präzisierung
                a, b = float(ts[i]), float(ts[i + 1])
                for _ in range(60):
                    mid = (a + b) / 2.0
                    zm = self.Z(mid)
                    if abs(zm) < 1e-12:
                        a = b = mid
                        break
                    if self.Z(a) * zm < 0:
                        b = mid
                    else:
                        a = mid
                zeros.append((a + b) / 2.0)

        return zeros


# ===========================================================================
# KLASSE: GramBlocks
# ===========================================================================

class GramBlocks:
    """
    Analyse der Gram-Blöcke und Rosser-Regel.

    **Gram-Block** (Bredberg, Odlyzko):
        Ein Gram-Block [gₙ, gₙ₊ₖ] ist ein maximales Intervall von k ≥ 1
        aufeinanderfolgenden Gram-Intervallen [gₙ₊ⱼ, gₙ₊ⱼ₊₁], in dem kein
        Gram-Intervall die Bedingung des Gram-Gesetzes erfüllt — bis auf
        das erste und letzte.

    **Rosser-Regel** (Rosser 1941):
        Wenn Z(gₙ) > 0 und Z(gₙ₊₂) > 0 (beide das richtige Vorzeichen),
        dann liegen GENAU eine Nullstelle in [gₙ, gₙ₊₁] und [gₙ₊₁, gₙ₊₂].
        → Verallgemeinerung sichert die Zählbarkeit aller Nullstellen.

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self, rs_extended: Optional[RiemannSiegelExtended] = None):
        """
        Initialisiert die GramBlocks-Analyse.

        @param rs_extended: RiemannSiegelExtended-Instanz (wird neu erstellt wenn None)
        @lastModified: 2026-03-12
        """
        self.rs = rs_extended or RiemannSiegelExtended()

    def identify_gram_blocks(self, n_start: int = 0, n_end: int = 30) -> List[Dict]:
        """
        Identifiziert Gram-Blöcke im Bereich [n_start, n_end].

        Gram-Block der Breite k: k aufeinanderfolgende Gram-Intervalle,
        die zusammen k Nullstellen enthalten (aber einzeln nicht).

        @param n_start: Start-Gram-Index
        @param n_end: End-Gram-Index
        @return: Liste von Gram-Block-Deskriptoren
        @lastModified: 2026-03-12
        """
        blocks = []
        signs = []

        # Berechne Z(gₙ) für alle Gram-Punkte im Bereich
        for n in range(n_start, n_end + 1):
            try:
                gn = self.rs.find_gram_point(n)
                z_val = self.rs.Z(gn)
                expected = self.rs.gram_sign(n)
                signs.append({
                    "n": n,
                    "gram_point": gn,
                    "Z_value": z_val,
                    "expected_sign": expected,
                    "gram_law_holds": (z_val * expected > 0),
                })
            except Exception:
                continue

        # Identifiziere Gram-Blöcke: zusammenhängende Bereiche mit Gram-Failures
        i = 0
        while i < len(signs):
            if signs[i]["gram_law_holds"]:
                # Einzelnes gutes Gram-Intervall (Block der Breite 1)
                blocks.append({
                    "type": "good",
                    "width": 1,
                    "start_n": signs[i]["n"],
                    "end_n": signs[i]["n"],
                })
                i += 1
            else:
                # Beginn eines Gram-Blocks: suche Ende
                start_i = i
                while i < len(signs) and not signs[i]["gram_law_holds"]:
                    i += 1
                # Block der Breite (i - start_i + 1)
                blocks.append({
                    "type": "gram_block",
                    "width": i - start_i + 1,
                    "start_n": signs[start_i]["n"],
                    "end_n": signs[min(i, len(signs) - 1)]["n"],
                })

        return blocks

    def rosser_rule_holds(self, n: int) -> bool:
        """
        Prüft ob die Rosser-Regel für den n-ten Gram-Punkt gilt.

        Rosser-Regel: Z(gₙ)·Z(gₙ₊₁) hat das erwartete kombinierte Vorzeichen.
        Genauer: Wenn Z(gₙ) > 0 (gerade n), dann hat [gₙ, gₙ₊₁] mindestens
        eine Nullstelle von Z (und damit von ζ auf der kritischen Geraden).

        @param n: Gram-Index
        @return: True wenn Rosser-Regel gilt
        @lastModified: 2026-03-12
        """
        try:
            _, _, z_n = self.rs.check_gram_law(n)
            _, _, z_n1 = self.rs.check_gram_law(n + 1)
            # Rosser-Regel: Z(gₙ)·(−1)ⁿ > 0 AND Z(gₙ₊₁)·(−1)ⁿ⁺¹ > 0
            return z_n * self.rs.gram_sign(n) > 0 and z_n1 * self.rs.gram_sign(n + 1) > 0
        except Exception:
            return False

    def gram_block_statistics(self, n_max: int = 100) -> Dict:
        """
        Statistische Auswertung der Gram-Blöcke bis n_max.

        Zählt: Blocks verschiedener Breite, Gram-Law-Anteil, Rosser-Anteil.

        @param n_max: Maximaler Gram-Index
        @return: Statistik-Dictionary
        @lastModified: 2026-03-12
        """
        blocks = self.identify_gram_blocks(0, n_max)
        width_dist: Dict[int, int] = {}
        good_count = 0
        for b in blocks:
            w = b["width"]
            width_dist[w] = width_dist.get(w, 0) + 1
            if b["type"] == "good":
                good_count += 1

        return {
            "total_blocks": len(blocks),
            "good_gram_intervals": good_count,
            "gram_block_count": len(blocks) - good_count,
            "width_distribution": width_dist,
            "gram_law_rate": good_count / len(blocks) if blocks else 0.0,
        }
