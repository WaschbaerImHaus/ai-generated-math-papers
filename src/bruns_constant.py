"""
@file bruns_constant.py
@brief Hochpräzise Berechnung und Analyse der Bruns Konstante B₂.
@description
    Dieses Modul berechnet und analysiert die Brun'sche Konstante B₂,
    die Summe der reziproken Zwillingsprimzahlen.

    **Bruns Konstante** (Viggo Brun, 1919):
        B₂ = Σ_{(p,p+2) Zwillingsprim} (1/p + 1/(p+2))
           = (1/3 + 1/5) + (1/5 + 1/7) + (1/11 + 1/13) + ...
           ≈ 1.9021605831...

    Die Konvergenz dieser Reihe (Bruns Theorem, 1919) ist nicht-trivial:
    Obwohl unklar ist, ob unendlich viele Zwillingsprimzahlpaare existieren,
    konvergiert die Summe ihrer Reziproken.

    **Hardy-Littlewood-Vermutung** (Conjecture B, 1923):
        π₂(x) ~ 2·C₂·x / (ln x)²   mit C₂ = Π_{p>2} p(p-2)/(p-1)² ≈ 0.6601618...

    **Bekannte Werte**:
        B₂ ≈ 1.9021605831 (Nicely, 1995; Sebah & Gourdon, 2002)
        Berechnet bis p < 10^{14}

    **Konvergenzverhalten**:
        Der Fehler bei Abschneiden bei x beträgt O(1 / (ln x)²) gemäß
        der Hardy-Littlewood-Heuristik für π₂(x).

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import time
from typing import Dict, List, Optional, Tuple

import mpmath
from mpmath import mp, mpf, log, nsum, inf


class BrunsKonstante:
    """
    @brief Berechnung und Analyse der Bruns Konstante B₂.
    @description
        Bietet Methoden zur Berechnung von B₂ mit variabler Präzision,
        Konvergenzanalyse, Primzahlzählung π₂(x) und Vergleich mit
        der Hardy-Littlewood-Vermutung.

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    # Hardy-Littlewood Konstante C₂ (Zwillingsprimzahlkonstante)
    # C₂ = Π_{p≥3 prim} p(p-2)/(p-1)²
    # Numerischer Wert auf 10 Stellen: 0.6601618158...
    HARDY_LITTLEWOOD_C2 = 0.6601618158468695739278121

    def berechne_zwillingsprimes(self, grenze: int) -> List[Tuple[int, int]]:
        """
        @brief Findet alle Zwillingsprimzahlpaare (p, p+2) mit p ≤ grenze.
        @description
            Zwillingsprimzahlen sind Paare aufeinanderfolgender Primzahlen
            mit Abstand 2:  (p, p+2) mit p und p+2 beide prim.

            Bekannte Paare: (3,5), (5,7), (11,13), (17,19), (29,31), ...

            Algorithmus: Sieb des Eratosthenes bis grenze+2, dann Paare filtern.

        @param grenze  Obere Grenze für den kleineren Partner p
        @return        Liste von Tupeln (p, p+2)
        @lastModified 2026-03-12
        """
        if grenze < 3:
            return []

        # Sieb des Eratosthenes bis grenze+2
        n = grenze + 2
        ist_prim = bytearray([1]) * (n + 1)
        ist_prim[0] = ist_prim[1] = 0

        for i in range(2, int(n ** 0.5) + 1):
            if ist_prim[i]:
                # Alle Vielfachen von i ab i² als nicht-prim markieren
                ist_prim[i * i::i] = bytearray(len(ist_prim[i * i::i]))

        # Zwillingsprimpaare sammeln
        paare = []
        for p in range(3, grenze + 1):
            if ist_prim[p] and ist_prim[p + 2]:
                paare.append((p, p + 2))

        return paare

    def bruns_summe(self, grenze: int) -> float:
        """
        @brief Berechnet die partielle Bruns Summe B₂(x) = Σ_{p≤x} (1/p + 1/(p+2)).
        @description
            B₂(x) = Σ_{(p,p+2) Zwillingsprim, p≤x} (1/p + 1/(p+2))

            Diese partielle Summe konvergiert gegen B₂ für x → ∞.

            Laufzeit: O(x / ln(x) · Siebaufwand)

        @param grenze  Obere Schranke für p
        @return        Partielle Bruns-Summe als Float
        @lastModified 2026-03-12
        """
        paare = self.berechne_zwillingsprimes(grenze)
        summe = 0.0
        for p, p2 in paare:
            summe += 1.0 / p + 1.0 / p2
        return summe

    def konvergenzanalyse(self, schranken_list: Optional[List[int]] = None) -> Dict:
        """
        @brief Zeigt die Konvergenz von B₂(x) für verschiedene obere Schranken.
        @description
            Berechnet B₂(x) für x ∈ schranken_list und zeigt den
            Konvergenzprozess.

            Erwartetes Verhalten:
            - B₂(x) ist monoton wachsend (neue Paare addieren positive Terme)
            - Inkremente: ΔB₂(x) ≈ 2·π₂'(x)/x (sehr langsam abfallend)
            - Grenzwert: B₂ ≈ 1.9021605...

        @param schranken_list  Liste von oberen Schranken (Standard: [100, 1000, 10^4, 10^5])
        @return                Dictionary mit B₂-Werten und Differenzen
        @lastModified 2026-03-12
        """
        if schranken_list is None:
            schranken_list = [100, 1_000, 10_000, 100_000]

        ergebnisse = []
        vorher = 0.0

        for x in sorted(schranken_list):
            b2_x = self.bruns_summe(x)
            anzahl_paare = len(self.berechne_zwillingsprimes(x))
            delta = b2_x - vorher
            ergebnisse.append({
                'x': x,
                'b2_x': b2_x,
                'anzahl_paare': anzahl_paare,
                'delta_zu_vorher': delta,
            })
            vorher = b2_x

        return {
            'schranken': schranken_list,
            'ergebnisse': ergebnisse,
            'bekannter_grenzwert': 1.9021605831,
            'erklärung': (
                'B₂(x) konvergiert sehr langsam gegen B₂ ≈ 1.9021605831. '
                'Der Fehler bei Abschneiden bei x ist ~O(1/(ln x)²).'
            )
        }

    def hochpraezise_berechnung(self,
                                grenze: int = 1_000_000,
                                dezimalstellen: int = 30) -> Dict:
        """
        @brief Hochpräzise Berechnung der Bruns Summe mittels mpmath.
        @description
            Berechnet B₂(x) mit `dezimalstellen` gültigen Dezimalstellen
            durch mpmath's hochpräzise Arithmetik.

            Für eine vollständig konvergierte Konstante wäre x → ∞ nötig.
            Mit Extrapolation kann B₂ auf ~12 Dezimalstellen geschätzt werden.

        @param grenze         Obere Grenze für Primzahlen (Standard: 10^6)
        @param dezimalstellen Gewünschte Präzision (Standard: 30)
        @return               Dictionary mit hochpräzisem B₂-Wert
        @lastModified 2026-03-12
        """
        # mpmath-Präzision setzen
        mp.dps = dezimalstellen + 10  # Etwas mehr als nötig für Rundungsfehler

        start = time.time()
        paare = self.berechne_zwillingsprimes(grenze)

        # Hochpräzise Summation
        summe = mpf('0')
        for p, p2 in paare:
            summe += mpf(1) / mpf(p) + mpf(1) / mpf(p2)

        laufzeit = time.time() - start

        # Zurück auf Standardpräzision
        mp.dps = 15

        return {
            'grenze': grenze,
            'dezimalstellen': dezimalstellen,
            'b2_hochpräzise': str(summe),
            'b2_float': float(summe),
            'anzahl_paare': len(paare),
            'laufzeit_sek': round(laufzeit, 3),
            'bekannter_wert': '1.9021605831040309591714',
            'erklärung': (
                f'Hochpräzise Berechnung bis p={grenze} mit {dezimalstellen} Dezimalstellen.'
            )
        }

    def meissel_lehmer_abschaetzung(self, x: float) -> Dict:
        """
        @brief Schätzung der Anzahl π₂(x) der Zwillingsprimzahlpaare bis x.
        @description
            Einfache Abschätzung für π₂(x):
                π₂(x) ≈ 2·C₂·x / (ln x)²

            Diese folgt aus der Hardy-Littlewood-Vermutung (Conjecture B).

            Bekannte Werte:
                π₂(10^4) = 205
                π₂(10^6) = 8169
                π₂(10^8) = 440312
                π₂(10^{10}) = 27412679

            Die Formel unterschätzt leicht für kleine x; für x → ∞ gilt Gleichheit
            (unter Annahme der Riemann-Hypothese).

        @param x  Obere Grenze
        @return   Dictionary mit Abschätzung und numerischen Referenzwerten
        @lastModified 2026-03-12
        """
        if x <= 2:
            return {'x': x, 'abschaetzung': 0, 'fehler': 'x zu klein'}

        ln_x = math.log(x)
        if ln_x < 0.01:
            return {'x': x, 'abschaetzung': 0, 'fehler': 'ln(x) zu klein'}

        # Hardy-Littlewood-Schätzung
        pi2_approx = 2 * self.HARDY_LITTLEWOOD_C2 * x / (ln_x ** 2)

        # Verbesserte Schätzung mit Logarithmuskorrekturen
        # (Brent/Harborth/Caldwell-Typ)
        pi2_verbessert = 2 * self.HARDY_LITTLEWOOD_C2 * x / (ln_x ** 2) * (
            1 + 2 / ln_x + (6 - 4 * math.log(math.log(x)) if x > 5 else 1) / (ln_x ** 2)
        ) if x > math.e else 0

        # Bekannte exakte Werte von π₂(x) (für Vergleich)
        bekannte_werte = {
            10**4: 205,
            10**5: 1224,
            10**6: 8169,
            10**8: 440312,
            10**10: 27412679,
        }

        return {
            'x': x,
            'hl_abschaetzung': pi2_approx,
            'verbesserte_abschaetzung': pi2_verbessert,
            'bekannte_werte': bekannte_werte,
            'c2': self.HARDY_LITTLEWOOD_C2,
            'formel': 'π₂(x) ≈ 2·C₂·x / (ln x)²',
        }

    def hardylittlewood_vorhersage(self, x: float) -> Dict:
        """
        @brief Berechnet die Hardy-Littlewood-Vorhersage π₂(x) ~ 2·C₂·x/(ln x)².
        @description
            Hardy-Littlewood Conjecture B (1923):
                π₂(x) ~ 2·C₂·x / (ln x)²  für x → ∞

            mit der Zwillingsprimzahlkonstante:
                C₂ = Π_{p>2 prim} p(p-2)/(p-1)²

            Die Faktoren p(p-2)/(p-1)² sind Korrekturen gegenüber dem Sieb-
            Faktor 1/2 für Dichteabschätzungen in arithmetischen Progressionen.

            Für p=3: 3·1/(4) = 3/4
            Für p=5: 5·3/(16) = 15/16
            Für p=7: 7·5/(36) = 35/36
            ...

            Berechnung von C₂ als Teilprodukt über die ersten 1000 Primzahlen.

        @param x  Obere Grenze für die Vorhersage
        @return   Dictionary mit Hardy-Littlewood-Vorhersage und C₂-Details
        @lastModified 2026-03-12
        """
        # C₂ berechnen als Teilprodukt (konvergiert schnell)
        from sympy import primerange as pr
        c2 = 1.0
        for p in pr(3, 1000):
            c2 *= p * (p - 2) / ((p - 1) ** 2)

        ln_x = math.log(x) if x > 1 else 1.0
        vorhersage = 2 * c2 * x / (ln_x ** 2)

        # Erste Faktoren des Produkts
        erste_faktoren = []
        produkt_laufend = 1.0
        for p in list(pr(3, 30)):
            faktor = p * (p - 2) / ((p - 1) ** 2)
            produkt_laufend *= faktor
            erste_faktoren.append({
                'p': p,
                'faktor': round(faktor, 8),
                'teilprodukt': round(produkt_laufend, 8),
            })

        return {
            'x': x,
            'c2_berechnet': c2,
            'c2_bekannt': self.HARDY_LITTLEWOOD_C2,
            'vorhersage_pi2_x': vorhersage,
            'erste_faktoren': erste_faktoren,
            'formel': 'π₂(x) ~ 2·C₂·x/(ln x)², C₂ = Π_{p>2} p(p-2)/(p-1)²',
            'status': 'VERMUTUNG — nicht bewiesen',
        }

    def vergleich_bruns_konstante_approximation(self) -> Dict:
        """
        @brief Vergleicht bekannte Schätzwerte und Approximationen für B₂.
        @description
            B₂ ist bekannt auf ~12 Dezimalstellen.

            Methoden zur Approximation:
            1. Direkte Summation (langsame Konvergenz)
            2. Euler-Maclaurin-Extrapolation
            3. Richardson-Extrapolation aus B₂(x₁), B₂(x₂), B₂(x₃)
            4. Asymptotische Korrektur via Hardy-Littlewood

            Die asymptotische Korrektur:
                B₂ ≈ B₂(x) + Σ_{p>x} (1/p + 1/(p+2))
                ≈ B₂(x) + 4·C₂ · ∫_x^∞ 1/(t · (ln t)²) dt
                ≈ B₂(x) + 4·C₂ / ln(x)    [grobe Näherung]

        @return   Dictionary mit verschiedenen Approximationsmethoden
        @lastModified 2026-03-12
        """
        import mpmath as mpm

        # Direkte Summation für verschiedene x
        grenzen = [10_000, 50_000, 100_000]
        summenwerte = {}
        for g in grenzen:
            summenwerte[g] = self.bruns_summe(g)

        # Richardson-Extrapolation: nutze B₂(x₁), B₂(x₂)
        # B₂ ≈ B₂(x) + C/ln(x) (asymptotisches Modell)
        # Aus B₂(x₁) = B₂ - C/ln(x₁) und B₂(x₂) = B₂ - C/ln(x₂):
        # B₂ ≈ (B₂(x₂)·ln(x₂) - B₂(x₁)·ln(x₁)) / (ln(x₂) - ln(x₁))
        x1, x2 = 50_000, 100_000
        b1, b2 = summenwerte[x1], summenwerte[x2]
        ln1, ln2 = math.log(x1), math.log(x2)
        richardson = (b2 * ln2 - b1 * ln1) / (ln2 - ln1)

        # Asymptotische Korrektur
        # ∫_x^∞ 1/(t·(ln t)²) dt = 1/ln(x) + ...
        g = 100_000
        korrektur = 4 * self.HARDY_LITTLEWOOD_C2 / math.log(g)
        asymptotisch = summenwerte[g] + korrektur

        return {
            'direkte_summation': summenwerte,
            'richardson_extrapolation': {
                'wert': richardson,
                'x1': x1, 'x2': x2,
            },
            'asymptotische_korrektur': {
                'b2_x': summenwerte[g],
                'korrektur': korrektur,
                'b2_approx': asymptotisch,
            },
            'bekannter_wert': 1.9021605831,
            'beste_schätzung': richardson,
            'erklärung': (
                'B₂ ≈ 1.9021605831 (Nicely 1995, Sebah-Gourdon 2002). '
                'Richardson-Extrapolation verbessert die direkte Summation erheblich.'
            )
        }


# ===========================================================================
# HAUPTPROGRAMM (Demo)
# ===========================================================================

if __name__ == "__main__":
    """
    Demonstriert die BrunsKonstante-Klasse.
    """
    bk = BrunsKonstante()

    print("=" * 65)
    print("BRUNS KONSTANTE B₂ — Analyse und Berechnung")
    print("=" * 65)

    # 1. Erste Zwillingsprimpaare
    print("\n1. Erste Zwillingsprimzahlpaare (bis 100):")
    paare = bk.berechne_zwillingsprimes(100)
    print(f"   {paare}")
    print(f"   Anzahl: {len(paare)}")

    # 2. Konvergenzanalyse
    print("\n2. Konvergenzanalyse von B₂(x):")
    ka = bk.konvergenzanalyse([100, 1_000, 10_000, 100_000])
    for r in ka['ergebnisse']:
        print(f"   x={r['x']:8d}: B₂(x)={r['b2_x']:.10f}  "
              f"(Paare: {r['anzahl_paare']:6d}, Δ={r['delta_zu_vorher']:.8f})")
    print(f"   Bekannter Grenzwert: {ka['bekannter_grenzwert']}")

    # 3. Hochpräzise Berechnung
    print("\n3. Hochpräzise Berechnung bis p=10^6 (30 Dezimalstellen):")
    hp = bk.hochpraezise_berechnung(grenze=1_000_000, dezimalstellen=30)
    print(f"   B₂(10^6) = {hp['b2_hochpräzise']}")
    print(f"   Bekannt:   {hp['bekannter_wert']}")
    print(f"   Laufzeit:  {hp['laufzeit_sek']} Sek.")

    # 4. B₂ für x = 10^4, 10^5, 10^6, 10^7 (partielle Summen)
    print("\n4. B₂(x) für x = 10^4, 10^5, 10^6:")
    for exp in [4, 5, 6]:
        x = 10 ** exp
        b2_x = bk.bruns_summe(x)
        paare_x = len(bk.berechne_zwillingsprimes(x))
        print(f"   B₂(10^{exp}) = {b2_x:.12f}  (Paare: {paare_x})")

    # 5. Hardy-Littlewood-Vorhersage
    print("\n5. Hardy-Littlewood-Vorhersage π₂(x):")
    for exp in [4, 6, 8, 10]:
        x = 10 ** exp
        hl = bk.meissel_lehmer_abschaetzung(float(x))
        print(f"   π₂(10^{exp}) ≈ {hl['hl_abschaetzung']:.0f}")

    # 6. C₂-Berechnung
    print("\n6. Hardy-Littlewood-Konstante C₂:")
    hl_detail = bk.hardylittlewood_vorhersage(1e8)
    print(f"   C₂ (berechnet, erste 1000 Primes) = {hl_detail['c2_berechnet']:.10f}")
    print(f"   C₂ (bekannt)                      = {hl_detail['c2_bekannt']:.10f}")
    print(f"   Status: {hl_detail['status']}")

    # 7. Approximationsvergleich
    print("\n7. Approximationsvergleich für B₂:")
    approx = bk.vergleich_bruns_konstante_approximation()
    print(f"   Direkte Summation (10^5): {approx['direkte_summation'][100_000]:.10f}")
    print(f"   Richardson-Extrapolation:  {approx['richardson_extrapolation']['wert']:.10f}")
    print(f"   Asymptotische Korrektur:   {approx['asymptotische_korrektur']['b2_approx']:.10f}")
    print(f"   Bekannter Wert:            {approx['bekannter_wert']:.10f}")
