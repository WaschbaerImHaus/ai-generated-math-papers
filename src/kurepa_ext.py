"""
@file kurepa_ext.py
@brief Erweiterte Analyse der Kurepa-Vermutung (1950): !p ≢ 0 (mod p).
@description
    Dieses Modul untersucht die Kurepa-Vermutung vertieft.

    **Kurepa-Vermutung** (Đuro Kurepa, 1950):
        Die Linksfakultät !p = 0! + 1! + 2! + ... + (p-1)! ist für keine
        Primzahl p durch p teilbar.

    Formal: ∀ Primzahl p gilt: !p ≢ 0 (mod p)

    Die Vermutung wurde numerisch verifiziert für alle Primzahlen p < 10^6.
    Ein Beweis ist bislang nicht bekannt.

    Inhalte:
    - Berechnung der Linksfakultät !n und !p mod p
    - Numerische Verifikation bis p ≤ p_max
    - p-adische Bewertungsanalyse
    - Restklassenstruktur nach verschiedenen Moduln
    - Verbindung zur Wilson-Bedingung (p-1)! ≡ -1 (mod p)
    - Rekursive Struktur: !(n) = !(n-1) + (n-1)!
    - Kandidatensuche: Welche p haben !p "nahe" bei 0 mod p?

    Mathematischer Hintergrund:
    - Beziehung zu Wilson: (p-1)! ≡ -1 (mod p) für Primzahlen p
    - Rekursion: !p = !(p-1) + (p-1)! = !(p-1) - 1 (mod p) nach Wilson
    - Daraus: !p ≡ !2 + Σ_{k=2}^{p-1} [Δ-Korrektur] (mod p)
    - p-adische Bewertung v_p(!p) wäre ≥ 1, falls die Vermutung falsch wäre

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import time
from functools import lru_cache
from typing import Dict, List, Tuple, Optional

import sympy
from sympy import isprime, nextprime, primerange


class KurepaExt:
    """
    @brief Erweiterte Analyse der Kurepa-Vermutung.
    @description
        Bietet numerische Verifikation, p-adische Analyse, Restklassenstruktur
        und Verbindungen zu anderen bekannten Ergebnissen (Wilson, Rekursion).

        Die Kurepa-Vermutung besagt, dass die Linksfakultät
            !p = Σ_{i=0}^{p-1} i!
        niemals durch p teilbar ist, wenn p eine Primzahl ist.

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def berechne_leftfakultaet(self, n: int) -> int:
        """
        @brief Berechnet die Linksfakultät !n = 0! + 1! + 2! + ... + (n-1)!
        @description
            Die Linksfakultät ist definiert als:
                !n = Σ_{i=0}^{n-1} i!

            Werte: !1=1, !2=2, !3=4, !4=10, !5=34, !6=154, !7=874, ...

            Berechnung erfolgt iterativ ohne erneute Auswertung von i! durch
            Akkumulation: akt_fak = (i-1)! → i! = i · (i-1)!

        @param n  Obere Grenze (nicht inklusive): Summe läuft von 0! bis (n-1)!
        @return   Wert der Linksfakultät !n als ganzzahligen Wert
        @lastModified 2026-03-12
        """
        if n <= 0:
            # Leere Summe = 0
            return 0
        if n == 1:
            # !1 = 0! = 1
            return 1

        # Iterative Berechnung: Akkumuliere !n und aktuellen Fakultätswert
        gesamt = 1      # Startet mit 0! = 1
        akt_fak = 1     # 0! = 1

        for i in range(1, n):
            akt_fak *= i       # i! = i · (i-1)!
            gesamt += akt_fak  # !n += i!

        return gesamt

    def kurepa_restklasse(self, p: int) -> int:
        """
        @brief Berechnet !p mod p für eine gegebene (Prim-)Zahl p.
        @description
            Gibt den Rest der Linksfakultät !p bei Division durch p zurück.

            Für die Kurepa-Vermutung interessiert uns: Ist dieser Rest = 0?

            Berechnung modular, um sehr große Zahlen zu vermeiden:
                !p ≡ Σ_{i=0}^{p-1} i! (mod p)

            Da i! ≡ 0 (mod p) für alle i ≥ p, reichen die Terme bis i = p-1.
            Ab i = p gilt i! = i · (i-1)! ≡ 0 · (i-1)! = 0 (mod p).

        @param p  Eine Primzahl (oder beliebige natürliche Zahl ≥ 1)
        @return   !p mod p (ganzzahlig in {0, ..., p-1})
        @lastModified 2026-03-12
        """
        if p <= 0:
            raise ValueError(f"p muss positiv sein, erhalten: {p}")

        # Modular rechnen
        rest = 0
        akt_fak_mod = 1  # 0! mod p = 1

        for i in range(p):
            if i == 0:
                # 0! = 1
                rest = (rest + 1) % p
            else:
                # i! = i * (i-1)!
                akt_fak_mod = (akt_fak_mod * i) % p
                rest = (rest + akt_fak_mod) % p

        return rest

    def numerische_verifikation(self, p_max: int) -> Dict:
        """
        @brief Überprüft die Kurepa-Vermutung für alle Primzahlen p ≤ p_max.
        @description
            Für jede Primzahl p ≤ p_max wird !p mod p berechnet.
            Kein Ergebnis = 0 bedeutet: Vermutung gilt in diesem Bereich.

            Laufzeit: O(p_max · π(p_max)) ≈ O(p_max² / ln(p_max))
            Für p_max = 50000 ca. einige Sekunden.

        @param p_max  Obere Grenze für die Verifikation
        @return       Dictionary mit:
                      - 'p_max': verwendete Grenze
                      - 'anzahl_geprüft': Anzahl geprüfter Primzahlen
                      - 'gegenbeispiele': Liste von Primzahlen mit !p ≡ 0 (mod p)
                      - 'verifiziert': True falls kein Gegenbeispiel gefunden
                      - 'laufzeit_sek': Berechnungszeit
                      - 'erste_reste': Liste der ersten 20 (!p mod p)-Werte
        @lastModified 2026-03-12
        """
        start_zeit = time.time()
        gegenbeispiele: List[int] = []
        erste_reste: List[Tuple[int, int]] = []

        anzahl = 0
        for p in primerange(2, p_max + 1):
            rest = self.kurepa_restklasse(p)
            anzahl += 1
            if rest == 0:
                gegenbeispiele.append(p)
            if len(erste_reste) < 20:
                erste_reste.append((p, rest))

        laufzeit = time.time() - start_zeit

        return {
            'p_max': p_max,
            'anzahl_geprüft': anzahl,
            'gegenbeispiele': gegenbeispiele,
            'verifiziert': len(gegenbeispiele) == 0,
            'laufzeit_sek': round(laufzeit, 4),
            'erste_reste': erste_reste,
            'schluss': (
                f'Kurepa-Vermutung verifiziert für alle {anzahl} Primzahlen p ≤ {p_max}.'
                if not gegenbeispiele
                else f'GEGENBEISPIEL GEFUNDEN: {gegenbeispiele}'
            )
        }

    def p_adische_bewertung_analyse(self, p: int) -> Dict:
        """
        @brief Analysiert die p-adische Bewertung v_p(!p) der Linksfakultät.
        @description
            Die p-adische Bewertung v_p(n) ist der größte Exponent k mit p^k | n.

            Kurepa-Vermutung ⟺ v_p(!p) = 0 für alle Primzahlen p.

            Zusätzlich analysieren wir:
            - v_p(k!) für k = 0, ..., p-1 via Legendre-Formel
            - Beitrag jedes Summanden zur p-adischen Bewertung
            - Die Term-Struktur der Summe

            Legendres Formel: v_p(k!) = Σ_{j=1}^∞ ⌊k / p^j⌋

        @param p  Eine Primzahl
        @return   Dictionary mit p-adischen Bewertungen der Summanden und Gesamtwert
        @lastModified 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        def v_p(n: int) -> int:
            """p-adische Bewertung von n (mit Sonderfällen n=0)."""
            if n == 0:
                return float('inf')  # type: ignore[return-value]
            count = 0
            while n % p == 0:
                n //= p
                count += 1
            return count

        def legendre_fak(k: int) -> int:
            """v_p(k!) via Legendres Formel."""
            result = 0
            pk = p
            while pk <= k:
                result += k // pk
                pk *= p
            return result

        # Bewertungen jedes Summanden i! in !p
        summanden_bewertung = []
        akt_fak = 1
        for i in range(p):
            if i == 0:
                summanden_bewertung.append({'i': 0, 'i_fak': 1, 'v_p': 0})
            else:
                akt_fak *= i
                bew = legendre_fak(i)
                summanden_bewertung.append({
                    'i': i,
                    'i_fak_vp': bew,
                    'v_p_legendre': bew
                })

        # !p tatsächlich berechnen (nur für kleine p)
        if p <= 200:
            left_fak = self.berechne_leftfakultaet(p)
            vp_total = v_p(left_fak)
        else:
            # Modular-Berechnung für v_p
            rest = self.kurepa_restklasse(p)
            vp_total = 0 if rest != 0 else 'unknown (≥1)'
            left_fak = None

        # Minimale Bewertung unter den Summanden — bestimmt v_p der Summe (falls eindeutig)
        min_bew = min(s['v_p_legendre'] for s in summanden_bewertung[1:]) if p > 1 else 0

        return {
            'p': p,
            'vp_left_fak': vp_total,
            'kurepa_gilt': vp_total == 0 if isinstance(vp_total, int) else None,
            'min_summanden_bewertung': min_bew,
            'summanden_detail': summanden_bewertung[:min(15, p)],
            'legendre_p_minus_1_fak': legendre_fak(p - 1),
            'erklärung': (
                'v_p(!p) = 0 bedeutet: p teilt !p nicht → Kurepa gilt für dieses p. '
                'Nach Wilson: v_p((p-1)!) = 0, da (p-1)! ≡ -1 (mod p) für p prim.'
            )
        }

    def restklassen_struktur(self) -> Dict:
        """
        @brief Analysiert !p mod p für p in verschiedenen Restklassen mod k.
        @description
            Untersucht, ob die Restklasse von p mod k das Verhalten von !p mod p
            beeinflusst. Dies könnte Muster aufdecken, die einen Beweis unterstützen.

            Analysiert werden:
            - p mod 3, p mod 4, p mod 6, p mod 8, p mod 12, p mod 24

        @return   Dictionary mit statistischer Auswertung nach Restklassen
        @lastModified 2026-03-12
        """
        moduln = [3, 4, 6, 8, 12, 24]
        primes = list(primerange(3, 500))  # Erste 97 ungerade Primzahlen bis 499

        ergebnis = {}
        for m in moduln:
            # Gruppiere Primzahlen nach ihrer Restklasse mod m
            gruppen: Dict[int, List[Tuple[int, int]]] = {}
            for p in primes:
                rest_p_mod_m = p % m
                if rest_p_mod_m not in gruppen:
                    gruppen[rest_p_mod_m] = []
                kurepa_rest = self.kurepa_restklasse(p)
                gruppen[rest_p_mod_m].append((p, kurepa_rest))

            # Statistik pro Restklasse
            statistik = {}
            for rk, liste in sorted(gruppen.items()):
                reste = [lr for _, lr in liste]
                statistik[rk] = {
                    'anzahl_primzahlen': len(liste),
                    'kurepa_reste': reste[:10],  # Erste 10 zeigen
                    'durchschnitt': sum(reste) / len(reste) if reste else 0,
                    'nullen': reste.count(0),    # Falls != 0 → Gegenbeispiel
                }
            ergebnis[m] = statistik

        return {
            'beschreibung': 'Analyse von !p mod p nach Restklassen von p',
            'primzahlen_bis': 499,
            'ergebnis_nach_modul': ergebnis,
        }

    def wilsons_verbindung(self, p: int) -> Dict:
        """
        @brief Verbindung zwischen Kurepa-Linksfakultät und dem Satz von Wilson.
        @description
            **Satz von Wilson**: Für eine Primzahl p gilt (p-1)! ≡ -1 (mod p).

            **Rekursionsstruktur**:
                !p = !(p-1) + (p-1)!
                ≡ !(p-1) + (-1)   (mod p, nach Wilson)
                ≡ !(p-1) - 1      (mod p)

            Also gilt die Rekursion:
                !p ≡ !(p-1) - 1 (mod p)

            und allgemein für Primzahlen p:
                !p ≡ !(q) - (p-q) ... (komplexer für nicht-konsekutive Primzahlen)

            Gezeigt wird auch:
                !p = Σ_{i=0}^{p-1} i!
                = (0! + 1! + ... + (p-2)!) + (p-1)!
                = !(p-1) + (p-1)!

            Für p prim: (p-1)! ≡ -1 (mod p)
            Daher:      !p ≡ !(p-1) - 1 (mod p)

        @param p  Eine Primzahl
        @return   Dictionary mit Verbindungsanalyse
        @lastModified 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        # Direkte Berechnung
        kurepa_p = self.kurepa_restklasse(p)
        kurepa_p_minus_1 = self.kurepa_restklasse(p)  # Wir brauchen !(p-1) mod p

        # !(p-1) mod p berechnen
        left_fak_p_minus_1 = 0
        akt_fak_mod = 1  # 0!
        for i in range(p - 1):
            if i == 0:
                left_fak_p_minus_1 = (left_fak_p_minus_1 + 1) % p
            else:
                akt_fak_mod = (akt_fak_mod * i) % p
                left_fak_p_minus_1 = (left_fak_p_minus_1 + akt_fak_mod) % p

        # (p-1)! mod p — Wilson: sollte -1 ≡ p-1 (mod p) sein
        wilson_wert = math.factorial(p - 1) % p

        # Verifikation der Rekursion
        rekursion_stimmt = (kurepa_p == (left_fak_p_minus_1 + wilson_wert) % p)

        return {
            'p': p,
            'kurepa_p_mod_p': kurepa_p,
            'kurepa_p_minus_1_mod_p': left_fak_p_minus_1,
            'wilson_p_minus_1_fak_mod_p': wilson_wert,
            'wilson_gilt': wilson_wert == p - 1,
            'rekursion_verifiziert': rekursion_stimmt,
            'formel': f'!{p} ≡ !{p-1} + ({p}-1)! ≡ !{p-1} + ({wilson_wert}) ≡ {kurepa_p} (mod {p})',
            'erklärung': (
                f'Nach Wilson: ({p}-1)! ≡ -1 ≡ {p-1} (mod {p}). '
                f'Rekursion: !{p} = !{p-1} + ({p}-1)! '
                f'≡ {left_fak_p_minus_1} + {wilson_wert} = {(left_fak_p_minus_1 + wilson_wert) % p} (mod {p}).'
            )
        }

    def rekursive_analyse(self, n: int) -> Dict:
        """
        @brief Nutzt die Rekursion !n = !(n-1) + (n-1)! zur Strukturanalyse.
        @description
            Die Linksfakultät genügt der Rekursion:
                !(n) = !(n-1) + (n-1)!

            Dies erlaubt eine schrittweise Analyse des Wertverlaufs.

            Für Primzahlen q < p gibt (q-1)! ≡ -1 (mod q) nach Wilson
            eine direkte Beziehung zwischen aufeinanderfolgenden Primzahlen.

            Wir berechnen !k für k = 1, ..., n und zeigen die Struktur.

        @param n  Obergrenze der Analyse
        @return   Dictionary mit Schrittfolge und Analyse
        @lastModified 2026-03-12
        """
        schritte = []
        left_fak = 0   # !(0) = 0 (leere Summe)
        akt_fak = 1    # 0! = 1

        for k in range(1, n + 1):
            if k == 1:
                # !(1) = 0! = 1
                left_fak = 1
                schritte.append({
                    'k': k,
                    'left_fak': left_fak,
                    'k_minus_1_fak': 1,
                    'ist_prim': isprime(k),
                    'kurepa_rest_wenn_prim': left_fak % k if isprime(k) else None,
                })
            else:
                # !k = !(k-1) + (k-1)!
                # (k-1)! = (k-1) * (k-2)!
                akt_fak *= (k - 1)  # Jetzt (k-1)!
                left_fak += akt_fak
                ist_prim = isprime(k)
                kurepa_rest = left_fak % k if ist_prim else None

                schritte.append({
                    'k': k,
                    'left_fak': left_fak,
                    'k_minus_1_fak': akt_fak,
                    'ist_prim': ist_prim,
                    'kurepa_rest_wenn_prim': kurepa_rest,
                    'kurepa_verletzt': kurepa_rest == 0 if ist_prim else None,
                })

        # Zusammenfassung
        primschritte = [s for s in schritte if s['ist_prim']]
        verletzungen = [s for s in primschritte if s.get('kurepa_verletzt', False)]

        return {
            'n': n,
            'schritte': schritte[:50],  # Erste 50 Schritte
            'prim_schritte': primschritte[:20],
            'kurepa_verletzungen': verletzungen,
            'verifiziert': len(verletzungen) == 0,
        }

    def suche_kandidaten_Wilson(self, p_max: int) -> Dict:
        """
        @brief Sucht Primzahlen, bei denen !p mod p "nahe bei 0" liegt.
        @description
            Findet Primzahlen p ≤ p_max mit kleinem relativen Rest:
                r = (!p mod p) / p ≈ 0

            Diese sind keine Gegenbeispiele, zeigen aber, wo die Vermutung
            "knapp" gilt. Solche p könnten Muster offenbaren.

            Aufgelistet werden Primzahlen p mit !p mod p ≤ δ · p für δ = 0.01
            (d.h. die ersten 1% sind null-nah).

        @param p_max  Obere Grenze
        @return       Dictionary mit sortierten "Kandidaten" nach relativem Rest
        @lastModified 2026-03-12
        """
        # Schwellwert: relativer Rest < 1% oder > 99% (nahe 0 oder nahe p)
        kandidaten_nahe_null: List[Dict] = []
        alle_reste: List[Tuple[int, int]] = []

        for p in primerange(2, p_max + 1):
            rest = self.kurepa_restklasse(p)
            alle_reste.append((p, rest))
            relativer_rest = rest / p
            if relativer_rest < 0.01:
                # Sehr nahe an 0 mod p (aber ≠ 0, sonst Gegenbeispiel)
                kandidaten_nahe_null.append({
                    'p': p,
                    'rest': rest,
                    'relativer_rest': relativer_rest,
                    'ist_gegenbeispiel': rest == 0,
                })

        # Sortiere nach relativem Rest aufsteigend
        kandidaten_nahe_null.sort(key=lambda x: x['relativer_rest'])

        # Minimaler Rest absolut (normalisiert auf [0, p/2])
        reste_normiert = [(p, min(r, p - r)) for p, r in alle_reste]
        reste_normiert.sort(key=lambda x: x[1])

        return {
            'p_max': p_max,
            'kandidaten_nahe_null': kandidaten_nahe_null[:20],
            'global_minimaler_rest': reste_normiert[:10],
            'anzahl_kandidaten': len(kandidaten_nahe_null),
            'hinweis': (
                'Kein Gegenbeispiel gefunden. Kandidaten mit sehr kleinem '
                'relativem Rest könnten Ausgangspunkte für theoretische Analyse sein.'
            )
        }


# ===========================================================================
# HAUPTPROGRAMM (Demo)
# ===========================================================================

if __name__ == "__main__":
    """
    Demonstriert die Funktionen der KurepaExt-Klasse.
    """
    k = KurepaExt()

    print("=" * 60)
    print("KUREPA-VERMUTUNG — Erweiterte Analyse")
    print("=" * 60)

    # 1. Erste Linksfakultäten
    print("\n1. Linksfakultäten !n für n = 1..10:")
    for n in range(1, 11):
        print(f"   !{n} = {k.berechne_leftfakultaet(n)}")

    # 2. !p mod p für die ersten Primzahlen
    print("\n2. !p mod p für kleine Primzahlen:")
    primes_demo = list(primerange(2, 30))
    for p in primes_demo:
        rest = k.kurepa_restklasse(p)
        print(f"   p={p}: !p mod p = {rest}  {'← GEGENBEISPIEL!' if rest == 0 else 'OK'}")

    # 3. Wilson-Verbindung für p=7
    print("\n3. Wilson-Verbindung für p=7:")
    w = k.wilsons_verbindung(7)
    print(f"   {w['formel']}")
    print(f"   Wilson gilt: {w['wilson_gilt']}, Rekursion OK: {w['rekursion_verifiziert']}")

    # 4. Rekursive Analyse bis n=15
    print("\n4. Rekursive Schrittanalyse bis n=15:")
    ra = k.rekursive_analyse(15)
    for s in ra['prim_schritte']:
        print(f"   p={s['k']}: !p={s['left_fak']}, !p mod p={s['kurepa_rest_wenn_prim']}")

    # 5. Numerische Verifikation bis p=50000
    print("\n5. Numerische Verifikation bis p=50000...")
    erg = k.numerische_verifikation(50000)
    print(f"   Geprüft: {erg['anzahl_geprüft']} Primzahlen")
    print(f"   Verifiziert: {erg['verifiziert']}")
    print(f"   Laufzeit: {erg['laufzeit_sek']} Sek.")
    print(f"   Gegenbeispiele: {erg['gegenbeispiele']}")

    # 6. p-adische Analyse für p=7
    print("\n6. p-adische Bewertungsanalyse für p=7:")
    pa = k.p_adische_bewertung_analyse(7)
    print(f"   v_7(!7) = {pa['vp_left_fak']}, Kurepa gilt: {pa['kurepa_gilt']}")

    # 7. Kandidatensuche
    print("\n7. Kandidaten mit kleinem !p mod p (bis p=1000):")
    kand = k.suche_kandidaten_Wilson(1000)
    print(f"   Kandidaten (rel. Rest < 1%): {kand['anzahl_kandidaten']}")
    for kk in kand['kandidaten_nahe_null'][:5]:
        print(f"   p={kk['p']}: rest={kk['rest']}, rel={kk['relativer_rest']:.6f}")
