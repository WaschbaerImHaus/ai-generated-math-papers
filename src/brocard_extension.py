"""
@file brocard_extension.py
@brief Erweiterte numerische und modulare Analyse der Brocard-Ramanujan-Vermutung.
@description
    Dieses Modul untersucht die Brocard-Ramanujan-Vermutung (1876/1913):

    **Vermutung**: Die Gleichung n! + 1 = m² hat nur die Lösungen
        (n, m) ∈ {(4, 5), (5, 11), (7, 71)}.

    **Bekannte Fakten**:
    - Verifiziert bis n = 10^9 (Berndt und Galway, 2000)
    - Die Lösungen bei n=4,5,7 sind seit Brocard (1876) und Ramanujan (1913) bekannt
    - Für n ≥ 8 gibt es verschiedene modulare Ausschlussargumente

    **Inhalt dieses Moduls**:
    1. Modulare Ausschlussanalyse (welche Moduli n!+1 als Nicht-Quadrat ausweisen)
    2. Numerische Verifikation bis n = 10000 (mit exakter Ganzzahlarithmetik)
    3. Schrankenanalyse für m in Abhängigkeit von n
    4. p-adische Bewertungsanalyse
    5. Suche nach Ausschluss-Kandidaten für einzelne n

    **Wichtige Einschränkung**:
    Für kleines n (n < prime) schlagen viele modulare Ausschlüsse fehl,
    da n! viele Nullen enthält. Die Methode wird effektiver für große n.

    **Mathematische Grundlage** (Quadratreste):
    m² mod p kann nur Werte aus der Menge der quadratischen Reste mod p annehmen.
    Für eine Primzahl p gibt es genau (p+1)/2 quadratische Reste mod p
    (inkl. 0). Wenn n!+1 mod p kein quadratischer Rest ist, dann ist
    n!+1 ≠ m² für alle ganzen Zahlen m.

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import time
from typing import Dict, List, Optional, Set, Tuple

import sympy
from sympy import factorint, isprime


# ===========================================================
# Hilfsfunktionen
# ===========================================================

def _quadratische_reste(m: int) -> Set[int]:
    """
    @brief Berechnet die Menge aller quadratischen Reste mod m.
    @description
        Quadratische Reste mod m: {k² mod m | k = 0, 1, ..., m-1}.
    @param m Modulus (positive ganze Zahl ≥ 2).
    @return Menge der quadratischen Reste mod m.
    @lastModified 2026-03-12
    """
    return {(k * k) % m for k in range(m)}


def _ist_quadrat(n: int) -> bool:
    """
    @brief Prüft ob n eine perfekte Quadratzahl ist (exakt, ohne Gleitkomma).
    @description
        Nutzt math.isqrt() für exakte Ganzzahl-Wurzel.
        Für sehr große n wird sympy.sqrt genutzt.
    @param n Ganze Zahl ≥ 0.
    @return True wenn n = k² für ein k ∈ ℕ₀.
    @lastModified 2026-03-12
    """
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def _fakultaet_mod(n: int, m: int) -> int:
    """
    @brief Berechnet n! mod m effizient.
    @description
        Für n ≥ m gilt n! ≡ 0 (mod m) wenn m eine zusammengesetzte Zahl ist
        und m ≤ n, da dann alle Faktoren von m in der Fakultät enthalten sind.
        Für Primzahlmoduln p > n muss man alle Faktoren berechnen.
    @param n Argument der Fakultät.
    @param m Modulus.
    @return n! mod m.
    @lastModified 2026-03-12
    """
    if m == 1:
        return 0
    # Für n ≥ m ist n! durch m teilbar (wenn m ≤ n, da m als Faktor in n! auftritt
    # oder m zusammengesetzt und seine Primfaktoren in n! enthalten sind)
    # Für Primzahl p: p | n! für n ≥ p
    ergebnis = 1
    for k in range(1, n + 1):
        ergebnis = (ergebnis * k) % m
        if ergebnis == 0:
            # Alle weiteren Faktoren bringen nichts mehr
            return 0
    return ergebnis


# ===========================================================
# Hauptklasse: BrocardExtension
# ===========================================================

class BrocardExtension:
    """
    @brief Erweiterte Analyse der Brocard-Ramanujan-Gleichung n! + 1 = m².
    @description
        Implementiert modulare Ausschlussanalyse, numerische Suche und
        theoretische Schranken für die Brocard-Ramanujan-Vermutung.

        **Bekannte Lösungen**: (n, m) ∈ {(4, 5), (5, 11), (7, 71)}.

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    # Bekannte Lösungen (Brocard 1876, Ramanujan 1913)
    BEKANNTE_LOESUNGEN: Dict[int, int] = {4: 5, 5: 11, 7: 71}

    # -------------------------------------------------------
    # Modulare Ausschlussanalyse
    # -------------------------------------------------------

    def modular_ausschluss(self, n: int, moduli_list: List[int]) -> Dict:
        """
        @brief Schließt n! + 1 = m² durch modulare Argumente aus.
        @description
            Für jeden Modulus m aus moduli_list:
            1. Berechne r = n! + 1 (mod m) — exakt via _fakultaet_mod.
            2. Bestimme QR(m) = Menge der quadratischen Reste mod m.
            3. Falls r ∉ QR(m): n ist als Lösung ausgeschlossen.

            **Wichtig**: Ein Ausschluss durch einen Modulus genügt.
            Falls kein Modulus ausschließt, ist das KEIN Beweis dass n Lösung ist.

        @param n Kandidat für die Brocard-Gleichung.
        @param moduli_list Liste von Moduli zum Testen.
        @return Dictionary mit Ausschlussinformation pro Modulus.
        @lastModified 2026-03-12
        """
        ergebnisse = {
            'n': n,
            'moduli': {},
            'ausgeschlossen': False,
            'ausschluss_modulus': None,
            'ausschluss_rest': None,
        }

        for m in moduli_list:
            if m < 2:
                continue

            # Berechne n! mod m
            fak_mod = _fakultaet_mod(n, m)
            # n! + 1 mod m
            rest = (fak_mod + 1) % m
            # Quadratische Reste mod m
            qr = _quadratische_reste(m)

            ist_qr = rest in qr
            ergebnisse['moduli'][m] = {
                'n_fak_mod': fak_mod,
                'n_fak_plus1_mod': rest,
                'quadratische_reste': sorted(qr),
                'rest_ist_quadrat': ist_qr,
                'ausschluss': not ist_qr,
            }

            if not ist_qr and not ergebnisse['ausgeschlossen']:
                ergebnisse['ausgeschlossen'] = True
                ergebnisse['ausschluss_modulus'] = m
                ergebnisse['ausschluss_rest'] = rest

        return ergebnisse

    def analysiere_restklassen(self, max_mod: int) -> Dict:
        """
        @brief Findet Moduli, die n!+1 ≢ □ (mod m) für viele n zeigen.
        @description
            Analysiert für Moduli m von 2 bis max_mod:
            - Welche Werte nimmt n! + 1 mod m an (für n ≥ 8)?
            - Welche dieser Werte sind quadratische Nicht-Reste?

            Ein Modulus m ist "nützlich", wenn für viele n ≥ 8 gilt:
            n! + 1 ≢ □ (mod m).

            **Beobachtung**: Für n ≥ p (p Primzahl) gilt p | n!, also
            n! ≡ 0 (mod p) und n! + 1 ≡ 1 (mod p).
            Da 1 immer ein quadratischer Rest ist, kann p kein Ausschluss
            für n ≥ p liefern.

            **Folgerung**: Nützliche Moduli sind p² (Primzahlquadrate)
            oder zusammengesetzte Zahlen, die nur Primteiler > n haben —
            also muss man n und m gemeinsam betrachten.

        @param max_mod Obere Schranke für Moduli.
        @return Statistik über nützliche Moduli.
        @lastModified 2026-03-12
        """
        ergebnis = {
            'max_mod': max_mod,
            'moduli_analyse': {},
            'nuetzliche_moduli': [],
        }

        for m in range(2, max_mod + 1):
            qr = _quadratische_reste(m)
            nicht_qr = set(range(m)) - qr

            # Für welche n (8 bis 50) schließt m aus?
            ausschluss_fuer = []
            for n in range(8, 51):
                r = (_fakultaet_mod(n, m) + 1) % m
                if r in nicht_qr:
                    ausschluss_fuer.append(n)

            ergebnis['moduli_analyse'][m] = {
                'quadratische_reste': sorted(qr),
                'nicht_qr_count': len(nicht_qr),
                'ausschluss_fuer_n_8_50': ausschluss_fuer,
            }

            if ausschluss_fuer:
                ergebnis['nuetzliche_moduli'].append({
                    'modulus': m,
                    'n_ausgeschlossen': ausschluss_fuer,
                    'anzahl': len(ausschluss_fuer),
                })

        return ergebnis

    def finde_ausschlusskandidaten(self, n: int) -> Dict:
        """
        @brief Findet möglichst viele Moduli, die n!+1 als Nicht-Quadrat ausweisen.
        @description
            Strategie: Suche Moduli m mit m > n (sonst gilt n! ≡ 0 mod m für Primes p|m,p≤n).

            Für Primzahl p mit p > n:
            - p ∤ n!, also n! mod p ≠ 0 (i.A.)
            - n! + 1 mod p hängt von n! mod p ab
            - Falls n!+1 mod p ein Nicht-Rest: Ausschluss!

            Für p² mit Primzahl p ≤ n:
            - p | n! aber p² ∤ n! (wenn p ≤ n < 2p oder p² > n)
            - n! mod p² kann variieren

            Zusätzlich: Suche in {2,...,200} ∪ Primzahlen bis 500.

        @param n Kandidat für die Brocard-Gleichung.
        @return Dictionary mit gefundenen Ausschluss-Moduli.
        @lastModified 2026-03-12
        """
        ausschluss_moduli = []
        kein_ausschluss = []

        # Sammle Kandidaten-Moduli
        kandidaten = list(range(2, 50))
        # Primzahlen größer als n (besonders interessant)
        for p in sympy.primerange(n + 1, n + 200):
            kandidaten.append(int(p))
        # Primzahlquadrate p² mit p ≤ n (aber p² > n)
        for p in sympy.primerange(2, n + 1):
            p2 = int(p) ** 2
            if p2 > n:
                kandidaten.append(p2)

        # Entferne Duplikate und sortiere
        kandidaten = sorted(set(kandidaten))

        qr_cache: Dict[int, Set[int]] = {}

        for m in kandidaten:
            if m not in qr_cache:
                qr_cache[m] = _quadratische_reste(m)
            qr = qr_cache[m]

            fak_mod = _fakultaet_mod(n, m)
            rest = (fak_mod + 1) % m

            if rest not in qr:
                ausschluss_moduli.append({
                    'modulus': m,
                    'rest': rest,
                    'n_fak_mod': fak_mod,
                    'quadratische_reste': sorted(qr),
                })
            else:
                kein_ausschluss.append(m)

        return {
            'n': n,
            'ausschluss_moduli': ausschluss_moduli,
            'kein_ausschluss': kein_ausschluss,
            'ausgeschlossen': len(ausschluss_moduli) > 0,
            'stärkstes_argument': ausschluss_moduli[0] if ausschluss_moduli else None,
        }

    # -------------------------------------------------------
    # Numerische Suche
    # -------------------------------------------------------

    def numerische_suche(self, n_max: int, verbose: bool = False) -> Dict:
        """
        @brief Sucht n!+1 = m² für alle n von 1 bis n_max.
        @description
            **Algorithmus**:
            Für jedes n:
            1. Berechne fak = n! (exakte Ganzzahl, mit math.factorial)
            2. Berechne val = fak + 1
            3. Prüfe ob val ein Quadrat ist:
               - Für n ≤ 100: direkt math.isqrt(val)² == val
               - Für n > 100: Gleiche Methode, aber Zahlen werden sehr groß

            **Performance**:
            n! wächst super-exponentiell. Für n = 1000 hat n! etwa 2568 Stellen.
            math.isqrt ist in Python für arbiträre Präzision implementiert.

            Laut Berndt-Galway (2000) wurde bis n = 10^9 verifiziert.
            Dieses Modul verifiziert praktisch bis n_max (typisch 1000–10000).

        @param n_max Obere Schranke für n.
        @param verbose Wenn True, zeige Fortschritt alle 100 Schritte.
        @return Dictionary mit gefundenen Lösungen und Laufzeit.
        @lastModified 2026-03-12
        """
        start_zeit = time.time()
        loesungen = []
        laufzeit_pro_100 = []

        fak = 1  # Initialisierung: 0! = 1
        t100 = time.time()

        for n in range(1, n_max + 1):
            # n! = (n-1)! * n — inkrementell berechnen
            fak *= n
            val = fak + 1

            # Exakte Quadratwurzelprüfung
            r = math.isqrt(val)
            if r * r == val:
                loesungen.append({'n': n, 'm': r, 'n_fak_plus1': val})
                if verbose:
                    print(f"  LÖSUNG GEFUNDEN: n={n}, m={r}")

            if verbose and n % 100 == 0:
                jetzt = time.time()
                laufzeit_pro_100.append({'n': n, 'sekunden': round(jetzt - t100, 4)})
                t100 = jetzt
                stellen = len(str(fak))
                print(f"  n={n}: n! hat {stellen} Stellen, {jetzt - start_zeit:.2f}s")

        gesamt_zeit = time.time() - start_zeit

        # Verifikation gegen bekannte Lösungen
        gefundene_n = {s['n'] for s in loesungen}
        bekannte_n = set(self.BEKANNTE_LOESUNGEN.keys())
        vermisst = bekannte_n - gefundene_n

        if n_max >= max(bekannte_n):
            korrekt = (gefundene_n == bekannte_n & gefundene_n)
        else:
            # Nur Lösungen bis n_max vergleichen
            erwartete_im_bereich = {n for n in bekannte_n if n <= n_max}
            korrekt = erwartete_im_bereich.issubset(gefundene_n)

        return {
            'n_max': n_max,
            'loesungen': loesungen,
            'anzahl_loesungen': len(loesungen),
            'bekannte_loesungen': self.BEKANNTE_LOESUNGEN,
            'alle_bekannten_gefunden': not bool(vermisst & set(range(1, n_max + 1))),
            'unerwartete_loesungen': [s for s in loesungen if s['n'] not in bekannte_n],
            'laufzeit_sekunden': round(gesamt_zeit, 4),
            'laufzeit_pro_100': laufzeit_pro_100,
            'vermutung_status': (
                'BESTÄTIGT (im geprüften Bereich)'
                if not [s for s in loesungen if s['n'] not in bekannte_n]
                else 'GEGENBEISPIEL GEFUNDEN!'
            ),
        }

    # -------------------------------------------------------
    # Theoretische Schranken
    # -------------------------------------------------------

    def schranken_analyse(self) -> Dict:
        """
        @brief Untere und obere Schranken für m in Abhängigkeit von n.
        @description
            Falls n! + 1 = m²:

            **Untere Schranke** für m:
            m = √(n!+1) > √(n!) = √n! ≈ √(n/e)^n · √(2πn) (Stirling)

            Für große n: m ≈ √(n!) ≈ (n/e)^(n/2) · (2πn)^(1/4)

            **Obere Schranke**:
            m = √(n!+1) < √(n!+n!) = √(2n!) für n! ≥ 1.
            Präziser: m = √(n!+1) < √(n!)+1 (binomische Näherung für große n!).

            **Wachstumsrate**:
            log(m) ≈ (1/2) · log(n!) ≈ (n/2)·log(n) - n/2 (Stirling)

            **Folgerung für Ausschluss**:
            Für jede natürliche Zahl m muss m² ≡ n!+1 exakt gelten.
            Da m ≈ √(n!), hat m selbst schon n/2·log₁₀(n) Stellen für große n.

        @return Dictionary mit Schrankenformeln und numerischen Beispielen.
        @lastModified 2026-03-12
        """
        beispiele = []
        for n in [4, 5, 7, 10, 20, 50, 100]:
            fak = math.factorial(n)
            m_exakt = math.isqrt(fak + 1)
            # Stirling-Näherung für log(n!)
            if n > 0:
                log_fak_stirling = n * math.log(n) - n + 0.5 * math.log(2 * math.pi * n)
                m_schranke_unten = math.exp(log_fak_stirling / 2)
            else:
                m_schranke_unten = 1.0

            stellen_m = len(str(m_exakt)) if m_exakt > 0 else 1
            stellen_fak = len(str(fak))

            beispiele.append({
                'n': n,
                'n_fak_stellen': stellen_fak,
                'm_exakte_wurzel': m_exakt,
                'm_stellen': stellen_m,
                'ist_loesung': m_exakt * m_exakt == fak + 1,
                'm_stirling_schranke': round(m_schranke_unten, 2),
            })

        return {
            'formel_untere_schranke': 'Stirling: m > (n/e)^(n/2) · (2πn)^(1/4)',
            'formel_obere_schranke': 'm = √(n!+1) < √(n!) + 1 (für n! >> 1)',
            'wachstumsrate': 'log(m) ≈ (n/2)·log(n) − n/2',
            'folgerung': 'm hat für n=100 ca. 80 Stellen, für n=1000 ca. 1285 Stellen',
            'beispiele': beispiele,
            'bekannte_loesungen': self.BEKANNTE_LOESUNGEN,
        }

    # -------------------------------------------------------
    # p-adische Analyse
    # -------------------------------------------------------

    def p_adische_analyse(self, p: int) -> Dict:
        """
        @brief p-adische Bewertung von n!+1 und ihre Implikationen.
        @description
            Die p-adische Bewertung v_p(k) gibt die höchste Potenz von p an,
            die k teilt: v_p(k) = max{j ∈ ℕ₀ | p^j | k}.

            **Legendres Formel** für v_p(n!):
                v_p(n!) = Σ_{k≥1} ⌊n/p^k⌋

            **Für n! + 1 = m²** müsste gelten:
                v_p(n! + 1) muss gerade sein (Bedingung für Quadrat).

            **Analyse von v_p(n!+1)**:
            - Falls p ∤ n! (d.h. p > n): v_p(n!+1) = v_p(n!+1) direkt aus n!+1.
            - Falls p | n! (d.h. p ≤ n): n! ≡ 0 (mod p), also n!+1 ≡ 1 (mod p).
              Damit v_p(n!+1) = 0 (da p ∤ n!+1), was gerade ist.
              → Kein Ausschluss durch p | n! allein!

            **Für p² | n! mit p ≤ n < p²**:
            Dann v_p(n!) = 1 (genau einmal p als Faktor).
            n! = p·k mit p ∤ k.
            n!+1 = p·k+1 ≡ 1 (mod p). Also v_p(n!+1) = 0.

            **Legendre-Analyse**:
            v_p(n!) = ⌊n/p⌋ + ⌊n/p²⌋ + ⌊n/p³⌋ + ...

            Das Ergebnis ist gerade ⟺ n! ist ein p-adisches Quadrat.
            n!+1 ist schwerer zu analysieren.

        @param p Primzahl für die Analyse.
        @return Dictionary mit p-adischer Analyse.
        @lastModified 2026-03-12
        """
        if not isprime(p):
            return {'fehler': f'{p} ist keine Primzahl'}

        beispiele = []
        for n in list(range(1, 30)) + [50, 100]:
            # Legendres Formel für v_p(n!)
            vp_fak = 0
            pk = p
            while pk <= n:
                vp_fak += n // pk
                pk *= p

            # n! + 1 mod p^(vp_fak+2)
            fak = math.factorial(n)
            val = fak + 1

            # Exakte p-adische Bewertung von n!+1
            vp_val = 0
            temp = val
            while temp % p == 0:
                vp_val += 1
                temp //= p

            ist_loesung = math.isqrt(val) ** 2 == val

            beispiele.append({
                'n': n,
                'v_p(n!)': vp_fak,
                'v_p(n!+1)': vp_val,
                'v_p_gerade': vp_val % 2 == 0,
                'ausschluss': vp_val % 2 == 1,  # ungerade v_p → kein Quadrat
                'ist_loesung': ist_loesung,
            })

        # Zähle Ausschlüsse
        ausschluesse = [e for e in beispiele if e['ausschluss']]

        # Bekannte Lösungen prüfen
        for n_loes, m_loes in self.BEKANNTE_LOESUNGEN.items():
            eintrag = next((e for e in beispiele if e['n'] == n_loes), None)
            if eintrag:
                # Konsistenzcheck: bekannte Lösung darf nicht ausgeschlossen sein
                assert not eintrag['ausschluss'], (
                    f"FEHLER: Bekannte Lösung n={n_loes} fälschlich ausgeschlossen durch p={p}"
                )

        return {
            'p': p,
            'legendre_formel': f'v_{p}(n!) = Σ ⌊n/{p}^k⌋',
            'bedingung_fuer_quadrat': f'v_{p}(n!+1) muss gerade sein',
            'erklaerung': (
                f'Für n ≥ {p}: {p} | n!, also n!+1 ≡ 1 (mod {p}), '
                f'also v_{p}(n!+1) = 0 (gerade). Kein Ausschluss durch {p} allein.'
            ),
            'beispiele': beispiele,
            'ausschluss_faelle': ausschluesse,
            'fazit': (
                f'p={p} liefert {len(ausschluesse)} Ausschlüsse für n ∈ [1,100].'
            ),
        }

    # -------------------------------------------------------
    # Vollständige Analyse-Hilfsmethode
    # -------------------------------------------------------

    def vollstaendige_modular_analyse(
        self, n_min: int = 8, n_max: int = 100
    ) -> Dict:
        """
        @brief Vollständige modulare Ausschlussanalyse für n_min bis n_max.
        @description
            Für jedes n im Bereich wird versucht, es durch modulare Argumente
            als Lösung auszuschließen. Verwendet eine feste Liste von Moduli.
        @param n_min Untere Grenze (Standard: 8, da 4,5,7 bekannte Lösungen).
        @param n_max Obere Grenze.
        @return Zusammenfassung der Ausschlüsse.
        @lastModified 2026-03-12
        """
        # Standardmoduli-Liste (bewusst verschiedene Typen)
        standard_moduli = [
            5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
            # Quadrate von Primzahlen
            4, 9, 25, 49, 121,
            # Zusammengesetzte Zahlen
            8, 12, 16, 24, 40,
        ]

        gesamt = {
            'n_min': n_min,
            'n_max': n_max,
            'ausgeschlossen': [],
            'nicht_ausgeschlossen': [],
            'anteil_ausgeschlossen': 0.0,
        }

        for n in range(n_min, n_max + 1):
            res = self.modular_ausschluss(n, standard_moduli)
            if res['ausgeschlossen']:
                gesamt['ausgeschlossen'].append({
                    'n': n,
                    'durch_modulus': res['ausschluss_modulus'],
                    'rest': res['ausschluss_rest'],
                })
            else:
                gesamt['nicht_ausgeschlossen'].append(n)

        total = n_max - n_min + 1
        gesamt['anteil_ausgeschlossen'] = (
            round(len(gesamt['ausgeschlossen']) / total * 100, 2) if total > 0 else 0.0
        )
        return gesamt


# ===========================================================
# Demonstration
# ===========================================================

def _demo():
    """
    @brief Demonstriert die Klasse BrocardExtension.
    @lastModified 2026-03-12
    """
    ext = BrocardExtension()

    print("=" * 60)
    print("BROCARD-RAMANUJAN-VERMUTUNG: Erweiterte Analyse")
    print("=" * 60)

    # Numerische Suche bis n = 1000
    print("\n--- Numerische Suche bis n = 1000 ---")
    start = time.time()
    res = ext.numerische_suche(1000, verbose=False)
    dauer = time.time() - start
    print(f"  Gefundene Lösungen: {[(s['n'], s['m']) for s in res['loesungen']]}")
    print(f"  Vermutungsstatus: {res['vermutung_status']}")
    print(f"  Laufzeit: {dauer:.3f}s")

    # Modulare Ausschlussanalyse für einzelne n
    print("\n--- Modulare Ausschlussanalyse ---")
    for n in [8, 9, 10, 15, 20]:
        kandidaten = ext.finde_ausschlusskandidaten(n)
        status = "AUSGESCHLOSSEN" if kandidaten['ausgeschlossen'] else "nicht ausgeschlossen"
        stärkstes = kandidaten.get('stärkstes_argument')
        mod_info = f" durch m={stärkstes['modulus']}" if stärkstes else ""
        print(f"  n={n:3d}: {status}{mod_info}")

    # Schrankenanalyse
    print("\n--- Schrankenanalyse ---")
    schranken = ext.schranken_analyse()
    for b in schranken['beispiele'][:6]:
        loes_marker = " ← LÖSUNG" if b['ist_loesung'] else ""
        print(f"  n={b['n']:4d}: n! hat {b['n_fak_stellen']:4d} Stellen, "
              f"m≈{b['m_exakte_wurzel']:>15}{loes_marker}")

    # p-adische Analyse für p=5
    print("\n--- p-adische Analyse (p=5) ---")
    pad = ext.p_adische_analyse(5)
    print(f"  {pad['erklaerung']}")
    ausschluss_faelle = [e['n'] for e in pad['ausschluss_faelle'][:10]]
    print(f"  Ausschlussfälle (erste 10): {ausschluss_faelle}")

    # Restklassenanalyse mod 8
    print("\n--- Restklassenanalyse mod 8 ---")
    analyse_mod8 = ext.modular_ausschluss(10, [8])
    r8 = analyse_mod8['moduli'][8]
    print(f"  n=10: 10!+1 ≡ {r8['n_fak_plus1_mod']} (mod 8)")
    print(f"  QR(8) = {r8['quadratische_reste']}")
    print(f"  Ausschluss: {r8['ausschluss']}")

    # Vollständige modulare Analyse
    print("\n--- Vollständige modulare Analyse n ∈ [8, 50] ---")
    voll = ext.vollstaendige_modular_analyse(8, 50)
    print(f"  Ausgeschlossen: {len(voll['ausgeschlossen'])}/{50-8+1} ({voll['anteil_ausgeschlossen']}%)")
    nicht = voll['nicht_ausgeschlossen']
    print(f"  Nicht ausgeschlossen: {nicht[:20]}")


if __name__ == "__main__":
    _demo()
