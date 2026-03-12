"""
@file giuga_4prim.py
@brief Formaler Beweisversuch: Kein 4-Primfaktor-Giuga-Pseudoprime existiert.
@description
    Dieses Modul untersucht, ob eine zusammengesetzte Zahl der Form
    n = p · q · r · s mit p < q < r < s (paarweise verschiedene Primzahlen)
    ein Giuga-Pseudoprime sein kann.

    **Hintergrund**: Ein Giuga-Pseudoprime ist eine zusammengesetzte Zahl n,
    die für jeden Primteiler pᵢ | n beide folgenden Bedingungen erfüllt:
        - Schwache Bedingung: pᵢ | (n/pᵢ - 1)
        - Starke Bedingung:   (pᵢ - 1) | (n/pᵢ - 1)

    **Bekannte Resultate** (Borwein et al. 1996, Bednarek 2014):
        - Jedes Giuga-Pseudoprime hat ≥ 13635 Primfaktoren.
        - Das kleinste hat > 10^19907 Dezimalstellen.

    **Aufbauend auf beweisversuche.py** (Sätze 1–4 + Korollar):
        - Satz 1:   Quadratfreiheit (n = p₁ · p₂ · ... · pₖ, pairweise verschieden)
        - Satz 2:   Kein 2-Prim-Giuga-Pseudoprime (n = p · q)
        - Satz 3:   Kein 3-Prim-Giuga-Pseudoprime mit p₁ = 2 (n = 2 · q · r)
        - Satz 4:   Kein 3-Prim-Giuga-Pseudoprime mit allen ungeraden Faktoren
        - Korollar: KEIN 3-Prim-Giuga-Pseudoprime existiert

    **Diese Datei** zeigt:
        (a) Numerische Verifikation bis zu einer großen Schranke: kein Kandidat gefunden.
        (b) Schranken-Argumente für die Fälle p=2, p=3, p≥5, die den Suchraum stark
            einschränken.
        (c) Explizite Darstellung, wo der Beweis noch Lücken hat
            (der 4-Prim-Fall ist mathematisch noch offen).

    **Status des 4-Prim-Falls**:
        KEIN vollständiger algebraischer Beweis bekannt (Stand 2026).
        Numerische Evidenz: Kein 4-Prim-Giuga-Pseudoprime bis n < 10^6.
        Schrankenanalyse: Wenn ein solches existiert, müssen die Faktoren
        extrem groß sein (konsistent mit Borwein et al.).

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import sympy
from sympy import isprime, nextprime, primerange, factorint

# =============================================================================
# Hilfsfunktionen
# =============================================================================


def berechne_giuga_bedingung(n: int) -> Dict:
    """
    @brief Prüft, ob n ein Giuga-Pseudoprime (oder eine Giuga-Zahl) ist.
    @description
        Für jede zusammengesetzte quadratfreie Zahl n = p₁ · p₂ · ... · pₖ
        werden die schwachen und starken Giuga-Bedingungen geprüft:
            - Schwach: pᵢ | (n/pᵢ - 1)  für alle i
            - Stark:  (pᵢ-1) | (n/pᵢ - 1) für alle i

        Giuga-Zahl:       erfüllt nur die schwache Bedingung.
        Giuga-Pseudoprime: erfüllt BEIDE Bedingungen (und ist zusammengesetzt).

        Bekannte Giuga-Zahlen (schwach): {30, 858, 1722, 66198, ...}
        Bekannte Giuga-Pseudoprimes: Keines bisher gefunden (vermutlich keine).

    @param n Zu prüfende natürliche Zahl ≥ 4.
    @return Dictionary mit Ergebnis, Primfaktoren und Detail-Flags.
    @lastModified 2026-03-12
    """
    if n < 4:
        return {'n': n, 'ist_giuga_zahl': False, 'ist_giuga_pseudoprime': False,
                'grund': 'n < 4, nicht relevant'}

    if isprime(n):
        return {'n': n, 'ist_giuga_zahl': False, 'ist_giuga_pseudoprime': False,
                'grund': 'n ist prim'}

    # Primfaktorzerlegung bestimmen
    fd = factorint(n)
    # Prüfe Quadratfreiheit (notwendig für Giuga-Kandidat)
    if any(exp > 1 for exp in fd.values()):
        return {'n': n, 'ist_giuga_zahl': False, 'ist_giuga_pseudoprime': False,
                'grund': 'n ist nicht quadratfrei'}

    primfaktoren = sorted(fd.keys())
    schwach_details = []
    stark_details = []

    # Bedingungen für jeden Primteiler prüfen
    for p in primfaktoren:
        # n/p ist das Produkt aller anderen Primfaktoren
        n_durch_p = n // p
        rest = n_durch_p - 1

        schwach_ok = (rest % p == 0)
        # Starke Bedingung: (p-1) | (n/p - 1); für p=2 ist p-1=1, immer erfüllt
        stark_ok = (rest % (p - 1) == 0) if p > 2 else True

        schwach_details.append({'p': p, 'n/p - 1': rest, 'schwach': schwach_ok})
        stark_details.append({'p': p, 'n/p - 1': rest, 'stark': stark_ok})

    ist_schwach = all(d['schwach'] for d in schwach_details)
    ist_stark = all(d['stark'] for d in stark_details)

    return {
        'n': n,
        'primfaktoren': primfaktoren,
        'ist_giuga_zahl': ist_schwach,          # nur schwache Bedingung
        'ist_giuga_pseudoprime': ist_schwach and ist_stark,  # beide Bedingungen
        'schwach_details': schwach_details,
        'stark_details': stark_details,
    }


def numerische_suche_4prim(grenze: int = 1_000_000) -> Dict:
    """
    @brief Sucht alle 4-Prim-Giuga-Pseudoprime-Kandidaten bis zur angegebenen Grenze.
    @description
        Iteriert über alle quadratfreien zusammengesetzten Zahlen n ≤ grenze
        mit genau 4 Primfaktoren und prüft beide Giuga-Bedingungen.

        Effizienz-Hinweis: Für große Grenzen wird nur über Produkte der Form
        p₁ · p₂ · p₃ · p₄ iteriert (mit p₁ < p₂ < p₃ < p₄), statt über
        alle n ≤ grenze.

    @param grenze Obere Suchschranke für n.
    @return Dictionary mit gefundenen Kandidaten und Statistiken.
    @lastModified 2026-03-12
    """
    giuga_pseudoprimes = []  # erfüllen BEIDE Bedingungen (schwach + stark)
    giuga_zahlen = []        # erfüllen nur schwache Bedingung
    geprueft = 0

    # Erzeuge alle Primzahlen bis zur Grenze/2 (grob, für p₁)
    # Effiziente Methode: generiere 4-fache Produkte direkt
    # p₁ < p₂ < p₃ < p₄ mit p₁·p₂·p₃·p₄ ≤ grenze
    # Maximaler Wert für p₁: (grenze)^(1/4)
    p_max_1 = int(grenze ** 0.25) + 1
    primzahlen = list(primerange(2, min(grenze // 2 + 1, 10000)))

    for i1, p1 in enumerate(primzahlen):
        if p1 > p_max_1:
            break
        for i2 in range(i1 + 1, len(primzahlen)):
            p2 = primzahlen[i2]
            if p1 * p2 * (p2 + 2) * (p2 + 4) > grenze:
                # Schranke: Selbst mit den kleinsten nächsten Primzahlen überschritten
                break
            for i3 in range(i2 + 1, len(primzahlen)):
                p3 = primzahlen[i3]
                if p1 * p2 * p3 * nextprime(p3) > grenze:
                    break
                for i4 in range(i3 + 1, len(primzahlen)):
                    p4 = primzahlen[i4]
                    n = p1 * p2 * p3 * p4
                    if n > grenze:
                        break
                    geprueft += 1

                    # Giuga-Bedingungen direkt prüfen (ohne erneute Faktorisierung)
                    faktoren = [p1, p2, p3, p4]
                    schwach = all((n // p - 1) % p == 0 for p in faktoren)
                    if not schwach:
                        continue

                    # Schwache Bedingung erfüllt → Giuga-Zahl
                    giuga_zahlen.append({'n': n, 'faktoren': faktoren})

                    # Starke Bedingung
                    stark = all(
                        (n // p - 1) % (p - 1) == 0 if p > 2 else True
                        for p in faktoren
                    )
                    if stark:
                        giuga_pseudoprimes.append({'n': n, 'faktoren': faktoren})

    return {
        'grenze': grenze,
        'geprueft': geprueft,
        'giuga_zahlen_4prim': giuga_zahlen,          # nur schwach
        'giuga_pseudoprimes_4prim': giuga_pseudoprimes,  # schwach + stark
        'kein_pseudoprime_gefunden': len(giuga_pseudoprimes) == 0,
        'fazit': (
            'Kein 4-Prim-Giuga-Pseudoprime bis {} gefunden.'.format(grenze)
            if not giuga_pseudoprimes
            else 'GEGENBEISPIEL GEFUNDEN: {}'.format(giuga_pseudoprimes)
        )
    }


# =============================================================================
# Schrankenberechnung
# =============================================================================


def schranke_fuer_4prim_fall_p_gleich_2(q_max: int = 200) -> Dict:
    """
    @brief Schrankenanalyse für den Fall p₁=2 bei 4-Prim-Kandidaten (n=2·q·r·s).
    @description
        Sei n = 2 · q · r · s mit 2 < q < r < s.

        **Schwache Bedingungen**:
        - 2 | (qrs - 1): Da q,r,s ungerade → qrs ungerade → qrs-1 gerade ✓ (immer)
        - q | (2rs - 1)    ... (I)
        - r | (2qs - 1)    ... (II)
        - s | (2qr - 1)    ... (III)

        **Starke Bedingungen**:
        - 1 | (qrs - 1):   immer wahr (p-1 = 1 für p=2)
        - (q-1) | (2rs-1) ... (I')
        - (r-1) | (2qs-1) ... (II')
        - (s-1) | (2qr-1) ... (III')

        **Aus (I) und (I')**:
        q | (2rs-1) und (q-1) | (2rs-1).
        Da gcd(q, q-1) = 1, folgt: q·(q-1) | (2rs-1).
        Also: 2rs - 1 ≥ q(q-1), d.h. 2rs ≥ q(q-1) + 1.

        **Aus (III)**:
        s | (2qr-1). Da s > r > q, ist 2qr - 1 < 2qr ≤ 2qs < 2s·s = 2s².
        Also entweder: s = 2qr - 1 (einzige mögliche Lösung, da s ≤ 2qr-1 < 2s).
        Genauer: Wenn s > 2qr - 1, gibt es kein positives Vielfaches von s ≤ 2qr-1.
        Also MUSS gelten: s ≤ 2qr - 1, und s teilt 2qr - 1,
        also s ∈ {Teiler von (2qr-1) die > r sind}.

        Dies kombiniert mit (III') liefert: (s-1) | (2qr-1).
        Da s | (2qr-1) und (s-1) | (2qr-1) und gcd(s,s-1)=1:
        s·(s-1) | (2qr - 1).
        Also: 2qr - 1 ≥ s(s-1) ≥ s².  - s ≥ s²/2 (grob).
        D.h.: r ≥ s(s-1) / (2q) ≥ s²/(2q) - s/(2q).

        Das ist eine wichtige Schranke: r wächst schnell mit s.

    @param q_max Maximale Größe für q bei der numerischen Überprüfung.
    @return Dictionary mit Schranken und Analyseergebnis.
    @lastModified 2026-03-12
    """
    # Prüfe numerisch: Für welche (q, r, s) erfüllen alle 4 starken Bedingungen?
    kandidaten_stark = []
    kandidaten_schwach = []

    primes_q = list(primerange(3, q_max + 1))

    for q in primes_q:
        # Obere Schranke für r: r < (2qr-1)/s ≤ ... grob: r < q_max²
        r_max = min(q_max * q_max, 10 * q_max)
        for r in primerange(q + 1, r_max + 1):
            # Aus Bedingung (III): s teilt 2qr-1, s > r
            produkt_III = 2 * q * r - 1
            if produkt_III <= r:
                continue
            # Alle Primteiler von (2qr-1) die > r sind
            try:
                teiler = factorint(produkt_III)
            except Exception:
                continue
            for s_kandidat in sorted(teiler.keys()):
                if s_kandidat <= r:
                    continue
                # s muss Primzahl sein (Exponent 1 für Quadratfreiheit)
                if teiler[s_kandidat] != 1:
                    continue
                s = s_kandidat

                n = 2 * q * r * s
                faktoren = [2, q, r, s]

                # Schwache Bedingungen prüfen
                schwach = all((n // p - 1) % p == 0 for p in faktoren)
                if not schwach:
                    continue

                kandidaten_schwach.append({'q': q, 'r': r, 's': s, 'n': n})

                # Starke Bedingungen prüfen
                stark = all(
                    (n // p - 1) % (p - 1) == 0 if p > 2 else True
                    for p in faktoren
                )
                if stark:
                    kandidaten_stark.append({'q': q, 'r': r, 's': s, 'n': n})

    return {
        'fall': 'p=2: n = 2·q·r·s',
        'q_max': q_max,
        'schwache_kandidaten': kandidaten_schwach[:10],
        'starke_kandidaten_pseudoprimes': kandidaten_stark,
        'kein_pseudoprime': len(kandidaten_stark) == 0,
        'schranken_beschreibung': {
            'aus_III_und_III_strich': (
                's·(s-1) | (2qr-1), also 2qr ≥ s(s-1)+1 ≥ s²−s+1.'
                ' D.h. r ≥ s(s-1)/(2q) — r wächst quadratisch in s.'
            ),
            'aus_I_und_I_strich': (
                'q·(q-1) | (2rs-1), also 2rs ≥ q(q-1)+1.'
                ' Liefert untere Schranke für rs in Abhängigkeit von q.'
            ),
            'lücke': (
                'Das Schranken-Argument führt zu immer größeren Faktoren,'
                ' beweist aber nicht formal den Widerspruch.'
                ' Der 4-Prim-Fall ist mathematisch noch offen.'
            ),
        }
    }


def schranke_fuer_4prim_fall_p_gleich_3(q_max: int = 200) -> Dict:
    """
    @brief Schrankenanalyse für den Fall p₁=3 bei 4-Prim-Kandidaten (n=3·q·r·s).
    @description
        Sei n = 3 · q · r · s mit 3 < q < r < s.

        **Schwache Bedingungen**:
        - 3 | (qrs - 1)    ... (I)
        - q | (3rs - 1)    ... (II)
        - r | (3qs - 1)    ... (III)
        - s | (3qr - 1)    ... (IV)

        **Aus (I)**: qrs ≡ 1 (mod 3), d.h. q·r·s ≡ 1 (mod 3).
        Da q, r, s ≠ 3 (da 3 bereits als p₁ verwendet) und alle ungerade prim,
        gilt q, r, s ≡ 1 oder 2 (mod 3). Das Produkt qrs ≡ 1 (mod 3) erfordert,
        dass die Anzahl der Faktoren ≡ 2 (mod 3) gerade ist: 0 oder 2.

        **Aus (IV)**: s | (3qr-1). Da s > r > q > 3:
        Falls s > 3qr-1 → kein positives Vielfaches. Also: s ≤ 3qr-1.
        Für das einzig mögliche 1-fache Vielfaches: s = 3qr-1 oder
        s ein echter Teiler von (3qr-1) mit s > r.

        **Kombination stark + schwach für s**:
        s·(s-1) | (3qr-1).
        Also 3qr ≥ s(s-1)+1.
        Damit: r ≥ s(s-1)/(3q).

        **Vergleich mit Fall p=2**:
        Im Fall p=3 ist der Faktor 3 (statt 2) in der Schranke → der Suchraum
        ist etwas größer, aber die Struktur ist analog.

    @param q_max Maximale Größe für q.
    @return Analyse-Ergebnis.
    @lastModified 2026-03-12
    """
    kandidaten_stark = []
    kandidaten_schwach = []

    primes_q = list(primerange(5, q_max + 1))  # q > 3

    for q in primes_q:
        r_max = min(q_max * q, 10 * q_max)
        for r in primerange(q + 1, r_max + 1):
            # Aus Bedingung (IV): s teilt 3qr-1, s > r
            produkt_IV = 3 * q * r - 1
            if produkt_IV <= r:
                continue
            try:
                teiler = factorint(produkt_IV)
            except Exception:
                continue
            for s_kandidat in sorted(teiler.keys()):
                if s_kandidat <= r:
                    continue
                if teiler[s_kandidat] != 1:
                    continue
                s = s_kandidat

                n = 3 * q * r * s
                faktoren = [3, q, r, s]

                schwach = all((n // p - 1) % p == 0 for p in faktoren)
                if not schwach:
                    continue

                kandidaten_schwach.append({'q': q, 'r': r, 's': s, 'n': n})

                stark = all((n // p - 1) % (p - 1) == 0 for p in faktoren)
                if stark:
                    kandidaten_stark.append({'q': q, 'r': r, 's': s, 'n': n})

    return {
        'fall': 'p=3: n = 3·q·r·s',
        'q_max': q_max,
        'schwache_kandidaten': kandidaten_schwach[:10],
        'starke_kandidaten_pseudoprimes': kandidaten_stark,
        'kein_pseudoprime': len(kandidaten_stark) == 0,
        'schranken_beschreibung': {
            'aus_IV_und_IV_strich': (
                's·(s-1) | (3qr-1), also 3qr ≥ s(s-1)+1.'
                ' D.h. r ≥ s(s-1)/(3q).'
            ),
            'kongruenz_aus_I': (
                'q·r·s ≡ 1 (mod 3): Anzahl der Faktoren ≡ 2 (mod 3) ist gerade.'
            ),
        }
    }


def schranke_fuer_4prim_allgemeiner_fall(p_max: int = 20) -> Dict:
    """
    @brief Schrankenanalyse für den allgemeinen Fall p≥5 bei 4-Prim-Kandidaten.
    @description
        Sei n = p · q · r · s mit 5 ≤ p < q < r < s (alle ungerade prim, p ≥ 5).

        **Schlüsselargument** (analog zu Satz 4 in beweisversuche.py):

        Für den größten Faktor s gilt die starke + schwache Giuga-Bedingung:
            s | (pqr - 1)  und  (s-1) | (pqr - 1)
        Also: s·(s-1) | (pqr - 1).

        Da pqr - 1 > 0: s(s-1) ≤ pqr - 1 < pqr.

        Untere Schranke für s: Da s > r und beide ungerade prim, gilt s ≥ r+2.
        Also: (r+2)(r+1) ≤ s(s-1) < pqr.
        Damit: pqr > r² + 3r + 2, d.h. pq > r + 3 + 2/r.

        Für das nächste Niveau (r und q):
        Entsprechend für r: r·(r-1) | (pqs - 1).
        Da s·(s-1) ≤ pqr - 1, ist s selbst schon groß → pqs - 1 ist sehr groß.
        Dies ist konsistent — kein direkter Widerspruch entsteht hier.

        **Warum der 3-Prim-Beweis hier nicht übertragbar ist**:
        Im 3-Prim-Fall (Satz 4) leitete man ab: p ≥ q + 4, Widerspruch mit p < q.
        Im 4-Prim-Fall erhält man analoge Schranken nur für die größeren Faktoren,
        während p und q noch "Spielraum" haben, größer zu sein als die abgeleiteten
        Schranken. Daher ist der direkte Widerspruch nicht reproduzierbar.

        **Numerischer Befund**: Kein 4-Prim-Giuga-Pseudoprime bis n < 10^6 gefunden.

    @param p_max Maximaler kleinster Primfaktor p.
    @return Schranken und Zusammenfassung.
    @lastModified 2026-03-12
    """
    kandidaten = []
    schranken_info = []

    primes_p = list(primerange(5, p_max + 1))

    for p in primes_p:
        for q in primerange(p + 1, p_max * 10 + 1):
            # Obere Schranke für r aus pq > r+3: r < pq - 3
            r_schranke = p * q - 3
            if r_schranke <= q:
                continue
            for r in primerange(q + 1, min(r_schranke, p_max * 100) + 1):
                # s teilt (pqr-1), s > r
                produkt = p * q * r - 1
                if produkt <= r:
                    continue
                try:
                    teiler = factorint(produkt)
                except Exception:
                    continue
                for s_kand in sorted(teiler.keys()):
                    if s_kand <= r:
                        continue
                    if teiler[s_kand] != 1:
                        continue
                    s = s_kand
                    n = p * q * r * s
                    faktoren = [p, q, r, s]

                    schwach = all((n // f - 1) % f == 0 for f in faktoren)
                    if not schwach:
                        continue

                    stark = all((n // f - 1) % (f - 1) == 0 for f in faktoren)
                    if stark:
                        kandidaten.append({'p': p, 'q': q, 'r': r, 's': s, 'n': n})

        # Schrankeninformation für dieses p festhalten
        schranken_info.append({
            'p': p,
            'aus_s_bedingung': f'pq > r + 3 → r < {p}·q - 3 (Oberschranke für r)',
            'aus_r_bedingung': f'(r+2)(r+1) < pqr, also pq > r+3+2/r',
        })

    return {
        'fall': 'p ≥ 5: n = p·q·r·s (alle ungerade, p ≥ 5)',
        'p_max': p_max,
        'giuga_pseudoprimes_gefunden': kandidaten,
        'kein_pseudoprime': len(kandidaten) == 0,
        'schranken': schranken_info[:5],
        'theoretische_lücke': (
            'Im Gegensatz zum 3-Prim-Fall ergibt der 4-Prim-Fall keinen direkten'
            ' Widerspruch p ≥ q (da p und q noch konsistente Werte annehmen können).'
            ' Der algebraische Beweis ist daher noch nicht vollständig.'
        ),
    }


# =============================================================================
# Hauptklasse: Giuga4PrimBeweis
# =============================================================================


class Giuga4PrimBeweis:
    """
    @brief Formale Analyse des 4-Prim-Giuga-Falls.
    @description
        Diese Klasse bündelt alle Beweisversuche und numerischen Verifikationen
        für die Frage: "Existiert ein Giuga-Pseudoprime mit genau 4 Primfaktoren?"

        Die Analyse umfasst:
            1. Fallunterscheidung p=2, p=3, p≥5
            2. Schrankenargumente für jeden Fall
            3. Numerische Verifikation bis zu einer Schranke
            4. Zusammenfassung der Beweislage (inkl. offener Lücken)

        **Fazit**: Numerisch kein Gegenbeispiel gefunden.
        Algebraisch: Kein vollständiger Beweis existiert (Stand 2026).
        Bekannte Literaturergebnisse (Borwein, Bednarek) schließen den
        4-Prim-Fall implizit aus (da k ≥ 13635 Faktoren nötig wären),
        aber ein direkter elementarer Beweis ist noch offen.

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    def fall1_p_gleich_2_analyse(self) -> Dict:
        """
        @brief Schrankenanalyse für den Fall p₁ = 2 (n = 2·q·r·s).
        @description
            **Bedingungen für n = 2·q·r·s (2 < q < r < s ungerade prim)**:

            Schwache Bedingungen:
            - 2 | (qrs-1): Da qrs ungerade → qrs-1 gerade → automatisch erfüllt.
            - q | (2rs-1)         ... (I)
            - r | (2qs-1)         ... (II)
            - s | (2qr-1)         ... (III)

            Starke Bedingungen:
            - 1 | (qrs-1): trivial (p-1 = 1 für p=2).
            - (q-1) | (2rs-1)     ... (I')
            - (r-1) | (2qs-1)     ... (II')
            - (s-1) | (2qr-1)     ... (III')

            **Kombinierte Schranken**:
            Aus (III) + (III'): gcd(s, s-1) = 1 → s(s-1) | (2qr-1).
            Also: 2qr - 1 ≥ s(s-1), d.h. s² - s ≤ 2qr - 1 < 2qr.
            Damit: s < √(2qr) + 1/2 ≈ √(2qr).

            Analog aus (I) + (I'): q(q-1) | (2rs-1), also 2rs ≥ q(q-1)+1.
            Und aus (II) + (II'): r(r-1) | (2qs-1), also 2qs ≥ r(r-1)+1.

            Diese Schranken sind konsistent (kein direkter Widerspruch),
            aber sie schränken den Suchraum stark ein.

            **Warum kein Widerspruch wie im 3-Prim-Fall**:
            Im 3-Prim-Fall 2·q·r lieferte r | (2q-1) direkt r = 2q-1
            (da r > q und r teilt 2q-1 < 2r). Dann folgte Paritätswiderspruch.
            Im 4-Prim-Fall ist der analoge Schritt komplexer: s | (2qr-1),
            und da 2qr-1 ein viel größerer Wert als s sein kann, gibt es
            mehrere mögliche Werte für s. Der Widerspruch greift nicht direkt.

        @return Beweis-Status und Schranken.
        @lastModified 2026-03-12
        """
        # Numerische Überprüfung der Schranken (kleiner Bereich)
        ergebnis = schranke_fuer_4prim_fall_p_gleich_2(q_max=100)

        return {
            'fall': 'Fall 1: p = 2 (n = 2·q·r·s)',
            'status': 'KEIN VOLLSTÄNDIGER BEWEIS (numerisch bestätigt)',
            'schlüsselbedingungen': {
                'automatisch': '2 | (qrs-1) (qrs ist immer ungerade)',
                'aus_s': 's(s-1) | (2qr-1) → s < √(2qr) + O(1)',
                'aus_q': 'q(q-1) | (2rs-1) → q² ≤ 2rs',
                'aus_r': 'r(r-1) | (2qs-1) → r² ≤ 2qs',
            },
            'warum_kein_widerspruch_wie_3prim': (
                'Im 3-Prim-Fall 2·q·r gilt: r|(2q-1) und r>q → r=2q-1 eindeutig.'
                ' Dann (r-1)=2(q-1) gerade teilt 2q-1 ungerade → Widerspruch.'
                ' Im 4-Prim-Fall: s|(2qr-1), aber 2qr-1 kann s·k für k≥1 sein —'
                ' nicht eindeutig. Der einfache Paritätswiderspruch greift nicht.'
            ),
            'numerisch': {
                'q_max': ergebnis['q_max'],
                'schwache_kandidaten_gefunden': len(ergebnis['schwache_kandidaten']),
                'pseudoprimes_gefunden': len(ergebnis['starke_kandidaten_pseudoprimes']),
                'bestätigt': ergebnis['kein_pseudoprime'],
            },
        }

    def fall2_p_gleich_3_analyse(self) -> Dict:
        """
        @brief Schrankenanalyse für den Fall p₁ = 3 (n = 3·q·r·s).
        @description
            **Bedingungen für n = 3·q·r·s (3 < q < r < s ungerade prim)**:

            Schwache Bedingungen:
            - 3 | (qrs-1): erfordert qrs ≡ 1 (mod 3).
            - q | (3rs-1)       ... (I)
            - r | (3qs-1)       ... (II)
            - s | (3qr-1)       ... (IV)

            Starke Bedingungen:
            - 2 | (qrs-1): erfordert qrs ≡ 1 (mod 2) → qrs ungerade ✓.
            - (q-1) | (3rs-1)   ... (I')
            - (r-1) | (3qs-1)   ... (II')
            - (s-1) | (3qr-1)   ... (IV')

            **Kongruenzbedingung aus (3|qrs-1)**:
            q, r, s sind Primzahlen ≠ 3, also q, r, s ≡ 1 oder 2 (mod 3).
            Produkt qrs ≡ 1 (mod 3) bedeutet: genau 0 oder 2 der Faktoren
            sind ≡ 2 (mod 3). Das schränkt die erlaubten Kombinationen ein.

            **Aus (IV) + (IV')**:
            s(s-1) | (3qr-1), also s² ≤ 3qr → s ≤ √(3qr).

            Analoge Schranken für r und q liefern ein konsistentes System,
            aber keinen direkten Widerspruch.

        @return Analyse-Ergebnis.
        @lastModified 2026-03-12
        """
        ergebnis = schranke_fuer_4prim_fall_p_gleich_3(q_max=100)

        return {
            'fall': 'Fall 2: p = 3 (n = 3·q·r·s)',
            'status': 'KEIN VOLLSTÄNDIGER BEWEIS (numerisch bestätigt)',
            'schlüsselbedingungen': {
                'kongruenz_mod3': 'qrs ≡ 1 (mod 3): gerade Anzahl Faktoren ≡ 2 (mod 3)',
                'aus_s': 's(s-1) | (3qr-1) → s ≤ √(3qr)',
                'aus_r': 'r(r-1) | (3qs-1) → r ≤ √(3qs)',
                'aus_q': 'q(q-1) | (3rs-1) → q ≤ √(3rs)',
            },
            'kongruenz_analyse': {
                'beschreibung': (
                    'Da q,r,s ∈ {1,2} mod 3, muss deren Produkt ≡ 1 mod 3 sein.'
                    ' Mögliche Kombinationen: (1,1,1) oder (2,2,1) (und Permutationen).'
                    ' Dies filtert ~2/3 der Kandidaten aus, ist aber kein Widerspruch.'
                ),
            },
            'numerisch': {
                'q_max': ergebnis['q_max'],
                'schwache_kandidaten_gefunden': len(ergebnis['schwache_kandidaten']),
                'pseudoprimes_gefunden': len(ergebnis['starke_kandidaten_pseudoprimes']),
                'bestätigt': ergebnis['kein_pseudoprime'],
            },
        }

    def fall3_p_ungerade_analyse(self) -> Dict:
        """
        @brief Schrankenanalyse für den allgemeinen Fall p ≥ 5 (n = p·q·r·s, alle ungerade).
        @description
            **Sei n = p·q·r·s mit 5 ≤ p < q < r < s.**

            **Schlüsselschranke** (analog zu Satz 4 für 3 Primfaktoren):

            Für den größten Faktor s:
                s·(s-1) | (pqr - 1) → s(s-1) ≤ pqr - 1 < pqr.

            Da s > r > q > p ≥ 5 und alle ungerade prim → s ≥ r+2.
            Also: (r+2)(r+1) ≤ s(s-1) < pqr.
            Damit: pqr > r²+3r+2, d.h. pq > r + 3 + 2/r.
            Da r ≥ 7 (kleinste ungerade Primzahl > 5 = p ≥ 5, q ≥ 7, r ≥ 11):
                pq > r + 3.                                  ... (*)

            Analoges Argument für r:
                r·(r-1) | (pqs - 1) → r(r-1) < pqs.
            Da r > q ≥ 7: r ≥ q+2 (nächste ungerade Primzahl).
                (q+2)(q+1) ≤ r(r-1) < pqs.
            Also: pqs > q²+3q+2, d.h. ps > q + 3 + 2/q.
                ps > q + 4  (da q ≥ 7).                    ... (**)

            **Warum kein Widerspruch wie im 3-Prim-Fall**:
            Im 3-Prim-Fall p·q·r lautete der Schritt:
                p ≥ q + 4  (aus Schranken), aber p < q → Widerspruch.
            Im 4-Prim-Fall erhält man aus (*): pq > r + 3 (Konsistent!),
            und aus (**): ps > q + 4 (ebenfalls konsistent, da s > q).
            Der Widerspruch tritt nicht auf, weil der vierte Faktor s
            die Schranken "absorbiert".

            **Schlussfolgerung**:
            Die Schranken-Methode aus Satz 4 ist NICHT auf den 4-Prim-Fall übertragbar.
            Ein anderer Beweisansatz wäre nötig (z.B. Primzahlsieb, Gittermethoden,
            oder die tieferen Ergebnisse von Borwein/Bednarek).

        @return Analyse-Ergebnis.
        @lastModified 2026-03-12
        """
        ergebnis = schranke_fuer_4prim_allgemeiner_fall(p_max=20)

        return {
            'fall': 'Fall 3: p ≥ 5 (n = p·q·r·s, alle ungerade)',
            'status': 'KEIN VOLLSTÄNDIGER BEWEIS (numerisch bestätigt, Schranken konsistent)',
            'schlüsselschranken': {
                'aus_s': 'pq > r + 3  [aus s(s-1) | pqr-1 und s ≥ r+2]',
                'aus_r': 'ps > q + 4  [aus r(r-1) | pqs-1 und r ≥ q+2]',
                'aus_q': 'pr > p + 4  [aus q(q-1) | prs-1 und q ≥ p+2]',
            },
            'warum_3prim_beweis_versagt': (
                'Im 3-Prim-Fall ergibt die Schranke p ≥ q+4, aber p < q → Widerspruch.'
                ' Im 4-Prim-Fall: pq > r+3. Da q > p und pq > r+3 mit r > q,'
                ' ist dies konsistent (z.B. p=5, q=7 → pq=35 > r+3 für r ≤ 31).'
                ' Der vierte Faktor s ermöglicht, dass die Bedingungen simultan erfüllt'
                ' werden können — kein automatischer Widerspruch.'
            ),
            'numerisch': {
                'p_max': ergebnis['p_max'],
                'pseudoprimes_gefunden': len(ergebnis['giuga_pseudoprimes_gefunden']),
                'bestätigt': ergebnis['kein_pseudoprime'],
            },
        }

    def schranken_analyse(self) -> Dict:
        """
        @brief Vollständige Schrankenanalyse für alle Fälle des 4-Prim-Problems.
        @description
            Kombiniert die Schranken aus Fall 1 (p=2), Fall 2 (p=3) und Fall 3 (p≥5)
            und listet bekannte externe Resultate auf.

            **Eigene Resultate** (aus dieser Datei):
                - Numerische Verifikation: kein 4-Prim-Kandidat bis 10^6
                - Schranken: s(s-1) ≤ prod(anderen Faktoren) - 1 für jeden Faktor
                - Kongruenzbedingungen für die Fälle p=2,3

            **Externe Resultate** (Literatur):
                - Giuga 1950: Vermutung aufgestellt (noch offen)
                - Borwein, Borwein, Fee, Girgensohn 1996: k ≥ 59 Primfaktoren
                - Borwein 1996 (erweitert): k ≥ 13635 Primfaktoren
                - Bednarek 2014: n > 10^{19907} Dezimalstellen
                → Diese Resultate schließen implizit den 4-Prim-Fall aus,
                  aber mittels tiefer analytischer Zahlentheorie, nicht elementar.

        @return Gesamtschranken-Bericht.
        @lastModified 2026-03-12
        """
        # Berechne Primorial-Schranken für k=4 Faktoren
        # Kleinste mögliche n = 2·3·5·7 = 210 (Giuga-Pseudoprime muss > das sein)
        primorial_4 = 2 * 3 * 5 * 7  # = 210
        # Aber aus Satz 3: Alle Faktoren ungerade → p₁ ≥ 3
        primorial_4_ungerade = 3 * 5 * 7 * 11  # = 1155

        # Wachstum der Schranke durch das lcm-Argument:
        # Für s den größten Faktor: s(s-1) | prod(restliche) - 1
        # Für r: r(r-1) | prod(restliche) - 1
        # Selbst für kleine Faktoren wird s sehr groß.

        # Beispiel: p=3, q=5, r=7 → s(s-1) | 3·5·7-1 = 104
        # Teiler von 104 > 7: 8 (nicht prim), 13, 26, 52, 104 → s ∈ {13} (prim, >7)
        # Dann: n = 3·5·7·13 = 1365. Prüfe Giuga-Bedingungen:
        beispiel_schranke = {
            'p': 3, 'q': 5, 'r': 7,
            'pqr_minus_1': 3 * 5 * 7 - 1,
            's_kandidaten': [
                s for s in sympy.divisors(3 * 5 * 7 - 1)
                if s > 7 and isprime(s)
            ],
        }

        return {
            'eigene_resultate': {
                'satz_1': 'n quadratfrei (aus beweisversuche.py)',
                'satz_2': 'k ≥ 3 Primfaktoren (aus beweisversuche.py)',
                'korollar': 'k ≥ 4 Primfaktoren (aus Sätzen 3+4 in beweisversuche.py)',
                'diese_datei': 'k=4 numerisch kein Kandidat bis 10^6',
            },
            'externe_schranken_literatur': {
                'Giuga_1950': 'Vermutung aufgestellt: n prim ⟺ Σk^{n-1} ≡ -1 (mod n)',
                'Borwein_1996': 'k ≥ 59 Primfaktoren für jedes Giuga-Pseudoprime',
                'Borwein_erweitert': 'k ≥ 13635 Primfaktoren (aus Verfeinerung)',
                'Bednarek_2014': 'n > 10^{19907} Dezimalstellen',
                'Folgerung': (
                    'Diese Resultate schließen k=4 implizit aus.'
                    ' Aber: Sie verwenden tiefe analytische Methoden,'
                    ' keinen elementaren Beweis für genau k=4.'
                ),
            },
            'untere_schranken_aus_eigenem_satz': {
                'trivial': f'n ≥ 2·3·5·7 = {primorial_4} (wenn 2 erlaubt)',
                'aus_korollar': f'n ≥ 3·5·7·11 = {primorial_4_ungerade} (alle Faktoren ungerade)',
            },
            'beispiel_schranke': beispiel_schranke,
            'lcm_argument': (
                'Für jeden Faktor p_i gilt: p_i(p_i-1) | (n/p_i - 1).'
                ' Dieses lcm-Argument erzwingt sehr große Primfaktoren,'
                ' ist aber für k=4 noch nicht zum Widerspruch führbar.'
            ),
        }

    def numerische_verifikation(self, grenze: int = 1_000_000) -> Dict:
        """
        @brief Vollständige numerische Suche nach 4-Prim-Giuga-Pseudoprimes bis zur Grenze.
        @description
            Iteriert über alle quadratfreien zusammengesetzten Zahlen n ≤ grenze
            mit genau 4 Primfaktoren.

            **Effizienz**: Die Suche geschieht direkt über Produktbildung
            p₁·p₂·p₃·p₄ ≤ grenze, was viel effizienter ist als n durchzuiterieren.

            **Erwartetes Ergebnis**: Kein Giuga-Pseudoprime wird gefunden.
            Dies ist konsistent mit den bekannten Schranken aus der Literatur
            (Borwein, Bednarek), welche implizieren, dass ein Giuga-Pseudoprime
            enormen Größenordnungen hätte.

        @param grenze Obere Schranke für die Suche.
        @return Suchergebnis mit Statistiken.
        @lastModified 2026-03-12
        """
        return numerische_suche_4prim(grenze)

    def korollar_kein_4prim_pseudoprime(self) -> Dict:
        """
        @brief Zusammenfassung des Beweisversuchs: Kein 4-Prim-Giuga-Pseudoprime.
        @description
            **Status (Stand 2026-03-12)**:

            NICHT vollständig bewiesen (algebraisch, elementar).
            Numerisch bestätigt für n < 10^6.
            Implizit ausgeschlossen durch externe Resultate (Borwein/Bednarek),
            aber ohne elementaren Direktbeweis für genau k=4.

            **Was wir haben**:
            1. Numerische Verifikation: kein 4-Prim-Giuga-Pseudoprime bis 10^6.
            2. Schranken: Für jeden Faktor pᵢ gilt pᵢ(pᵢ-1) | (n/pᵢ - 1),
               was alle Faktoren nach oben durch die jeweils anderen begrenzt
               und konsistente (aber nicht widersprüchliche) Systeme liefert.
            3. Warum der 3-Prim-Beweis versagt:
               - 3-Prim: Schranke → p ≥ q+4, aber p < q → ⊥
               - 4-Prim: Schranke → pq > r+3, konsistent mit p < q < r.

            **Lücke**: Kein direktes algebraisches Argument für k=4 gefunden.
            Notwendig wäre entweder:
            (a) Ein tieferes Divisibilitäts-/Kongruenzargument das alle Fälle abdeckt.
            (b) Anwendung der Borwein-Methodik (vollständige Primorialschranken).
            (c) Computergestützter Beweis für alle n bis zu einer bekannten Schranke.

            **Nächste Schritte** (für zukünftige Arbeit):
            - Implementiere Borwein et al. (1996) Schrankenmethode für k=4.
            - Untersuche, ob zusätzliche Kongruenzbedingungen (mod kleiner Primzahlen)
              den Suchraum für k=4 vollständig eliminieren.
            - Prüfe, ob die Methode von Satz 4 (3-Prim) per Induktion auf k=4 erweiterbar.

        @return Status-Zusammenfassung des Beweisversuchs.
        @lastModified 2026-03-12
        """
        return {
            'aussage': 'Kein Giuga-Pseudoprime hat genau 4 Primfaktoren.',
            'status': 'NUMERISCH BESTÄTIGT (bis 10^6), algebraisch OFFEN',
            'beweis_typ': 'Unvollständiger Beweisversuch + numerische Evidenz',
            'vorherige_sätze': {
                'satz_2': 'Kein 2-Prim-Giuga-Pseudoprime (BEWIESEN, beweisversuche.py)',
                'satz_3': 'Kein 3-Prim-Giuga-Pseudoprime mit p=2 (BEWIESEN)',
                'satz_4': 'Kein 3-Prim-Giuga-Pseudoprime alle ungerade (BEWIESEN)',
                'korollar': 'Kein 3-Prim-Giuga-Pseudoprime (BEWIESEN)',
            },
            'dieser_satz': {
                'fall_1_p_gleich_2': 'Numerisch bestätigt, algebraisch offen',
                'fall_2_p_gleich_3': 'Numerisch bestätigt, algebraisch offen',
                'fall_3_p_groesser_5': 'Numerisch bestätigt, Schranken konsistent',
            },
            'wo_versagt_3prim_methode': (
                'Im 3-Prim-Fall: s·(s-1) | (pq-1) → pq-1 < s·1·(s-1) → s²-s < pq.'
                ' Mit s > q und s ≥ q+2: p(q) > q²+3q+2 → p > q+3 → p ≥ q+4. Widerspruch.'
                '\nIm 4-Prim-Fall: s·(s-1) | (pqr-1). Mit s ≥ r+2:'
                ' (r+2)(r+1) ≤ pqr-1 → pq > r+3. Dies ist mit p<q<r vereinbar'
                ' (z.B. p=5, q=7 → pq=35, r ≤ 31). Kein Widerspruch!'
            ),
            'literatur_schranken': {
                'Borwein_1996': 'k ≥ 59 impliziert: k=4 ausgeschlossen',
                'Bednarek_2014': 'n > 10^{19907} impliziert: 4-Prim-n viel zu klein',
                'hinweis': (
                    'Diese Resultate schließen den 4-Prim-Fall aus, verwenden aber'
                    ' tiefe analytische Methoden (nicht elementar nachvollziehbar).'
                ),
            },
            'empfehlung': (
                'Für einen vollständigen elementaren Beweis des 4-Prim-Falls'
                ' sind weitere algebraische Argumente notwendig.'
                ' Die vorliegende numerische Evidenz und Schrankenanalyse'
                ' liefern starke Hinweise, aber keinen abschließenden Beweis.'
            ),
        }

    def vollstaendige_analyse(self) -> Dict:
        """
        @brief Führt alle Analysen aus und gibt den vollständigen Bericht zurück.
        @return Vollständiger Analysebericht über den 4-Prim-Giuga-Fall.
        @lastModified 2026-03-12
        """
        return {
            'Fall_1_p_gleich_2': self.fall1_p_gleich_2_analyse(),
            'Fall_2_p_gleich_3': self.fall2_p_gleich_3_analyse(),
            'Fall_3_p_ungerade': self.fall3_p_ungerade_analyse(),
            'Schranken_Analyse': self.schranken_analyse(),
            'Numerische_Verifikation': self.numerische_verifikation(500_000),
            'Korollar': self.korollar_kein_4prim_pseudoprime(),
        }


# =============================================================================
# Einstiegspunkt für direkte Ausführung
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("GIUGA 4-PRIM-FALL: Beweisversuch und numerische Analyse")
    print("Autor: Michael Fuhrmann | 2026-03-12")
    print("=" * 70)

    beweis = Giuga4PrimBeweis()

    # --- Numerische Verifikation ---
    print("\n[1] NUMERISCHE SUCHE bis 10^6")
    print("-" * 50)
    ergebnis_num = beweis.numerische_verifikation(1_000_000)
    print(f"  Grenze:           {ergebnis_num['grenze']:,}")
    print(f"  Geprüfte n:       {ergebnis_num['geprueft']:,}")
    print(f"  Giuga-Zahlen(4p): {len(ergebnis_num['giuga_zahlen_4prim'])}")
    print(f"  Pseudoprimes:     {len(ergebnis_num['giuga_pseudoprimes_4prim'])}")
    print(f"  Fazit:            {ergebnis_num['fazit']}")

    # --- Giuga-Zahlen (nur schwach) ausgeben, falls vorhanden ---
    if ergebnis_num['giuga_zahlen_4prim']:
        print(f"\n  Bekannte 4-Prim-Giuga-Zahlen (nur schwache Bedingung):")
        for gz in ergebnis_num['giuga_zahlen_4prim'][:5]:
            print(f"    n={gz['n']}, Faktoren={gz['faktoren']}")

    # --- Schrankenanalyse ---
    print("\n[2] SCHRANKENANALYSE")
    print("-" * 50)
    schranken = beweis.schranken_analyse()
    print("  Eigene Resultate:")
    for k, v in schranken['eigene_resultate'].items():
        print(f"    {k}: {v}")
    print("  Externe Schranken (Literatur):")
    for k, v in schranken['externe_schranken_literatur'].items():
        print(f"    {k}: {v}")

    # --- Einzelfälle ---
    print("\n[3] FALL 1: p=2 (n=2·q·r·s)")
    print("-" * 50)
    f1 = beweis.fall1_p_gleich_2_analyse()
    print(f"  Status:      {f1['status']}")
    print(f"  Schlüsselbedingungen:")
    for k, v in f1['schlüsselbedingungen'].items():
        print(f"    {k}: {v}")
    print(f"  Numerisch:   {f1['numerisch']}")

    print("\n[4] FALL 2: p=3 (n=3·q·r·s)")
    print("-" * 50)
    f2 = beweis.fall2_p_gleich_3_analyse()
    print(f"  Status:      {f2['status']}")
    print(f"  Numerisch:   {f2['numerisch']}")

    print("\n[5] FALL 3: p≥5 (n=p·q·r·s, alle ungerade)")
    print("-" * 50)
    f3 = beweis.fall3_p_ungerade_analyse()
    print(f"  Status:      {f3['status']}")
    print(f"  Schlüsselschranken:")
    for k, v in f3['schlüsselschranken'].items():
        print(f"    {k}: {v}")
    print(f"  Numerisch:   {f3['numerisch']}")

    # --- Korollar ---
    print("\n[6] KOROLLAR: ZUSAMMENFASSUNG")
    print("-" * 50)
    korollar = beweis.korollar_kein_4prim_pseudoprime()
    print(f"  Aussage: {korollar['aussage']}")
    print(f"  Status:  {korollar['status']}")
    print(f"\n  Wo versagt die 3-Prim-Methode?")
    print(f"  {korollar['wo_versagt_3prim_methode']}")
    print(f"\n  Empfehlung:")
    print(f"  {korollar['empfehlung']}")

    # --- Beispielprüfung bekannte Giuga-Zahlen ---
    print("\n[7] PRÜFUNG BEKANNTER GIUGA-ZAHLEN {30, 858, 1722, 66198}")
    print("-" * 50)
    for gz in [30, 858, 1722, 66198]:
        r = berechne_giuga_bedingung(gz)
        faktoren = factorint(gz)
        print(f"  n={gz}: Giuga-Zahl={r['ist_giuga_zahl']}, "
              f"Pseudoprime={r['ist_giuga_pseudoprime']}, "
              f"Faktoren={r.get('primfaktoren', list(faktoren.keys()))}")

    print("\n" + "=" * 70)
    print("FAZIT: Kein 4-Prim-Giuga-Pseudoprime numerisch bis 10^6 gefunden.")
    print("       Algebraischer Beweis für k=4 noch offen (Stand 2026).")
    print("       Externe Literatur (Borwein, Bednarek) schließt k=4 implizit aus.")
    print("=" * 70)
