"""
@file lehmer_tau.py
@brief Analyse der Lehmer-Tau-Vermutung: τ(n) ≠ 0 für alle n ≥ 1.
@description
    Dieses Modul untersucht die Lehmer-Tau-Vermutung (1947).

    **Lehmer-Tau-Vermutung** (Derrick Henry Lehmer, 1947):
        Die Ramanujan-τ-Funktion erfüllt τ(n) ≠ 0 für alle n ≥ 1.

    **Ramanujan-τ-Funktion**:
        Definiert durch die Erzeugende Funktion:
            Δ(q) = q · Π_{n=1}^∞ (1 - q^n)^24 = Σ_{n=1}^∞ τ(n) q^n

        Dies ist die eindeutige normalisierte Spitzenform des Gewichts 12
        bezüglich der vollen Modulgruppe SL(2,ℤ).

    **Bekannte Werte**:
        τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830,
        τ(6)=-6048, τ(7)=-16744, τ(8)=84480, τ(9)=-113643, τ(10)=-115920

    **Wichtige Eigenschaften**:
        - Multiplikativität: τ(mn) = τ(m)·τ(n) für ggT(m,n)=1 (Hecke-Relationen)
        - Deligne-Schranke: |τ(p)| ≤ 2p^{11/2} für Primzahlen p (Beweis 1974)
        - Viele Kongruenzbedingungen (Swinnerton-Dyer u.a.)

    **Verifizierungsstand**:
        Die Vermutung gilt für alle n ≤ 10^22 (numerisch verifiziert).

    **Kongruenzbedingungen** (bekannte Notwendigkeiten für τ(p) = 0):
        Falls τ(p) = 0 für eine Primzahl p, so gelten simultane Kongruenzen
        bezüglich der Moduln 2, 3, 5, 7, 23, 691 — bislang nie alle erfüllt.

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import sympy
from sympy import isprime, primerange


# ===========================================================================
# Tau-Funktion: Koeffizientenextraktion aus Δ(q)
# ===========================================================================

def _berechne_tau_liste(n_max: int) -> List[int]:
    """
    @brief Berechnet τ(1), ..., τ(n_max) via Polynommultiplikation.
    @description
        Δ(q) = q · Π_{n=1}^∞ (1-q^n)^24

        Berechnet die ersten n_max Koeffizienten durch schrittweise
        Multiplikation mit (1-q^m)^24 für m=1,...,n_max.

        (1-q^m)^24 = Σ_{k=0}^{24} C(24,k)·(-1)^k · q^{k·m}

        Die Koeffizienten sind exakte ganze Zahlen.

    @param n_max  Maximales n für die Berechnung
    @return       Liste [τ(1), τ(2), ..., τ(n_max)]
    @lastModified 2026-03-12
    """
    N = n_max

    # Binomialkoeffizienten C(24,k) für k=0..24
    binom24 = [1] * 25
    for j in range(1, 25):
        binom24[j] = binom24[j - 1] * (25 - j) // j  # C(24,j)

    # Initialisierung: p(q) = 1
    coeffs = [0] * (N + 1)
    coeffs[0] = 1

    # Multipliziere sukzessiv mit (1-q^m)^24
    for m in range(1, N + 1):
        new_coeffs = [0] * (N + 1)
        for k in range(25):
            shift = k * m
            if shift > N:
                break
            sign_binom = ((-1) ** k) * binom24[k]
            for i in range(N + 1 - shift):
                new_coeffs[i + shift] += sign_binom * coeffs[i]
        coeffs = new_coeffs

    # τ(n) = Koeffizient von q^{n-1} in p(q)  [Δ = q·p(q)]
    tau_liste = []
    for n in range(1, n_max + 1):
        if n - 1 <= N:
            tau_liste.append(coeffs[n - 1])
        else:
            tau_liste.append(0)

    return tau_liste


class LehmerTauAnalyse:
    """
    @brief Analyse der Lehmer-Tau-Vermutung τ(n) ≠ 0.
    @description
        Stellt Methoden zur Berechnung der Ramanujan-τ-Funktion,
        numerischer Verifikation der Lehmer-Vermutung sowie
        zur Analyse von Kongruenzbedingungen bereit.

        Die Ramanujan-τ-Funktion ist die Fourier-Entwicklung der
        einzigartigen normalisierten Spitzenform Δ ∈ S_12(SL₂(ℤ)).

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    # Cache für berechnete τ-Werte (verhindert Doppelberechnungen)
    _tau_cache: Dict[int, List[int]] = {}

    def _hole_tau_liste(self, n_max: int) -> List[int]:
        """
        @brief Gibt τ-Liste zurück, verwendet Cache wenn möglich.
        @param n_max  Benötigter Index
        @return       Liste [τ(1),...,τ(n_max)]
        @lastModified 2026-03-12
        """
        # Bestehenden Cache verwenden falls vorhanden und ausreichend
        for cached_max, cached_liste in LehmerTauAnalyse._tau_cache.items():
            if cached_max >= n_max:
                return cached_liste[:n_max]

        # Neu berechnen
        liste = _berechne_tau_liste(n_max)
        LehmerTauAnalyse._tau_cache[n_max] = liste
        return liste

    def berechne_tau(self, n: int) -> int:
        """
        @brief Berechnet τ(n), die n-te Ramanujan-Tau-Funktion.
        @description
            τ(n) ist der n-te Fourier-Koeffizient der Diskriminant-Modulform
                Δ(τ) = q · Π_{k=1}^∞ (1-q^k)^24,   q = e^{2πiτ}

            Berechnung über Polynommultiplikation (exakt, für beliebiges n).

            Bekannte Werte:
                τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830

        @param n  Index ≥ 1
        @return   τ(n) als ganze Zahl
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        liste = self._hole_tau_liste(n)
        return liste[n - 1]

    def verifiziere_lehmer(self, n_max: int) -> Dict:
        """
        @brief Überprüft τ(n) ≠ 0 für alle n = 1, ..., n_max.
        @description
            Berechnet alle τ(n) für n ≤ n_max und sucht nach Nullstellen.
            Eine Nullstelle τ(n₀) = 0 wäre ein Gegenbeispiel zur Lehmer-Vermutung.

            Da τ multiplikativ ist, genügt es, Primzahlpotenzen zu prüfen:
                Falls τ(p) ≠ 0 für alle Primzahlen p ≤ n_max und
                     τ(p^k) ≠ 0 für p^k ≤ n_max,
                dann gilt τ(n) ≠ 0 für alle n ≤ n_max.

        @param n_max  Obere Grenze
        @return       Dictionary mit Verifikationsergebnis und Nullstellenliste
        @lastModified 2026-03-12
        """
        tau_liste = self._hole_tau_liste(n_max)
        nullstellen = []
        erste_werte = []

        for n in range(1, n_max + 1):
            tau_n = tau_liste[n - 1]
            if n <= 20:
                erste_werte.append((n, tau_n))
            if tau_n == 0:
                nullstellen.append({'n': n, 'tau': tau_n, 'ist_prim': isprime(n)})

        return {
            'n_max': n_max,
            'verifiziert': len(nullstellen) == 0,
            'nullstellen': nullstellen,
            'erste_tau_werte': erste_werte,
            'schluss': (
                f'Lehmer-Vermutung verifiziert: τ(n) ≠ 0 für alle n ≤ {n_max}.'
                if not nullstellen
                else f'GEGENBEISPIEL: {nullstellen}'
            )
        }

    def kongruenz_analyse(self, p: int, mod_list: Optional[List[int]] = None) -> Dict:
        """
        @brief Analysiert τ(p) mod m für verschiedene Moduli m.
        @description
            Für eine Primzahl p werden die Reste τ(p) mod m berechnet.

            Bekannte notwendige Bedingungen für τ(p) = 0 (Swinnerton-Dyer):
            - τ(p) ≡ 0 (mod 691)
            - τ(p) ≡ 0 (mod 3)
            - τ(p) ≡ 0 (mod 5) oder τ(p) ≡ ±p^2 (mod 5)
            - τ(p) ≡ 0 (mod 7) oder τ(p) ≡ ±p·? (mod 7)
            - τ(p) ≡ 0 (mod 2^8) (Serre)
            - Weitere...

            Falls τ(p) = 0, müssen ALLE diese Kongruenzen gelten.
            Bisher ist kein p mit all diesen Bedingungen bekannt.

        @param p         Primzahl
        @param mod_list  Liste von Moduli (Standard: [2, 3, 5, 7, 23, 691])
        @return          Dictionary mit Kongruenz-Resten
        @lastModified 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        if mod_list is None:
            mod_list = [2, 3, 5, 7, 11, 13, 23, 691]

        tau_p = self.berechne_tau(p)
        kongruenzen = {}
        for m in mod_list:
            rest = tau_p % m
            kongruenzen[m] = {
                'tau_mod_m': rest,
                'null': rest == 0,
            }

        # Bekannte notwendige Bedingungen für τ(p) = 0
        notwendig_fuer_null = {
            691: kongruenzen.get(691, {}).get('null', False),
            3: kongruenzen.get(3, {}).get('null', False),
            5: kongruenzen.get(5, {}).get('null', False),
        }
        alle_notwendig_erfüllt = all(notwendig_fuer_null.values())

        return {
            'p': p,
            'tau_p': tau_p,
            'kongruenzen': kongruenzen,
            'notwendige_bedingungen_erfüllt': alle_notwendig_erfüllt,
            'könnte_null_sein': alle_notwendig_erfüllt,
            'erklärung': (
                f'τ({p}) = {tau_p}. '
                f'Für τ(p)=0 nötig: τ≡0 mod 691, mod 3, mod 5 (und weitere). '
                f'Alle erfüllt: {alle_notwendig_erfüllt}.'
            )
        }

    def hecke_eigenvalue_schranke(self, p: int) -> Dict:
        """
        @brief Gibt die Deligne-Schranke |τ(p)| ≤ 2·p^{11/2} und prüft sie.
        @description
            **Ramanujan-Vermutung** (1916), bewiesen von **Deligne** (1974):
                |τ(p)| ≤ 2·p^{11/2}  für alle Primzahlen p

            Deligne bewies dies als Spezialfall der Weil-Vermutungen
            (Fields-Medaille 1978).

            Diese Schranke ist scharf: τ(p) ~ 2·p^{11/2}·cos(θ_p) mit θ_p ∈ [0,π].

            Falls τ(p) = 0, dann wäre |τ(p)|/p^{11/2} = 0 ≤ 2, also erlaubt.
            Die Schranke schließt τ(p) = 0 NICHT aus.

        @param p  Primzahl
        @return   Dictionary mit Schrankenangabe und Verifikation
        @lastModified 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        tau_p = self.berechne_tau(p)
        schranke = 2 * (p ** (11 / 2))
        abs_tau = abs(tau_p)
        verhältnis = abs_tau / schranke

        return {
            'p': p,
            'tau_p': tau_p,
            'abs_tau_p': abs_tau,
            'deligne_schranke': schranke,
            'schranke_gilt': abs_tau <= schranke,
            'verhältnis_abs_tau_zu_schranke': round(verhältnis, 6),
            'erklärung': (
                f'|τ({p})| = {abs_tau} ≤ 2·{p}^{{11/2}} ≈ {schranke:.2f}. '
                f'Verhältnis: {verhältnis:.4f}.'
            )
        }

    def ausschluss_durch_kongruenzen(self, p: int) -> Dict:
        """
        @brief Versucht τ(p) = 0 durch Kongruenzbedingungen auszuschließen.
        @description
            Bekannte Kongruenzbedingungen (Swinnerton-Dyer, Serre, u.a.)
            liefern notwendige Bedingungen für τ(p) = 0. Falls auch nur
            eine Bedingung verletzt ist, gilt τ(p) ≠ 0.

            Verwendete Kongruenzen:
            (i)   τ(p) ≡ 1+p^{11} (mod 691)  [Ramanujan-Kongruenz]
                  → τ(p) = 0 erfordert: p^{11} ≡ -1 (mod 691)
            (ii)  τ(p) ≡ 0 (mod 2^3)         [Hecke-Symmetrie mod 2]
            (iii) τ(p) ≡ p^2+p^9 (mod 5)     [mod-5-Kongruenz]
            (iv)  Weitere (mod 7, mod 23, ...)

            Falls eine Bedingung nicht erfüllt ist, kann τ(p) = 0 ausgeschlossen
            werden. Sind alle bekannten Bedingungen erfüllt, bleibt τ(p) = 0
            theoretisch möglich.

        @param p  Primzahl
        @return   Dictionary mit Ausschlussanalyse
        @lastModified 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        tau_p = self.berechne_tau(p)
        ausschluss_gruende = []
        bedingungen = {}

        # (i) mod 691: τ(p) ≡ 1 + p^11 (mod 691)
        ramanujan_691 = (1 + pow(p, 11, 691)) % 691
        tau_mod_691 = tau_p % 691
        bed_691 = tau_mod_691 == ramanujan_691
        bedingungen['mod_691'] = {
            'erwarteter_rest': ramanujan_691,
            'tatsächlicher_rest': tau_mod_691,
            'erfüllt': bed_691,
        }
        if not bed_691:
            ausschluss_gruende.append(f'τ(p) mod 691 = {tau_mod_691} ≠ {ramanujan_691} (erwartet)')

        # (ii) mod 3: τ(p) ≡ 1 + p^{11} (mod 3)?
        ramanujan_3 = (1 + pow(p, 11, 3)) % 3
        tau_mod_3 = tau_p % 3
        bed_3 = tau_mod_3 == ramanujan_3
        bedingungen['mod_3'] = {
            'tau_mod_3': tau_mod_3,
            'erwarteter_rest': ramanujan_3,
            'erfüllt': bed_3,
        }

        # (iii) Für τ(p)=0 muss gelten: tau_mod_691 = 0 → p^{11} ≡ -1 (mod 691)
        falls_null_notwendig_691 = (pow(p, 11, 691) == 690)  # d.h. p^{11} ≡ -1 ≡ 690 (mod 691)
        bedingungen['null_notwendig_691'] = {
            'p_11_mod_691': pow(p, 11, 691),
            'erfordert_p_11_equiv_minus1_mod_691': falls_null_notwendig_691,
        }

        # Kann τ(p) = 0 ausgeschlossen werden?
        ausgeschlossen = tau_p != 0  # Direkte Prüfung
        if tau_p != 0:
            ausschluss_gruende.insert(0, f'Direkte Berechnung: τ({p}) = {tau_p} ≠ 0.')

        return {
            'p': p,
            'tau_p': tau_p,
            'tau_p_gleich_null': tau_p == 0,
            'ausgeschlossen': ausgeschlossen,
            'ausschluss_gruende': ausschluss_gruende,
            'kongruenz_bedingungen': bedingungen,
            'erklärung': (
                f'τ({p}) = {tau_p}. '
                'Durch direkte Berechnung und Kongruenzbedingungen analysiert.'
            )
        }

    def bekannte_kongruenzen_table(self) -> Dict:
        """
        @brief Gibt eine tabellarische Übersicht bekannter Kongruenzbedingungen für τ.
        @description
            Bekannte Kongruenzen für die Ramanujan-τ-Funktion:

            | Modul  | Kongruenz                        | Quelle                  |
            |--------|----------------------------------|-------------------------|
            | 2      | τ(n) ≡ σ_{11}(n) (mod 2^3)      | Ramanujan (1916)        |
            | 3      | τ(n) ≡ σ_{11}(n) (mod 3)        | Ramanujan (1916)        |
            | 5      | τ(n) ≡ σ_{11}(n) (mod 5)        | Ramanujan (1916)        |
            | 7      | τ(n) ≡ σ_{11}(n) (mod 7)        | Ramanujan (1916)        |
            | 23     | τ(p) ≡ 0 oder ±2 (mod 23)        | Swinnerton-Dyer (1973)  |
            | 691    | τ(n) ≡ σ_{11}(n) (mod 691)      | Ramanujan (1916)        |

            σ_{11}(n) = Σ_{d|n} d^{11}  (Teilersumme der 11. Potenzen)

            Hinweis: τ(n) ≡ σ₁₁(n) (mod 3) gilt NICHT für alle n.
            Korrekte mod-3-Kongruenz: τ(p) ≡ 1 + p^{11} (mod 3) für Primzahlen p.
            Die in der Tabelle gelistete Kongruenz mod 3 gilt nur für spezielle n-Klassen.

            Die 691-Kongruenz ist besonders fundamental: 691 teilt den Zähler
            von B₁₂ (Bernoulli-Zahl), was die Kongruenz erklärt.

        @return   Dictionary mit Kongruenztabelle und Erläuterungen
        @lastModified 2026-03-12
        """
        # σ_{11}(n) = Summe der 11. Potenzen der Teiler von n
        def sigma_11(n: int) -> int:
            return sum(d ** 11 for d in sympy.divisors(n))

        # Verifiziere die Kongruenzen für n=1..20
        tau_liste = self._hole_tau_liste(20)
        verifikation = []
        for n in range(1, 21):
            tau_n = tau_liste[n - 1]
            sig_n = sigma_11(n)
            # mod-3 korrekte Formel: τ(p) ≡ 1 + p^{11} (mod 3) gilt nur für Primzahlen p
            # Für nicht-Primzahlen gibt es komplexere Relationen
            ist_prim_n = sympy.isprime(n)
            mod3_korrekt = (tau_n % 3 == (1 + pow(n, 11, 3)) % 3) if ist_prim_n else None
            verifikation.append({
                'n': n,
                'tau_n': tau_n,
                'sigma_11_n': sig_n,
                'mod_691_stimmt': tau_n % 691 == sig_n % 691,
                'mod_3_stimmt': tau_n % 3 == sig_n % 3,  # Gilt nicht für alle n!
                'mod_3_primformel_stimmt': mod3_korrekt,  # Gilt für Primzahlen
                'mod_7_stimmt': tau_n % 7 == sig_n % 7,
            })

        kongruenzen_tabelle = [
            {'modul': 2, 'kongruenz': 'τ(n) ≡ σ₁₁(n) (mod 8)',
             'quelle': 'Ramanujan 1916', 'für_null_notwendig': 'τ ≡ 0 (mod 8)'},
            {'modul': 3, 'kongruenz': 'τ(n) ≡ σ₁₁(n) (mod 3)',
             'quelle': 'Ramanujan 1916', 'für_null_notwendig': 'σ₁₁(n) ≡ 0 (mod 3)'},
            {'modul': 5, 'kongruenz': 'τ(n) ≡ σ₁₁(n) (mod 5)',
             'quelle': 'Ramanujan 1916', 'für_null_notwendig': 'σ₁₁(n) ≡ 0 (mod 5)'},
            {'modul': 7, 'kongruenz': 'τ(n) ≡ σ₁₁(n) (mod 7)',
             'quelle': 'Ramanujan 1916', 'für_null_notwendig': 'σ₁₁(n) ≡ 0 (mod 7)'},
            {'modul': 23, 'kongruenz': 'τ(p) ≡ 0 oder ±2 (mod 23)',
             'quelle': 'Swinnerton-Dyer 1973', 'für_null_notwendig': 'τ(p) ≡ 0 (mod 23)'},
            {'modul': 691, 'kongruenz': 'τ(n) ≡ σ₁₁(n) (mod 691)',
             'quelle': 'Ramanujan 1916', 'für_null_notwendig': 'σ₁₁(n) ≡ 0 (mod 691)'},
        ]

        return {
            'kongruenzen_tabelle': kongruenzen_tabelle,
            'verifikation_erste_n': verifikation,
            'erklärung': (
                'Diese Kongruenzen sind notwendige Bedingungen für τ(p) = 0. '
                'Alle müssen gleichzeitig gelten. Bislang ist kein p bekannt, '
                'das alle Bedingungen erfüllt. Verifiziert bis n ≤ 10^22.'
            ),
            'bernoulli_691': (
                'B₁₂ = -691/2730, der Zähler 691 erklärt die 691-Kongruenz '
                'durch die Eisenstein-Reihe E₁₂ ≡ 1 + (65520/691)·q + ... (mod 691).'
            )
        }


# ===========================================================================
# HAUPTPROGRAMM (Demo)
# ===========================================================================

if __name__ == "__main__":
    """
    Demonstriert die LehmerTauAnalyse-Klasse.
    """
    ana = LehmerTauAnalyse()

    print("=" * 60)
    print("LEHMER-TAU-VERMUTUNG — Analyse")
    print("=" * 60)

    # 1. τ-Werte für kleine n
    print("\n1. Ramanujan-τ(n) für n = 1..15:")
    for n in range(1, 16):
        tau = ana.berechne_tau(n)
        print(f"   τ({n:2d}) = {tau}")

    # 2. Verifikation bis n = 1000
    print("\n2. Lehmer-Verifikation bis n = 1000...")
    erg = ana.verifiziere_lehmer(1000)
    print(f"   Verifiziert: {erg['verifiziert']}")
    print(f"   Nullstellen: {erg['nullstellen']}")
    print(f"   {erg['schluss']}")

    # 3. Kongruenzanalyse für p=2,3,5,7
    print("\n3. Kongruenzanalyse:")
    for p in [2, 3, 5, 7, 11, 13]:
        k = ana.kongruenz_analyse(p)
        print(f"   p={p}: τ(p)={k['tau_p']}, könnte null sein: {k['könnte_null_sein']}")

    # 4. Deligne-Schranke für erste Primzahlen
    print("\n4. Deligne-Schranken |τ(p)| ≤ 2·p^{11/2}:")
    for p in [2, 3, 5, 7, 11, 13]:
        hs = ana.hecke_eigenvalue_schranke(p)
        print(f"   p={p}: |τ(p)|={hs['abs_tau_p']}, Schranke≈{hs['deligne_schranke']:.0f}, "
              f"gilt: {hs['schranke_gilt']}")

    # 5. Kongruenztabelle
    print("\n5. Bekannte Kongruenzbedingungen:")
    tab = ana.bekannte_kongruenzen_table()
    print("   Verifikation τ(n) ≡ σ₁₁(n) (mod 691):")
    for v in tab['verifikation_erste_n'][:5]:
        print(f"   n={v['n']}: τ={v['tau_n']}, σ₁₁={v['sigma_11_n']}, "
              f"691-Kongr: {v['mod_691_stimmt']}")

    # 6. Ausschlussanalyse
    print("\n6. Ausschlussanalyse für p=5:")
    aus = ana.ausschluss_durch_kongruenzen(5)
    print(f"   τ(5) = {aus['tau_p']}")
    print(f"   Ausgeschlossen: {aus['ausgeschlossen']}")
    for grund in aus['ausschluss_gruende']:
        print(f"   Grund: {grund}")
