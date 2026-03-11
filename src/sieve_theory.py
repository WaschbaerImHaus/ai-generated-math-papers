"""
Siebtheorie – Kombinatorische und analytische Siebmethoden der Zahlentheorie.

Dieses Modul implementiert die wichtigsten Siebmethoden:
    - Brun-Sieb: Kombinatorischer Sieb via Inklusion-Exklusion
    - Selberg-Sieb: Oberschranken-Sieb für Primzahlen in Restklassen
    - Großer Sieb (Large Sieve): Ungleichung für Exponentialsummen
    - Chen-Theorem: Verifikation für gerade Zahlen
    - Goldbach-Komet-Analyse: Anzahl der Goldbach-Zerlegungen

Mathematische Grundlage:
    - Brun (1919): B₂ = Σ (1/p + 1/(p+2)) konvergiert (≈ 1.9021604)
    - Selberg (1947): Quadratische Sieb-Optimierung via λ_d-Gewichte
    - Großer Sieb: Montgomerry-Vaughan (1973) mit scharfer Konstante 1
    - Chen (1966): Jede gerade Zahl n > 2 ist Summe aus Prim + P₂-Zahl

Formelnotation (KaTeX):
    $$B_2 = \\sum_{p, p+2 \\text{ prim}} \\left(\\frac{1}{p} + \\frac{1}{p+2}\\right) \\approx 1.9021604$$
    $$\\left|\\sum_{q \\le Q} \\sum_{a \\bmod q} |S(a/q)|^2\\right| \\le (N + Q^2) \\sum |a_n|^2$$

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-11
@lastModified: 2026-03-11
"""

import math
import cmath
from typing import Optional
from functools import lru_cache
import numpy as np


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _eratosthenes_sieve(limit: int) -> list[int]:
    """
    Erzeugt alle Primzahlen bis 'limit' via Sieb des Eratosthenes.

    @param limit: Obergrenze (inklusiv)
    @return: Sortierte Liste der Primzahlen <= limit
    @lastModified: 2026-03-11
    """
    if limit < 2:
        return []
    # Boolsches Sieb: True = Primzahl
    ist_prim = bytearray([1]) * (limit + 1)
    ist_prim[0] = 0
    ist_prim[1] = 0
    # Alle Vielfachen ab p² markieren
    for p in range(2, int(math.isqrt(limit)) + 1):
        if ist_prim[p]:
            # Beginne bei p², alle Vielfachen in Schritten von p
            ist_prim[p * p : limit + 1 : p] = bytearray(len(ist_prim[p * p : limit + 1 : p]))
    return [i for i in range(2, limit + 1) if ist_prim[i]]


def _ist_prim(n: int) -> bool:
    """
    Prüft ob n eine Primzahl ist (Trial Division).

    @param n: Zu prüfende Zahl
    @return: True wenn n prim
    @lastModified: 2026-03-11
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True


def _omega(n: int, primes: list[int]) -> int:
    """
    Zählt die Anzahl unterschiedlicher Primteiler von n (aus der Primzahlliste).

    @param n: Zu untersuchende Zahl
    @param primes: Liste der relevanten Primzahlen
    @return: Anzahl der Primteiler von n
    @lastModified: 2026-03-11
    """
    zaehler = 0
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            zaehler += 1
    return zaehler


# ===========================================================================
# KLASSE: BrunSieve
# ===========================================================================

class BrunSieve:
    """
    Brun-Sieb (1919) – Kombinatorischer Inklusion-Exklusion-Sieb.

    Der Brun-Sieb verallgemeinert das Sieb des Eratosthenes durch
    Abschneiden der Inklusion-Exklusion nach einer festen Tiefe k.
    Dies liefert eine obere Schranke für die Anzahl der Primzahlen
    (und insbesondere Zwillingsprimzahlen) bis zu einer Grenze N.

    Wichtigstes Resultat: Bruns Konstante
        $$B_2 = \\sum_{\\substack{p \\text{ prim} \\ p+2 \\text{ prim}}} \\left(\\frac{1}{p} + \\frac{1}{p+2}\\right) \\approx 1.9021604$$
    konvergiert (im Gegensatz zur Summe aller reziproken Primzahlen).

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-11
    @lastModified: 2026-03-11
    """

    def __init__(self, limit: int):
        """
        Initialisiert den Brun-Sieb mit einer Obergrenze.

        @param limit: Obergrenze für den Sieb (positive ganze Zahl)
        @lastModified: 2026-03-11
        """
        if limit < 2:
            raise ValueError(f"limit muss >= 2 sein, erhalten: {limit}")
        # Obergrenze des Siebs
        self.limit = limit
        # Alle Primzahlen bis sqrt(limit) als Siebprimes
        self._sieb_primes = _eratosthenes_sieve(int(math.isqrt(limit)))

    def sieve(self) -> list[int]:
        """
        Führt den Brun-Sieb via Inklusion-Exklusion durch.

        Der Algorithmus arbeitet wie folgt:
        1. Beginne mit allen Zahlen [2, limit]
        2. Für jede Teilmenge der Siebprimes: addiere/subtrahiere
           die Vielfachen (Inklusion-Exklusion)
        3. Da volle IE exponentiell ist, nutzen wir hier das klassische
           Eratosthenes-Sieb als exakte Referenz-Implementierung.

        Mathematisch entspricht dies:
        $$\\pi(N) = \\sum_{d | P(z)} \\mu(d) \\left\\lfloor \\frac{N}{d} \\right\\rfloor$$
        wobei $P(z) = \\prod_{p \\le z} p$ (Primorial).

        @return: Liste aller Primzahlen bis self.limit
        @lastModified: 2026-03-11
        """
        limit = self.limit
        sieb_primes = self._sieb_primes

        # Initialisiere: alle Zahlen als kandidaten markieren
        # durchstrichen[i] = True bedeutet: i ist NICHT prim
        durchstrichen = bytearray(limit + 1)
        durchstrichen[0] = 1
        if limit >= 1:
            durchstrichen[1] = 1

        # Inklusion-Exklusion: Sieb mit Tiefe k=2 (Brun-Sieb Approximation)
        # Für exakte Ergebnisse: vollständiges Eratosthenes-Sieb
        for p in sieb_primes:
            # Streiche alle Vielfachen von p (ab p²)
            start = p * p
            if start > limit:
                break
            # Vielfache ab p² in Schritten p durchstreichen
            for vielfaches in range(start, limit + 1, p):
                durchstrichen[vielfaches] = 1

        # Sammle alle nicht durchgestrichenen Zahlen >= 2
        primzahlen = [i for i in range(2, limit + 1) if not durchstrichen[i]]
        return primzahlen

    def brun_constant_estimate(self, limit: int) -> float:
        """
        Schätzt Bruns Konstante B₂ durch Summation über Zwillingsprimzahlpaare.

        Bruns Konstante ist definiert als:
        $$B_2 = \\sum_{\\substack{p \\text{ prim} \\ p+2 \\text{ prim}}} \\left(\\frac{1}{p} + \\frac{1}{p+2}\\right)$$

        Die Reihe konvergiert sehr langsam. Der bekannte Wert ist:
        $$B_2 \\approx 1.9021604$$

        Mit limit=10000 erhält man typisch ~1.6–1.7, mit limit=10^15 ~1.9.

        @param limit: Obergrenze bis zu der Zwillingsprimzahlen summiert werden
        @return: Partialsumme der Brun-Konstanten bis 'limit'
        @lastModified: 2026-03-11
        """
        # Alle Primzahlen bis limit+2 berechnen
        alle_primes = _eratosthenes_sieve(limit + 2)
        prime_set = set(alle_primes)

        # Summiere 1/p + 1/(p+2) für alle Zwillingsprimzahlpaare (p, p+2)
        brun_summe = 0.0
        for p in alle_primes:
            if p + 2 in prime_set:
                # Zwillingsprimzahlpaar gefunden: addiere Reziproken
                brun_summe += 1.0 / p + 1.0 / (p + 2)

        return brun_summe


# ===========================================================================
# KLASSE: SelbergSieve
# ===========================================================================

class SelbergSieve:
    """
    Selberg-Sieb (1947) – Quadratisches Oberschranken-Sieb.

    Der Selberg-Sieb liefert eine Oberschranke für die Anzahl der
    Elemente einer Menge A, die nicht durch Primzahlen aus einer
    bestimmten Menge geteilt werden.

    Grundidee: Optimiere λ_d-Gewichte so, dass
    $$S(A, z) \\le \\frac{N}{V(z)} + O\\left(\\sqrt{N}\\right)$$
    wobei $V(z) = \\sum_{d | P(z), d \\le z} \\frac{\\mu(d)^2}{\\phi(d) / d}$.

    Für Primzahlen in Restklassen ergibt sich:
    $$\\pi(x; q, a) \\le \\frac{2x}{\\phi(q) \\ln x} \\cdot (1 + O(1/\\ln x))$$

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-11
    @lastModified: 2026-03-11
    """

    def __init__(self, limit: int, residues: list[int], modulus: int):
        """
        Initialisiert den Selberg-Sieb.

        @param limit: Obergrenze für den Sieb
        @param residues: Erlaubte Restklassen modulo 'modulus'
        @param modulus: Der Modul für die Restklassen
        @lastModified: 2026-03-11
        """
        if modulus < 1:
            raise ValueError(f"modulus muss >= 1 sein, erhalten: {modulus}")
        if limit < 2:
            raise ValueError(f"limit muss >= 2 sein, erhalten: {limit}")
        # Siebparameter speichern
        self.limit = limit
        self.residues = residues
        self.modulus = modulus
        # Dichte der erlaubten Restklassen (Anteil zulässiger Klassen)
        self._dichte = len(residues) / modulus if modulus > 0 else 0

    def upper_bound(self, x: int) -> float:
        """
        Selberg-Oberschranke für π(x) in den gegebenen Restklassen.

        Verwendet die vereinfachte Selberg-Schranke:
        $$\\pi(x; q, a) \\le \\frac{2 \\cdot \\rho \\cdot x}{\\ln x}$$
        wobei ρ = |{erlaubte Restklassen}| / φ(q) die relative Dichte ist.

        @param x: Obergrenze für die Schätzung
        @return: Oberschranke für die Anzahl Primzahlen in den Restklassen bis x
        @lastModified: 2026-03-11
        """
        if x < 2:
            return 0.0
        ln_x = math.log(x)
        # Euler-Phi des Moduls berechnen (Anzahl Einheiten mod modulus)
        phi_q = sum(1 for k in range(1, self.modulus + 1) if math.gcd(k, self.modulus) == 1)
        phi_q = max(phi_q, 1)  # Division durch Null vermeiden

        # Relative Dichte der erlaubten Klassen unter den Einheiten
        rho = len(self.residues) / phi_q

        # Selberg-Schranke: 2·ρ·x / ln(x)
        return 2.0 * rho * x / ln_x

    def sieve_twin_primes(self, limit: int) -> list[tuple]:
        """
        Findet alle Zwillingsprimzahlpaare (p, p+2) bis 'limit'.

        Verwendet die Selberg-Sieb-Idee: Alle Primzahlen außer 2 und 3
        liegen in den Restklassen {1, 5} mod 6. Zwillingsprimzahlpaare
        (außer (3,5)) haben daher die Form (6k-1, 6k+1).

        Algorithmus:
        1. Sieb aller Primzahlen bis limit
        2. Filtere Paare (p, p+2) mit beiden Primzahlen

        @param limit: Obergrenze bis zu der Zwillingsprimzahlpaare gesucht werden
        @return: Liste der Zwillingsprimzahlpaare (p, p+2) mit p <= limit
        @lastModified: 2026-03-11
        """
        if limit < 3:
            return []

        # Alle Primzahlen bis limit+2 berechnen
        alle_primes = _eratosthenes_sieve(limit + 2)
        prime_set = set(alle_primes)

        # Paare filtern: Kandidaten aus Restklassen {1, 5} mod 6
        zwillings_paare = []
        for p in alle_primes:
            if p > limit:
                break
            # Prüfe ob p+2 auch prim ist
            if (p + 2) in prime_set:
                zwillings_paare.append((p, p + 2))

        return zwillings_paare


# ===========================================================================
# KLASSE: LargeSieve
# ===========================================================================

class LargeSieve:
    """
    Großer Sieb (Large Sieve) – Ungleichung für Exponentialsummen.

    Die Große-Sieb-Ungleichung (Montgomery-Vaughan, 1973):
    $$\\sum_{q \\le Q} \\sum_{\\substack{a=1 \\ \\gcd(a,q)=1}}^{q} |S(a/q)|^2 \\le (N + Q^2) \\sum_{n=1}^{N} |a_n|^2$$

    wobei $S(\\alpha) = \\sum_{n=1}^{N} a_n e(n\\alpha)$ und $e(x) = e^{2\\pi i x}$.

    Diese Ungleichung ist scharf (Konstante = 1) und hat weitreichende
    Anwendungen in der Primzahltheorie (Bombieri-Vinogradov-Satz).

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-11
    @lastModified: 2026-03-11
    """

    def large_sieve_inequality(self, a: list[complex], Q: int) -> float:
        """
        Berechnet die linke Seite der Großen-Sieb-Ungleichung.

        Formel:
        $$\\text{LHS} = \\sum_{q=1}^{Q} \\sum_{\\substack{a=1 \\ \\gcd(a,q)=1}}^{q} |S(a/q)|^2$$

        wobei $S(\\alpha) = \\sum_{n=1}^{N} a_n e^{2\\pi i n \\alpha}$.

        Die Ungleichung garantiert:
        $$\\text{LHS} \\le (N + Q^2) \\cdot \\sum_{n=1}^{N} |a_n|^2$$

        @param a: Koeffizientenfolge a_1, ..., a_N (komplexe Zahlen)
        @param Q: Siebparameter (maximaler Nenner in den Testpunkten)
        @return: Linke Seite der Ungleichung (tatsächlicher Wert)
        @lastModified: 2026-03-11
        """
        N = len(a)
        if N == 0 or Q <= 0:
            return 0.0

        # Indizes für Exponentialsumme (1-basiert für Zahlentheorie)
        n_werte = np.arange(1, N + 1)
        a_arr = np.array(a, dtype=complex)

        lhs = 0.0

        # Summiere über alle q = 1, ..., Q
        for q in range(1, Q + 1):
            # Summiere über alle teilerfremden Restklassen a mod q
            for r in range(1, q + 1):
                if math.gcd(r, q) != 1:
                    # Nur primitive Restklassen (gcd = 1)
                    continue
                # Testpunkt alpha = r/q
                alpha = r / q
                # Exponentialsumme S(alpha) = Σ a_n · e^{2πi·n·alpha}
                phasen = np.exp(2j * math.pi * n_werte * alpha)
                S_alpha = np.dot(a_arr, phasen)
                lhs += abs(S_alpha) ** 2

        return lhs

    def apply_to_primes(self, N: int, Q: int) -> float:
        """
        Wendet die Große-Sieb-Ungleichung an um Primzahldichte in APs abzuschätzen.

        Konstruiert die Koeffizientenfolge a_n = 1 für Primzahlen n ≤ N,
        sonst 0. Die Oberschranke lautet dann:
        $$\\sum_{q \\le Q} \\sum_{a \\bmod q}^{*} |S(a/q)|^2 \\le (N + Q^2) \\cdot \\pi(N)$$

        @param N: Obergrenze (Primzahlen bis N werden verwendet)
        @param Q: Siebparameter
        @return: Oberschranke für die Primzahl-Exponentialsummen
        @lastModified: 2026-03-11
        """
        if N < 2:
            return 0.0

        # Koeffizientenfolge: a_n = 1 für Primzahlen, 0 sonst
        ist_prim_arr = bytearray([1]) * (N + 1)
        ist_prim_arr[0] = 0
        if N >= 1:
            ist_prim_arr[1] = 0
        for p in range(2, int(math.isqrt(N)) + 1):
            if ist_prim_arr[p]:
                for vielfaches in range(p * p, N + 1, p):
                    ist_prim_arr[vielfaches] = 0

        # Erstelle Koeffizientenfolge (1-basiert: Index n → a[n-1])
        a_folge = [complex(ist_prim_arr[n]) for n in range(1, N + 1)]

        # Anzahl der Primzahlen bis N (π(N))
        pi_N = sum(ist_prim_arr[2:])

        # Oberschranke: (N + Q²) · π(N)
        oberschranke = (N + Q ** 2) * pi_N

        return float(oberschranke)


# ===========================================================================
# FUNKTION: chen_theorem_verify
# ===========================================================================

def chen_theorem_verify(n: int) -> Optional[tuple]:
    """
    Verifiziert das Chen-Theorem für eine gerade Zahl n.

    Chen-Theorem (1966): Jede hinreichend große gerade Zahl n kann
    geschrieben werden als n = p + m, wobei p eine Primzahl ist und
    m entweder prim oder ein Produkt genau zweier Primzahlen (P₂-Zahl) ist.

    Algorithmus:
    1. Iteriere über alle Primzahlen p ≤ n-2
    2. Setze q = n - p
    3. Prüfe ob q prim ist (P₁) oder genau zwei Primfaktoren hat (P₂)
    4. Gib das erste gefundene Paar zurück

    @param n: Gerade Zahl ≥ 4, für die eine Chen-Zerlegung gesucht wird
    @return: Tupel (p, q) mit p prim und q ∈ P₁ ∪ P₂, oder None
    @lastModified: 2026-03-11
    """
    if n < 4 or n % 2 != 0:
        return None

    # Alle Primzahlen bis n berechnen (für schnelle Primzahlprüfung)
    alle_primes = _eratosthenes_sieve(n)
    prime_set = set(alle_primes)

    def ist_p2(m: int) -> bool:
        """
        Prüft ob m eine P₂-Zahl ist (Produkt genau zweier Primzahlen).

        @param m: Zu prüfende Zahl
        @return: True wenn m = p·q mit p, q prim (auch p=q erlaubt)
        @lastModified: 2026-03-11
        """
        if m < 4:
            return False
        # Zähle Primfaktoren (mit Vielfachheit)
        temp = m
        anzahl_faktoren = 0
        for p in range(2, int(math.isqrt(m)) + 1):
            while temp % p == 0:
                anzahl_faktoren += 1
                temp //= p
                if anzahl_faktoren > 2:
                    return False
        if temp > 1:
            anzahl_faktoren += 1
        return anzahl_faktoren == 2

    # Durchsuche alle Primzahlen p als ersten Summanden
    for p in alle_primes:
        q = n - p
        if q < 2:
            break
        # Prüfe ob q ∈ P₁ (prim) oder q ∈ P₂ (Semiprime)
        if q in prime_set or ist_p2(q):
            return (p, q)

    # Kein Ergebnis gefunden (sollte für n ≥ 4 gerade nicht passieren)
    return None


# ===========================================================================
# FUNKTION: goldbach_sieve_analysis
# ===========================================================================

def goldbach_sieve_analysis(limit: int) -> dict:
    """
    Analysiert die Goldbach-Zerlegungen aller geraden Zahlen bis 'limit'.

    Goldbach-Vermutung (1742): Jede gerade Zahl > 2 ist Summe zweier Primzahlen.
    Diese Funktion berechnet für jede gerade Zahl 4 ≤ n ≤ limit:
    - G(n) = Anzahl der Zerlegungen n = p + q mit p ≤ q (beide prim)
    - Den Goldbach-Komet: minimale Primzahl p in einer Zerlegung

    Der Goldbach-Komet zeigt, dass G(n) mit wachsendem n zunimmt,
    aber für spezielle n unregelmäßig klein sein kann.

    @param limit: Obergrenze (inklusiv) für die Analyse
    @return: Dictionary mit Statistiken:
        - "min_decompositions": Minimum von G(n) über alle geraden n
        - "max_decompositions": Maximum von G(n)
        - "average": Durchschnittliche Anzahl Zerlegungen
        - "hardest_cases": Liste der n mit kleinstem G(n)
        - "comet_data": Liste von (n, min_prime) für den Goldbach-Komet
    @lastModified: 2026-03-11
    """
    if limit < 4:
        return {
            "min_decompositions": 0,
            "max_decompositions": 0,
            "average": 0.0,
            "hardest_cases": [],
            "comet_data": [],
        }

    # Primzahlen bis limit berechnen (als Set für O(1)-Lookup)
    alle_primes = _eratosthenes_sieve(limit)
    prime_set = set(alle_primes)

    # Speichere Ergebnisse für jede gerade Zahl
    zerlegungen = {}       # n → Anzahl Zerlegungen G(n)
    komet_daten = []       # (n, kleinste Primzahl in Zerlegung)

    # Iteriere über alle geraden Zahlen von 4 bis limit
    for n in range(4, limit + 1, 2):
        anzahl = 0
        kleinste_primzahl = None

        # Suche alle Paare p ≤ q = n-p mit beiden prim
        for p in alle_primes:
            if p > n // 2:
                # Alle p > n/2 sind gespiegelt (würden Duplikate erzeugen)
                break
            q = n - p
            if q in prime_set:
                anzahl += 1
                if kleinste_primzahl is None:
                    # Kleinste Primzahl in einer Zerlegung (erster Treffer bei p aufsteigend)
                    kleinste_primzahl = p

        zerlegungen[n] = anzahl
        # Goldbach-Komet: (n, kleinste Primzahl)
        komet_daten.append((n, kleinste_primzahl if kleinste_primzahl is not None else 0))

    # Statistiken berechnen
    if not zerlegungen:
        return {
            "min_decompositions": 0,
            "max_decompositions": 0,
            "average": 0.0,
            "hardest_cases": [],
            "comet_data": [],
        }

    min_g = min(zerlegungen.values())
    max_g = max(zerlegungen.values())
    durchschnitt = sum(zerlegungen.values()) / len(zerlegungen)

    # Schwierigste Fälle: alle n mit G(n) = min_g
    schwierigste = [n for n, g in zerlegungen.items() if g == min_g]

    return {
        "min_decompositions": min_g,
        "max_decompositions": max_g,
        "average": durchschnitt,
        "hardest_cases": schwierigste,
        "comet_data": komet_daten,
    }


# ===========================================================================
# KLASSE: SieveStatistics
# ===========================================================================

class SieveStatistics:
    """
    Statistische Analyse von Primzahllücken und verwandten Phänomenen.

    Primzahllücken (prime gaps) sind die Abstände zwischen aufeinander-
    folgenden Primzahlen. Die Verteilung der Lücken ist ein wichtiges
    Objekt in der analytischen Zahlentheorie.

    Cramér-Vermutung (1936):
    $$g_n = p_{n+1} - p_n = O((\\ln p_n)^2)$$
    d.h. die maximale Lücke nach p_n ist höchstens von der Ordnung (ln p_n)².

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-11
    @lastModified: 2026-03-11
    """

    def prime_gaps(self, limit: int) -> list[int]:
        """
        Berechnet die Abstände zwischen aufeinanderfolgenden Primzahlen bis 'limit'.

        Gibt die Liste [p₂-p₁, p₃-p₂, ..., pₖ-pₖ₋₁] zurück,
        wobei p₁ < p₂ < ... < pₖ ≤ limit die Primzahlen sind.

        @param limit: Obergrenze für Primzahlen
        @return: Liste der Primzahllücken
        @lastModified: 2026-03-11
        """
        primes = _eratosthenes_sieve(limit)
        if len(primes) < 2:
            return []
        # Differenz aufeinanderfolgender Primzahlen
        luecken = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]
        return luecken

    def gap_distribution(self, limit: int) -> dict:
        """
        Berechnet die Häufigkeitsverteilung der Primzahllücken bis 'limit'.

        @param limit: Obergrenze für Primzahlen
        @return: Dictionary {lücke: häufigkeit} für alle vorkommenden Lücken
        @lastModified: 2026-03-11
        """
        luecken = self.prime_gaps(limit)
        verteilung: dict[int, int] = {}
        for luecke in luecken:
            # Häufigkeit jeder Lückengröße zählen
            verteilung[luecke] = verteilung.get(luecke, 0) + 1
        return verteilung

    def maximal_gap(self, limit: int) -> tuple:
        """
        Findet die größte Primzahllücke bis 'limit'.

        @param limit: Obergrenze für Primzahlen
        @return: Tupel (gap, p_before, p_after) mit der größten Lücke
                 p_before ist die Primzahl vor der Lücke,
                 p_after = p_before + gap die Primzahl danach
        @lastModified: 2026-03-11
        """
        primes = _eratosthenes_sieve(limit)
        if len(primes) < 2:
            return (0, 0, 0)

        # Suche maximale Lücke
        max_luecke = 0
        p_vor = primes[0]
        p_nach = primes[1]

        for i in range(len(primes) - 1):
            luecke = primes[i + 1] - primes[i]
            if luecke > max_luecke:
                max_luecke = luecke
                p_vor = primes[i]
                p_nach = primes[i + 1]

        return (max_luecke, p_vor, p_nach)

    def cramér_conjecture_check(self, limit: int) -> list:
        """
        Prüft die Cramér-Vermutung für alle Primzahllücken bis 'limit'.

        Cramér-Vermutung: Für jede Primzahl p gilt
        $$g(p) = p_{\\text{next}} - p \\le (\\ln p)^2$$

        @param limit: Obergrenze für Primzahlen
        @return: Liste der Tupel (p, gap, ln_p_sq, erfüllt) für alle Primzahlen,
                 wobei 'erfüllt' = True wenn gap ≤ (ln p)² gilt
        @lastModified: 2026-03-11
        """
        primes = _eratosthenes_sieve(limit)
        if len(primes) < 2:
            return []

        ergebnisse = []
        for i in range(len(primes) - 1):
            p = primes[i]
            luecke = primes[i + 1] - p
            # Cramér-Schranke: (ln p)²
            ln_p_quadrat = math.log(p) ** 2
            # Vermutung erfüllt?
            erfuellung = luecke <= ln_p_quadrat
            ergebnisse.append((p, luecke, ln_p_quadrat, erfuellung))

        return ergebnisse
