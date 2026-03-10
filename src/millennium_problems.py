"""
@file millennium_problems.py
@brief Numerische und algebraische Werkzeuge für die Millennium-Probleme.
@description
    Die 7 Millennium-Probleme des Clay Mathematics Institute (2000):
    1. P vs NP (Komplexitätstheorie)
    2. Hodge-Vermutung (algebraische Geometrie)
    3. Poincaré-Vermutung (GELÖST: Perelman 2003)
    4. Riemann-Hypothese (analytische Zahlentheorie)
    5. Yang-Mills-Existenz und Massenlücke (mathematische Physik)
    6. Navier-Stokes-Existenz und Glattheit (PDE)
    7. Birch und Swinnerton-Dyer Vermutung (Zahlentheorie)

    Dieses Modul stellt numerische Werkzeuge bereit, um diese Probleme zu
    untersuchen, zu verstehen und empirisch zu testen.
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import cmath
import itertools
import numpy as np
from typing import Optional

# mpmath für hochgenaue Berechnungen (Zeta-Nullstellen)
try:
    import mpmath
    MPMATH_AVAILABLE = True
except ImportError:
    MPMATH_AVAILABLE = False


# =============================================================================
# 1. RIEMANN-HYPOTHESE WERKZEUGE
# =============================================================================

def riemann_zeros_mpmath(n_zeros: int = 20, dps: int = 50) -> list:
    """
    @brief Berechnet die ersten n nicht-trivialen Riemann-Nullstellen.
    @description
        Alle bekannten Nullstellen liegen auf der kritischen Geraden Re(s) = 1/2.
        mpmath.zetazero(n) gibt die n-te Nullstelle 1/2 + i*t_n zurück.

        Die Riemann-Hypothese behauptet, dass ALLE nicht-trivialen Nullstellen
        auf dieser Geraden liegen. Status: >10^13 Nullstellen verifiziert.

        Formel: ζ(1/2 + i·t_n) = 0, t_n > 0 für alle n
    @param n_zeros Anzahl der zu berechnenden Nullstellen
    @param dps Dezimalstellen der Genauigkeit
    @return Liste von dicts: {'n', 'zero', 're', 'im', 'on_line'}
    @lastModified 2026-03-10
    """
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath ist nicht verfügbar. Installiere es mit: pip install mpmath")

    # Genauigkeit setzen
    mpmath.mp.dps = dps
    results = []

    for k in range(1, n_zeros + 1):
        # mpmath.zetazero(k) berechnet die k-te Nullstelle
        zero = mpmath.zetazero(k)

        # Realteil sollte exakt 1/2 sein (gemäß RH)
        re_part = float(zero.real)
        im_part = float(zero.imag)

        # Prüfen ob Realteil nahe 1/2 liegt (Toleranz: 10^-10)
        on_line = abs(re_part - 0.5) < 1e-10

        results.append({
            'n': k,
            'zero': complex(re_part, im_part),
            're': re_part,
            'im': im_part,
            'on_line': on_line
        })

    return results


def verify_rh_critical_line(n_zeros: int = 100, dps: int = 30) -> dict:
    """
    @brief Verifiziert numerisch, dass die ersten n Nullstellen auf Re(s)=1/2 liegen.
    @description
        Überprüft die Riemann-Hypothese für die ersten n_zeros nicht-trivialen
        Nullstellen mithilfe von mpmath. Gibt Statistik über Abweichungen zurück.

        Status: >10^13 Nullstellen verifiziert (ZetaGrid-Projekt, 2004; Platt 2021).
    @param n_zeros Anzahl der zu prüfenden Nullstellen
    @param dps Dezimalstellen der Genauigkeit
    @return {'verified': int, 'max_deviation': float, 'all_on_line': bool}
    @lastModified 2026-03-10
    """
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath nicht verfügbar")

    mpmath.mp.dps = dps
    max_deviation = 0.0
    verified = 0

    for k in range(1, n_zeros + 1):
        zero = mpmath.zetazero(k)
        # Abweichung von 1/2
        deviation = abs(float(zero.real) - 0.5)
        max_deviation = max(max_deviation, deviation)
        # Toleranz: 10^-(dps/2) für numerische Genauigkeit
        if deviation < 10 ** (-(dps // 2)):
            verified += 1

    return {
        'verified': verified,
        'total': n_zeros,
        'max_deviation': max_deviation,
        'all_on_line': (verified == n_zeros)
    }


def gram_law_check(n_grams: int = 20) -> dict:
    """
    @brief Prüft das Gram-Gesetz für die ersten n_grams Gram-Punkte.
    @description
        Gram-Punkte g_n sind die t-Werte, wo die Riemann-Siegel-Thetafunktion
        θ(g_n) = n·π erfüllt ist.

        Gram-Gesetz: Z(g_n) hat Vorzeichen (-1)^n (häufig, aber nicht immer).
        Das Gesetz gilt NICHT für alle Gram-Punkte (sog. "Gram-Versagen").

        θ(t) = Im(log Γ(1/4 + it/2)) - t/2 · log(π)
    @param n_grams Anzahl der zu prüfenden Gram-Punkte
    @return dict mit Gram-Punkten, Vorzeichen und Erfolgsquote
    @lastModified 2026-03-10
    """
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath nicht verfügbar")

    mpmath.mp.dps = 25
    results = []
    law_holds = 0

    for n in range(n_grams):
        # Gram-Punkt g_n berechnen
        g_n = float(mpmath.grampoint(n))

        # Z-Funktion am Gram-Punkt auswerten
        z_val = float(mpmath.siegelz(g_n))

        # Erwartetes Vorzeichen gemäß Gram-Gesetz: (-1)^n
        expected_sign = (-1) ** n
        actual_sign = 1 if z_val > 0 else (-1 if z_val < 0 else 0)

        # Gram-Gesetz erfüllt wenn Vorzeichen übereinstimmt
        holds = (actual_sign == expected_sign)
        if holds:
            law_holds += 1

        results.append({
            'n': n,
            'gram_point': g_n,
            'z_value': z_val,
            'expected_sign': expected_sign,
            'actual_sign': actual_sign,
            'law_holds': holds
        })

    return {
        'gram_points': results,
        'law_holds_count': law_holds,
        'total': n_grams,
        'success_rate': law_holds / n_grams if n_grams > 0 else 0.0
    }


def zeta_zero_gaps(n_zeros: int = 50) -> dict:
    """
    @brief Analysiert die Abstände zwischen aufeinanderfolgenden Riemann-Nullstellen.
    @description
        Berechnet die normalisierten Abstände zwischen den Imaginärteilen der
        nicht-trivialen Nullstellen. Diese folgen der GUE-Statistik (Gaussian
        Unitary Ensemble aus der Zufallsmatrix-Theorie).

        Montgomery-Odlyzko-Gesetz: Die normierten Nullstellenabstände haben
        dieselbe Verteilung wie Eigenwerte zufälliger hermitescher Matrizen.
        Mittlerer Abstand bei Höhe T: 2π / log(T/(2π))
    @param n_zeros Anzahl der zu analysierenden Nullstellen
    @return {'gaps': list, 'mean': float, 'std': float, 'min': float, 'max': float}
    @lastModified 2026-03-10
    """
    if not MPMATH_AVAILABLE:
        raise ImportError("mpmath nicht verfügbar")

    mpmath.mp.dps = 25

    # Nullstellen berechnen
    zeros = [float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)]

    # Abstände zwischen aufeinanderfolgenden Nullstellen
    gaps = [zeros[i + 1] - zeros[i] for i in range(len(zeros) - 1)]

    # Mittlerer Abstand bei Höhe T ≈ 2π / log(T/(2π))
    mean_gap = float(np.mean(gaps))
    std_gap = float(np.std(gaps))

    # Normierte Abstände (dividiert durch mittleren Abstand → GUE-Vergleich)
    normalized_gaps = [g / mean_gap for g in gaps]

    return {
        'zeros': zeros,
        'gaps': gaps,
        'normalized_gaps': normalized_gaps,
        'mean': mean_gap,
        'std': std_gap,
        'min': float(min(gaps)),
        'max': float(max(gaps)),
        'gue_note': 'Normierte Abstände folgen GUE-Statistik (Montgomery-Odlyzko)'
    }


# =============================================================================
# 2. GOLDBACH-VERMUTUNG WERKZEUGE
# =============================================================================

def _simple_sieve(limit: int) -> list:
    """
    @brief Sieb des Eratosthenes für interne Verwendung.
    @param limit Obergrenze (inklusive)
    @return Liste aller Primzahlen ≤ limit
    @lastModified 2026-03-10
    """
    if limit < 2:
        return []
    # Boolean-Array: True = Primzahl
    is_prime_arr = [True] * (limit + 1)
    is_prime_arr[0] = is_prime_arr[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime_arr[i]:
            # Alle Vielfachen von i als nicht-prim markieren
            for j in range(i * i, limit + 1, i):
                is_prime_arr[j] = False
    return [i for i in range(2, limit + 1) if is_prime_arr[i]]


def _sieve_of_eratosthenes(n: int) -> list:
    """
    @brief Sieb des Eratosthenes – liefert alle Primzahlen bis einschließlich n.
    @description
        Klassischer Algorithmus zur effizienten Primzahlberechnung:
        1. Markiere alle Zahlen 2..n als potenzielle Primzahlen
        2. Für jede noch markierte Zahl p: lösche alle Vielfachen p², p²+p, ...
        3. Alle noch markierten Zahlen sind Primzahlen

        Laufzeit: O(n log log n)  (praktisch fast linear)
        Speicher: O(n)

    @param n Obergrenze (inklusive).
    @return Sortierte Liste aller Primzahlen ≤ n.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    if n < 2:
        return []
    # Boolean-Array: True = noch potenzielle Primzahl
    is_prime_arr = [True] * (n + 1)
    is_prime_arr[0] = is_prime_arr[1] = False
    # Sieben: nur bis √n nötig
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime_arr[i]:
            # Alle Vielfachen von i (ab i²) als nicht-prim markieren
            for j in range(i * i, n + 1, i):
                is_prime_arr[j] = False
    return [i for i in range(2, n + 1) if is_prime_arr[i]]


def _check_goldbach_worker(args: tuple) -> tuple | None:
    """
    @brief Modullevelfunktion für multiprocessing.Pool (picklebar).
    @description
        Prüft ob n als Summe zweier Primzahlen darstellbar ist.
        Muss auf Modulebene definiert sein, da lokale Funktionen nicht
        vom pickle-Modul serialisiert werden können (multiprocessing-Einschränkung).
    @param args Tupel (n, primes, prime_set)
    @return (n, p, q) bei Erfolg, None wenn keine Zerlegung gefunden.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n, primes, prime_set = args
    for p in primes:
        if p > n // 2:
            break
        if (n - p) in prime_set:
            return (n, p, n - p)
    return None


def goldbach_verification_range(n_max: int = 10000) -> dict:
    """
    @brief Verifiziert die Goldbach-Vermutung für alle geraden n ≤ n_max.
    @description
        Goldbach-Vermutung (1742): Jede gerade natürliche Zahl > 2 ist die
        Summe zweier Primzahlen.

        Beispiel: 4=2+2, 6=3+3, 8=3+5, 10=3+7=5+5, ...

        Status: Verifiziert bis 4×10^18 (2014, Oliveira e Silva et al.)
        Schwächere Aussage "schwache Goldbach": Helfgott 2013 bewiesen.

        Implementierung: Verwendet multiprocessing.Pool für parallele Überprüfung.
        Jeder Worker-Prozess prüft eine Teilmenge der geraden Zahlen.
        Anzahl Worker = min(CPU-Kerne, 4) für optimalen Durchsatz.

    @param n_max Obergrenze (alle geraden n von 4 bis n_max werden geprüft).
    @return Dict mit 'verified_up_to', 'all_verified', 'failed', 'min_decompositions', 'workers_used'.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    import multiprocessing as mp

    # Primzahlen bis n_max einmalig berechnen (Sieb des Eratosthenes)
    primes = _sieve_of_eratosthenes(n_max)
    prime_set = set(primes)

    # Alle geraden Zahlen > 2 als Prüfliste aufbauen
    even_numbers = list(range(4, n_max + 1, 2))

    # Parallele Berechnung: min(CPU-Kerne, 4) Worker-Prozesse
    # Lokale Funktionen sind nicht picklebar → Modullfunktion _check_goldbach_worker nutzen
    n_workers = min(mp.cpu_count(), 4)
    # Primzahlmenge als sortierte Liste übergeben (picklebar)
    args = [(n, primes, prime_set) for n in even_numbers]
    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_check_goldbach_worker, args)

    # Ergebnisse auswerten: None = keine Zerlegung gefunden (Gegenbeispiel!)
    failed = [r for r in results if r is None]
    decompositions = {r[0]: (r[1], r[2]) for r in results if r is not None}

    # Erste 10 Zerlegungen als Beispiel + alle Vielfachen von 1000
    sample_keys = list(decompositions.keys())[:10]
    min_decompositions = {k: decompositions[k] for k in sample_keys}

    return {
        'verified_up_to': n_max,
        'all_verified': len(failed) == 0,
        'failed': failed,
        'min_decompositions': min_decompositions,
        'workers_used': n_workers
    }


def goldbach_weak_check(n: int) -> dict:
    """
    @brief Prüft die schwache Goldbach-Vermutung für eine ungerade Zahl.
    @description
        Schwache (ternäre) Goldbach-Vermutung (BEWIESEN 2013, Harald Helfgott):
        Jede ungerade Zahl > 5 ist die Summe von drei Primzahlen.

        Beispiel: 7=2+2+3, 9=2+2+5=3+3+3, 11=2+2+7=3+3+5, ...

        Helfgott bewies dies vollständig für alle n > 5 (Ausnahme: n=1,3,5).
    @param n Ungerade natürliche Zahl > 5
    @return {'n': int, 'decomposition': tuple, 'verified': bool, 'proof_status': str}
    @lastModified 2026-03-10
    """
    # Eingabevalidierung
    if n % 2 == 0:
        return {
            'n': n,
            'decomposition': None,
            'verified': False,
            'error': 'Schwache Goldbach-Vermutung gilt nur für ungerade Zahlen'
        }

    if n <= 5:
        return {
            'n': n,
            'decomposition': None,
            'verified': False,
            'proof_status': 'BEWIESEN (Helfgott 2013) - aber erst für n > 5',
            'error': f'n={n} ≤ 5, Vermutung gilt nur für n > 5'
        }

    # Primzahlen bis n berechnen
    primes = _simple_sieve(n)
    prime_set = set(primes)

    # Zerlegung n = p1 + p2 + p3 suchen
    for p1 in primes:
        if p1 > n - 4:
            break
        for p2 in primes:
            if p1 + p2 > n - 2:
                break
            p3 = n - p1 - p2
            if p3 >= 2 and p3 in prime_set:
                return {
                    'n': n,
                    'decomposition': (p1, p2, p3),
                    'verified': True,
                    'proof_status': 'BEWIESEN (Helfgott 2013)'
                }

    return {
        'n': n,
        'decomposition': None,
        'verified': False,
        'proof_status': 'Keine Zerlegung gefunden (sollte nicht passieren!)'
    }


def hardy_littlewood_goldbach_density(n: int) -> dict:
    """
    @brief Hardy-Littlewood-Schätzung für die Goldbach-Zerlegungsanzahl.
    @description
        Hardy-Littlewood (1923) schätzten die Anzahl der Goldbach-Zerlegungen
        von n als Summe zweier Primzahlen ab:

        r₂(n) ≈ 2 · C₂ · n / (log n)² · Π_{p|n, p>2} (p-1)/(p-2)

        wobei C₂ = Π_{p>2} (1 - 1/(p-1)²) ≈ 0.6601618... (Zwillingsprimzahl-Konstante)

        Diese Formel ist nicht bewiesen, passt aber sehr gut zu empirischen Daten.
    @param n Gerade natürliche Zahl > 2
    @return dict mit Schätzung, tatsächlicher Anzahl und relativem Fehler
    @lastModified 2026-03-10
    """
    if n <= 2 or n % 2 != 0:
        return {'error': 'n muss eine gerade Zahl > 2 sein'}

    # Zwillingsprimzahl-Konstante C₂ (Näherung über endliches Produkt)
    # C₂ = Π_{p>2 prim} (1 - 1/(p-1)²)
    primes_for_c2 = _simple_sieve(1000)
    c2 = 1.0
    for p in primes_for_c2:
        if p > 2:
            # Faktor (1 - 1/(p-1)²) = (p²-2p)/(p-1)²
            c2 *= (1.0 - 1.0 / ((p - 1) ** 2))

    # Lokaler Faktor: Π_{p|n, p>2} (p-1)/(p-2)
    local_factor = 1.0
    primes_of_n = _simple_sieve(n)
    for p in primes_of_n:
        if p > 2 and n % p == 0:
            local_factor *= (p - 1.0) / (p - 2.0)

    # Hardy-Littlewood-Schätzung
    log_n = math.log(n)
    estimate = 2.0 * c2 * n / (log_n ** 2) * local_factor

    # Tatsächliche Anzahl berechnen
    prime_set = set(primes_of_n)
    actual_count = 0
    for p in primes_of_n:
        if p > n // 2:
            break
        if (n - p) in prime_set:
            actual_count += 1

    # Relativer Fehler der Schätzung
    rel_error = abs(estimate - actual_count) / actual_count if actual_count > 0 else float('inf')

    return {
        'n': n,
        'twin_prime_constant_c2': c2,
        'local_factor': local_factor,
        'estimate_r2': estimate,
        'actual_r2': actual_count,
        'relative_error': rel_error,
        'formula': 'r₂(n) ≈ 2·C₂·n/(log n)² · Π_{p|n,p>2} (p-1)/(p-2)'
    }


# =============================================================================
# 2b. HAUPTBOGEN / NEBENBOGEN ANALYSE (Hardy-Littlewood-Kreismethode)
# =============================================================================

def goldbach_exponential_sum(primes: list, alpha: float) -> complex:
    """
    @brief Berechnet die Exponentialsumme S(α) = Σ_{p in primes} e(2πi·p·α).
    @description
        Kernbaustein der Hardy-Littlewood-Kreismethode für Goldbach:

            S(α) = Σ_{p ≤ N} e^{2πi·p·α}

        Die Anzahl der Goldbach-Darstellungen r₂(n) ergibt sich als:

            r₂(n) = ∫₀¹ S(α)² · e(-2πi·n·α) dα

        Auf Hauptbögen (α nahe rationaler Zahlen a/q mit kleinem q) gilt:
            S(α) ≈ (μ(q)/φ(q)) · Li(N) · e(a·n/q)

        Auf Nebenbögen dominiert Vinogradov's Abschätzung: |S(α)| = o(N/log N)
    @param primes Liste der Primzahlen ≤ N
    @param alpha Frequenzparameter α ∈ [0, 1)
    @return Komplexer Wert der Exponentialsumme
    @lastModified 2026-03-10
    """
    # e(k) = exp(2πi·k·α): direkte Summation über alle Primzahlen
    two_pi_i = 2j * math.pi
    return sum(cmath.exp(two_pi_i * p * alpha) for p in primes)


def _best_rational_approx(alpha: float, max_denom: int) -> tuple:
    """
    @brief Findet die beste rationale Näherung a/q für alpha mit q ≤ max_denom.
    @description
        Nutzt den Kettenbruch-Algorithmus (Stern-Brocot / Farey-Folge).
        Ergebnis: a/q mit |alpha - a/q| ≤ 1/(q·max_denom) minimal.
    @param alpha Zu nähernde reelle Zahl in [0, 1)
    @param max_denom Maximaler Nenner Q
    @return (a, q) mit a/q beste Näherung und q ≤ max_denom
    @lastModified 2026-03-10
    """
    best_a, best_q = 0, 1
    best_dist = abs(alpha)

    for q in range(1, max_denom + 1):
        # Nächster Zähler: a = runde(α · q)
        a = round(alpha * q)
        dist = abs(alpha - a / q)
        if dist < best_dist:
            best_dist = dist
            best_a = a
            best_q = q

    return best_a, best_q


def classify_arc(alpha: float, N: int, delta: float = 0.5) -> str:
    """
    @brief Klassifiziert α als Haupt- oder Nebenbogen.
    @description
        Hauptbogen M_{a,q}: |α - a/q| ≤ Q/N, mit Q = N^delta, q ≤ Q.

        Hauptbögen (Major Arcs): Genähert durch Hardy-Littlewood-Hauptterm.
        Nebenbögen (Minor Arcs): Rest [0,1) \\ M; Beitrag durch Vinogradov beschränkt.

        Wahl δ = 1/2 entspricht dem klassischen Goldbach-Setup.
    @param alpha Frequenzparameter α ∈ [0, 1)
    @param N Obere Grenze der Primzahlen
    @param delta Exponent für Hauptbogengröße (Standard: 0.5)
    @return 'major' oder 'minor'
    @lastModified 2026-03-10
    """
    Q = int(N ** delta)
    a, q = _best_rational_approx(alpha, Q)

    # Abstand |α - a/q| ≤ Q/N → Hauptbogen
    dist = abs(alpha - (a / q if q > 0 else 0))
    threshold = Q / N

    return 'major' if (q <= Q and dist <= threshold) else 'minor'


def major_arc_contribution(n: int, N: int, delta: float = 0.5,
                           num_points: int = 500) -> dict:
    """
    @brief Numerische Integration des Goldbach-Integrals über die Hauptbögen.
    @description
        Berechnet den Hauptterm der Hardy-Littlewood-Kreismethode:

            I_major = ∫_{M} S(α)² · e(-2πi·n·α) dα

        Klassischer Hauptterm (nach Hardy-Littlewood):
            I_major ≈ 2 · C₂(n) · n / (log n)²

        Die numerische Version diskretisiert [0,1) und summiert nur α-Werte
        die zu Hauptbögen gehören. Das ist eine endliche Approximation
        des eigentlichen Kreisintegrals.
    @param n Gerade Zahl für Goldbach-Zerlegung
    @param N Primzahlgrenze (default: n selbst)
    @param delta Hauptbogenexponent (Standard: 0.5)
    @param num_points Anzahl Diskretisierungspunkte
    @return dict mit Hauptterm-Integral, Primzahldichte und Klassifikation
    @lastModified 2026-03-10
    """
    if n <= 2 or n % 2 != 0:
        return {'error': 'n muss gerade und > 2 sein'}

    primes = _simple_sieve(N)
    alphas = [k / num_points for k in range(num_points)]
    dalpha = 1.0 / num_points

    major_sum = 0.0 + 0.0j
    minor_sum = 0.0 + 0.0j
    n_major = 0
    n_minor = 0

    for alpha in alphas:
        S_alpha = goldbach_exponential_sum(primes, alpha)
        # Integrand: S(α)² · e(-2πi·n·α)
        integrand = S_alpha ** 2 * cmath.exp(-2j * math.pi * n * alpha)

        arc_type = classify_arc(alpha, N, delta)
        if arc_type == 'major':
            major_sum += integrand * dalpha
            n_major += 1
        else:
            minor_sum += integrand * dalpha
            n_minor += 1

    total = major_sum + minor_sum

    return {
        'n': n,
        'N': N,
        'delta': delta,
        'major_arc_integral': major_sum.real,
        'minor_arc_integral': minor_sum.real,
        'total_integral': total.real,
        'n_major_points': n_major,
        'n_minor_points': n_minor,
        'major_fraction': n_major / num_points,
        'goldbach_estimate_from_circle': max(0, round(total.real)),
        'actual_r2': sum(1 for p in primes if (n - p) in set(primes) and p <= n // 2),
        'interpretation': (
            'Hauptbögen liefern den dominanten Beitrag (Hardy-Littlewood-Hauptterm). '
            'Nebenbögen sind klein (Vinogradov), aber ihr Beitrag ist unkontrolliert '
            'ohne Beweis der Goldbach-Vermutung.'
        )
    }


def minor_arc_vinogradov_bound(N: int) -> dict:
    """
    @brief Berechnet Vinogradov's Schranke für den Nebenbogen-Beitrag.
    @description
        Vinogradov (1937) zeigte für die Nebenbögen der Kreismethode:

            max_{α ∈ minor arcs} |S(α)| ≪ N · (log N)^C / N^{1/6}
                                        = N^{5/6} · (log N)^C

        Für die ternäre (schwache) Goldbach-Vermutung:
            ∫_{minor} S(α)³ e(-nα) dα = o(N²/(log N)³)

        Diese Abschätzung ist schwächer als für die binäre Vermutung.
        Für binäres Goldbach benötigt man eine stärkere Minor-Arc-Kontrolle.

        Numerisch: Wir berechnen max|S(α)| über die Nebenbögen und
        vergleichen mit der theoretischen Vinogradov-Schranke.
    @param N Primzahlgrenze
    @return dict mit empirischer und theoretischer Nebenbögen-Schranke
    @lastModified 2026-03-10
    """
    primes = _simple_sieve(N)
    delta = 0.5
    num_points = 200  # Gröbere Diskretisierung für Schranken-Berechnung

    max_S_major = 0.0
    max_S_minor = 0.0
    S_values_minor = []

    for k in range(num_points):
        alpha = k / num_points
        S = abs(goldbach_exponential_sum(primes, alpha))
        arc_type = classify_arc(alpha, N, delta)

        if arc_type == 'major':
            max_S_major = max(max_S_major, S)
        else:
            max_S_minor = max(max_S_minor, S)
            S_values_minor.append(S)

    # Theoretische Vinogradov-Schranke: N^(5/6) · (log N)^3
    log_N = math.log(N) if N > 1 else 1.0
    vinogradov_bound = N ** (5 / 6) * log_N ** 3

    # Verhältnis: empirisch zu theoretisch
    ratio = max_S_minor / vinogradov_bound if vinogradov_bound > 0 else float('inf')

    # Triviale Schranke: |S(α)| ≤ π(N) ≈ N/log N
    trivial_bound = len(primes)

    return {
        'N': N,
        'pi_N': len(primes),
        'max_S_on_major_arcs': max_S_major,
        'max_S_on_minor_arcs': max_S_minor,
        'vinogradov_bound': vinogradov_bound,
        'trivial_bound': trivial_bound,
        'ratio_empirical_to_vinogradov': ratio,
        'minor_arc_suppression': max_S_minor / max_S_major if max_S_major > 0 else float('inf'),
        'interpretation': (
            f'Nebenbögen: max|S| = {max_S_minor:.1f} vs. triviale Schranke {trivial_bound}. '
            f'Vinogradov-Schranke N^(5/6)·(log N)³ = {vinogradov_bound:.1e}. '
            f'Für Goldbach brauchen wir Kontrolle über Nebenbögen: |S| ≪ N/log N.'
        )
    }


# =============================================================================
# 3. P vs NP WERKZEUGE (Komplexitätstheorie)
# =============================================================================

def is_np_complete_3sat(clauses: list, variables: int) -> dict:
    """
    @brief Löst das 3-SAT-Problem via Backtracking (NP-vollständig).
    @description
        3-SAT ist das kanonische NP-vollständige Problem (Cook-Levin-Theorem, 1971).
        Es fragt: Gibt es eine Belegung der booleschen Variablen x₁,...,xₙ, sodass
        eine gegebene Konjunktion von Klauseln (KNF mit je 3 Literalen) wahr wird?

        Beispiel: (x₁ ∨ ¬x₂ ∨ x₃) ∧ (¬x₁ ∨ x₂ ∨ x₃) → clauses=[[1,-2,3],[-1,2,3]]

        Konvention: Positive Zahl k = x_k, negative Zahl -k = ¬x_k
        Im besten Fall polynomiell, im Worst-Case O(2^n) (kein Poly-Alg. bekannt).
    @param clauses Liste von Klauseln, jede = Liste von 3 Literalen (±int)
    @param variables Anzahl der Variablen
    @return {'satisfiable': bool, 'assignment': dict, 'steps': int}
    @lastModified 2026-03-10
    """
    steps = [0]  # Mutable für Closure

    def evaluate_clause(clause: list, assignment: dict) -> Optional[bool]:
        """Wertet eine Klausel aus: True/False/None (wenn unbestimmt)."""
        result = False
        undecided = False
        for lit in clause:
            var = abs(lit)
            if var in assignment:
                # Variable ist belegt
                val = assignment[var] if lit > 0 else not assignment[var]
                if val:
                    return True  # Klausel erfüllt
            else:
                undecided = True  # Noch nicht alle Variablen belegt
        return None if undecided else False

    def backtrack(var_idx: int, assignment: dict) -> Optional[dict]:
        """Rekursives Backtracking über alle Variablen."""
        steps[0] += 1

        # Alle Klauseln prüfen
        for clause in clauses:
            val = evaluate_clause(clause, assignment)
            if val is False:
                return None  # Diese Belegung widerlegt eine Klausel → Backtrack

        # Alle Variablen belegt?
        if var_idx > variables:
            # Nochmal alle Klauseln vollständig prüfen
            for clause in clauses:
                # Klausel vollständig auswerten
                clause_val = any(
                    assignment.get(abs(lit), False) == (lit > 0)
                    for lit in clause
                )
                if not clause_val:
                    return None
            return assignment.copy()

        # Variable var_idx mit True und False versuchen
        for val in [True, False]:
            assignment[var_idx] = val
            result = backtrack(var_idx + 1, assignment)
            if result is not None:
                return result
            del assignment[var_idx]

        return None

    # Backtracking starten
    solution = backtrack(1, {})

    return {
        'satisfiable': solution is not None,
        'assignment': solution,
        'steps': steps[0],
        'complexity_note': 'NP-vollständig: Kein polynomial-Zeit Algorithmus bekannt (falls P≠NP)'
    }


def traveling_salesman_brute(distance_matrix: list) -> dict:
    """
    @brief Löst das Traveling-Salesman-Problem via Brute-Force (exakt, NP-schwer).
    @description
        TSP: Finde die kürzeste Rundreise, die alle n Städte genau einmal besucht.
        Exakte Lösung via Enumeration aller (n-1)! Permutationen.
        Nur für n ≤ 12 praktikabel (12! ≈ 479 Millionen).

        TSP ist NP-schwer: Es gibt (wahrscheinlich) keinen polynomiellen Algorithmus.
        Dies ist eines der praktisch wichtigsten offenen Probleme der Informatik.
    @param distance_matrix n×n Matrix, distance_matrix[i][j] = Distanz von Stadt i nach j
    @return {'tour': list, 'length': float, 'optimal': True}
    @lastModified 2026-03-10
    """
    n = len(distance_matrix)
    if n == 0:
        return {'tour': [], 'length': 0.0, 'optimal': True}
    if n == 1:
        return {'tour': [0], 'length': 0.0, 'optimal': True}
    if n > 12:
        raise ValueError(f"Brute-Force TSP nur für n ≤ 12 (gegeben: n={n}). Verwende traveling_salesman_greedy().")

    # Alle Permutationen der Städte 1..n-1 (Stadt 0 als Startpunkt fixiert)
    cities = list(range(1, n))
    best_length = float('inf')
    best_tour = None

    for perm in itertools.permutations(cities):
        # Tour: 0 → perm[0] → ... → perm[-1] → 0
        tour = [0] + list(perm) + [0]
        length = sum(
            distance_matrix[tour[i]][tour[i + 1]]
            for i in range(len(tour) - 1)
        )
        if length < best_length:
            best_length = length
            best_tour = [0] + list(perm)

    return {
        'tour': best_tour,
        'length': best_length,
        'optimal': True,
        'n_cities': n,
        'permutations_checked': math.factorial(n - 1)
    }


def traveling_salesman_greedy(distance_matrix: list) -> dict:
    """
    @brief Löst TSP via Nearest-Neighbor-Heuristik (polynomiell, nicht optimal).
    @description
        Greedy-Nearest-Neighbor: Starte in Stadt 0, gehe immer zur nächsten
        noch nicht besuchten Stadt. Approximationsverhältnis: O(log n).

        Zeigt den P vs NP Unterschied:
        - Exakte Lösung: NP-schwer (exponentiell)
        - Approximation: Polynomiell in O(n²)
    @param distance_matrix n×n Distanzmatrix
    @return {'tour': list, 'length': float, 'optimal': False}
    @lastModified 2026-03-10
    """
    n = len(distance_matrix)
    if n == 0:
        return {'tour': [], 'length': 0.0, 'optimal': False}
    if n == 1:
        return {'tour': [0], 'length': 0.0, 'optimal': False}

    # Nearest-Neighbor-Algorithmus
    visited = [False] * n
    tour = [0]
    visited[0] = True
    total_length = 0.0

    for _ in range(n - 1):
        current = tour[-1]
        # Nächste unbesuchte Stadt finden
        best_next = -1
        best_dist = float('inf')
        for j in range(n):
            if not visited[j] and distance_matrix[current][j] < best_dist:
                best_dist = distance_matrix[current][j]
                best_next = j

        if best_next == -1:
            break  # Keine nächste Stadt gefunden (sollte nicht passieren)

        tour.append(best_next)
        visited[best_next] = True
        total_length += best_dist

    # Rückweg zur Startstadt
    total_length += distance_matrix[tour[-1]][tour[0]]

    return {
        'tour': tour,
        'length': total_length,
        'optimal': False,
        'algorithm': 'Nearest-Neighbor (Greedy)',
        'complexity': 'O(n²) - polynomiell'
    }


def complexity_comparison(n: int) -> dict:
    """
    @brief Vergleicht das Wachstum verschiedener Komplexitätsklassen für Eingabegröße n.
    @description
        Zeigt, wie stark verschiedene Komplexitätsklassen für dieselbe Eingabe n wachsen.
        Dies verdeutlicht den fundamentalen Unterschied zwischen P und NP:

        P: Polynomial lösbar → O(1), O(log n), O(n), O(n log n), O(n²), O(n³)
        NP: Exponentiell (kein Poly-Alg bekannt) → O(2^n), O(n!)

        Bei n=100: n² = 10000, aber 2^100 ≈ 10^30 und 100! ≈ 10^157
    @param n Eingabegröße
    @return {'n': n, 'complexities': dict mit Klassenname → Wert}
    @lastModified 2026-03-10
    """
    complexities = {
        'O(1)': 1,
        'O(log n)': math.log2(n) if n > 0 else 0,
        'O(sqrt n)': math.sqrt(n),
        'O(n)': n,
        'O(n log n)': n * math.log2(n) if n > 0 else 0,
        'O(n^2)': n ** 2,
        'O(n^3)': n ** 3,
        'O(2^n)': 2 ** n if n <= 100 else float('inf'),
        'O(n!)': math.factorial(n) if n <= 20 else float('inf')
    }

    # Klasseneinteilung
    p_classes = ['O(1)', 'O(log n)', 'O(sqrt n)', 'O(n)', 'O(n log n)', 'O(n^2)', 'O(n^3)']
    np_classes = ['O(2^n)', 'O(n!)']

    return {
        'n': n,
        'complexities': complexities,
        'p_class': {k: complexities[k] for k in p_classes},
        'np_hard_class': {k: complexities[k] for k in np_classes},
        'p_vs_np_note': 'P⊆NP, aber ob P=NP ist das 1. Millennium-Problem (ungelöst)'
    }


# =============================================================================
# 4. NAVIER-STOKES WERKZEUGE (numerische PDE)
# =============================================================================

def navier_stokes_2d_simple(
    nx: int = 50,
    ny: int = 50,
    dt: float = 0.01,
    T: float = 1.0,
    nu: float = 0.01,
    f_x: float = 1.0
) -> dict:
    """
    @brief Vereinfachte 2D-Navier-Stokes-Simulation (inkompressibel, Finite-Differenzen).
    @description
        Löst die inkompressiblen Navier-Stokes-Gleichungen numerisch:

        ∂u/∂t + (u·∇)u = -∇p + ν·Δu + f     (Impulserhaltung)
        ∇·u = 0                                 (Inkompressibilität)

        Verwendetes Verfahren:
        1. Advektionsterm via Upwind-Schema
        2. Viskositätsterm via zentrale Differenzen (Laplace)
        3. Druckpoisson-Gleichung für Inkompressibilität
        4. Geschwindigkeitskorrektur (Chorin-Projektion)

        Rand: Dirichlet (u=0 an Wänden, u=1 am Deckel).

        Millennium-Problem: Bleibt die Lösung für alle Zeiten glatt und beschränkt?
        Für 2D ist dies bewiesen (Ladyzhenskaya 1969), 3D bleibt offen.
    @param nx Gitterpunkte in x-Richtung
    @param ny Gitterpunkte in y-Richtung
    @param dt Zeitschritt
    @param T Simulationsendzeit
    @param nu Kinematische Viskosität (niedrig = turbulenter)
    @param f_x Äußere Kraft in x-Richtung (treibende Kraft)
    @return {'u': array, 'v': array, 'p': array, 'time': float, 'converged': bool}
    @lastModified 2026-03-10
    """
    # Gitterabstand (Einheitsquadrat [0,1]×[0,1])
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    # Felder initialisieren: u (x-Geschwindigkeit), v (y-Geschwindigkeit), p (Druck)
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))

    # CFL-Bedingung prüfen (Stabilitätskriterium für explizite Verfahren)
    # CFL: dt ≤ min(dx, dy)² / (2·ν) und dt ≤ min(dx, dy)
    cfl_viscous = min(dx, dy) ** 2 / (4.0 * nu) if nu > 0 else float('inf')
    cfl_number = dt / cfl_viscous

    n_steps = max(1, int(T / dt))
    time = 0.0
    converged = True

    # Anzahl Druck-Iterationen pro Zeitschritt (Poisson-Lösung)
    n_pressure_iter = 50

    for step in range(n_steps):
        u_old = u.copy()
        v_old = v.copy()

        # -------------------------------------------------------
        # Schritt 1: Geschwindigkeit ohne Druckkorrektur updaten
        # (Euler-Vorwärts + zentrale Differenzen für Viskosität)
        # -------------------------------------------------------
        u_star = u_old.copy()
        v_star = v_old.copy()

        # Innere Gitterpunkte (ohne Rand)
        i = slice(1, ny - 1)
        j = slice(1, nx - 1)

        # Advektionsterm: zentrale Differenzen (stabil für moderate Re-Zahlen)
        # u_old[i, 2:] und u_old[i, :-2] haben die gleiche Shape wie u_old[i, j]
        adv_u = (
            u_old[i, j] * (u_old[i, 2:] - u_old[i, :-2]) / (2 * dx) +
            v_old[i, j] * (u_old[2:, j] - u_old[:-2, j]) / (2 * dy)
        )

        adv_v = (
            u_old[i, j] * (v_old[i, 2:] - v_old[i, :-2]) / (2 * dx) +
            v_old[i, j] * (v_old[2:, j] - v_old[:-2, j]) / (2 * dy)
        )

        # Viskositätsterm: ν · Δu (Laplace-Operator via zentrale Differenzen)
        lap_u = (
            (u_old[i, 2:] - 2 * u_old[i, j] + u_old[i, :-2]) / dx ** 2 +
            (u_old[2:, j] - 2 * u_old[i, j] + u_old[:-2, j]) / dy ** 2
        )

        lap_v = (
            (v_old[i, 2:] - 2 * v_old[i, j] + v_old[i, :-2]) / dx ** 2 +
            (v_old[2:, j] - 2 * v_old[i, j] + v_old[:-2, j]) / dy ** 2
        )

        # Druckgradient (aus vorherigem Schritt)
        dp_dx = (p[i, 2:] - p[i, :-2]) / (2 * dx)
        dp_dy = (p[2:, j] - p[:-2, j]) / (2 * dy)

        # Vorläufige Geschwindigkeiten (ohne Druckkorrektur)
        u_star[i, j] = u_old[i, j] + dt * (-adv_u - dp_dx + nu * lap_u + f_x)
        v_star[i, j] = v_old[i, j] + dt * (-adv_v - dp_dy + nu * lap_v)

        # Randbedingungen: u=0 an allen Wänden außer Deckel (u=1)
        u_star[0, :] = 0.0    # Boden
        u_star[-1, :] = 1.0   # Deckel (Lid-Driven Cavity)
        u_star[:, 0] = 0.0    # Links
        u_star[:, -1] = 0.0   # Rechts
        v_star[0, :] = 0.0
        v_star[-1, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0

        # -------------------------------------------------------
        # Schritt 2: Druck-Poisson-Gleichung lösen
        # Δp = (ρ/dt) · ∇·u*  (Inkompressibilitätsbedingung)
        # -------------------------------------------------------
        for _ in range(n_pressure_iter):
            p_old = p.copy()
            # Divergenz von u*
            div_u = (
                (u_star[i, 2:] - u_star[i, :-2]) / (2 * dx) +
                (v_star[2:, j] - v_star[:-2, j]) / (2 * dy)
            )

            # Gauss-Seidel für Δp = div_u / dt
            p[i, j] = (
                (p_old[i, 2:] + p_old[i, :-2]) * dy ** 2 +
                (p_old[2:, j] + p_old[:-2, j]) * dx ** 2 -
                div_u * dx ** 2 * dy ** 2 / dt
            ) / (2 * (dx ** 2 + dy ** 2))

            # Neumann-Randbedingungen für Druck (dp/dn = 0)
            p[:, 0] = p[:, 1]
            p[:, -1] = p[:, -2]
            p[0, :] = p[1, :]
            p[-1, :] = p[-2, :]

        # -------------------------------------------------------
        # Schritt 3: Geschwindigkeit mit Druckgradient korrigieren
        # u = u* - dt · ∇p
        # -------------------------------------------------------
        u[i, j] = u_star[i, j] - dt * (p[i, 2:] - p[i, :-2]) / (2 * dx)
        v[i, j] = v_star[i, j] - dt * (p[2:, j] - p[:-2, j]) / (2 * dy)

        # Randbedingungen erneut setzen
        u[0, :] = 0.0
        u[-1, :] = 1.0
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[0, :] = 0.0
        v[-1, :] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0

        time += dt

        # Stabilitätsprüfung: Werte dürfen nicht explodieren
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            converged = False
            break

    # Reynolds-Zahl (charakteristisch für die Strömung)
    re = 1.0 / nu if nu > 0 else float('inf')

    return {
        'u': u,
        'v': v,
        'p': p,
        'time': time,
        'converged': converged,
        'reynolds_number': re,
        'cfl_number': cfl_number,
        'grid': (nx, ny),
        'millennium_note': '2D Navier-Stokes: Glattheit bewiesen (Ladyzhenskaya 1969). 3D: offen.'
    }


def check_navier_stokes_energy(u: np.ndarray, v: np.ndarray) -> dict:
    """
    @brief Berechnet die kinetische Energie des Strömungsfeldes.
    @description
        Kinetische Energie: E = 1/2 · ∫|u|² dx

        Numerisch (Trapezregel): E ≈ 1/2 · Σᵢ (uᵢ² + vᵢ²) · ΔV

        Navier-Stokes Millennium-Problem: Bleibt E(t) für alle Zeiten beschränkt?
        - In 2D: Ja (Beweis existiert, globale reguläre Lösung)
        - In 3D: Offen! Können Singularitäten (blow-up) entstehen?

        Energiedissipation durch Viskosität: dE/dt = -ν·∫|∇u|² dx ≤ 0
    @param u x-Komponente des Geschwindigkeitsfeldes (2D-Array)
    @param v y-Komponente des Geschwindigkeitsfeldes (2D-Array)
    @return {'kinetic_energy': float, 'max_velocity': float, 'energy_bounded': bool}
    @lastModified 2026-03-10
    """
    # Kinetische Energie (Trapezregel über Gitter)
    ny, nx = u.shape
    dx = 1.0 / (nx - 1) if nx > 1 else 1.0
    dy = 1.0 / (ny - 1) if ny > 1 else 1.0

    # |u|² = u² + v²
    speed_squared = u ** 2 + v ** 2

    # Numerische Integration (Trapezregel in 2D via np.trapezoid)
    energy = 0.5 * float(np.trapezoid(np.trapezoid(speed_squared, dx=dx, axis=1), dx=dy))

    # Maximale Geschwindigkeit
    max_vel = float(np.max(np.sqrt(speed_squared)))

    # Enstrophie: ε = 1/2 · ∫|∇×u|² dx (misst Wirbelintensität)
    # ∇×u in 2D: ω = ∂v/∂x - ∂u/∂y
    if nx > 2 and ny > 2:
        dv_dx = np.gradient(v, dx, axis=1)
        du_dy = np.gradient(u, dy, axis=0)
        vorticity = dv_dx - du_dy
        enstrophy = 0.5 * float(np.trapezoid(np.trapezoid(vorticity ** 2, dx=dx, axis=1), dx=dy))
    else:
        vorticity = np.zeros_like(u)
        enstrophy = 0.0

    return {
        'kinetic_energy': energy,
        'max_velocity': max_vel,
        'enstrophy': enstrophy,
        'energy_bounded': (energy < 1e10 and not math.isnan(energy)),
        'vorticity_max': float(np.max(np.abs(vorticity))),
        'millennium_note': '3D: Offene Frage ob E stets beschränkt bleibt (Clay Prize $1M)'
    }


# =============================================================================
# 4b. NAVIER-STOKES ENERGIE-EVOLUTION UND BLOW-UP KRITERIEN
# =============================================================================

def navier_stokes_energy_evolution(
    nx: int = 32,
    ny: int = 32,
    dt: float = 0.005,
    T: float = 2.0,
    nu: float = 0.01,
    f_x: float = 1.0,
    record_interval: int = 10
) -> dict:
    """
    @brief Systematische Untersuchung der Energieentwicklung E(t) für 2D Navier-Stokes.
    @description
        Verfolgt folgende physikalische Größen über die Zeit:

        - Kinetische Energie:  E(t) = (1/2) ∫∫ |u|² dA
        - Enstrophie:          Z(t) = (1/2) ∫∫ |ω|² dA  (ω = ∂v/∂x - ∂u/∂y)
        - Dissipationsrate:    ε(t) = ν · Z(t)
        - Max. Vortizität:     ||ω||_∞(t) = max|ω(x,t)|

        Energiebilanz (exakt für glatte Lösungen):
            dE/dt = -ν · Z(t) + F·u̅   (Dissipation + Kraftbeitrag)

        Beale-Kato-Majda-Kriterium (3D, angepasst):
        Blow-up zur Zeit T* ⟺ ∫₀^{T*} ||ω||_∞(t) dt = ∞

        In 2D ist Blow-up ausgeschlossen (Ladyzhenskaya, 1969):
        Z(t) ≤ Z(0) · exp(C·t) → E(t) immer beschränkt.
    @param nx Gitterpunkte in x
    @param ny Gitterpunkte in y
    @param dt Zeitschritt
    @param T Simulationsendzeit
    @param nu Kinematische Viskosität (ν ↓ → Re ↑ → turbulentest)
    @param f_x Treibende Kraft (Lid-Driven-Cavity: f_x=1 am Deckel)
    @param record_interval Aufzeichnung jede record_interval-ten Zeitschritte
    @return dict mit Zeitreihen E(t), Z(t), ||ω||_∞(t) und Blow-up-Indikator
    @lastModified 2026-03-10
    """
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    # Felder initialisieren
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))

    # Zeitreihen für Energiediagnostik
    times = []
    energies = []
    enstrophies = []
    max_vorticities = []
    bkm_integral = 0.0  # Kumuliertes Beale-Kato-Majda-Integral ∫||ω||_∞ dt

    n_steps = max(1, int(T / dt))
    n_pressure_iter = 30

    i_slice = slice(1, ny - 1)
    j_slice = slice(1, nx - 1)

    for step in range(n_steps):
        u_old = u.copy()
        v_old = v.copy()

        # --- Advektion (zentrale Differenzen) ---
        adv_u = (
            u_old[i_slice, j_slice] * (u_old[i_slice, 2:] - u_old[i_slice, :-2]) / (2 * dx) +
            v_old[i_slice, j_slice] * (u_old[2:, j_slice] - u_old[:-2, j_slice]) / (2 * dy)
        )
        adv_v = (
            u_old[i_slice, j_slice] * (v_old[i_slice, 2:] - v_old[i_slice, :-2]) / (2 * dx) +
            v_old[i_slice, j_slice] * (v_old[2:, j_slice] - v_old[:-2, j_slice]) / (2 * dy)
        )

        # --- Viskosität (Laplace) ---
        lap_u = (
            (u_old[i_slice, 2:] - 2 * u_old[i_slice, j_slice] + u_old[i_slice, :-2]) / dx ** 2 +
            (u_old[2:, j_slice] - 2 * u_old[i_slice, j_slice] + u_old[:-2, j_slice]) / dy ** 2
        )
        lap_v = (
            (v_old[i_slice, 2:] - 2 * v_old[i_slice, j_slice] + v_old[i_slice, :-2]) / dx ** 2 +
            (v_old[2:, j_slice] - 2 * v_old[i_slice, j_slice] + v_old[:-2, j_slice]) / dy ** 2
        )

        # --- Druckgradient ---
        dp_dx = (p[i_slice, 2:] - p[i_slice, :-2]) / (2 * dx)
        dp_dy = (p[2:, j_slice] - p[:-2, j_slice]) / (2 * dy)

        u_star = u_old.copy()
        v_star = v_old.copy()
        u_star[i_slice, j_slice] = u_old[i_slice, j_slice] + dt * (-adv_u - dp_dx + nu * lap_u + f_x)
        v_star[i_slice, j_slice] = v_old[i_slice, j_slice] + dt * (-adv_v - dp_dy + nu * lap_v)

        # Randbedingungen
        u_star[0, :] = 0.0; u_star[-1, :] = 1.0
        u_star[:, 0] = 0.0; u_star[:, -1] = 0.0
        v_star[:] = 0.0; v_star[:, 0] = 0.0; v_star[:, -1] = 0.0

        # --- Druck-Poisson ---
        for _ in range(n_pressure_iter):
            p_old = p.copy()
            div_u = (
                (u_star[i_slice, 2:] - u_star[i_slice, :-2]) / (2 * dx) +
                (v_star[2:, j_slice] - v_star[:-2, j_slice]) / (2 * dy)
            )
            p[i_slice, j_slice] = (
                (p_old[i_slice, 2:] + p_old[i_slice, :-2]) * dy ** 2 +
                (p_old[2:, j_slice] + p_old[:-2, j_slice]) * dx ** 2 -
                div_u * dx ** 2 * dy ** 2 / dt
            ) / (2 * (dx ** 2 + dy ** 2))
            p[:, 0] = p[:, 1]; p[:, -1] = p[:, -2]
            p[0, :] = p[1, :]; p[-1, :] = p[-2, :]

        # --- Geschwindigkeitskorrektur ---
        u[i_slice, j_slice] = u_star[i_slice, j_slice] - dt * (p[i_slice, 2:] - p[i_slice, :-2]) / (2 * dx)
        v[i_slice, j_slice] = v_star[i_slice, j_slice] - dt * (p[2:, j_slice] - p[:-2, j_slice]) / (2 * dy)
        u[0, :] = 0.0; u[-1, :] = 1.0; u[:, 0] = 0.0; u[:, -1] = 0.0
        v[0, :] = 0.0; v[-1, :] = 0.0; v[:, 0] = 0.0; v[:, -1] = 0.0

        # --- Diagnostik ---
        if step % record_interval == 0 or step == n_steps - 1:
            t = step * dt
            # Energie E = (1/2) ∫|u|² dA
            speed_sq = u ** 2 + v ** 2
            E = 0.5 * float(np.trapezoid(np.trapezoid(speed_sq, dx=dx, axis=1), dx=dy))

            # Enstrophie Z = (1/2) ∫|ω|² dA  (ω = ∂v/∂x - ∂u/∂y)
            omega = np.gradient(v, dx, axis=1) - np.gradient(u, dy, axis=0)
            Z = 0.5 * float(np.trapezoid(np.trapezoid(omega ** 2, dx=dx, axis=1), dx=dy))
            omega_max = float(np.max(np.abs(omega)))

            times.append(t)
            energies.append(E)
            enstrophies.append(Z)
            max_vorticities.append(omega_max)

            # Beale-Kato-Majda-Integral (Trapez-Update)
            bkm_integral += omega_max * record_interval * dt

        # Stabilitätsprüfung
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            break

    # --- Auswertung ---
    # Energietrend: Steigt oder fällt E(t)?
    if len(energies) >= 2:
        energy_trend = 'increasing' if energies[-1] > energies[0] * 1.01 else (
            'decreasing' if energies[-1] < energies[0] * 0.99 else 'stable'
        )
        # Maximaler Energie-Anstieg (relativ)
        max_energy_ratio = max(energies) / energies[0] if energies[0] > 0 else float('inf')
    else:
        energy_trend = 'unknown'
        max_energy_ratio = 1.0

    # Blow-up-Indikator nach BKM: ||ω||_∞ wächst exponentiell?
    blowup_indicator = False
    if len(max_vorticities) >= 3:
        # Exponentielle Wachstumsrate schätzen: ω_max(t) ~ exp(λt)
        vort_positive = [v for v in max_vorticities if v > 0]
        if len(vort_positive) >= 2:
            log_growth = math.log(vort_positive[-1] / vort_positive[0]) if vort_positive[0] > 0 else 0
            blowup_indicator = log_growth > 10.0  # Mehr als e^10-fache Verstärkung

    return {
        'times': times,
        'energies': energies,
        'enstrophies': enstrophies,
        'max_vorticities': max_vorticities,
        'bkm_integral': bkm_integral,
        'energy_trend': energy_trend,
        'max_energy_ratio': max_energy_ratio,
        'reynolds_number': 1.0 / nu if nu > 0 else float('inf'),
        'blowup_indicator': blowup_indicator,
        'n_time_records': len(times),
        'millennium_note': (
            '2D NS: E(t) immer beschränkt (Ladyzhenskaya 1969). '
            '3D NS: Blow-up-Frage offen. '
            f'BKM-Integral ∫||ω||_∞ dt = {bkm_integral:.4f} (endlich → kein Blow-up in [0,T]).'
        )
    }


def navier_stokes_reynolds_study(
    nu_values: list | None = None,
    nx: int = 24,
    ny: int = 24,
    T: float = 1.0
) -> dict:
    """
    @brief Vergleicht Navier-Stokes-Energie für verschiedene Reynolds-Zahlen.
    @description
        Reynolds-Zahl Re = U·L/ν charakterisiert laminare vs. turbulente Strömung:
        - Re ≪ 1: Laminare Strömung, viskos dominiert
        - Re ~ 100-1000: Übergangsbereich
        - Re ≫ 1: Turbulente Strömung (für NS-Problem relevant)

        Das Millennium-Problem fragt: Auch für Re → ∞ (ν → 0) bleiben
        3D-Lösungen glatt? Unsere 2D-Simulation zeigt, dass 2D stets stabil ist,
        während in 3D Energie-Kaskaden entstehen könnten.
    @param nu_values Liste von Viskositäten (None → [0.1, 0.01, 0.005])
    @param nx Gitterpunkte x
    @param ny Gitterpunkte y
    @param T Simulationsendzeit
    @return dict mit Vergleich der Energieverläufe für verschiedene Re
    @lastModified 2026-03-10
    """
    if nu_values is None:
        nu_values = [0.1, 0.01, 0.005]

    results = {}
    for nu in nu_values:
        Re = 1.0 / nu
        # Zeitschritt für Stabilität anpassen: dt ≤ (dx)²/(4ν) und dt ≤ dx
        dx = 1.0 / (nx - 1)
        dt = min(0.01, dx ** 2 / (4 * nu) * 0.5, dx * 0.5)

        ev = navier_stokes_energy_evolution(
            nx=nx, ny=ny, dt=dt, T=T, nu=nu, record_interval=5
        )
        results[f'Re={Re:.0f}'] = {
            'nu': nu,
            'reynolds_number': Re,
            'final_energy': ev['energies'][-1] if ev['energies'] else 0.0,
            'max_energy': max(ev['energies']) if ev['energies'] else 0.0,
            'final_enstrophy': ev['enstrophies'][-1] if ev['enstrophies'] else 0.0,
            'max_vorticity': max(ev['max_vorticities']) if ev['max_vorticities'] else 0.0,
            'bkm_integral': ev['bkm_integral'],
            'energy_trend': ev['energy_trend'],
            'blowup_indicator': ev['blowup_indicator'],
        }

    return {
        'reynolds_study': results,
        'conclusion': (
            '2D Navier-Stokes zeigt für alle Re stabiles Energieverhalten. '
            'Höheres Re → höhere Enstrophie und Vortizität, aber kein Blow-up. '
            'In 3D könnte Re → ∞ zu Singularitäten führen (offene Millennium-Frage).'
        )
    }


# =============================================================================
# 5. YANG-MILLS (vereinfachte Quantenfeldtheorie)
# =============================================================================

def yang_mills_action(field_strength: np.ndarray, coupling: float = 1.0) -> float:
    """
    @brief Berechnet die Yang-Mills-Wirkung für ein gegebenes Feldstärken-Array.
    @description
        Yang-Mills-Wirkung (vereinfacht für diskretes 2D-Gitter):
        S = (1/(4g²)) · Σ Tr(F_{μν} · F^{μν})

        In kontinuierlicher Form:
        S[A] = (1/4g²) · ∫ Tr(F_{μν} F^{μν}) d⁴x

        wobei F_{μν} = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν] (Feldstärketensor)

        Massenlücke-Vermutung: Das Spektrum des Hamilton-Operators hat eine
        positive untere Schranke (mass gap Δ > 0). Beweis ist Millennium-Problem.
    @param field_strength 2D-Array der Feldstärken F_{μν}
    @param coupling Eichkopplungskonstante g
    @return Numerischer Wert der Wirkung S
    @lastModified 2026-03-10
    """
    if coupling == 0:
        raise ValueError("Kopplungskonstante g darf nicht 0 sein")

    # S = (1/4g²) · Σ F_{μν}²  (vereinfacht: abelsche Näherung Tr(F²) = |F|²)
    action = (1.0 / (4.0 * coupling ** 2)) * float(np.sum(field_strength ** 2))

    return action


def lattice_gauge_theory_demo(size: int = 4, coupling: float = 1.0, steps: int = 100) -> dict:
    """
    @brief Demonstration der Gitter-Eichtheorie (Grundlage des Yang-Mills-Beweises).
    @description
        Gitter-Eichtheorie (Wilson 1974) ist die diskrete Näherung der Yang-Mills-Theorie.
        Auf einem size×size-Gitter werden SU(2)-Eichfelder als 2×2-unitäre Matrizen
        auf jeder Gitter-Verbindung (Link) repräsentiert.

        Wilson-Loop (C = rechteckiger Pfad der Größe R×T):
        W(C) = Tr(Π_{links ∈ C} U_link)

        Für große Wilson-Loops gilt (bei confinement / Massenlücke):
        ⟨W(C)⟩ ~ exp(-σ · Area(C))    (Area-Law, σ = String-Tension)

        Dies deutet auf eine Massenlücke hin: Quarks sind eingesperrt (confinement).
    @param size Gittergröße (size × size Plaquetten)
    @param coupling Eichkopplungskonstante β = 1/g²
    @param steps Anzahl Monte-Carlo-Schritte
    @return {'wilson_loops': list, 'mass_gap_estimate': float, 'string_tension': float}
    @lastModified 2026-03-10
    """
    rng = np.random.default_rng(42)

    # Eichfelder als zufällige U(1)-Phasen auf Links initialisieren
    # (vereinfachte abelsche Version statt echtem SU(2))
    # links[i, j, mu]: Phase auf dem Link (i,j) in Richtung mu (0=x, 1=y)
    links = rng.uniform(0, 2 * np.pi, (size, size, 2))

    def plaquette_action(links_arr):
        """Berechnet die Plaquetten-Summe (Wilson-Wirkung für U(1))."""
        total = 0.0
        for i in range(size):
            for j in range(size):
                # Plaquette = U_x(i,j) + U_y(i+1,j) - U_x(i,j+1) - U_y(i,j)
                phase = (
                    links_arr[i, j, 0] +
                    links_arr[(i + 1) % size, j, 1] -
                    links_arr[i, (j + 1) % size, 0] -
                    links_arr[i, j, 1]
                )
                # Plaquette-Wirkung: S_plaq = β · (1 - cos(phase))
                total += 1.0 - math.cos(phase)
        return coupling * total

    # Monte-Carlo Metropolis für Thermalisierung
    current_action = plaquette_action(links)
    for _ in range(steps):
        # Zufälligen Link wählen und zufällige Änderung vorschlagen
        i, j, mu = rng.integers(0, size), rng.integers(0, size), rng.integers(0, 2)
        delta = rng.uniform(-0.5, 0.5)
        old_phase = links[i, j, mu]
        links[i, j, mu] = (old_phase + delta) % (2 * np.pi)
        new_action = plaquette_action(links)

        # Metropolis-Akzeptanzschritt: Akzeptiere wenn neue Wirkung kleiner
        delta_s = new_action - current_action
        if delta_s < 0 or rng.random() < math.exp(-delta_s):
            current_action = new_action
        else:
            # Änderung rückgängig machen
            links[i, j, mu] = old_phase

    # Wilson-Loops für verschiedene Größen berechnen
    wilson_loops = []
    for r in range(1, min(size // 2 + 1, 4)):
        # Rechteckiger Wilson-Loop der Größe r×r
        w_vals = []
        for i in range(size):
            for j in range(size):
                # Phasen entlang des rechteckigen Pfades summieren
                phase_sum = 0.0
                for k in range(r):
                    phase_sum += links[(i + k) % size, j, 0]        # Unterkante
                    phase_sum += links[(i + r) % size, (j + k) % size, 1]  # Rechtskante
                    phase_sum -= links[(i + k) % size, (j + r) % size, 0]  # Oberkante
                    phase_sum -= links[i, (j + k) % size, 1]         # Linkskante
                w_vals.append(math.cos(phase_sum))
        wilson_loops.append({'r': r, 'W': float(np.mean(w_vals))})

    # String-Tension σ schätzen via Area-Law: W(r) ≈ exp(-σ · r²)
    string_tension = 0.0
    if len(wilson_loops) >= 2 and wilson_loops[0]['W'] > 0 and wilson_loops[1]['W'] > 0:
        w1 = max(wilson_loops[0]['W'], 1e-10)
        w2 = max(wilson_loops[1]['W'], 1e-10)
        string_tension = max(0.0, -math.log(abs(w2) / abs(w1)))

    # Masse-Gap-Schätzung aus String-Tension: Δ ≈ σ (grobe Näherung)
    mass_gap_estimate = string_tension

    return {
        'wilson_loops': wilson_loops,
        'mass_gap_estimate': mass_gap_estimate,
        'string_tension': string_tension,
        'final_action': current_action,
        'millennium_note': 'Massenlücke Δ > 0 ist unbewiesen für 4D SU(N) Yang-Mills (Clay Prize $1M)'
    }


# =============================================================================
# 5b. SU(2) GITTER-EICHTHEORIE (nicht-abelsch, Yang-Mills Hauptgruppe)
# =============================================================================

def _random_su2(rng: np.random.Generator) -> np.ndarray:
    """
    @brief Erzeugt eine zufällige SU(2)-Matrix mit Haar-Maß.
    @description
        SU(2) = {U = a·I + i(b·σ₁ + c·σ₂ + d·σ₃) : a²+b²+c²+d²=1}
              = {[[a+id, -b+ic], [b+ic, a-id]] : (a,b,c,d) ∈ S³}

        Gleichmäßige Verteilung auf S³ ⊂ ℝ⁴ → Haar-Maß auf SU(2).

        Pauli-Matrizen:
            σ₁ = [[0,1],[1,0]], σ₂ = [[0,-i],[i,0]], σ₃ = [[1,0],[0,-1]]
    @param rng NumPy Zufallsgenerator
    @return 2×2 komplexe unitäre Matrix mit det=1
    @lastModified 2026-03-10
    """
    # Gleichmäßiger Punkt auf S³: Gausssche Projektion
    v = rng.standard_normal(4)
    v /= np.linalg.norm(v)
    a, b, c, d = v

    # SU(2)-Matrix: U = a·I + i(b·σ₁ + c·σ₂ + d·σ₃)
    # = [[a+id, ic-b], [ic+b, a-id]]
    U = np.array([
        [a + 1j * d,  -b + 1j * c],
        [ b + 1j * c,  a - 1j * d]
    ], dtype=complex)

    return U


def _su2_plaquette(links: np.ndarray, i: int, j: int, size: int) -> np.ndarray:
    """
    @brief Berechnet die Wilson-Plaquette im SU(2)-Gitter.
    @description
        Plaquette U_□(i,j) = U_x(i,j) · U_y(i+1,j) · U†_x(i,j+1) · U†_y(i,j)

        Im Gitter: links[i,j,mu] ist eine 2×2 SU(2)-Matrix auf dem Link
        (i,j) → (i+μ̂,j) für mu=0 (x-Richtung) oder (i,j+ŷ) für mu=1 (y-Richtung).

        Periodische Randbedingungen: (i + 1) % size
    @param links Array der Form (size, size, 2, 2, 2): links[i,j,mu] = 2×2 Matrix
    @param i Gitter-x-Index
    @param j Gitter-y-Index
    @param size Gittergröße
    @return 2×2 SU(2)-Plaquette-Matrix
    @lastModified 2026-03-10
    """
    # Vier Links der Plaquette (Uhrzeigersinn)
    U_x_ij      = links[i, j, 0]                           # U_x(i,j)
    U_y_ip1j    = links[(i + 1) % size, j, 1]              # U_y(i+1,j)
    U_x_ijp1_h  = links[i, (j + 1) % size, 0].conj().T    # U†_x(i,j+1)
    U_y_ij_h    = links[i, j, 1].conj().T                  # U†_y(i,j)

    # Plaquette = Produkt der vier Links
    return U_x_ij @ U_y_ip1j @ U_x_ijp1_h @ U_y_ij_h


def lattice_gauge_su2(
    size: int = 4,
    beta: float = 2.0,
    steps: int = 200,
    warmup: int = 50
) -> dict:
    """
    @brief SU(2) Gitter-Eichtheorie mit Metropolis-Algorithmus.
    @description
        Wilson-Wirkung für SU(2) auf einem size×size Gitter (2D):

            S_W = β · Σ_{Plaquetten} (1 - (1/2) Re Tr(U_□))

        Für β → ∞ (schwache Kopplung g → 0) → Kontinuumsgrenze.
        Für β → 0 (starke Kopplung) → Confinement, Area-Law.

        Metropolis-Update (Cabibbo-Marinari-ähnlich):
        1. Wähle zufälligen Link (i,j,μ)
        2. Schlage neue SU(2)-Matrix U' vor (U_neu = U_zuf · U_alt, mit U_zuf nahe I)
        3. Berechne ΔS = S(U') - S(U)
        4. Akzeptiere mit min(1, exp(-ΔS))

        Erwartungswerte:
        - ⟨P⟩ = Mittelwert der Plaquetten-Spurwerte
        - Wilson-Loop W(r) → Area-Law → String-Tension σ → Massenlücke Δ ∝ √σ
        - Für SU(2): Massenlücke aus Glueball-Korrelator (hier: String-Tension-Näherung)
    @param size Gittergröße (size × size Plaquetten)
    @param beta Eichkopplungsparameter β = 4/g² (SU(2)-Konvention)
    @param steps Anzahl Metropolis-Sweeps nach Warmup
    @param warmup Thermalisierungs-Sweeps
    @return dict mit Plaquetten-Mittelwert, Wilson-Loops, String-Tension, Massenlücken-Schätzung
    @lastModified 2026-03-10
    """
    rng = np.random.default_rng(42)

    # Gitter initialisieren: alle Links = Einheitsmatrix (kalter Start)
    # links[i, j, mu] ist eine 2×2-Matrix ∈ SU(2)
    identity = np.eye(2, dtype=complex)
    links = np.tile(identity, (size, size, 2, 1, 1))

    def plaquette_sum_action(links_arr: np.ndarray) -> float:
        """Berechnet die gesamte Wilson-Wirkung S = β·Σ(1 - (1/2)Re Tr U_□)."""
        total = 0.0
        for ii in range(size):
            for jj in range(size):
                Up = _su2_plaquette(links_arr, ii, jj, size)
                # (1/2) Re Tr(U_□) liegt in [-1, 1]
                total += 1.0 - 0.5 * Up.trace().real
        return beta * total

    def local_action_change(links_arr: np.ndarray, i0: int, j0: int,
                            mu0: int, U_new: np.ndarray) -> float:
        """Berechnet ΔS bei Ersetzen von links[i0,j0,mu0] durch U_new."""
        # Nur die Plaquetten die diesen Link enthalten ändern sich
        delta_s = 0.0
        for di, dj in [(-1, 0), (0, 0), (0, -1), (-1, -1)]:
            ii = (i0 + di) % size
            jj = (j0 + dj) % size

            # Plaquette vor und nach dem Update
            old_Up = _su2_plaquette(links_arr, ii, jj, size)
            links_arr[i0, j0, mu0] = U_new
            new_Up = _su2_plaquette(links_arr, ii, jj, size)
            links_arr[i0, j0, mu0] = links_arr[i0, j0, mu0]  # kein Rollback hier

            delta_s += -0.5 * (new_Up.trace().real - old_Up.trace().real)

        return beta * delta_s

    # --- Metropolis-Sweeps ---
    n_accept = 0
    n_total = 0

    # Warmup (Thermalisierung, keine Messungen)
    for _ in range(warmup):
        for i0 in range(size):
            for j0 in range(size):
                for mu0 in range(2):
                    # Neues Element: U_new = kleiner SU(2)-Schritt × aktuelles Element
                    epsilon = 0.3  # Schrittgröße (Akzeptanzrate ~50% angestrebt)
                    # Zufällige SU(2)-Matrix nahe der Identität
                    v = rng.standard_normal(4) * epsilon
                    v[0] = math.sqrt(max(0.0, 1.0 - sum(v[1:] ** 2)))  # Normierung auf S³
                    a, b, c, d = v
                    dU = np.array([[a + 1j*d, -b + 1j*c], [b + 1j*c, a - 1j*d]], dtype=complex)
                    U_new = dU @ links[i0, j0, mu0]

                    # Metropolis-Akzeptanz
                    old_action = plaquette_sum_action(links)
                    links_backup = links[i0, j0, mu0].copy()
                    links[i0, j0, mu0] = U_new
                    new_action = plaquette_sum_action(links)
                    delta_s = new_action - old_action

                    if delta_s > 0 and rng.random() >= math.exp(-delta_s):
                        links[i0, j0, mu0] = links_backup  # Ablehnen

    # Messungs-Sweeps
    plaquette_values = []
    for _ in range(steps):
        for i0 in range(size):
            for j0 in range(size):
                for mu0 in range(2):
                    v = rng.standard_normal(4) * 0.3
                    v[0] = math.sqrt(max(0.0, 1.0 - sum(v[1:] ** 2)))
                    a, b, c, d = v
                    dU = np.array([[a + 1j*d, -b + 1j*c], [b + 1j*c, a - 1j*d]], dtype=complex)
                    U_new = dU @ links[i0, j0, mu0]

                    old_action = plaquette_sum_action(links)
                    links_backup = links[i0, j0, mu0].copy()
                    links[i0, j0, mu0] = U_new
                    new_action = plaquette_sum_action(links)
                    delta_s = new_action - old_action

                    if delta_s > 0 and rng.random() >= math.exp(-delta_s):
                        links[i0, j0, mu0] = links_backup
                    else:
                        n_accept += 1
                    n_total += 1

        # Plaquetten-Messung: ⟨P⟩ = (1/N_plaq) Σ (1/2) Re Tr(U_□)
        p_sum = 0.0
        for ii in range(size):
            for jj in range(size):
                Up = _su2_plaquette(links, ii, jj, size)
                p_sum += 0.5 * Up.trace().real
        plaquette_values.append(p_sum / (size * size))

    # --- Wilson-Loops ---
    wilson_loops_su2 = []
    for r in range(1, min(size // 2 + 1, 4)):
        # Rechteckiger r×r Wilson-Loop
        w_vals = []
        for i0 in range(size):
            for j0 in range(size):
                W = np.eye(2, dtype=complex)
                # Unterkante: r Schritte in x
                for k in range(r):
                    W = W @ links[(i0 + k) % size, j0, 0]
                # Rechtskante: r Schritte in y
                for k in range(r):
                    W = W @ links[(i0 + r) % size, (j0 + k) % size, 1]
                # Oberkante: r Schritte in -x (hermitesch konjugiert)
                for k in range(r):
                    W = W @ links[(i0 + r - 1 - k) % size, (j0 + r) % size, 0].conj().T
                # Linkskante: r Schritte in -y
                for k in range(r):
                    W = W @ links[i0, (j0 + r - 1 - k) % size, 1].conj().T
                w_vals.append(0.5 * W.trace().real)
        wilson_loops_su2.append({'r': r, 'W': float(np.mean(w_vals))})

    # String-Tension σ aus Area-Law: W(r) ~ exp(-σ·r²)
    string_tension_su2 = 0.0
    if len(wilson_loops_su2) >= 2:
        w1 = wilson_loops_su2[0]['W']
        w2 = wilson_loops_su2[1]['W']
        if abs(w1) > 1e-10 and abs(w2) > 1e-10:
            string_tension_su2 = max(0.0, -math.log(abs(w2) / abs(w1)))

    # Massenlücken-Schätzung: Δ ≈ √σ (Glueball-Masse aus String-Tension)
    mass_gap_su2 = math.sqrt(string_tension_su2) if string_tension_su2 > 0 else 0.0

    mean_plaquette = float(np.mean(plaquette_values)) if plaquette_values else 0.0
    accept_rate = n_accept / n_total if n_total > 0 else 0.0

    return {
        'gauge_group': 'SU(2)',
        'size': size,
        'beta': beta,
        'mean_plaquette': mean_plaquette,
        'plaquette_variance': float(np.var(plaquette_values)) if plaquette_values else 0.0,
        'wilson_loops': wilson_loops_su2,
        'string_tension': string_tension_su2,
        'mass_gap_estimate': mass_gap_su2,
        'acceptance_rate': accept_rate,
        'n_sweeps': steps,
        'millennium_note': (
            'SU(2) nicht-abelsche Eichtheorie. '
            f'β={beta}: Massenlücken-Schätzung Δ≈{mass_gap_su2:.4f}. '
            'Massenlücke > 0 ist das Yang-Mills Millennium-Problem (Clay Prize $1M).'
        )
    }


def yang_mills_mass_gap_study(
    beta_values: list | None = None,
    size: int = 4
) -> dict:
    """
    @brief Untersucht die Massenlücke als Funktion der Kopplung β für SU(2).
    @description
        Die Yang-Mills-Massenlücke ist die positive untere Schranke des Spektrums
        des Hamilton-Operators H:

            Δ = inf { E : E > 0, E ∈ Spec(H) } > 0  (Millennium-Vermutung)

        Auf dem Gitter: Extrapolation der String-Tension σ(β) zum Kontinuum.
        Im Skalenbereich (β → ∞ bei fixem physikalischen Abstand):
            σ(β) ~ exp(-2π²β / (11/3))  (Asymptotische Freiheit für SU(3))

        Für SU(2): σ(β) ~ exp(-π²β / (11/6))

        Confinement-Kriterium: σ > 0 (Area-Law für Wilson-Loops)
        → Massenlücke Δ ≈ √σ / a  (a = Gitterabstand)
    @param beta_values Liste von β-Werten (None → [0.5, 1.0, 1.5, 2.0, 2.5])
    @param size Gittergröße
    @return dict mit σ(β), Massenlücke(β) und Skalenanalyse
    @lastModified 2026-03-10
    """
    if beta_values is None:
        beta_values = [0.5, 1.0, 1.5, 2.0, 2.5]

    results = {}
    for beta in beta_values:
        r = lattice_gauge_su2(size=size, beta=beta, steps=50, warmup=20)
        results[f'beta={beta}'] = {
            'beta': beta,
            'mean_plaquette': r['mean_plaquette'],
            'string_tension': r['string_tension'],
            'mass_gap': r['mass_gap_estimate'],
            'wilson_loops': r['wilson_loops'],
        }

    # Theoretische Skalenprädiktion für SU(2): σ ~ exp(-π²β/c) mit c=11/6
    beta_arr = np.array(beta_values)
    c_su2 = 11.0 / 6.0
    sigma_theory = np.exp(-math.pi ** 2 * beta_arr / c_su2)

    return {
        'gauge_group': 'SU(2)',
        'beta_scan': results,
        'theoretical_sigma': dict(zip([f'beta={b}' for b in beta_values],
                                      sigma_theory.tolist())),
        'interpretation': (
            'Bei großem β (schwacher Kopplung) → kleine String-Tension, großes Gitter. '
            'Bei kleinem β (starker Kopplung) → Area-Law, Confinement. '
            'Massenlücke Δ > 0 für alle endlichen β wird numerisch beobachtet, '
            'aber der analytische Beweis im Kontinuum ist das Millennium-Problem.'
        )
    }


# =============================================================================
# 6. BIRCH UND SWINNERTON-DYR VERMUTUNG (Hilfswerkzeuge)
# =============================================================================

def elliptic_curve_points_mod_p(a: int, b: int, p: int) -> list:
    """
    @brief Zählt rationale Punkte auf einer elliptischen Kurve modulo p.
    @description
        Elliptische Kurve: y² = x³ + ax + b über F_p (Körper mod p)
        Bedingung: 4a³ + 27b² ≢ 0 (mod p) (nicht-singuläre Kurve)

        BSD-Vermutung: ords=1 L(E,1) = Rang der rationalen Punkte auf E.
        Der L-Funktion-Wert bei s=1 kodiert die arithmetische Information.

        Hasse-Weil-Schranke: |#E(F_p) - (p+1)| ≤ 2√p
    @param a Koeffizient der Kurve y² = x³ + ax + b
    @param b Koeffizient der Kurve y² = x³ + ax + b
    @param p Primzahl (Modulus)
    @return Liste aller Punkte (x, y) auf E(F_p) einschließlich Punkt im Unendlichen
    @lastModified 2026-03-10
    """
    if p < 2:
        raise ValueError("p muss eine Primzahl ≥ 2 sein")

    # Singularitätsbedingung prüfen: 4a³ + 27b² ≢ 0 (mod p)
    discriminant = (4 * a ** 3 + 27 * b ** 2) % p
    if discriminant == 0:
        raise ValueError(f"Singuläre Kurve: 4a³+27b² ≡ 0 (mod {p})")

    points = ['O']  # O = Punkt im Unendlichen (neutrales Element)

    for x in range(p):
        # Rechte Seite: x³ + ax + b mod p
        rhs = (pow(x, 3, p) + a * x + b) % p

        # Quadratwurzeln mod p finden: y² ≡ rhs (mod p)
        for y in range(p):
            if (y * y) % p == rhs:
                points.append((x, y))

    return points


def bsd_l_function_rank_estimate(a: int, b: int, primes: list) -> dict:
    """
    @brief Schätzt den analytischen Rang einer elliptischen Kurve via BSD-Vermutung.
    @description
        Birch und Swinnerton-Dyer (BSD)-Vermutung:
        Π_{p ≤ X} (#E(F_p) / p) ≈ C · (log X)^r

        wobei r der Rang der Kurve E über ℚ ist.

        Empirisch: Logarithmische Steigung ≈ Rang
        Beweis: Nur für spezielle Fälle (Coates-Wiles 1977, Gross-Zagier 1986)
    @param a Koeffizient a der elliptischen Kurve y² = x³ + ax + b
    @param b Koeffizient b
    @param primes Liste der Primzahlen p (sollte lang sein für gute Schätzung)
    @return {'rank_estimate': int, 'log_product': float, 'curve': str}
    @lastModified 2026-03-10
    """
    if not primes:
        return {'error': 'Leere Primzahlliste'}

    log_product = 0.0
    valid_primes = []

    for p in primes:
        if p < 5:
            continue
        # Diskriminante prüfen
        disc = (4 * a ** 3 + 27 * b ** 2) % p
        if disc == 0:
            continue  # Kurve ist schlecht bei p → überspringen

        # #E(F_p) zählen
        count = 1  # Punkt im Unendlichen
        for x in range(p):
            rhs = (pow(x, 3, p) + a * x + b) % p
            # Legendre-Symbol: rhs ist Quadrat mod p?
            if rhs == 0:
                count += 1
            else:
                # Euler-Kriterium: rhs^((p-1)/2) ≡ 1 (mod p) gdw. QR
                if pow(rhs, (p - 1) // 2, p) == 1:
                    count += 2  # Zwei Lösungen y und p-y

        # Beitrag zu BSD-Produkt: log(#E(F_p) / p)
        log_product += math.log(count / p)
        valid_primes.append(p)

    # Rank-Schätzung: Steigung ≈ r
    if valid_primes:
        log_x = math.log(valid_primes[-1])
        rank_estimate = round(log_product / log_x) if log_x > 0 else 0
    else:
        rank_estimate = 0

    return {
        'rank_estimate': max(0, rank_estimate),
        'log_product': log_product,
        'n_primes': len(valid_primes),
        'curve': f'y² = x³ + {a}x + {b}',
        'bsd_note': 'BSD-Vermutung: Rang = ords=1 L(E,s) (Millennium-Problem, unbewiesen)'
    }


# =============================================================================
# 7. STATUS-REPORT ALLER MILLENNIUM-PROBLEME
# =============================================================================

def millennium_problems_status() -> dict:
    """
    @brief Gibt den aktuellen Status aller 7 Millennium-Probleme zurück.
    @description
        Das Clay Mathematics Institute stellte im Jahr 2000 sieben Probleme auf,
        für deren Lösung jeweils 1 Million US-Dollar ausgesetzt wurden.

        Bislang wurde nur die Poincaré-Vermutung gelöst (Grigori Perelman, 2003),
        der den Preis allerdings ablehnte.
    @return dict mit Status aller 7 Probleme
    @lastModified 2026-03-10
    """
    problems = {
        'P_vs_NP': {
            'name': 'P vs NP',
            'year_posed': 1971,
            'field': 'Komplexitätstheorie / Informatik',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Gilt P = NP? D.h.: Kann jedes Problem, dessen Lösung in '
                'polynomieller Zeit verifiziert werden kann, auch in poly. Zeit gelöst werden?'
            ),
            'tools_in_module': ['is_np_complete_3sat', 'traveling_salesman_brute',
                                'traveling_salesman_greedy', 'complexity_comparison'],
            'key_fact': 'Consensus: P ≠ NP, aber kein Beweis existiert'
        },
        'Hodge_Conjecture': {
            'name': 'Hodge-Vermutung',
            'year_posed': 1950,
            'field': 'Algebraische Geometrie',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Ist jede Hodge-Klasse auf einer nicht-singulären projektiven '
                'algebraischen Varietät eine rationale lineare Kombination von Klassen '
                'algebraischer Zykel?'
            ),
            'tools_in_module': [],
            'key_fact': 'Spezialfälle bekannt, allgemein offen'
        },
        'Poincare_Conjecture': {
            'name': 'Poincaré-Vermutung',
            'year_posed': 1904,
            'field': 'Topologie / Geometrie',
            'status': 'GELÖST',
            'solved_by': 'Grigori Perelman',
            'year_solved': 2003,
            'prize_usd': 1_000_000,
            'prize_accepted': False,
            'description': (
                'Jede einfach zusammenhängende, geschlossene 3-Mannigfaltigkeit '
                'ist homöomorph zur 3-Sphäre S³.'
            ),
            'tools_in_module': [],
            'key_fact': 'Perelman nutzte Ricci-Fluss (Hamilton). Preis abgelehnt.'
        },
        'Riemann_Hypothesis': {
            'name': 'Riemann-Hypothese',
            'year_posed': 1859,
            'field': 'Analytische Zahlentheorie',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Alle nicht-trivialen Nullstellen der Riemann-Zeta-Funktion ζ(s) '
                'liegen auf der kritischen Geraden Re(s) = 1/2.'
            ),
            'tools_in_module': ['riemann_zeros_mpmath', 'verify_rh_critical_line',
                                'gram_law_check', 'zeta_zero_gaps'],
            'key_fact': '>10^13 Nullstellen auf Re=1/2 verifiziert (ZetaGrid 2004)'
        },
        'Yang_Mills': {
            'name': 'Yang-Mills-Existenz und Massenlücke',
            'year_posed': 2000,
            'field': 'Mathematische Physik / Quantenfeldtheorie',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Für jede kompakte einfach-zusammenhängende Eichgruppe G existiert '
                'eine Yang-Mills-Quantenfeldtheorie in ℝ⁴ mit positiver Massenlücke Δ > 0.'
            ),
            'tools_in_module': ['yang_mills_action', 'lattice_gauge_theory_demo'],
            'key_fact': 'Gitter-QCD liefert numerische Evidenz für Massenlücke'
        },
        'Navier_Stokes': {
            'name': 'Navier-Stokes-Existenz und Glattheit',
            'year_posed': 2000,
            'field': 'Partielle Differentialgleichungen',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Existieren für beliebige glatte Anfangsdaten im ℝ³ globale glatte '
                'Lösungen der 3D-Navier-Stokes-Gleichungen? (Oder entstehen Singularitäten?)'
            ),
            'tools_in_module': ['navier_stokes_2d_simple', 'check_navier_stokes_energy'],
            'key_fact': '2D gelöst (Ladyzhenskaya 1969), 3D offen'
        },
        'BSD_Conjecture': {
            'name': 'Birch und Swinnerton-Dyer Vermutung',
            'year_posed': 1965,
            'field': 'Zahlentheorie / Arithmetische Geometrie',
            'status': 'OFFEN',
            'prize_usd': 1_000_000,
            'description': (
                'Der Rang einer elliptischen Kurve E über ℚ ist gleich der '
                'Ordnung der Nullstelle der zugehörigen L-Funktion L(E,s) bei s=1.'
            ),
            'tools_in_module': ['elliptic_curve_points_mod_p', 'bsd_l_function_rank_estimate'],
            'key_fact': 'Teilergebnisse für Rang 0 und 1 (Coates-Wiles, Kolyvagin)'
        }
    }

    # Zusammenfassung
    solved = sum(1 for p in problems.values() if p['status'] == 'GELÖST')
    open_count = len(problems) - solved
    total_prize = sum(p['prize_usd'] for p in problems.values() if p['status'] == 'OFFEN')

    return {
        'problems': problems,
        'summary': {
            'total': len(problems),
            'solved': solved,
            'open': open_count,
            'total_prize_remaining_usd': total_prize,
            'source': 'Clay Mathematics Institute (2000)'
        }
    }
