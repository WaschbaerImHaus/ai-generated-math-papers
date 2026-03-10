"""
@file information_theory.py
@brief Informationstheorie nach Shannon – Entropie, Kanalkapazität, Codierung, Fehlerkorrektur.
@description
    Dieses Modul implementiert die fundamentalen Konzepte der Informationstheorie,
    die von Claude Shannon in seinem bahnbrechenden Werk
    "A Mathematical Theory of Communication" (1948) begründet wurden.

    Enthaltene Klassen:
    - ShannonEntropy        : Shannon-Entropie und verwandte Größen
    - KLDivergence          : Kullback-Leibler-Divergenz und Informationsmaße
    - ChannelCapacity       : Kanalkapazität, Blahut-Arimoto-Algorithmus
    - SourceCoding          : Quellencodierung (Huffman, Arithmetik, LZ)
    - ErrorCorrection       : Fehlerkorrektur-Grundlagen (Hamming, Singleton)
    - DifferentialEntropy   : Differentielle Entropie für kontinuierliche Zufallsvariablen

    Standalone-Funktionen:
    - surprisal()                         : Informationsgehalt eines Ereignisses
    - entropy_rate_markov()               : Entropierate einer Markov-Kette
    - data_processing_inequality_demo()   : Demonstration der Datenverarbeitungs-Ungleichung

    Mathematische Grundlagen:
        Shannon-Entropie (diskret):
            H(X) = -Σ p_i · log₂(p_i)  [in Bit]
            H(X) = -Σ p_i · ln(p_i)    [in Nat]

        Maximale Entropie (Gleichverteilung):
            H_max = log₂(n) bei n gleichwahrscheinlichen Ereignissen

        Kanalkapazität (Shannon):
            C = max_{p(x)} I(X;Y)
            C_BSC = 1 - H(p)        (Binärer symmetrischer Kanal)
            C_AWGN = log₂(1 + SNR)  (Additiver weißer Gaußscher Rauschkanal)

        Quellencodierungssatz (Shannon):
            H(X) ≤ L̄ < H(X) + 1   (L̄ = mittlere Codelänge in Bit)

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import math
import heapq
import numpy as np
from typing import Optional


# ===========================================================================
# KLASSE: SHANNON-ENTROPIE
# ===========================================================================

class ShannonEntropy:
    """
    Klasse für Shannon-Entropie-Berechnungen.

    Die Shannon-Entropie misst die durchschnittliche Menge an Information
    (Überraschung) einer diskreten Zufallsvariable.

    Definition:
        H(X) = -Σ_{i} p_i · log_b(p_i)

    wobei b die Basis des Logarithmus ist:
        b = 2   → Einheit: Bit (gebräuchlichste Form)
        b = e   → Einheit: Nat
        b = 10  → Einheit: Hartley (Digit)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def entropy(probabilities: list[float], base: float = 2) -> float:
        """
        Berechnet die Shannon-Entropie H(X) einer diskreten Verteilung.

        Formel:
            H(X) = -Σ p_i · log_b(p_i)

        Konvention: 0 · log(0) := 0 (Grenzwert für p → 0).

        @param probabilities: Liste von Wahrscheinlichkeiten p_i (müssen sich zu 1 summieren)
        @param base: Logarithmusbasis (Standard: 2 → Bit)
        @return: Entropie H(X) in der gewählten Einheit
        @raises ValueError: wenn Wahrscheinlichkeiten negativ oder Summe ≠ 1
        @lastModified 2026-03-10
        """
        probs = np.asarray(probabilities, dtype=float)

        # Negative Wahrscheinlichkeiten sind physikalisch sinnlos
        if np.any(probs < 0):
            raise ValueError("Wahrscheinlichkeiten dürfen nicht negativ sein.")

        # Normalisierungsprüfung mit Toleranz für Rundungsfehler
        total = np.sum(probs)
        if not math.isclose(total, 1.0, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"Wahrscheinlichkeiten müssen sich zu 1 summieren, erhalten: {total}"
            )

        # Entropieberechnung: p=0 wird ausgelassen (0·log(0) = 0 per Konvention)
        h = 0.0
        for p in probs:
            if p > 0:
                h -= p * math.log(p, base)
        return h

    @staticmethod
    def joint_entropy(joint_probs: list[list[float]], base: float = 2) -> float:
        """
        Berechnet die gemeinsame Entropie H(X,Y) für zwei Zufallsvariablen.

        Formel:
            H(X,Y) = -Σ_{x,y} p(x,y) · log p(x,y)

        @param joint_probs: 2D-Wahrscheinlichkeitsmatrix P[i][j] = P(X=i, Y=j)
        @param base: Logarithmusbasis (Standard: 2)
        @return: gemeinsame Entropie H(X,Y)
        @lastModified 2026-03-10
        """
        matrix = np.asarray(joint_probs, dtype=float)

        # Normalisierungsprüfung der gesamten Matrix
        total = np.sum(matrix)
        if not math.isclose(total, 1.0, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"Verbundwahrscheinlichkeiten müssen sich zu 1 summieren, erhalten: {total}"
            )

        # Summiere -p·log(p) über alle Einträge (flache Iteration)
        h = 0.0
        for p in matrix.flatten():
            if p > 0:
                h -= p * math.log(p, base)
        return h

    @staticmethod
    def conditional_entropy(joint_probs: list[list[float]], base: float = 2) -> float:
        """
        Berechnet die bedingte Entropie H(Y|X).

        Formel:
            H(Y|X) = H(X,Y) - H(X)
                   = -Σ_{x,y} p(x,y) · log p(y|x)

        Die bedingte Entropie misst die verbleibende Unsicherheit über Y,
        wenn X bekannt ist.

        @param joint_probs: 2D-Verbundwahrscheinlichkeitsmatrix P[i][j] = P(X=i, Y=j)
        @param base: Logarithmusbasis
        @return: bedingte Entropie H(Y|X)
        @lastModified 2026-03-10
        """
        matrix = np.asarray(joint_probs, dtype=float)

        # Randverteilung von X (Zeilensummen)
        marginal_x = matrix.sum(axis=1)

        h_cond = 0.0
        for i, px in enumerate(marginal_x):
            if px > 0:
                for j in range(matrix.shape[1]):
                    p_xy = matrix[i, j]
                    if p_xy > 0:
                        # p(y|x) = p(x,y) / p(x)
                        p_y_given_x = p_xy / px
                        h_cond -= p_xy * math.log(p_y_given_x, base)
        return h_cond

    @staticmethod
    def mutual_information(joint_probs: list[list[float]], base: float = 2) -> float:
        """
        Berechnet die gegenseitige Information I(X;Y).

        Formel:
            I(X;Y) = H(X) + H(Y) - H(X,Y)
                   = H(X) - H(X|Y)
                   = H(Y) - H(Y|X)
                   = Σ_{x,y} p(x,y) · log [ p(x,y) / (p(x)·p(y)) ]

        Die gegenseitige Information ist immer ≥ 0 und misst die
        Abhängigkeit zwischen X und Y.

        @param joint_probs: 2D-Verbundwahrscheinlichkeitsmatrix
        @param base: Logarithmusbasis
        @return: gegenseitige Information I(X;Y) ≥ 0
        @lastModified 2026-03-10
        """
        matrix = np.asarray(joint_probs, dtype=float)

        # Randverteilungen
        marginal_x = matrix.sum(axis=1)
        marginal_y = matrix.sum(axis=0)

        mi = 0.0
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                p_xy = matrix[i, j]
                px = marginal_x[i]
                py = marginal_y[j]
                if p_xy > 0 and px > 0 and py > 0:
                    # I(X;Y) = Σ p(x,y) · log [ p(x,y) / (p(x)·p(y)) ]
                    mi += p_xy * math.log(p_xy / (px * py), base)
        return max(0.0, mi)  # numerische Fehler können minimal negativ sein

    @staticmethod
    def entropy_maximum(n: int, base: float = 2) -> float:
        """
        Berechnet die maximale Entropie für n gleichwahrscheinliche Ereignisse.

        Formel:
            H_max = log_b(n)

        Diese wird durch die Gleichverteilung erreicht (Maximum-Entropie-Prinzip).

        @param n: Anzahl der möglichen Ereignisse (n ≥ 1)
        @param base: Logarithmusbasis
        @return: maximale Entropie log_b(n)
        @raises ValueError: wenn n < 1
        @lastModified 2026-03-10
        """
        if n < 1:
            raise ValueError("Anzahl der Ereignisse muss mindestens 1 sein.")
        if n == 1:
            return 0.0  # Deterministisches Ereignis: keine Unsicherheit
        return math.log(n, base)

    @staticmethod
    def entropy_binary(p: float, base: float = 2) -> float:
        """
        Berechnet die binäre Entropiefunktion H(p) für Bernoulli-Verteilung.

        Formel:
            H(p) = -p·log_b(p) - (1-p)·log_b(1-p)

        Diese Funktion ist symmetrisch um p=0.5 und erreicht ihr Maximum
        H(0.5) = 1 Bit (bei b=2).

        @param p: Wahrscheinlichkeit für Ereignis A (0 ≤ p ≤ 1)
        @param base: Logarithmusbasis
        @return: binäre Entropie H(p) in der gewählten Einheit
        @raises ValueError: wenn p außerhalb [0,1]
        @lastModified 2026-03-10
        """
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p muss in [0,1] liegen, erhalten: {p}")

        # Randwerte: H(0) = H(1) = 0 (deterministisch)
        if p == 0.0 or p == 1.0:
            return 0.0

        return -p * math.log(p, base) - (1 - p) * math.log(1 - p, base)


# ===========================================================================
# KLASSE: KULLBACK-LEIBLER-DIVERGENZ
# ===========================================================================

class KLDivergence:
    """
    Kullback-Leibler-Divergenz und verwandte Informationsmaße.

    Die KL-Divergenz (auch relative Entropie) misst, wie sehr eine
    Wahrscheinlichkeitsverteilung P von einer Referenzverteilung Q abweicht.

    Formel:
        D_KL(P||Q) = Σ p_i · log(p_i / q_i)

    Wichtige Eigenschaften:
        - D_KL(P||Q) ≥ 0 (Gibb'sche Ungleichung)
        - D_KL(P||Q) = 0 ⟺ P = Q
        - Nicht symmetrisch: D_KL(P||Q) ≠ D_KL(Q||P)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def kl_divergence(p: list[float], q: list[float], base: float = 2) -> float:
        """
        Berechnet die Kullback-Leibler-Divergenz D_KL(P||Q).

        Formel:
            D_KL(P||Q) = Σ p_i · log_b(p_i / q_i)

        Konvention:
            - 0 · log(0/q) := 0
            - p · log(p/0) := +∞ (wenn p > 0 und q = 0)

        @param p: Verteilung P (Zähler)
        @param q: Verteilung Q (Referenz/Nenner)
        @param base: Logarithmusbasis
        @return: D_KL(P||Q) ≥ 0
        @raises ValueError: wenn q_i = 0 aber p_i > 0 (nicht im Träger)
        @lastModified 2026-03-10
        """
        p_arr = np.asarray(p, dtype=float)
        q_arr = np.asarray(q, dtype=float)

        if len(p_arr) != len(q_arr):
            raise ValueError("P und Q müssen die gleiche Länge haben.")

        kl = 0.0
        for pi, qi in zip(p_arr, q_arr):
            if pi > 0:
                if qi <= 0:
                    raise ValueError(
                        "Q darf keine Null haben, wo P > 0 ist (Träger-Bedingung verletzt)."
                    )
                kl += pi * math.log(pi / qi, base)
        return max(0.0, kl)

    @staticmethod
    def js_divergence(p: list[float], q: list[float], base: float = 2) -> float:
        """
        Berechnet die Jensen-Shannon-Divergenz JSD(P||Q).

        Formel:
            M = (P + Q) / 2
            JSD(P||Q) = (D_KL(P||M) + D_KL(Q||M)) / 2

        Die JSD ist symmetrisch und beschränkt:
            0 ≤ JSD(P||Q) ≤ log_b(2)

        Sie ist eine geglättete, symmetrische Version der KL-Divergenz
        und hat eine Wurzel, die eine echte Metrik ist.

        @param p: Verteilung P
        @param q: Verteilung Q
        @param base: Logarithmusbasis
        @return: Jensen-Shannon-Divergenz in [0, log_b(2)]
        @lastModified 2026-03-10
        """
        p_arr = np.asarray(p, dtype=float)
        q_arr = np.asarray(q, dtype=float)

        # Mittelpunkt-Verteilung
        m = (p_arr + q_arr) / 2.0

        # JSD als Mittel der KL-Divergenzen zum Mittelpunkt
        kl_pm = KLDivergence.kl_divergence(p_arr, m, base)
        kl_qm = KLDivergence.kl_divergence(q_arr, m, base)

        return (kl_pm + kl_qm) / 2.0

    @staticmethod
    def cross_entropy(p: list[float], q: list[float], base: float = 2) -> float:
        """
        Berechnet die Kreuzentropie H(P, Q).

        Formel:
            H(P,Q) = -Σ p_i · log_b(q_i)
                   = H(P) + D_KL(P||Q)

        Die Kreuzentropie misst die durchschnittliche Anzahl von Bits,
        die benötigt werden, um Ereignisse aus P mit einem Code zu beschreiben,
        der für Q optimiert ist.

        @param p: wahre Verteilung P
        @param q: geschätzte Verteilung Q (für den Code)
        @param base: Logarithmusbasis
        @return: Kreuzentropie H(P,Q) ≥ H(P)
        @lastModified 2026-03-10
        """
        p_arr = np.asarray(p, dtype=float)
        q_arr = np.asarray(q, dtype=float)

        ce = 0.0
        for pi, qi in zip(p_arr, q_arr):
            if pi > 0:
                if qi <= 0:
                    raise ValueError("q_i muss > 0 sein, wo p_i > 0.")
                ce -= pi * math.log(qi, base)
        return ce

    @staticmethod
    def renyi_entropy(probabilities: list[float], alpha: float, base: float = 2) -> float:
        """
        Berechnet die Rényi-Entropie H_α(X) der Ordnung α.

        Formel:
            H_α(X) = (1 / (1-α)) · log_b(Σ p_i^α)   für α ≠ 1

        Spezialfälle:
            α → 0  : Hartley-Entropie (Anzahl möglicher Ereignisse)
            α → 1  : Shannon-Entropie (Grenzwert)
            α → ∞  : Min-Entropie = -log(max p_i)

        @param probabilities: Wahrscheinlichkeitsverteilung
        @param alpha: Ordnungsparameter α ≥ 0, α ≠ 1
        @param base: Logarithmusbasis
        @return: Rényi-Entropie H_α(X)
        @raises ValueError: wenn α = 1 (dort Shannon-Entropie verwenden)
        @lastModified 2026-03-10
        """
        probs = np.asarray(probabilities, dtype=float)

        if math.isclose(alpha, 1.0, rel_tol=1e-9):
            # Grenzwert für α→1 ist die Shannon-Entropie
            return ShannonEntropy.entropy(probs, base)

        if alpha < 0:
            raise ValueError("alpha muss ≥ 0 sein.")

        if math.isclose(alpha, 0.0, abs_tol=1e-12):
            # H_0 = log(|Träger|) = Hartley-Entropie
            support_size = np.sum(probs > 0)
            return math.log(support_size, base) if support_size > 1 else 0.0

        if math.isinf(alpha):
            # Min-Entropie: H_∞ = -log(max p_i)
            max_p = np.max(probs[probs > 0])
            return -math.log(max_p, base)

        # Allgemeiner Fall: H_α = log(Σ p^α) / (1 - α)
        sum_p_alpha = np.sum(probs[probs > 0] ** alpha)
        return math.log(sum_p_alpha, base) / (1.0 - alpha)

    @staticmethod
    def tsallis_entropy(probabilities: list[float], q: float, base: float = 2) -> float:
        """
        Berechnet die Tsallis-Entropie S_q(X).

        Formel:
            S_q(X) = (1 / (q-1)) · (1 - Σ p_i^q)   für q ≠ 1

        Die Tsallis-Entropie ist eine nicht-additive Verallgemeinerung der
        Shannon-Entropie und wird in der statistischen Mechanik verwendet.

        Grenzwert q → 1: Shannon-Entropie (in Nat)

        @param probabilities: Wahrscheinlichkeitsverteilung
        @param q: Entraxieparameter q > 0, q ≠ 1
        @param base: Logarithmusbasis (beeinflusst die Einheit)
        @return: Tsallis-Entropie S_q(X)
        @lastModified 2026-03-10
        """
        probs = np.asarray(probabilities, dtype=float)
        valid = probs[probs > 0]

        if math.isclose(q, 1.0, rel_tol=1e-9):
            # Grenzwert: Shannon-Entropie
            return ShannonEntropy.entropy(probs, base)

        if q <= 0:
            raise ValueError("q muss > 0 sein.")

        # S_q = (1 - Σ p^q) / (q - 1)
        sum_p_q = np.sum(valid ** q)
        return (1.0 - sum_p_q) / (q - 1.0)

    @staticmethod
    def total_variation_distance(p: list[float], q: list[float]) -> float:
        """
        Berechnet die Total-Variation-Distanz TV(P, Q).

        Formel:
            TV(P,Q) = (1/2) · Σ |p_i - q_i|

        Die TV-Distanz ist eine echte Metrik auf Wahrscheinlichkeitsmaßen
        mit Werten in [0, 1].

        Beziehung zu anderen Maßen:
            TV(P,Q) ≤ √(JSD(P||Q) / 2)   (Pinsker-ähnliche Ungleichung)

        @param p: Verteilung P
        @param q: Verteilung Q
        @return: TV-Distanz in [0, 1]
        @lastModified 2026-03-10
        """
        p_arr = np.asarray(p, dtype=float)
        q_arr = np.asarray(q, dtype=float)

        if len(p_arr) != len(q_arr):
            raise ValueError("P und Q müssen die gleiche Länge haben.")

        return 0.5 * np.sum(np.abs(p_arr - q_arr))


# ===========================================================================
# KLASSE: KANALKAPAZITÄT
# ===========================================================================

class ChannelCapacity:
    """
    Kanalkapazität nach Shannon und verwandte Kanalmodelle.

    Shannons fundamentales Theorem (1948):
        Für jeden Kanal mit Kapazität C > 0 und jede Rate R < C
        existiert ein Codierungsschema, das beliebig zuverlässige
        Übertragung ermöglicht. Für R > C ist das nicht möglich.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def binary_symmetric_channel_capacity(p: float) -> float:
        """
        Berechnet die Kapazität des binären symmetrischen Kanals (BSC).

        Modell: Das gesendete Bit wird mit Wahrscheinlichkeit p verfälscht.

        Formel:
            C_BSC(p) = 1 - H(p) = 1 + p·log₂(p) + (1-p)·log₂(1-p)

        Grenzwerte:
            C(0) = 1 Bit (perfekter Kanal)
            C(0.5) = 0 Bit (komplett verrauschter Kanal)
            C(1) = 1 Bit (inverser perfekter Kanal)

        @param p: Bitfehlerwahrscheinlichkeit (0 ≤ p ≤ 1)
        @return: Kanalkapazität in Bit pro Kanalnutzung
        @lastModified 2026-03-10
        """
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p muss in [0,1] liegen, erhalten: {p}")

        h_p = ShannonEntropy.entropy_binary(p)
        return 1.0 - h_p

    @staticmethod
    def binary_erasure_channel_capacity(epsilon: float) -> float:
        """
        Berechnet die Kapazität des binären Auslöschungskanals (BEC).

        Modell: Das gesendete Bit wird mit Wahrscheinlichkeit ε ausgelöscht
        (Symbol '?') und mit Wahrscheinlichkeit 1-ε korrekt übertragen.

        Formel:
            C_BEC(ε) = 1 - ε

        @param epsilon: Auslöschungswahrscheinlichkeit (0 ≤ ε ≤ 1)
        @return: Kanalkapazität in Bit pro Kanalnutzung
        @lastModified 2026-03-10
        """
        if not 0.0 <= epsilon <= 1.0:
            raise ValueError(f"epsilon muss in [0,1] liegen, erhalten: {epsilon}")

        return 1.0 - epsilon

    @staticmethod
    def awgn_channel_capacity(snr_db: float) -> float:
        """
        Berechnet die Kanalkapazität des AWGN-Kanals (Shannon-Hartley-Theorem).

        Formel (Shannon-Hartley):
            C = log₂(1 + SNR)  [Bit/s/Hz]

        Dabei ist SNR das Signal-Rausch-Verhältnis (linear, nicht in dB).

        @param snr_db: Signal-Rausch-Verhältnis in Dezibel (dB)
        @return: Kanalkapazität in Bit pro Sekunde pro Hz Bandbreite
        @lastModified 2026-03-10
        """
        # Umrechnung dB → linear: SNR_lin = 10^(SNR_dB / 10)
        snr_linear = 10.0 ** (snr_db / 10.0)
        return math.log2(1.0 + snr_linear)

    @staticmethod
    def channel_capacity_general(transition_matrix: list[list[float]],
                                  n_iter: int = 200) -> float:
        """
        Berechnet die Kanalkapazität einer allgemeinen diskreten Kanalmatrix.

        Verwendet den Blahut-Arimoto-Algorithmus zur Optimierung.

        @param transition_matrix: Übergangsmatrix P[i][j] = P(Y=j | X=i)
                                   Zeilen: Eingangssymbole, Spalten: Ausgangssymbole
        @param n_iter: Anzahl der Iterationen für Blahut-Arimoto
        @return: Kanalkapazität in Bit pro Kanalnutzung
        @lastModified 2026-03-10
        """
        capacity, _ = ChannelCapacity.blahut_arimoto(transition_matrix, n_iter)
        return capacity

    @staticmethod
    def blahut_arimoto(transition_matrix: list[list[float]],
                        n_iter: int = 100) -> tuple[float, list[float]]:
        """
        Blahut-Arimoto-Algorithmus zur Berechnung der Kanalkapazität.

        Der Algorithmus wechselt zwischen zwei Schritten:
        1. Berechne Q(x|y) aus p(x) und P(y|x)
        2. Aktualisiere p(x) aus Q(x|y)

        Algorithmus (Blahut 1972, Arimoto 1972):
            Initialisierung: p(x) = 1/|X|
            Iteration:
                q(x|y) = p(x)·P(y|x) / Σ_{x'} p(x')·P(y|x')
                c(x) = Σ_y P(y|x) · log[ P(y|x) / Σ_{x'} p(x')·P(y|x') ]
                p_new(x) ∝ exp(c(x))

        Konvergenz: Die gegenseitige Information steigt monoton.

        @param transition_matrix: Kanalmatrix P(Y|X)
        @param n_iter: maximale Iterationsanzahl
        @return: Tupel (Kapazität in Bit, optimale Eingangsverteilung p*)
        @lastModified 2026-03-10
        """
        q_matrix = np.asarray(transition_matrix, dtype=float)
        n_x, n_y = q_matrix.shape

        # Initialisierung mit Gleichverteilung
        p = np.ones(n_x) / n_x

        for _ in range(n_iter):
            # Schritt 1: Rückwärtskanal Q(x|y) berechnen
            # p(x)·P(y|x) für alle x,y
            joint = p[:, np.newaxis] * q_matrix  # shape: (n_x, n_y)

            # Normalisierung über x für festes y → Q(x|y)
            marginal_y = joint.sum(axis=0)  # P(y) = Σ_x p(x)P(y|x)
            # Vermeidung von Division durch Null
            marginal_y = np.maximum(marginal_y, 1e-300)

            # Schritt 2: c(x) = Σ_y P(y|x) · log[ P(y|x) / P(y) ]
            # Berechne log[ P(y|x) / P(y) ] für alle x,y
            ratio = np.zeros_like(q_matrix)
            for i in range(n_x):
                for j in range(n_y):
                    if q_matrix[i, j] > 0 and marginal_y[j] > 0:
                        ratio[i, j] = math.log2(q_matrix[i, j] / marginal_y[j])

            # c(x) = Σ_y P(y|x) · log[...]
            c = np.sum(q_matrix * ratio, axis=1)

            # Schritt 3: p_new(x) ∝ p(x) · exp(c(x))
            # Normalisierung für numerische Stabilität
            c_shifted = c - np.max(c)
            p_new = p * np.exp(c_shifted)
            p_new /= p_new.sum()
            p = p_new

        # Berechne die Kanalkapazität mit der optimalen Eingangsverteilung
        marginal_y = (p[:, np.newaxis] * q_matrix).sum(axis=0)
        marginal_y = np.maximum(marginal_y, 1e-300)

        # C = Σ_{x,y} p(x)·P(y|x)·log[ P(y|x) / P(y) ]
        capacity = 0.0
        for i in range(n_x):
            for j in range(n_y):
                if q_matrix[i, j] > 0 and marginal_y[j] > 0:
                    capacity += p[i] * q_matrix[i, j] * math.log2(
                        q_matrix[i, j] / marginal_y[j]
                    )

        return max(0.0, capacity), p.tolist()

    @staticmethod
    def noisy_channel_coding_theorem_demo(capacity: float, rate: float) -> dict:
        """
        Demonstriert Shannons Kanalcodierungstheorem.

        Das Theorem besagt:
            - Für R < C: Es existiert ein Code mit verschwindender Fehlerwahrscheinlichkeit
            - Für R > C: Keine zuverlässige Übertragung möglich
            - R = C ist die Grenze (Shannonsche Kanalkapazität)

        @param capacity: Kanalkapazität C in Bit pro Kanalnutzung
        @param rate: gewünschte Codierungsrate R in Bit pro Symbol
        @return: Dictionary mit Analyseergebnis
        @lastModified 2026-03-10
        """
        margin = capacity - rate
        achievable = rate < capacity

        return {
            "capacity": capacity,
            "rate": rate,
            "margin": margin,
            "achievable": achievable,
            "message": (
                f"Rate R={rate:.4f} < C={capacity:.4f}: Zuverlässige Übertragung möglich."
                if achievable else
                f"Rate R={rate:.4f} ≥ C={capacity:.4f}: Zuverlässige Übertragung NICHT möglich."
            )
        }


# ===========================================================================
# KLASSE: QUELLENCODIERUNG
# ===========================================================================

class SourceCoding:
    """
    Quellencodierung: Huffman-Code, Arithmetisches Codieren, Lempel-Ziv.

    Shannons Quellencodierungssatz (1948):
        Für eine diskrete Quelle mit Entropie H(X) gilt für die mittlere
        Codelänge L̄ eines eindeutig dekodierbaren Präfixcodes:
            H(X) ≤ L̄ < H(X) + 1

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def huffman_code(symbols: list, probabilities: list[float]) -> dict:
        """
        Erzeugt einen optimalen Huffman-Präfixcode.

        Algorithmus (David Huffman, 1952):
        1. Erstelle einen Blattknoten für jedes Symbol mit seiner Wahrscheinlichkeit
        2. Entnimm wiederholt die zwei Knoten mit kleinster Wahrscheinlichkeit
        3. Verbinde sie zu einem neuen Knoten mit Summen-Wahrscheinlichkeit
        4. Der entstehende Baum definiert den optimalen Code

        Der Huffman-Code erfüllt den Quellencodierungssatz und ist optimal
        unter allen eindeutig dekodierbaren Präfixcodes.

        @param symbols: Liste der zu codierenden Symbole
        @param probabilities: zugehörige Wahrscheinlichkeiten
        @return: Dictionary {Symbol: Binärcode-String}
        @raises ValueError: wenn Längen nicht übereinstimmen
        @lastModified 2026-03-10
        """
        if len(symbols) != len(probabilities):
            raise ValueError("symbols und probabilities müssen gleich lang sein.")
        if len(symbols) == 0:
            return {}

        # Sonderfall: nur ein Symbol → Code ist "0"
        if len(symbols) == 1:
            return {symbols[0]: "0"}

        # Knoten werden als Dictionaries dargestellt, um Tupel-Längen-Konflikte zu vermeiden.
        # Heap-Einträge: [Wahrscheinlichkeit, tie-breaker, Knoten-Dict]
        # Knoten-Dict: {"prob": float, "sym": symbol|None, "left": Knoten|None, "right": Knoten|None}

        heap = []
        for idx, (sym, prob) in enumerate(zip(symbols, probabilities)):
            # Blattknoten
            node = {"prob": prob, "sym": sym, "left": None, "right": None}
            heapq.heappush(heap, [prob, idx, node])

        counter = len(symbols)  # Eindeutiger Tie-Breaker-Zähler

        # Huffman-Baum iterativ aufbauen
        while len(heap) > 1:
            # Zwei Knoten mit kleinster Wahrscheinlichkeit entnehmen
            prob1, _, node1 = heapq.heappop(heap)
            prob2, _, node2 = heapq.heappop(heap)

            # Neuer innerer Knoten: kombiniert node1 (links, "0") und node2 (rechts, "1")
            inner = {
                "prob": prob1 + prob2,
                "sym": None,
                "left": node1,
                "right": node2,
            }
            heapq.heappush(heap, [prob1 + prob2, counter, inner])
            counter += 1

        # Baumwurzel
        root_node = heap[0][2]

        # Codes durch Tiefensuche bestimmen
        codes = {}

        def _traverse(node: dict, current_code: str):
            """Rekursive Tiefensuche durch den Huffman-Baum."""
            if node["left"] is None and node["right"] is None:
                # Blattknoten: Code zuweisen
                codes[node["sym"]] = current_code if current_code else "0"
                return
            if node["left"]:
                _traverse(node["left"], current_code + "0")
            if node["right"]:
                _traverse(node["right"], current_code + "1")

        _traverse(root_node, "")
        return codes

    @staticmethod
    def huffman_average_length(huffman_dict: dict, probabilities: list[float]) -> float:
        """
        Berechnet die mittlere Codelänge L̄ eines Huffman-Codes.

        Formel:
            L̄ = Σ p_i · |c_i|

        wobei |c_i| die Länge des Codewortes für Symbol i ist.

        Der Quellencodierungssatz garantiert:
            H(X) ≤ L̄ < H(X) + 1

        @param huffman_dict: Dictionary {Symbol: Binärcode-String}
        @param probabilities: Wahrscheinlichkeiten in der gleichen Reihenfolge
        @return: mittlere Codelänge L̄ in Bit pro Symbol
        @lastModified 2026-03-10
        """
        symbols = list(huffman_dict.keys())
        if len(symbols) != len(probabilities):
            raise ValueError("Anzahl Symbole und Wahrscheinlichkeiten stimmen nicht überein.")

        total = 0.0
        for sym, prob in zip(symbols, probabilities):
            code_length = len(huffman_dict[sym])
            total += prob * code_length
        return total

    @staticmethod
    def arithmetic_coding_demo(symbols: list, probabilities: list[float],
                                message: list) -> dict:
        """
        Demonstriert das arithmetische Codierungsverfahren.

        Das arithmetische Codieren codiert eine gesamte Nachricht als
        eine einzige rationale Zahl im Intervall [0, 1).

        Algorithmus:
        1. Starte mit [low, high) = [0, 1)
        2. Für jedes Symbol: Teile [low, high) proportional zu den Wahrscheinlichkeiten
        3. Wähle das Teilintervall des aktuellen Symbols

        Effizienz: Nähert sich H(X) Bit/Symbol für lange Nachrichten.

        @param symbols: Alphabetliste
        @param probabilities: zugehörige Wahrscheinlichkeiten
        @param message: zu codierende Nachricht (Liste von Symbolen)
        @return: Dictionary mit Codierungsergebnis (low, high, bits)
        @lastModified 2026-03-10
        """
        if not message:
            return {"low": 0.0, "high": 1.0, "bits": 0, "message": []}

        # Kumulierte Wahrscheinlichkeiten für Intervallaufteilung
        cum_probs = [0.0]
        for p in probabilities:
            cum_probs.append(cum_probs[-1] + p)

        sym_to_idx = {s: i for i, s in enumerate(symbols)}

        low = 0.0
        high = 1.0

        for sym in message:
            if sym not in sym_to_idx:
                raise ValueError(f"Symbol '{sym}' nicht im Alphabet.")
            idx = sym_to_idx[sym]
            width = high - low

            # Neues Intervall proportional zur Symbolwahrscheinlichkeit
            high = low + width * cum_probs[idx + 1]
            low = low + width * cum_probs[idx]

        # Benötigte Bits: -log₂(Intervallbreite)
        interval_width = high - low
        bits_needed = math.ceil(-math.log2(interval_width)) if interval_width > 0 else 0

        return {
            "low": low,
            "high": high,
            "midpoint": (low + high) / 2.0,
            "interval_width": interval_width,
            "bits": bits_needed,
            "message": message
        }

    @staticmethod
    def lempel_ziv_complexity(sequence: list) -> int:
        """
        Berechnet die Lempel-Ziv-Komplexität einer Sequenz.

        Die LZ-Komplexität (Lempel & Ziv, 1976) zählt die Anzahl
        verschiedener Teilstrings/Wörter im LZ76-Parsing-Schema.
        Sie ist ein Maß für die Zufälligkeit/Kompressibilität einer Sequenz.

        Algorithmus (LZ76):
        - Zerlege S in minimale Teilstrings, die noch nicht als Präfix erschienen sind.

        @param sequence: Eingangssequenz (Liste von Symbolen)
        @return: LZ76-Komplexität (Anzahl der Wörter im LZ-Parsing)
        @lastModified 2026-03-10
        """
        if not sequence:
            return 0

        # Konvertiere in Tuple für Hashbarkeit
        seq = tuple(sequence)
        n = len(seq)

        complexity = 0
        seen = set()  # bisher gesehene Teilstrings
        c = 0         # Startindex des aktuellen Wortes
        l = 1         # Länge des aktuellen Wortes

        while c + l <= n:
            current = seq[c:c + l]
            if current in seen:
                # Noch nicht minimal: Wort verlängern
                l += 1
            else:
                # Neues Wort gefunden
                seen.add(current)
                complexity += 1
                c += l   # Nächstes Wort beginnt nach diesem
                l = 1     # Länge zurücksetzen

        return complexity

    @staticmethod
    def kolmogorov_complexity_demo(n: int) -> dict:
        """
        Demonstration der Kolmogorov-Komplexität K(x).

        Die Kolmogorov-Komplexität K(x) ist die Länge des kürzesten
        Programms, das x ausgibt. Sie ist theoretisch nicht berechenbar
        (Alan Turing, 1936; Andrey Kolmogorov, 1963).

        Kompressionsbasierte Näherung: K(x) ≈ |compress(x)|

        @param n: Länge der Testsequenz
        @return: Dictionary mit Demonstrationsdaten
        @lastModified 2026-03-10
        """
        import zlib

        # Verschiedene Sequenzen zum Vergleich
        all_zeros = bytes([0] * n)
        alternating = bytes([i % 2 for i in range(n)])
        random_bytes = bytes([i % 256 for i in range(n)])

        def compressed_length(data: bytes) -> int:
            return len(zlib.compress(data, level=9))

        return {
            "n": n,
            "all_zeros": {
                "raw_bytes": n,
                "compressed_bytes": compressed_length(all_zeros),
                "description": "Konstante Folge – sehr niedrige Komplexität"
            },
            "alternating": {
                "raw_bytes": n,
                "compressed_bytes": compressed_length(alternating),
                "description": "Alternierende Folge – niedrige Komplexität"
            },
            "counting": {
                "raw_bytes": n,
                "compressed_bytes": compressed_length(random_bytes),
                "description": "Zählende Folge – mittlere Komplexität"
            },
            "note": "K(x) ist nicht berechenbar; Kompression ist eine obere Schranke."
        }

    @staticmethod
    def source_coding_theorem_bound(probabilities: list[float], base: float = 2) -> dict:
        """
        Berechnet die Schranken des Quellencodierungssatzes.

        Shannons Quellencodierungssatz besagt:
            H(X) ≤ L̄ < H(X) + 1

        Eine perfekte Codierung würde genau H(X) Bits benötigen.
        Huffman-Codes erfüllen diese Schranke.

        @param probabilities: Wahrscheinlichkeitsverteilung der Quelle
        @param base: Logarithmusbasis (Standard: 2)
        @return: Dictionary mit H(X), Schranken und Effizienz
        @lastModified 2026-03-10
        """
        h = ShannonEntropy.entropy(probabilities, base)
        return {
            "entropy": h,
            "lower_bound": h,           # L̄ ≥ H(X)
            "upper_bound": h + 1.0,     # L̄ < H(X) + 1
            "optimal_bits_per_symbol": h,
            "unit": "Bit" if base == 2 else f"log_{base}"
        }


# ===========================================================================
# KLASSE: FEHLERKORREKTUR
# ===========================================================================

class ErrorCorrection:
    """
    Grundlagen der Fehlerkorrektur-Codes (Channel Coding).

    Fehlerkorrektur-Codes ermöglichen die Erkennung und Korrektur von
    Übertragungsfehlern. Ein linearer (n,k)-Code überträgt k Informationsbits
    in n Bits (n > k), die d = n-k Redundanzbits enthalten.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def hamming_distance(codeword1: list[int], codeword2: list[int]) -> int:
        """
        Berechnet den Hamming-Abstand zwischen zwei Codewörtern.

        Formel:
            d_H(c1, c2) = |{i : c1[i] ≠ c2[i]}|

        Der Hamming-Abstand ist eine Metrik auf dem Vektorraum GF(2)^n.

        @param codeword1: erstes Binär-Codewort (Liste von 0/1)
        @param codeword2: zweites Binär-Codewort (gleiche Länge)
        @return: Hamming-Abstand (Anzahl unterschiedlicher Positionen)
        @raises ValueError: wenn Längen unterschiedlich
        @lastModified 2026-03-10
        """
        if len(codeword1) != len(codeword2):
            raise ValueError("Codewörter müssen gleich lang sein.")

        return sum(b1 != b2 for b1, b2 in zip(codeword1, codeword2))

    @staticmethod
    def hamming_weight(codeword: list[int]) -> int:
        """
        Berechnet das Hamming-Gewicht (Anzahl der Einsen) eines Codewortes.

        Formel:
            w_H(c) = |{i : c[i] = 1}| = d_H(c, 0...0)

        @param codeword: Binär-Codewort
        @return: Hamming-Gewicht (Anzahl der Einsen)
        @lastModified 2026-03-10
        """
        return sum(1 for b in codeword if b == 1)

    @staticmethod
    def minimum_distance(code: list[list[int]]) -> int:
        """
        Berechnet die Minimaldistanz d_min eines binären Codes.

        Formel:
            d_min = min_{c1 ≠ c2 ∈ C} d_H(c1, c2)

        Die Minimaldistanz bestimmt die Fehlerkorrekturkapazität:
            - Erkennung von t Fehlern: d_min ≥ t + 1
            - Korrektur von t Fehlern: d_min ≥ 2t + 1

        @param code: Liste aller Codewörter
        @return: Minimaldistanz d_min
        @raises ValueError: wenn der Code weniger als 2 Codewörter hat
        @lastModified 2026-03-10
        """
        if len(code) < 2:
            raise ValueError("Code muss mindestens 2 Codewörter enthalten.")

        min_dist = float('inf')
        for i in range(len(code)):
            for j in range(i + 1, len(code)):
                d = ErrorCorrection.hamming_distance(code[i], code[j])
                if d < min_dist:
                    min_dist = d

        return int(min_dist)

    @staticmethod
    def singleton_bound(n: int, k: int) -> int:
        """
        Berechnet die Singleton-Schranke für lineare (n,k)-Codes.

        Formel:
            d_max ≤ n - k + 1

        Codes, die diese Schranke erreichen, heißen MDS-Codes
        (Maximum Distance Separable), z.B. Reed-Solomon-Codes.

        @param n: Codewortlänge
        @param k: Anzahl der Informationsbits
        @return: maximale erreichbare Minimaldistanz d ≤ n-k+1
        @lastModified 2026-03-10
        """
        if k > n:
            raise ValueError("k darf nicht größer als n sein.")
        return n - k + 1

    @staticmethod
    def hamming_bound(n: int, t: int) -> int:
        """
        Berechnet die Hamming-Schranke (Kugelpackungsschranke).

        Für einen binären (n,k)-Code, der t Fehler korrigiert, gilt:
            2^n / 2^k ≤ Σ_{i=0}^{t} C(n,i)

        Codes, die die Gleichheit erreichen, heißen perfekte Codes.
        Das bekannteste Beispiel ist der (7,4)-Hamming-Code (t=1).

        @param n: Codewortlänge
        @param t: Fehlerkorrekturvermögen
        @return: Volumen der Hamming-Kugel Σ_{i=0}^{t} C(n,i)
        @lastModified 2026-03-10
        """
        total = 0
        for i in range(t + 1):
            total += math.comb(n, i)
        return total

    @staticmethod
    def hamming_code_74() -> dict:
        """
        Erzeugt den klassischen (7,4)-Hamming-Code.

        Der (7,4)-Hamming-Code:
            - n = 7 Bits Codewortlänge
            - k = 4 Bits Information
            - d_min = 3 (korrigiert 1 Fehler, erkennt 2)
            - 2^7 / 2^4 = 8 = C(7,0)+C(7,1) = 1+7 → perfekter Code!

        Generatormatrix G (4×7, systematische Form):
            Informationsbits + Paritätsbits

        Prüfmatrix H (3×7):
            H·c = 0 für alle gültigen Codewörter c

        @return: Dictionary mit G, H und Code-Parametern
        @lastModified 2026-03-10
        """
        # Generatormatrix G in systematischer Form [I_k | P]
        # Zeilen: Informationsbits, Spalten: Codeword-Positionen
        G = np.array([
            [1, 0, 0, 0,  1, 0, 1],  # d1
            [0, 1, 0, 0,  1, 1, 0],  # d2
            [0, 0, 1, 0,  1, 1, 1],  # d3
            [0, 0, 0, 1,  0, 1, 1],  # d4
        ], dtype=int)

        # Prüfmatrix H (3×7): H·G^T = 0 (mod 2)
        H = np.array([
            [1, 1, 1, 0,  1, 0, 0],  # p1
            [0, 1, 1, 1,  0, 1, 0],  # p2
            [1, 0, 1, 1,  0, 0, 1],  # p3
        ], dtype=int)

        return {
            "n": 7,
            "k": 4,
            "d_min": 3,
            "t_correctable": 1,
            "generator_matrix": G,
            "parity_check_matrix": H,
            "code_rate": 4 / 7,
            "is_perfect": True,
            "description": "(7,4)-Hamming-Code: perfekter 1-Fehler-korrigierender Code"
        }

    @staticmethod
    def parity_check_code(n: int) -> dict:
        """
        Erzeugt einen einfachen Paritätscode der Länge n.

        Der (n, n-1)-Paritätscode fügt ein einzelnes Paritätsbit hinzu:
            d_min = 2 (erkennt 1 Fehler, korrigiert 0)
            Coderate = (n-1)/n

        @param n: Gesamtlänge des Codewortes (n ≥ 2)
        @return: Dictionary mit Code-Parametern und Generatormatrix
        @lastModified 2026-03-10
        """
        if n < 2:
            raise ValueError("n muss mindestens 2 sein.")

        k = n - 1

        # Generatormatrix: Einheitsmatrix + Paritätsspalte (Zeilensummen mod 2)
        G = np.zeros((k, n), dtype=int)
        G[:, :k] = np.eye(k, dtype=int)
        # Letztes Bit: XOR aller Informationsbits
        G[:, -1] = 1  # Jedes Informationsbit trägt zum Paritätsbit bei

        # Prüfmatrix: [1, 1, ..., 1] (Summe aller Bits muss 0 mod 2 sein)
        H = np.ones((1, n), dtype=int)

        return {
            "n": n,
            "k": k,
            "d_min": 2,
            "t_correctable": 0,
            "t_detectable": 1,
            "generator_matrix": G,
            "parity_check_matrix": H,
            "code_rate": k / n
        }

    @staticmethod
    def repetition_code(message: list[int], k: int) -> list[int]:
        """
        Codiert eine Nachricht mit dem k-fachen Wiederholungscode.

        Jedes Bit wird k-mal wiederholt:
            m → m, m, ..., m  (k-mal)

        Eigenschaften:
            d_min = k, korrigiert ⌊(k-1)/2⌋ Fehler
            Coderate = 1/k

        Dekodierung durch Mehrheitsentscheid:
            Bit = 1 wenn mehr als k/2 Wiederholungen = 1

        @param message: Bitnachricht (Liste von 0/1)
        @param k: Wiederholungsfaktor (k ≥ 1)
        @return: codierte Nachricht der Länge k · |message|
        @lastModified 2026-03-10
        """
        if k < 1:
            raise ValueError("k muss mindestens 1 sein.")

        encoded = []
        for bit in message:
            encoded.extend([bit] * k)
        return encoded


# ===========================================================================
# KLASSE: DIFFERENTIELLE ENTROPIE
# ===========================================================================

class DifferentialEntropy:
    """
    Differentielle Entropie für kontinuierliche Zufallsvariablen.

    Die differentielle Entropie ist das kontinuierliche Analogon zur
    Shannon-Entropie für stetige Wahrscheinlichkeitsverteilungen:

        h(X) = -∫ f(x) · ln f(x) dx

    Im Gegensatz zur diskreten Entropie kann h(X) negativ sein.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def differential_entropy_normal(sigma: float, base: float = 2) -> float:
        """
        Berechnet die differentielle Entropie der Normalverteilung N(μ, σ²).

        Formel (unabhängig von μ):
            h(X) = (1/2) · log_b(2πe·σ²)
                 = (1/2) · log_b(2π) + (1/2)·log_b(e) + log_b(σ)

        Die Normalverteilung maximiert die differentielle Entropie
        unter allen Verteilungen mit fester Varianz σ².

        @param sigma: Standardabweichung σ > 0
        @param base: Logarithmusbasis (Standard: 2)
        @return: differentielle Entropie h(X) in der gewählten Einheit
        @raises ValueError: wenn σ ≤ 0
        @lastModified 2026-03-10
        """
        if sigma <= 0:
            raise ValueError(f"Standardabweichung muss > 0 sein, erhalten: {sigma}")

        # h = (1/2) · log(2πe·σ²) = (1/2) · log(2πe) + log(σ)
        h_nat = 0.5 * math.log(2 * math.pi * math.e * sigma ** 2)

        # Umrechnung in gewünschte Basis
        return h_nat / math.log(base)

    @staticmethod
    def differential_entropy_exponential(lambda_: float, base: float = 2) -> float:
        """
        Berechnet die differentielle Entropie der Exponentialverteilung Exp(λ).

        Formel:
            h(X) = 1 - log_b(λ)  [in Nat: h = 1 - ln(λ)]

        Die Exponentialverteilung maximiert die differentielle Entropie
        unter allen Verteilungen mit festem Erwartungswert E[X] = 1/λ.

        @param lambda_: Ratenparameter λ > 0
        @param base: Logarithmusbasis
        @return: differentielle Entropie h(X)
        @raises ValueError: wenn λ ≤ 0
        @lastModified 2026-03-10
        """
        if lambda_ <= 0:
            raise ValueError(f"Ratenparameter muss > 0 sein, erhalten: {lambda_}")

        # h = 1 - log_b(λ) = 1/ln(b) - ln(λ)/ln(b)
        h_nat = 1.0 - math.log(lambda_)
        return h_nat / math.log(base)

    @staticmethod
    def differential_entropy_uniform(a: float, b: float, base: float = 2) -> float:
        """
        Berechnet die differentielle Entropie der Gleichverteilung U(a,b).

        Formel:
            h(X) = log_b(b - a)

        Für b-a < 1 ist h(X) < 0 (kein Widerspruch – differentielle
        Entropie kann negativ sein).

        @param a: untere Grenze des Intervalls
        @param b: obere Grenze (a < b)
        @param base: Logarithmusbasis
        @return: differentielle Entropie h(X) = log_b(b-a)
        @raises ValueError: wenn a ≥ b
        @lastModified 2026-03-10
        """
        if a >= b:
            raise ValueError(f"Es muss a < b gelten, erhalten: a={a}, b={b}")

        return math.log(b - a, base)

    @staticmethod
    def mutual_information_gaussian(rho: float, base: float = 2) -> float:
        """
        Berechnet die gegenseitige Information einer bivariaten Normalverteilung.

        Für (X,Y) ~ N(0, [[1,ρ],[ρ,1]]) mit Korrelationskoeffizient ρ:

        Formel:
            I(X;Y) = -log_b(√(1 - ρ²))
                   = -(1/2) · log_b(1 - ρ²)

        @param rho: Korrelationskoeffizient ρ ∈ (-1, 1)
        @param base: Logarithmusbasis
        @return: gegenseitige Information I(X;Y) ≥ 0
        @raises ValueError: wenn |ρ| ≥ 1
        @lastModified 2026-03-10
        """
        if not -1.0 < rho < 1.0:
            raise ValueError(f"Korrelationskoeffizient muss in (-1,1) liegen, erhalten: {rho}")

        # I(X;Y) = -(1/2) · log_b(1 - ρ²)
        return -0.5 * math.log(1.0 - rho ** 2, base)


# ===========================================================================
# STANDALONE-FUNKTIONEN
# ===========================================================================

def surprisal(p: float, base: float = 2) -> float:
    """
    Berechnet den Informationsgehalt (Surprisal) eines Ereignisses.

    Formel:
        I(x) = -log_b(p(x))

    Interpretation:
        - p = 1 → I = 0 Bit (sicheres Ereignis, keine Information)
        - p = 0.5 → I = 1 Bit (Münzwurf)
        - p → 0 → I → ∞ (sehr seltenes Ereignis, viel Information)

    Die Shannon-Entropie ist der Erwartungswert des Surprisal:
        H(X) = E[I(X)] = Σ p_i · I(x_i)

    @param p: Wahrscheinlichkeit des Ereignisses (0 < p ≤ 1)
    @param base: Logarithmusbasis (Standard: 2 → Bit)
    @return: Informationsgehalt in Bit (oder Nat/Hartley je nach base)
    @raises ValueError: wenn p ≤ 0 oder p > 1
    @lastModified 2026-03-10
    """
    if p <= 0 or p > 1.0:
        raise ValueError(f"p muss in (0, 1] liegen, erhalten: {p}")

    return -math.log(p, base)


def entropy_rate_markov(transition_matrix: list[list[float]],
                         stationary_dist: list[float],
                         base: float = 2) -> float:
    """
    Berechnet die Entropierate H einer stationären Markov-Kette.

    Formel:
        H(X) = -Σ_{i} π_i · Σ_{j} P_{ij} · log P_{ij}

    wobei π die stationäre Verteilung und P die Übergangsmatrix ist.

    Die Entropierate gibt die durchschnittliche Entropie pro Schritt an.

    @param transition_matrix: Übergangsmatrix P[i][j] = P(X_{n+1}=j | X_n=i)
    @param stationary_dist: stationäre Verteilung π
    @param base: Logarithmusbasis
    @return: Entropierate H in Bit (oder gewählter Einheit) pro Schritt
    @lastModified 2026-03-10
    """
    P = np.asarray(transition_matrix, dtype=float)
    pi = np.asarray(stationary_dist, dtype=float)

    n = len(pi)
    h_rate = 0.0

    for i in range(n):
        # Bedingte Entropie des Zustands j gegeben Zustand i
        for j in range(n):
            if P[i, j] > 0:
                h_rate -= pi[i] * P[i, j] * math.log(P[i, j], base)

    return h_rate


def data_processing_inequality_demo() -> dict:
    """
    Demonstriert die Daten-Verarbeitungs-Ungleichung (DPI).

    Theorem (Data Processing Inequality):
        Für die Markov-Kette X → Y → Z gilt:
            I(X;Z) ≤ I(X;Y)

    Informationsgehalt kann durch Verarbeitung nicht zunehmen.
    Anwendung: Tiefe neuronale Netzwerke müssen Information über
    die Eingabe in jeder Schicht bewahren.

    @return: Dictionary mit Demonstrations-Daten und Verifikation
    @lastModified 2026-03-10
    """
    # Beispiel-Kanal X → Y (BSC mit p=0.1)
    p_xy = 0.1  # Bitfehlerrate X → Y
    channel_XY = [
        [1 - p_xy, p_xy],
        [p_xy, 1 - p_xy]
    ]

    # Weiterer Kanal Y → Z (BSC mit p=0.2)
    p_yz = 0.2  # Bitfehlerrate Y → Z
    channel_YZ = [
        [1 - p_yz, p_yz],
        [p_yz, 1 - p_yz]
    ]

    # Kombinierter Kanal X → Z (BSC mit p_xz = p_xy + p_yz - 2*p_xy*p_yz)
    p_xz = p_xy + p_yz - 2 * p_xy * p_yz
    channel_XZ = [
        [1 - p_xz, p_xz],
        [p_xz, 1 - p_xz]
    ]

    # Kanalkapazitäten
    c_xy = ChannelCapacity.binary_symmetric_channel_capacity(p_xy)
    c_xz = ChannelCapacity.binary_symmetric_channel_capacity(p_xz)

    # Gegenseitige Informationen mit Gleichverteilung p(X) = [0.5, 0.5]
    joint_xy = [[0.5 * (1 - p_xy), 0.5 * p_xy],
                 [0.5 * p_yz, 0.5 * (1 - p_yz)]]

    # Direkter Ansatz: I(X;Y) = C für optimale Eingangsverteilung bei BSC
    i_xy = c_xy
    i_xz = c_xz

    return {
        "theorem": "I(X;Z) ≤ I(X;Y) für Markov-Kette X → Y → Z",
        "channel_XY": {"error_prob": p_xy, "capacity": c_xy},
        "channel_XZ": {"error_prob": p_xz, "capacity": c_xz},
        "I_XY": i_xy,
        "I_XZ": i_xz,
        "dpi_satisfied": i_xz <= i_xy + 1e-10,
        "information_loss": i_xy - i_xz
    }
