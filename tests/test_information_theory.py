"""
@file test_information_theory.py
@brief Umfassende Tests für das information_theory-Modul.
@description
    Testet alle Klassen und Funktionen des Informationstheorie-Moduls:

    - ShannonEntropy      : Entropie, bedingte/gemeinsame Entropie, gegenseitige Information
    - KLDivergence        : KL-Divergenz, JS-Divergenz, Kreuzentropie, Rényi, Tsallis, TV
    - ChannelCapacity     : BSC, BEC, AWGN, Blahut-Arimoto, allgemeine Kanäle
    - SourceCoding        : Huffman-Code, Lempel-Ziv, Quellencodierungssatz
    - ErrorCorrection     : Hamming-Distanz, -Gewicht, -Schranke, (7,4)-Code
    - DifferentialEntropy : Normal, Exponential, Uniform, Gauß'sche MI
    - Standalone          : surprisal, entropy_rate_markov, DPI-Demo

    Edge Cases:
    - p = 0 (deterministisch)
    - p = 1 (sicheres Ereignis)
    - Gleichverteilung (maximale Entropie)
    - leere Sequenzen
    - Einzel-Symbol-Codierung

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import math
import sys
import os
import pytest
import numpy as np

# Sicherstellen, dass src/ im Suchpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from information_theory import (
    ShannonEntropy,
    KLDivergence,
    ChannelCapacity,
    SourceCoding,
    ErrorCorrection,
    DifferentialEntropy,
    surprisal,
    entropy_rate_markov,
    data_processing_inequality_demo,
)


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def approx_eq(a: float, b: float, tol: float = 1e-6) -> bool:
    """Prüft ob zwei Zahlen näherungsweise gleich sind."""
    return abs(a - b) <= tol


# ===========================================================================
# TESTS: SHANNON-ENTROPIE
# ===========================================================================

class TestShannonEntropy:
    """Tests für die Klasse ShannonEntropy."""

    # --- Grundlegende Entropie ---

    def test_entropy_uniform_two_symbols(self):
        """Gleichverteilung über 2 Symbole → H = 1 Bit."""
        h = ShannonEntropy.entropy([0.5, 0.5])
        assert approx_eq(h, 1.0), f"Erwartet 1.0, erhalten {h}"

    def test_entropy_uniform_four_symbols(self):
        """Gleichverteilung über 4 Symbole → H = 2 Bit."""
        h = ShannonEntropy.entropy([0.25, 0.25, 0.25, 0.25])
        assert approx_eq(h, 2.0), f"Erwartet 2.0, erhalten {h}"

    def test_entropy_uniform_eight_symbols(self):
        """Gleichverteilung über 8 Symbole → H = 3 Bit."""
        h = ShannonEntropy.entropy([1/8] * 8)
        assert approx_eq(h, 3.0)

    def test_entropy_deterministic_zero(self):
        """Deterministisches Ereignis p=1 → H = 0."""
        h = ShannonEntropy.entropy([1.0, 0.0])
        assert approx_eq(h, 0.0)

    def test_entropy_deterministic_only_one(self):
        """Einzelnes Ereignis p=1 → H = 0."""
        h = ShannonEntropy.entropy([1.0])
        assert approx_eq(h, 0.0)

    def test_entropy_with_zeros(self):
        """Null-Wahrscheinlichkeiten werden korrekt behandelt (0·log(0)=0)."""
        h = ShannonEntropy.entropy([0.5, 0.5, 0.0])
        assert approx_eq(h, 1.0)

    def test_entropy_non_uniform(self):
        """Nicht-uniforme Verteilung: H < log₂(n)."""
        h = ShannonEntropy.entropy([0.9, 0.1])
        expected = -0.9 * math.log2(0.9) - 0.1 * math.log2(0.1)
        assert approx_eq(h, expected)

    def test_entropy_base_nats(self):
        """Entropie in Nat (Basis e)."""
        h = ShannonEntropy.entropy([0.5, 0.5], base=math.e)
        assert approx_eq(h, math.log(2))

    def test_entropy_base_10(self):
        """Entropie in Hartley (Basis 10)."""
        h = ShannonEntropy.entropy([0.1] * 10, base=10)
        assert approx_eq(h, 1.0)

    def test_entropy_negative_prob_raises(self):
        """Negative Wahrscheinlichkeit löst ValueError aus."""
        with pytest.raises(ValueError):
            ShannonEntropy.entropy([-0.1, 1.1])

    def test_entropy_sum_not_one_raises(self):
        """Summe ≠ 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ShannonEntropy.entropy([0.3, 0.3, 0.3])

    def test_entropy_maximized_by_uniform(self):
        """Gleichverteilung maximiert die Entropie."""
        h_uniform = ShannonEntropy.entropy([1/4, 1/4, 1/4, 1/4])
        h_nonuniform = ShannonEntropy.entropy([0.7, 0.1, 0.1, 0.1])
        assert h_uniform >= h_nonuniform

    # --- Gemeinsame Entropie ---

    def test_joint_entropy_independent(self):
        """Unabhängige X,Y: H(X,Y) = H(X) + H(Y)."""
        # p(x,y) = p(x)·p(y) mit je p(x)=[0.5,0.5], p(y)=[0.5,0.5]
        joint = [[0.25, 0.25], [0.25, 0.25]]
        h_joint = ShannonEntropy.joint_entropy(joint)
        assert approx_eq(h_joint, 2.0)  # H(X) + H(Y) = 1 + 1

    def test_joint_entropy_deterministic(self):
        """Perfekte Abhängigkeit: H(X,Y) = H(X)."""
        joint = [[0.5, 0.0], [0.0, 0.5]]
        h_joint = ShannonEntropy.joint_entropy(joint)
        assert approx_eq(h_joint, 1.0)

    def test_joint_entropy_sum_not_one_raises(self):
        """Joint-Matrix mit Summe ≠ 1 löst Fehler aus."""
        with pytest.raises(ValueError):
            ShannonEntropy.joint_entropy([[0.3, 0.3], [0.3, 0.3]])

    # --- Bedingte Entropie ---

    def test_conditional_entropy_independent(self):
        """Unabhängige X,Y: H(Y|X) = H(Y)."""
        joint = [[0.25, 0.25], [0.25, 0.25]]
        h_cond = ShannonEntropy.conditional_entropy(joint)
        assert approx_eq(h_cond, 1.0)

    def test_conditional_entropy_determined(self):
        """Y vollständig durch X bestimmt: H(Y|X) = 0."""
        joint = [[0.5, 0.0], [0.0, 0.5]]
        h_cond = ShannonEntropy.conditional_entropy(joint)
        assert approx_eq(h_cond, 0.0)

    def test_conditional_entropy_chain_rule(self):
        """Kettenregel: H(Y|X) = H(X,Y) - H(X)."""
        joint = [[0.3, 0.1], [0.2, 0.4]]
        h_joint = ShannonEntropy.joint_entropy(joint)
        h_x = ShannonEntropy.entropy([0.4, 0.6])  # Randverteilung X
        h_cond = ShannonEntropy.conditional_entropy(joint)
        assert approx_eq(h_cond, h_joint - h_x, tol=1e-5)

    # --- Gegenseitige Information ---

    def test_mutual_information_independent(self):
        """Unabhängige X,Y: I(X;Y) = 0."""
        joint = [[0.25, 0.25], [0.25, 0.25]]
        mi = ShannonEntropy.mutual_information(joint)
        assert approx_eq(mi, 0.0)

    def test_mutual_information_perfect_dependence(self):
        """Perfekte Abhängigkeit: I(X;Y) = H(X)."""
        joint = [[0.5, 0.0], [0.0, 0.5]]
        mi = ShannonEntropy.mutual_information(joint)
        assert approx_eq(mi, 1.0)

    def test_mutual_information_nonnegative(self):
        """I(X;Y) ≥ 0 für alle Verbundverteilungen."""
        joint = [[0.1, 0.4], [0.3, 0.2]]
        mi = ShannonEntropy.mutual_information(joint)
        assert mi >= -1e-9

    def test_mutual_information_symmetry(self):
        """I(X;Y) = I(Y;X) (Symmetrie)."""
        joint = [[0.2, 0.3], [0.1, 0.4]]
        joint_T = [[joint[j][i] for j in range(2)] for i in range(2)]
        mi_xy = ShannonEntropy.mutual_information(joint)
        mi_yx = ShannonEntropy.mutual_information(joint_T)
        assert approx_eq(mi_xy, mi_yx)

    # --- Maximale Entropie ---

    def test_entropy_maximum_two(self):
        """Maximale Entropie für n=2: log₂(2) = 1."""
        assert approx_eq(ShannonEntropy.entropy_maximum(2), 1.0)

    def test_entropy_maximum_one(self):
        """Maximale Entropie für n=1: log₂(1) = 0."""
        assert approx_eq(ShannonEntropy.entropy_maximum(1), 0.0)

    def test_entropy_maximum_power_of_two(self):
        """Maximale Entropie für n=8: log₂(8) = 3."""
        assert approx_eq(ShannonEntropy.entropy_maximum(8), 3.0)

    def test_entropy_maximum_invalid_raises(self):
        """n < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ShannonEntropy.entropy_maximum(0)

    # --- Binäre Entropie ---

    def test_binary_entropy_half(self):
        """H(0.5) = 1 Bit (Maximum)."""
        assert approx_eq(ShannonEntropy.entropy_binary(0.5), 1.0)

    def test_binary_entropy_zero(self):
        """H(0) = 0 (deterministisch)."""
        assert approx_eq(ShannonEntropy.entropy_binary(0.0), 0.0)

    def test_binary_entropy_one(self):
        """H(1) = 0 (deterministisch)."""
        assert approx_eq(ShannonEntropy.entropy_binary(1.0), 0.0)

    def test_binary_entropy_symmetry(self):
        """H(p) = H(1-p) (Symmetrie um 0.5)."""
        for p in [0.1, 0.2, 0.3, 0.4]:
            assert approx_eq(
                ShannonEntropy.entropy_binary(p),
                ShannonEntropy.entropy_binary(1 - p)
            )

    def test_binary_entropy_invalid_raises(self):
        """p außerhalb [0,1] löst ValueError aus."""
        with pytest.raises(ValueError):
            ShannonEntropy.entropy_binary(1.5)
        with pytest.raises(ValueError):
            ShannonEntropy.entropy_binary(-0.1)

    def test_binary_entropy_matches_general_entropy(self):
        """Binäre Entropie stimmt mit allgemeiner Entropie überein."""
        p = 0.3
        h_binary = ShannonEntropy.entropy_binary(p)
        h_general = ShannonEntropy.entropy([p, 1 - p])
        assert approx_eq(h_binary, h_general)


# ===========================================================================
# TESTS: KL-DIVERGENZ UND INFORMATIONSMASSE
# ===========================================================================

class TestKLDivergence:
    """Tests für die Klasse KLDivergence."""

    def test_kl_zero_for_equal_distributions(self):
        """D_KL(P||P) = 0."""
        p = [0.3, 0.5, 0.2]
        assert approx_eq(KLDivergence.kl_divergence(p, p), 0.0)

    def test_kl_nonnegative(self):
        """D_KL(P||Q) ≥ 0 (Gibb'sche Ungleichung)."""
        p = [0.4, 0.6]
        q = [0.3, 0.7]
        assert KLDivergence.kl_divergence(p, q) >= 0

    def test_kl_asymmetric(self):
        """D_KL(P||Q) ≠ D_KL(Q||P) im Allgemeinen."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        kl_pq = KLDivergence.kl_divergence(p, q)
        kl_qp = KLDivergence.kl_divergence(q, p)
        assert not approx_eq(kl_pq, kl_qp)

    def test_kl_with_zero_p_allowed(self):
        """p_i = 0 ist erlaubt (0·log(0) = 0)."""
        p = [0.5, 0.5, 0.0]
        q = [0.4, 0.4, 0.2]
        result = KLDivergence.kl_divergence(p, q)
        assert result >= 0

    def test_kl_zero_q_with_positive_p_raises(self):
        """q_i = 0 mit p_i > 0 löst ValueError aus."""
        p = [0.5, 0.5]
        q = [0.5, 0.0]
        with pytest.raises(ValueError):
            KLDivergence.kl_divergence(p, q)

    def test_kl_different_lengths_raises(self):
        """Verschiedene Längen lösen ValueError aus."""
        with pytest.raises(ValueError):
            KLDivergence.kl_divergence([0.5, 0.5], [1.0])

    def test_js_divergence_symmetric(self):
        """JSD ist symmetrisch: JSD(P||Q) = JSD(Q||P)."""
        p = [0.3, 0.7]
        q = [0.6, 0.4]
        assert approx_eq(
            KLDivergence.js_divergence(p, q),
            KLDivergence.js_divergence(q, p)
        )

    def test_js_divergence_bounded(self):
        """JSD ∈ [0, 1] für Basis 2."""
        p = [0.5, 0.5]
        q = [1.0, 0.0]
        js = KLDivergence.js_divergence(p, q)
        assert 0 <= js <= 1.0 + 1e-9

    def test_js_divergence_zero_for_equal(self):
        """JSD(P||P) = 0."""
        p = [0.3, 0.4, 0.3]
        assert approx_eq(KLDivergence.js_divergence(p, p), 0.0)

    def test_cross_entropy_relation(self):
        """H(P,Q) = H(P) + D_KL(P||Q)."""
        p = [0.4, 0.6]
        q = [0.3, 0.7]
        h_p = ShannonEntropy.entropy(p)
        kl = KLDivergence.kl_divergence(p, q)
        ce = KLDivergence.cross_entropy(p, q)
        assert approx_eq(ce, h_p + kl, tol=1e-5)

    def test_cross_entropy_minimum_at_p_equals_q(self):
        """H(P,Q) ≥ H(P) mit Gleichheit bei P = Q."""
        p = [0.4, 0.6]
        ce_pp = KLDivergence.cross_entropy(p, p)
        h_p = ShannonEntropy.entropy(p)
        assert approx_eq(ce_pp, h_p)

    def test_renyi_entropy_alpha_one_equals_shannon(self):
        """H_1(X) = Shannon-Entropie (Grenzwert)."""
        probs = [0.3, 0.5, 0.2]
        h_shannon = ShannonEntropy.entropy(probs)
        h_renyi = KLDivergence.renyi_entropy(probs, alpha=1.0)
        assert approx_eq(h_renyi, h_shannon, tol=1e-5)

    def test_renyi_entropy_decreasing_in_alpha(self):
        """H_α(X) ist monoton fallend in α."""
        probs = [0.3, 0.5, 0.2]
        h_05 = KLDivergence.renyi_entropy(probs, alpha=0.5)
        h_2 = KLDivergence.renyi_entropy(probs, alpha=2.0)
        assert h_05 >= h_2

    def test_renyi_entropy_alpha_zero_hartley(self):
        """H_0(X) = Hartley-Entropie = log₂(|Träger|)."""
        probs = [0.3, 0.5, 0.2]
        h_0 = KLDivergence.renyi_entropy(probs, alpha=0.0)
        expected = math.log2(3)  # 3 Symbole mit p > 0
        assert approx_eq(h_0, expected)

    def test_tsallis_entropy_q_one_equals_shannon(self):
        """S_1(X) = Shannon-Entropie (Grenzwert q→1)."""
        probs = [0.4, 0.6]
        s_tsallis = KLDivergence.tsallis_entropy(probs, q=1.0)
        h_shannon = ShannonEntropy.entropy(probs)
        assert approx_eq(s_tsallis, h_shannon, tol=1e-5)

    def test_tsallis_entropy_nonnegative(self):
        """Tsallis-Entropie ist ≥ 0 für gültige Verteilungen."""
        probs = [0.3, 0.5, 0.2]
        for q in [0.5, 2.0, 3.0]:
            s = KLDivergence.tsallis_entropy(probs, q=q)
            assert s >= -1e-9

    def test_total_variation_zero_for_equal(self):
        """TV(P,P) = 0."""
        p = [0.3, 0.4, 0.3]
        assert approx_eq(KLDivergence.total_variation_distance(p, p), 0.0)

    def test_total_variation_bounded(self):
        """TV ∈ [0, 1]."""
        p = [1.0, 0.0]
        q = [0.0, 1.0]
        tv = KLDivergence.total_variation_distance(p, q)
        assert approx_eq(tv, 1.0)

    def test_total_variation_symmetric(self):
        """TV(P,Q) = TV(Q,P)."""
        p = [0.4, 0.6]
        q = [0.7, 0.3]
        assert approx_eq(
            KLDivergence.total_variation_distance(p, q),
            KLDivergence.total_variation_distance(q, p)
        )

    def test_total_variation_different_lengths_raises(self):
        """Verschiedene Längen lösen ValueError aus."""
        with pytest.raises(ValueError):
            KLDivergence.total_variation_distance([0.5, 0.5], [1.0])


# ===========================================================================
# TESTS: KANALKAPAZITÄT
# ===========================================================================

class TestChannelCapacity:
    """Tests für die Klasse ChannelCapacity."""

    def test_bsc_perfect_channel(self):
        """BSC mit p=0: Kapazität = 1 Bit."""
        c = ChannelCapacity.binary_symmetric_channel_capacity(0.0)
        assert approx_eq(c, 1.0)

    def test_bsc_complete_noise(self):
        """BSC mit p=0.5: Kapazität = 0 Bit."""
        c = ChannelCapacity.binary_symmetric_channel_capacity(0.5)
        assert approx_eq(c, 0.0)

    def test_bsc_inverse_channel(self):
        """BSC mit p=1: Kapazität = 1 Bit (inverser Kanal)."""
        c = ChannelCapacity.binary_symmetric_channel_capacity(1.0)
        assert approx_eq(c, 1.0)

    def test_bsc_symmetric_around_half(self):
        """BSC-Kapazität ist symmetrisch um p=0.5."""
        for p in [0.1, 0.2, 0.3]:
            assert approx_eq(
                ChannelCapacity.binary_symmetric_channel_capacity(p),
                ChannelCapacity.binary_symmetric_channel_capacity(1 - p)
            )

    def test_bsc_invalid_p_raises(self):
        """p außerhalb [0,1] löst ValueError aus."""
        with pytest.raises(ValueError):
            ChannelCapacity.binary_symmetric_channel_capacity(1.5)

    def test_bec_no_erasure(self):
        """BEC mit ε=0: Kapazität = 1 Bit."""
        c = ChannelCapacity.binary_erasure_channel_capacity(0.0)
        assert approx_eq(c, 1.0)

    def test_bec_full_erasure(self):
        """BEC mit ε=1: Kapazität = 0 Bit."""
        c = ChannelCapacity.binary_erasure_channel_capacity(1.0)
        assert approx_eq(c, 0.0)

    def test_bec_half_erasure(self):
        """BEC mit ε=0.5: Kapazität = 0.5 Bit."""
        c = ChannelCapacity.binary_erasure_channel_capacity(0.5)
        assert approx_eq(c, 0.5)

    def test_bec_linear(self):
        """BEC-Kapazität ist linear in ε."""
        for eps in [0.1, 0.3, 0.7]:
            c = ChannelCapacity.binary_erasure_channel_capacity(eps)
            assert approx_eq(c, 1 - eps)

    def test_awgn_zero_snr_db(self):
        """AWGN mit 0 dB SNR (SNR_lin=1): C = log₂(2) = 1 Bit."""
        c = ChannelCapacity.awgn_channel_capacity(0.0)
        assert approx_eq(c, 1.0)

    def test_awgn_high_snr(self):
        """AWGN mit hohem SNR nähert sich log₂(SNR)."""
        snr_db = 30  # SNR = 1000
        c = ChannelCapacity.awgn_channel_capacity(snr_db)
        snr_lin = 10 ** (snr_db / 10)
        expected = math.log2(1 + snr_lin)
        assert approx_eq(c, expected)

    def test_awgn_capacity_increases_with_snr(self):
        """Höheres SNR → größere Kapazität."""
        c1 = ChannelCapacity.awgn_channel_capacity(0.0)
        c2 = ChannelCapacity.awgn_channel_capacity(10.0)
        assert c2 > c1

    def test_blahut_arimoto_bsc(self):
        """Blahut-Arimoto auf BSC mit p=0.1 → C ≈ 1-H(0.1)."""
        p = 0.1
        P = [[1 - p, p], [p, 1 - p]]
        capacity, p_opt = ChannelCapacity.blahut_arimoto(P, n_iter=500)
        expected = ChannelCapacity.binary_symmetric_channel_capacity(p)
        assert approx_eq(capacity, expected, tol=1e-3)

    def test_blahut_arimoto_optimal_distribution(self):
        """Blahut-Arimoto findet Gleichverteilung als optimal für BSC."""
        p = 0.1
        P = [[1 - p, p], [p, 1 - p]]
        _, p_opt = ChannelCapacity.blahut_arimoto(P, n_iter=500)
        assert approx_eq(p_opt[0], p_opt[1], tol=1e-3)

    def test_channel_capacity_perfect_channel(self):
        """Perfekter Kanal (Einheitsmatrix) → C = log₂(n)."""
        P = [[1.0, 0.0], [0.0, 1.0]]
        c = ChannelCapacity.channel_capacity_general(P)
        assert approx_eq(c, 1.0, tol=1e-3)

    def test_noisy_channel_theorem_achievable(self):
        """Rate < Kapazität → erreichbar."""
        result = ChannelCapacity.noisy_channel_coding_theorem_demo(0.9, 0.5)
        assert result["achievable"] is True

    def test_noisy_channel_theorem_not_achievable(self):
        """Rate > Kapazität → nicht erreichbar."""
        result = ChannelCapacity.noisy_channel_coding_theorem_demo(0.5, 0.9)
        assert result["achievable"] is False

    def test_noisy_channel_theorem_at_capacity(self):
        """Rate = Kapazität → nicht erreichbar (strikt < nötig)."""
        result = ChannelCapacity.noisy_channel_coding_theorem_demo(0.7, 0.7)
        assert result["achievable"] is False


# ===========================================================================
# TESTS: QUELLENCODIERUNG
# ===========================================================================

class TestSourceCoding:
    """Tests für die Klasse SourceCoding."""

    def test_huffman_single_symbol(self):
        """Einzelnes Symbol → Code '0'."""
        codes = SourceCoding.huffman_code(['a'], [1.0])
        assert 'a' in codes
        assert codes['a'] == '0'

    def test_huffman_two_symbols_equal(self):
        """Zwei gleich wahrscheinliche Symbole → Codes der Länge 1."""
        codes = SourceCoding.huffman_code(['a', 'b'], [0.5, 0.5])
        assert len(codes['a']) == 1
        assert len(codes['b']) == 1

    def test_huffman_prefix_free(self):
        """Huffman-Code ist präfixfrei."""
        symbols = ['a', 'b', 'c', 'd']
        probs = [0.4, 0.3, 0.2, 0.1]
        codes = SourceCoding.huffman_code(symbols, probs)
        code_list = list(codes.values())
        for i, ci in enumerate(code_list):
            for j, cj in enumerate(code_list):
                if i != j:
                    assert not cj.startswith(ci), f"Präfix verletzt: {ci} ist Präfix von {cj}"

    def test_huffman_covers_all_symbols(self):
        """Jedes Symbol hat einen Code."""
        symbols = ['a', 'b', 'c', 'd', 'e']
        probs = [0.35, 0.25, 0.20, 0.15, 0.05]
        codes = SourceCoding.huffman_code(symbols, probs)
        assert set(codes.keys()) == set(symbols)

    def test_huffman_empty_returns_empty(self):
        """Leere Symbolliste → leeres Dictionary."""
        codes = SourceCoding.huffman_code([], [])
        assert codes == {}

    def test_huffman_length_mismatch_raises(self):
        """Verschiedene Längen → ValueError."""
        with pytest.raises(ValueError):
            SourceCoding.huffman_code(['a', 'b'], [0.5])

    def test_huffman_average_length_bounds(self):
        """Mittlere Codelänge erfüllt H(X) ≤ L̄ < H(X) + 1."""
        symbols = ['a', 'b', 'c', 'd']
        probs = [0.4, 0.3, 0.2, 0.1]
        codes = SourceCoding.huffman_code(symbols, probs)
        avg_len = SourceCoding.huffman_average_length(codes, probs)
        h = ShannonEntropy.entropy(probs)
        assert h <= avg_len + 1e-9
        assert avg_len < h + 1.0

    def test_huffman_optimal_for_power_of_two(self):
        """Für 2^n gleich wahrscheinliche Symbole: L̄ = H(X) = n."""
        probs = [0.25, 0.25, 0.25, 0.25]
        codes = SourceCoding.huffman_code(['a', 'b', 'c', 'd'], probs)
        avg_len = SourceCoding.huffman_average_length(codes, probs)
        h = ShannonEntropy.entropy(probs)
        assert approx_eq(avg_len, h)

    def test_arithmetic_coding_single_symbol(self):
        """Arithmetisches Codieren einer einzigen Nachricht."""
        result = SourceCoding.arithmetic_coding_demo(['a', 'b'], [0.5, 0.5], ['a'])
        assert 0.0 <= result['low'] < result['high'] <= 1.0

    def test_arithmetic_coding_empty_message(self):
        """Leere Nachricht → Intervall [0, 1)."""
        result = SourceCoding.arithmetic_coding_demo(['a'], [1.0], [])
        assert approx_eq(result['low'], 0.0)
        assert approx_eq(result['high'], 1.0)

    def test_arithmetic_coding_interval_shrinks(self):
        """Längere Nachricht → kleineres Intervall."""
        syms = ['a', 'b']
        probs = [0.7, 0.3]
        r1 = SourceCoding.arithmetic_coding_demo(syms, probs, ['a'])
        r2 = SourceCoding.arithmetic_coding_demo(syms, probs, ['a', 'b'])
        assert r2['interval_width'] < r1['interval_width']

    def test_lempel_ziv_empty_sequence(self):
        """Leere Sequenz → Komplexität 0."""
        assert SourceCoding.lempel_ziv_complexity([]) == 0

    def test_lempel_ziv_constant_low_complexity(self):
        """Konstante Sequenz hat niedrige Komplexität."""
        c_const = SourceCoding.lempel_ziv_complexity([0] * 100)
        c_varied = SourceCoding.lempel_ziv_complexity(list(range(50)) * 2)
        assert c_const <= c_varied

    def test_lempel_ziv_single_element(self):
        """Einzelnes Element → Komplexität 1."""
        assert SourceCoding.lempel_ziv_complexity([0]) == 1

    def test_source_coding_bound_structure(self):
        """Quellencodierungssatz gibt korrekte Struktur zurück."""
        probs = [0.4, 0.6]
        result = SourceCoding.source_coding_theorem_bound(probs)
        assert 'entropy' in result
        assert 'lower_bound' in result
        assert 'upper_bound' in result
        assert approx_eq(result['lower_bound'], result['entropy'])
        assert approx_eq(result['upper_bound'], result['entropy'] + 1.0)

    def test_kolmogorov_complexity_demo(self):
        """Kolmogorov-Demo gibt sinnvolle Struktur zurück."""
        result = SourceCoding.kolmogorov_complexity_demo(100)
        assert result['n'] == 100
        # Konstante Folge ist kompressibler
        assert result['all_zeros']['compressed_bytes'] < result['all_zeros']['raw_bytes']


# ===========================================================================
# TESTS: FEHLERKORREKTUR
# ===========================================================================

class TestErrorCorrection:
    """Tests für die Klasse ErrorCorrection."""

    def test_hamming_distance_identical(self):
        """Identische Codewörter → Distanz 0."""
        c = [1, 0, 1, 1, 0]
        assert ErrorCorrection.hamming_distance(c, c) == 0

    def test_hamming_distance_one_bit(self):
        """Ein Bit Unterschied → Distanz 1."""
        c1 = [1, 0, 1]
        c2 = [1, 1, 1]
        assert ErrorCorrection.hamming_distance(c1, c2) == 1

    def test_hamming_distance_all_different(self):
        """Alle Bits verschieden → Distanz = n."""
        c1 = [0, 0, 0, 0]
        c2 = [1, 1, 1, 1]
        assert ErrorCorrection.hamming_distance(c1, c2) == 4

    def test_hamming_distance_symmetric(self):
        """d(c1,c2) = d(c2,c1)."""
        c1 = [1, 0, 1, 0]
        c2 = [0, 1, 0, 1]
        assert (ErrorCorrection.hamming_distance(c1, c2) ==
                ErrorCorrection.hamming_distance(c2, c1))

    def test_hamming_distance_length_mismatch_raises(self):
        """Verschiedene Längen → ValueError."""
        with pytest.raises(ValueError):
            ErrorCorrection.hamming_distance([1, 0], [1, 0, 1])

    def test_hamming_weight_all_zeros(self):
        """Nullvektor → Gewicht 0."""
        assert ErrorCorrection.hamming_weight([0, 0, 0, 0]) == 0

    def test_hamming_weight_all_ones(self):
        """Einservektor → Gewicht = n."""
        assert ErrorCorrection.hamming_weight([1, 1, 1, 1]) == 4

    def test_hamming_weight_relation_to_distance(self):
        """w(c) = d(c, 0...0)."""
        c = [1, 0, 1, 1, 0]
        zeros = [0] * len(c)
        assert ErrorCorrection.hamming_weight(c) == ErrorCorrection.hamming_distance(c, zeros)

    def test_minimum_distance_two_codewords(self):
        """Minimaldistanz von zwei Codewörtern."""
        code = [[1, 0, 0], [0, 1, 1]]
        assert ErrorCorrection.minimum_distance(code) == ErrorCorrection.hamming_distance(
            [1, 0, 0], [0, 1, 1]
        )

    def test_minimum_distance_single_codeword_raises(self):
        """Weniger als 2 Codewörter → ValueError."""
        with pytest.raises(ValueError):
            ErrorCorrection.minimum_distance([[1, 0]])

    def test_singleton_bound_basic(self):
        """Singleton-Schranke: d ≤ n-k+1."""
        assert ErrorCorrection.singleton_bound(7, 4) == 4  # n=7, k=4 → d≤4

    def test_singleton_bound_mds(self):
        """Reed-Solomon: n=7, k=3 → d=5 (MDS)."""
        assert ErrorCorrection.singleton_bound(7, 3) == 5

    def test_singleton_bound_k_equals_n_raises(self):
        """k > n löst ValueError aus."""
        with pytest.raises(ValueError):
            ErrorCorrection.singleton_bound(3, 5)

    def test_hamming_bound_n7_t1(self):
        """Hamming-Schranke für (7,*,3)-Code: Σ C(7,i) i=0..1 = 1+7 = 8."""
        assert ErrorCorrection.hamming_bound(7, 1) == 8

    def test_hamming_bound_t0(self):
        """t=0 (kein Fehlerkorrektur): Kugelvolumen = 1."""
        assert ErrorCorrection.hamming_bound(5, 0) == 1

    def test_hamming_code_74_structure(self):
        """(7,4)-Hamming-Code hat korrekte Parameter."""
        result = ErrorCorrection.hamming_code_74()
        assert result['n'] == 7
        assert result['k'] == 4
        assert result['d_min'] == 3
        assert result['is_perfect'] is True

    def test_hamming_code_74_generator_shape(self):
        """Generatormatrix hat Form 4×7."""
        result = ErrorCorrection.hamming_code_74()
        G = result['generator_matrix']
        assert G.shape == (4, 7)

    def test_hamming_code_74_parity_check_shape(self):
        """Prüfmatrix hat Form 3×7."""
        result = ErrorCorrection.hamming_code_74()
        H = result['parity_check_matrix']
        assert H.shape == (3, 7)

    def test_hamming_code_74_orthogonality(self):
        """H · G^T = 0 (mod 2)."""
        result = ErrorCorrection.hamming_code_74()
        G = result['generator_matrix']
        H = result['parity_check_matrix']
        product = (H @ G.T) % 2
        assert np.all(product == 0)

    def test_parity_check_code_structure(self):
        """Paritätscode der Länge 4: n=4, k=3, d_min=2."""
        result = ErrorCorrection.parity_check_code(4)
        assert result['n'] == 4
        assert result['k'] == 3
        assert result['d_min'] == 2

    def test_parity_check_code_invalid_raises(self):
        """n < 2 löst ValueError aus."""
        with pytest.raises(ValueError):
            ErrorCorrection.parity_check_code(1)

    def test_repetition_code_length(self):
        """Wiederholungscode: Ausgabelänge = k · |Nachricht|."""
        msg = [1, 0, 1]
        encoded = ErrorCorrection.repetition_code(msg, 3)
        assert len(encoded) == 9

    def test_repetition_code_content(self):
        """Jedes Bit wird k-mal wiederholt."""
        msg = [1, 0]
        encoded = ErrorCorrection.repetition_code(msg, 3)
        assert encoded == [1, 1, 1, 0, 0, 0]

    def test_repetition_code_k1_identity(self):
        """k=1: kein Wiederholungscode → identisch."""
        msg = [1, 0, 1]
        assert ErrorCorrection.repetition_code(msg, 1) == msg

    def test_repetition_code_invalid_k_raises(self):
        """k < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ErrorCorrection.repetition_code([1], 0)


# ===========================================================================
# TESTS: DIFFERENTIELLE ENTROPIE
# ===========================================================================

class TestDifferentialEntropy:
    """Tests für die Klasse DifferentialEntropy."""

    def test_normal_entropy_sigma_one(self):
        """Normalverteilung σ=1: h = (1/2)·log₂(2πe)."""
        h = DifferentialEntropy.differential_entropy_normal(1.0)
        expected = 0.5 * math.log2(2 * math.pi * math.e)
        assert approx_eq(h, expected)

    def test_normal_entropy_increases_with_sigma(self):
        """Größere Varianz → größere Entropie."""
        h1 = DifferentialEntropy.differential_entropy_normal(1.0)
        h2 = DifferentialEntropy.differential_entropy_normal(2.0)
        assert h2 > h1

    def test_normal_entropy_invalid_sigma_raises(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            DifferentialEntropy.differential_entropy_normal(0.0)
        with pytest.raises(ValueError):
            DifferentialEntropy.differential_entropy_normal(-1.0)

    def test_exponential_entropy_lambda_one(self):
        """Exp(1): h = 1 - log₂(1) = 1 Bit."""
        h = DifferentialEntropy.differential_entropy_exponential(1.0)
        # h_nat = 1 - ln(1) = 1, h_bit = 1/ln(2)
        expected = 1.0 / math.log(2)
        assert approx_eq(h, expected)

    def test_exponential_entropy_decreases_with_lambda(self):
        """Größeres λ → kleinere Entropie (weniger Unsicherheit)."""
        h1 = DifferentialEntropy.differential_entropy_exponential(1.0)
        h2 = DifferentialEntropy.differential_entropy_exponential(2.0)
        assert h1 > h2

    def test_exponential_entropy_invalid_lambda_raises(self):
        """λ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            DifferentialEntropy.differential_entropy_exponential(0.0)

    def test_uniform_entropy_unit_interval(self):
        """U(0,1): h = log₂(1) = 0 Bit."""
        h = DifferentialEntropy.differential_entropy_uniform(0.0, 1.0)
        assert approx_eq(h, 0.0)

    def test_uniform_entropy_wider_interval(self):
        """U(0,4): h = log₂(4) = 2 Bit."""
        h = DifferentialEntropy.differential_entropy_uniform(0.0, 4.0)
        assert approx_eq(h, 2.0)

    def test_uniform_entropy_can_be_negative(self):
        """U(0,0.5): h = log₂(0.5) = -1 Bit (negativ ist erlaubt)."""
        h = DifferentialEntropy.differential_entropy_uniform(0.0, 0.5)
        assert approx_eq(h, -1.0)

    def test_uniform_entropy_invalid_raises(self):
        """a ≥ b löst ValueError aus."""
        with pytest.raises(ValueError):
            DifferentialEntropy.differential_entropy_uniform(1.0, 0.5)

    def test_gaussian_mi_zero_correlation(self):
        """Unkorreliertheit → I(X;Y) = 0."""
        mi = DifferentialEntropy.mutual_information_gaussian(0.0)
        assert approx_eq(mi, 0.0)

    def test_gaussian_mi_positive(self):
        """Positive Korrelation → I(X;Y) > 0."""
        mi = DifferentialEntropy.mutual_information_gaussian(0.5)
        assert mi > 0

    def test_gaussian_mi_increases_with_correlation(self):
        """Höhere Korrelation → mehr gegenseitige Information."""
        mi1 = DifferentialEntropy.mutual_information_gaussian(0.3)
        mi2 = DifferentialEntropy.mutual_information_gaussian(0.9)
        assert mi2 > mi1

    def test_gaussian_mi_symmetric(self):
        """I(X;Y) ist symmetrisch in ρ und -ρ."""
        assert approx_eq(
            DifferentialEntropy.mutual_information_gaussian(0.5),
            DifferentialEntropy.mutual_information_gaussian(-0.5)
        )

    def test_gaussian_mi_invalid_rho_raises(self):
        """|ρ| ≥ 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            DifferentialEntropy.mutual_information_gaussian(1.0)
        with pytest.raises(ValueError):
            DifferentialEntropy.mutual_information_gaussian(-1.0)


# ===========================================================================
# TESTS: STANDALONE-FUNKTIONEN
# ===========================================================================

class TestStandaloneFunctions:
    """Tests für die Standalone-Funktionen."""

    def test_surprisal_certain_event(self):
        """Sicheres Ereignis p=1: I = 0 Bit."""
        assert approx_eq(surprisal(1.0), 0.0)

    def test_surprisal_half(self):
        """p=0.5: I = 1 Bit."""
        assert approx_eq(surprisal(0.5), 1.0)

    def test_surprisal_quarter(self):
        """p=0.25: I = 2 Bit."""
        assert approx_eq(surprisal(0.25), 2.0)

    def test_surprisal_nonnegative(self):
        """Informationsgehalt ist immer ≥ 0."""
        for p in [0.1, 0.5, 0.9, 1.0]:
            assert surprisal(p) >= 0

    def test_surprisal_zero_prob_raises(self):
        """p = 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            surprisal(0.0)

    def test_surprisal_out_of_range_raises(self):
        """p > 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            surprisal(1.5)

    def test_surprisal_base_e(self):
        """I in Nat (Basis e)."""
        i = surprisal(math.e ** (-1), base=math.e)
        assert approx_eq(i, 1.0)

    def test_entropy_rate_markov_two_state(self):
        """Entropierate einer 2-Zustands-Kette mit symmetrischer Matrix."""
        P = [[0.5, 0.5], [0.5, 0.5]]
        pi = [0.5, 0.5]
        h = entropy_rate_markov(P, pi)
        # Jeder Übergang hat p=0.5 → H = 1 Bit pro Schritt
        assert approx_eq(h, 1.0)

    def test_entropy_rate_markov_deterministic(self):
        """Deterministische Kette → Entropierate 0."""
        P = [[1.0, 0.0], [0.0, 1.0]]
        pi = [0.5, 0.5]
        h = entropy_rate_markov(P, pi)
        assert approx_eq(h, 0.0)

    def test_entropy_rate_markov_nonnegative(self):
        """Entropierate ist immer ≥ 0."""
        P = [[0.7, 0.3], [0.4, 0.6]]
        pi = [0.4 / 0.7, 0.3 / 0.7]  # approximative stationäre Verteilung
        # Normalisieren
        total = sum(pi)
        pi = [x / total for x in pi]
        h = entropy_rate_markov(P, pi)
        assert h >= -1e-9

    def test_data_processing_inequality_demo_structure(self):
        """DPI-Demo gibt korrekte Struktur zurück."""
        result = data_processing_inequality_demo()
        assert 'I_XY' in result
        assert 'I_XZ' in result
        assert 'dpi_satisfied' in result

    def test_data_processing_inequality_satisfied(self):
        """DPI: I(X;Z) ≤ I(X;Y) muss gelten."""
        result = data_processing_inequality_demo()
        assert result['dpi_satisfied'] is True

    def test_data_processing_information_loss(self):
        """Information geht durch weitere Verarbeitung verloren."""
        result = data_processing_inequality_demo()
        assert result['information_loss'] >= -1e-9


# ===========================================================================
# ERWEITERTE EDGE-CASE-TESTS
# ===========================================================================

class TestEdgeCases:
    """Tests für Grenzfälle und spezielle Situationen."""

    def test_entropy_three_symbols_manual(self):
        """Manuelle Berechnung für 3-Symbol-Verteilung."""
        p = [0.5, 0.25, 0.25]
        h = ShannonEntropy.entropy(p)
        expected = -0.5 * math.log2(0.5) - 0.25 * math.log2(0.25) - 0.25 * math.log2(0.25)
        assert approx_eq(h, expected)

    def test_entropy_equals_entropy_maximum_for_uniform(self):
        """Gleichverteilung erreicht maximale Entropie."""
        for n in [2, 4, 8, 16]:
            probs = [1/n] * n
            h = ShannonEntropy.entropy(probs)
            h_max = ShannonEntropy.entropy_maximum(n)
            assert approx_eq(h, h_max)

    def test_kl_divergence_manual_calculation(self):
        """Manuelle KL-Berechnung für 2-Symbol-Fall."""
        p = [0.4, 0.6]
        q = [0.5, 0.5]
        expected = 0.4 * math.log2(0.4 / 0.5) + 0.6 * math.log2(0.6 / 0.5)
        result = KLDivergence.kl_divergence(p, q)
        assert approx_eq(result, expected)

    def test_huffman_skewed_distribution(self):
        """Stark schiefe Verteilung erzeugt kurzen Code für häufigstes Symbol."""
        symbols = ['a', 'b', 'c', 'd']
        probs = [0.9, 0.05, 0.03, 0.02]
        codes = SourceCoding.huffman_code(symbols, probs)
        # Das häufigste Symbol sollte den kürzesten Code haben
        len_a = len(codes['a'])
        for sym in ['b', 'c', 'd']:
            assert len_a <= len(codes[sym])

    def test_awgn_negative_snr_db(self):
        """Negatives SNR in dB → Kapazität < 1."""
        c = ChannelCapacity.awgn_channel_capacity(-10.0)
        assert c < 1.0
        assert c > 0.0

    def test_lempel_ziv_periodic_lower_than_random(self):
        """Periodische Sequenz hat niedrigere LZ-Komplexität als diverse."""
        periodic = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        diverse = list(range(12))
        c_periodic = SourceCoding.lempel_ziv_complexity(periodic)
        c_diverse = SourceCoding.lempel_ziv_complexity(diverse)
        assert c_periodic <= c_diverse

    def test_bsc_capacity_formula_matches(self):
        """BSC-Kapazität stimmt mit Formel C = 1 - H(p) überein."""
        for p in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
            c = ChannelCapacity.binary_symmetric_channel_capacity(p)
            h = ShannonEntropy.entropy_binary(p)
            assert approx_eq(c, 1.0 - h)

    def test_mutual_information_bounded_by_marginals(self):
        """I(X;Y) ≤ min(H(X), H(Y))."""
        joint = [[0.3, 0.1], [0.2, 0.4]]
        mi = ShannonEntropy.mutual_information(joint)
        h_x = ShannonEntropy.entropy([0.4, 0.6])
        h_y = ShannonEntropy.entropy([0.5, 0.5])
        assert mi <= min(h_x, h_y) + 1e-9

    def test_differential_entropy_normal_matches_formula(self):
        """h(N(0,σ²)) = (1/2)log₂(2πeσ²) für verschiedene σ."""
        for sigma in [0.5, 1.0, 2.0, 5.0]:
            h = DifferentialEntropy.differential_entropy_normal(sigma)
            expected = 0.5 * math.log2(2 * math.pi * math.e * sigma ** 2)
            assert approx_eq(h, expected)

    def test_surprisal_is_expectation_of_entropy(self):
        """Erwartungswert von I(X) = H(X)."""
        probs = [0.4, 0.3, 0.2, 0.1]
        expected_surprisal = sum(p * surprisal(p) for p in probs if p > 0)
        h = ShannonEntropy.entropy(probs)
        assert approx_eq(expected_surprisal, h)
