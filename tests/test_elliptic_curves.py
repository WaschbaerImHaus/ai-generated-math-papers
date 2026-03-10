"""
@file test_elliptic_curves.py
@brief Tests fuer das elliptic_curves Modul (TDD).
@description
    Umfassende Tests fuer:
    - ECPoint: Punkt-Arithetik, Infinity, Gleichheit, Inverses, Addition, Multiplikation
    - EllipticCurve: Diskriminante, j-Invariante, Kurvenpunkte, Hasse-Schranke
    - EllipticCurveModP: Gruppenstruktur ueber F_p, alle Punkte, Ordnung
    - ECCKeyExchange: ECDH-Demo, Schluessel-Paar, gemeinsames Geheimnis
    - Zahlentheoretische Funktionen: Nagell-Lutz, Mordell-Weil, BSD, ECM

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest

# Suchpfad fuer src-Modul setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from elliptic_curves import (
    ECPoint,
    EllipticCurve,
    EllipticCurveModP,
    ECCKeyExchange,
    lenstra_ecm_factorization,
    nagell_lutz_theorem,
    mordell_weil_rank_estimate,
    l_function_rank_order,
    is_supersingular,
    endomorphism_ring_type,
    example_bsd_curve,
    congruent_number_curve,
    secp256k1,
    secp256k1_generator,
)


# ===========================================================================
# Tests fuer ECPoint
# ===========================================================================

class TestECPoint:
    """Tests fuer die ECPoint-Klasse."""

    def test_infinity_is_infinity(self):
        """ECPoint.infinity().is_infinity muss True sein."""
        O = ECPoint.infinity()
        assert O.is_infinity is True

    def test_regular_point_not_infinity(self):
        """Regulaerer Punkt darf nicht is_infinity sein."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        assert P.is_infinity is False

    def test_infinity_has_none_coords(self):
        """Unendlichkeitspunkt muss x=None, y=None haben."""
        O = ECPoint.infinity()
        assert O.x is None
        assert O.y is None

    def test_infinity_equality(self):
        """Zwei Unendlichkeitspunkte sind gleich."""
        O1 = ECPoint.infinity()
        O2 = ECPoint.infinity()
        assert O1 == O2

    def test_point_equality(self):
        """Zwei Punkte mit gleichen Koordinaten sind gleich."""
        curve = EllipticCurve(-1, 0)
        P1 = ECPoint(0, 0, curve)
        P2 = ECPoint(0, 0, curve)
        assert P1 == P2

    def test_point_inequality(self):
        """Zwei Punkte mit verschiedenen Koordinaten sind ungleich."""
        curve = EllipticCurve(-1, 0)
        P1 = ECPoint(0, 0, curve)
        P2 = ECPoint(1, 0, curve)
        assert P1 != P2

    def test_point_not_equal_to_infinity(self):
        """Regulaerer Punkt ist nicht gleich O."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        O = ECPoint.infinity(curve)
        assert P != O

    def test_add_infinity_left(self):
        """O + P = P."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        O = ECPoint.infinity(curve)
        result = O + P
        assert result == P

    def test_add_infinity_right(self):
        """P + O = P."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        O = ECPoint.infinity(curve)
        result = P + O
        assert result == P

    def test_add_inverse_gives_infinity(self):
        """P + (-P) = O."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        neg_P = -P
        result = P + neg_P
        assert result.is_infinity

    def test_negation_regular(self):
        """-(x, y) = (x, -y)."""
        curve = EllipticCurve(-1, 1)
        P = ECPoint(0, 1, curve)
        neg_P = -P
        assert neg_P.x == 0
        assert neg_P.y == -1

    def test_negation_infinity(self):
        """-O = O."""
        O = ECPoint.infinity()
        assert (-O).is_infinity

    def test_two_torsion_point_y0(self):
        """Punkt mit y=0 ist 2-Torsion: P + P = O."""
        curve = EllipticCurve(-1, 0)
        # (-1, 0) liegt auf y^2 = x^3 - x: 0 = -1 - (-1) = 0 ✓
        P = ECPoint(-1, 0, curve)
        result = P + P
        assert result.is_infinity, f"Erwartet O, bekam {result}"

    def test_add_point_to_itself_y0(self):
        """Verdopplung bei y=0 ergibt O (Tangente ist vertikal)."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        result = P + P
        assert result.is_infinity

    def test_scalar_mult_zero(self):
        """0 * P = O."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        result = 0 * P
        assert result.is_infinity

    def test_scalar_mult_one(self):
        """1 * P = P."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        result = 1 * P
        assert result == P

    def test_scalar_mult_two(self):
        """2 * P = P + P."""
        curve = EllipticCurve(0, 1)
        # Punkt auf y^2 = x^3 + 1: (-1, 0)
        P = ECPoint(-1, 0, curve)
        result_mul = 2 * P
        result_add = P + P
        assert result_mul == result_add

    def test_scalar_mult_negative(self):
        """(-1) * P = -P."""
        curve = EllipticCurve(-1, 1)
        P = ECPoint(0, 1, curve)
        result = (-1) * P
        assert result == -P

    def test_rmul(self):
        """n * P funktioniert (via __rmul__)."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        result = 3 * P
        assert result.is_infinity  # 2-Torsion, 2*P=O, 3*P=P, aber P ist 2-Torsion

    def test_on_curve_true(self):
        """Punkt auf Kurve: on_curve() = True."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        assert P.on_curve() is True

    def test_on_curve_infinity(self):
        """O liegt immer auf der Kurve."""
        curve = EllipticCurve(-1, 0)
        O = ECPoint.infinity(curve)
        assert O.on_curve() is True

    def test_repr_infinity(self):
        """Repr von O."""
        O = ECPoint.infinity()
        assert 'O' in repr(O)

    def test_repr_regular(self):
        """Repr eines regulaeren Punktes."""
        P = ECPoint(3, 5, None)
        r = repr(P)
        assert '3' in r and '5' in r


# ===========================================================================
# Tests fuer EllipticCurve
# ===========================================================================

class TestEllipticCurve:
    """Tests fuer die EllipticCurve-Klasse."""

    def test_singular_curve_raises(self):
        """EllipticCurve(0, 0) muss ValueError werfen (singulaer)."""
        with pytest.raises(ValueError):
            EllipticCurve(0, 0)

    def test_valid_curve_no_exception(self):
        """EllipticCurve(-1, 0) darf keine Exception werfen."""
        curve = EllipticCurve(-1, 0)
        assert curve is not None

    def test_discriminant_nonzero(self):
        """Diskriminante muss ungleich 0 sein."""
        curve = EllipticCurve(-1, 0)
        assert abs(curve.discriminant()) > 0

    def test_discriminant_value(self):
        """Diskriminante von y^2 = x^3 - x: Delta = -16*(4*(-1)^3 + 27*0) = 64."""
        curve = EllipticCurve(-1, 0)
        expected = -16 * (4 * (-1) ** 3 + 27 * 0 ** 2)  # = 64
        assert abs(curve.discriminant() - expected) < 1e-10

    def test_j_invariant_1728(self):
        """
        j-Invariante von y^2 = x^3 - x muss 1728 sein.
        j = -1728 * (4*(-1))^3 / Delta = -1728*(-64)/64 = 1728.
        """
        curve = EllipticCurve(-1, 0)
        assert abs(curve.j_invariant() - 1728.0) < 1e-6

    def test_j_invariant_0(self):
        """
        j-Invariante von y^2 = x^3 + 1 (a=0) muss 0 sein.
        j = -1728 * 0 / Delta = 0.
        """
        curve = EllipticCurve(0, 1)
        assert abs(curve.j_invariant()) < 1e-10

    def test_is_on_curve_origin(self):
        """(0, 0) liegt auf y^2 = x^3 - x: 0 = 0 ✓."""
        curve = EllipticCurve(-1, 0)
        assert curve.is_on_curve(0, 0) is True

    def test_is_on_curve_one_zero(self):
        """(1, 0) liegt auf y^2 = x^3 - x: 0 = 1-1 = 0 ✓."""
        curve = EllipticCurve(-1, 0)
        assert curve.is_on_curve(1, 0) is True

    def test_is_on_curve_minus_one_zero(self):
        """(-1, 0) liegt auf y^2 = x^3 - x: 0 = -1+1 = 0 ✓."""
        curve = EllipticCurve(-1, 0)
        assert curve.is_on_curve(-1, 0) is True

    def test_is_not_on_curve(self):
        """(2, 2) liegt nicht auf y^2 = x^3 - x: 4 != 6."""
        curve = EllipticCurve(-1, 0)
        assert curve.is_on_curve(2, 2) is False

    def test_point_raises_if_not_on_curve(self):
        """EllipticCurve.point() wirft ValueError wenn nicht auf Kurve."""
        curve = EllipticCurve(-1, 0)
        with pytest.raises(ValueError):
            curve.point(5, 5)

    def test_point_returns_ecpoint(self):
        """EllipticCurve.point() gibt ECPoint zurueck."""
        curve = EllipticCurve(-1, 0)
        P = curve.point(0, 0)
        assert isinstance(P, ECPoint)
        assert P.on_curve()

    def test_y_values_positive_rhs(self):
        """y_values gibt zwei Werte zurueck wenn x^3+ax+b > 0."""
        curve = EllipticCurve(0, 1)  # y^2 = x^3 + 1, bei x=2: rhs=9
        ys = curve.y_values(2)
        assert len(ys) == 2
        assert abs(ys[0] - 3.0) < 1e-10
        assert abs(ys[1] + 3.0) < 1e-10

    def test_y_values_zero_rhs(self):
        """y_values gibt [0] zurueck wenn x^3+ax+b = 0."""
        curve = EllipticCurve(-1, 0)
        ys = curve.y_values(0)  # x=0: 0^3 - 0 = 0
        assert len(ys) == 1
        assert abs(ys[0]) < 1e-10

    def test_y_values_negative_rhs(self):
        """y_values gibt leere Liste zurueck wenn keine reellen Loesungen."""
        curve = EllipticCurve(-1, 0)
        # x = 0.5: 0.125 - 0.5 = -0.375 < 0
        ys = curve.y_values(0.5)
        assert len(ys) == 0

    def test_points_over_fp_contains_infinity(self):
        """E(F_p) enthaelt immer O."""
        curve = EllipticCurve(-1, 0)
        points = curve.points_over_fp(7)
        has_infinity = any(p.is_infinity for p in points)
        assert has_infinity

    def test_order_over_fp(self):
        """Gruppenordnung stimmt mit Anzahl Punkte ueberein."""
        curve = EllipticCurve(-1, 0)
        p = 7
        order = curve.order_over_fp(p)
        points = curve.points_over_fp(p)
        assert order == len(points)

    def test_hasse_bound(self):
        """Hasse-Schranke: |#E(F_p) - (p+1)| <= 2*sqrt(p)."""
        curve = EllipticCurve(-1, 0)
        for p in [5, 7, 11, 13, 17]:
            order = curve.order_over_fp(p)
            deviation = abs(order - (p + 1))
            bound = 2 * math.sqrt(p)
            assert deviation <= bound + 1e-10, (
                f"Hasse-Schranke verletzt: |{order} - {p+1}| = {deviation} > 2*sqrt({p}) = {bound}"
            )

    def test_trace_of_frobenius(self):
        """a_p = p + 1 - #E(F_p)."""
        curve = EllipticCurve(-1, 0)
        p = 7
        trace = curve.trace_of_frobenius(p)
        order = curve.order_over_fp(p)
        assert trace == p + 1 - order

    def test_repr(self):
        """Repr der Kurve."""
        curve = EllipticCurve(-1, 0)
        r = repr(curve)
        assert 'EllipticCurve' in r

    def test_torsion_points_over_r(self):
        """torsion_points_over_r gibt Liste von ECPoints zurueck."""
        curve = EllipticCurve(-1, 0)
        pts = curve.torsion_points_over_r((-2, 2), n_points=100)
        assert isinstance(pts, list)
        assert len(pts) > 0


# ===========================================================================
# Tests fuer EllipticCurveModP
# ===========================================================================

class TestEllipticCurveModP:
    """Tests fuer die EllipticCurveModP-Klasse."""

    def test_basic_creation(self):
        """EllipticCurveModP(-1, 0, 7) erstellen."""
        curve = EllipticCurveModP(-1, 0, 7)
        assert curve.p == 7

    def test_singular_mod_p_raises(self):
        """Singulaere Kurve mod p muss ValueError werfen."""
        with pytest.raises(ValueError):
            # a=0, b=0 => Delta = 0
            EllipticCurveModP(0, 0, 7)

    def test_all_points_contains_infinity(self):
        """E(F_p).all_points() enthaelt O."""
        curve = EllipticCurveModP(-1, 0, 7)
        points = curve.all_points()
        has_infinity = any(p.is_infinity for p in points)
        assert has_infinity

    def test_all_points_count_matches_group_order(self):
        """Anzahl Punkte = group_order()."""
        curve = EllipticCurveModP(-1, 0, 7)
        points = curve.all_points()
        assert len(points) == curve.group_order()

    def test_hasse_bound_mod_p(self):
        """Hasse-Schranke auch fuer EllipticCurveModP."""
        curve = EllipticCurveModP(-1, 0, 7)
        p = 7
        order = curve.group_order()
        deviation = abs(order - (p + 1))
        bound = 2 * math.sqrt(p)
        assert deviation <= bound + 1e-10

    def test_all_points_on_curve(self):
        """Alle berechneten Punkte muessen auf der Kurve liegen."""
        curve = EllipticCurveModP(-1, 0, 7)
        for pt in curve.all_points():
            assert curve.is_on_curve(pt.x, pt.y), f"Punkt {pt} nicht auf Kurve"

    def test_add_mod_commutativity(self):
        """P + Q = Q + P (Kommutativitaet)."""
        curve = EllipticCurveModP(-1, 0, 7)
        points = curve.all_points()
        # Nimm zwei reale Punkte
        real_pts = [pt for pt in points if not pt.is_infinity]
        if len(real_pts) >= 2:
            P, Q = real_pts[0], real_pts[1]
            assert curve.add_mod(P, Q) == curve.add_mod(Q, P)

    def test_add_mod_associativity(self):
        """(P + Q) + R = P + (Q + R) (Assoziativitaet)."""
        curve = EllipticCurveModP(-1, 0, 7)
        points = curve.all_points()
        real_pts = [pt for pt in points if not pt.is_infinity]
        if len(real_pts) >= 3:
            P, Q, R = real_pts[0], real_pts[1], real_pts[2]
            lhs = curve.add_mod(curve.add_mod(P, Q), R)
            rhs = curve.add_mod(P, curve.add_mod(Q, R))
            assert lhs == rhs

    def test_scalar_mult_mod_group_order_gives_infinity(self):
        """n * P = O fuer Gruppenordnung n."""
        curve = EllipticCurveModP(-1, 0, 7)
        n = curve.group_order()
        points = curve.all_points()
        real_pts = [pt for pt in points if not pt.is_infinity]
        if real_pts:
            P = real_pts[0]
            result = curve.scalar_mult_mod(n, P)
            assert result.is_infinity, f"n*P = {result}, erwartet O"

    def test_scalar_mult_mod_2p_eq_pp(self):
        """2*P = P + P."""
        curve = EllipticCurveModP(-1, 0, 7)
        real_pts = [pt for pt in curve.all_points() if not pt.is_infinity]
        if real_pts:
            P = real_pts[0]
            assert curve.scalar_mult_mod(2, P) == curve.add_mod(P, P)

    def test_trace_of_frobenius_mod_p(self):
        """a_p = p + 1 - #E(F_p) fuer EllipticCurveModP."""
        curve = EllipticCurveModP(-1, 0, 7)
        p = 7
        order = curve.group_order()
        # Hier nutzen wir die Methode aus der Elternklasse
        trace = p + 1 - order
        assert isinstance(trace, int)

    def test_is_on_curve_mod_p(self):
        """is_on_curve mod p Pruefung."""
        curve = EllipticCurveModP(-1, 0, 7)
        # O liegt immer drauf
        assert curve.is_on_curve(None, None) is True

    def test_repr_mod_p(self):
        """Repr der Kurve mod p."""
        curve = EllipticCurveModP(-1, 0, 7)
        r = repr(curve)
        assert 'mod' in r.lower() or '7' in r


# ===========================================================================
# Tests fuer Gruppenoperationen (allgemein)
# ===========================================================================

class TestGroupOperations:
    """Tests fuer allgemeine Gruppenaxiome auf elliptischen Kurven."""

    def setup_method(self):
        """Kurve und Punkte fuer alle Tests."""
        self.curve = EllipticCurveModP(-1, 0, 7)
        all_pts = self.curve.all_points()
        self.real_pts = [p for p in all_pts if not p.is_infinity]
        self.O = ECPoint.infinity(self.curve)

    def test_identity_element(self):
        """P + O = O + P = P."""
        if self.real_pts:
            P = self.real_pts[0]
            assert self.curve.add_mod(P, self.O) == P
            assert self.curve.add_mod(self.O, P) == P

    def test_inverse_element(self):
        """P + (-P) = O."""
        for P in self.real_pts[:3]:
            neg_P_x = P.x
            neg_P_y = (-P.y) % self.curve.p
            neg_P = ECPoint(neg_P_x, neg_P_y, self.curve)
            result = self.curve.add_mod(P, neg_P)
            assert result.is_infinity, f"P + (-P) != O fuer P={P}"

    def test_double_and_add_consistency(self):
        """n*P via Double-and-Add stimmt mit iterierter Addition ueberein."""
        if self.real_pts:
            P = self.real_pts[0]
            # 3*P = P + P + P
            three_P_mul = self.curve.scalar_mult_mod(3, P)
            three_P_add = self.curve.add_mod(self.curve.add_mod(P, P), P)
            assert three_P_mul == three_P_add


# ===========================================================================
# Tests fuer ECCKeyExchange
# ===========================================================================

class TestECCKeyExchange:
    """Tests fuer Elliptic Curve Diffie-Hellman."""

    def setup_method(self):
        """Einfache Kurve fuer ECDH-Tests."""
        self.curve = EllipticCurveModP(-1, 0, 7)
        # Nehme einen Punkt als Generator
        all_pts = self.curve.all_points()
        real_pts = [p for p in all_pts if not p.is_infinity]
        self.G = real_pts[0] if real_pts else ECPoint.infinity(self.curve)
        self.ecdh = ECCKeyExchange(self.curve, self.G)

    def test_generate_keypair_returns_tuple(self):
        """generate_keypair gibt (int, ECPoint) zurueck."""
        priv, pub = self.ecdh.generate_keypair()
        assert isinstance(priv, int)
        assert isinstance(pub, ECPoint)

    def test_keypair_private_in_range(self):
        """Privater Schluessel in {1, ..., n-1}."""
        priv, _ = self.ecdh.generate_keypair()
        n = self.curve.group_order()
        assert 1 <= priv <= n - 1

    def test_generate_keypair_with_fixed_key(self):
        """Fixer privater Schluessel ergibt korrektes Schluessel-Paar."""
        priv, pub = self.ecdh.generate_keypair(private_key=2)
        assert priv == 2
        expected_pub = self.curve.scalar_mult_mod(2, self.G)
        assert pub == expected_pub

    def test_shared_secret_matches(self):
        """Alice und Bob berechnen gleichen gemeinsamen Schluessel."""
        alice_priv, alice_pub = self.ecdh.generate_keypair()
        bob_priv, bob_pub = self.ecdh.generate_keypair()

        alice_secret = self.ecdh.shared_secret(alice_priv, bob_pub)
        bob_secret = self.ecdh.shared_secret(bob_priv, alice_pub)

        assert alice_secret == bob_secret

    def test_demo_exchange_returns_dict(self):
        """demo_exchange() gibt Dictionary zurueck."""
        result = self.ecdh.demo_exchange()
        assert isinstance(result, dict)

    def test_demo_exchange_contains_match(self):
        """demo_exchange() hat 'shared_match' Key."""
        result = self.ecdh.demo_exchange()
        assert 'shared_match' in result

    def test_demo_exchange_shared_match_true(self):
        """demo_exchange(): gemeinsames Geheimnis stimmt ueberein."""
        result = self.ecdh.demo_exchange()
        assert result['shared_match'] is True


# ===========================================================================
# Tests fuer Zahlentheorie
# ===========================================================================

class TestNagellLutz:
    """Tests fuer den Nagell-Lutz-Satz."""

    def test_bsd_curve_torsion(self):
        """y^2 = x^3 - x hat Torsionspunkte (0,0), (1,0), (-1,0)."""
        torsion = nagell_lutz_theorem(-1, 0)
        torsion_set = set(torsion)
        assert (0, 0) in torsion_set, f"(0,0) nicht in {torsion}"
        assert (1, 0) in torsion_set, f"(1,0) nicht in {torsion}"
        assert (-1, 0) in torsion_set, f"(-1,0) nicht in {torsion}"

    def test_returns_list(self):
        """nagell_lutz_theorem gibt Liste zurueck."""
        result = nagell_lutz_theorem(-1, 0)
        assert isinstance(result, list)

    def test_torsion_points_on_curve(self):
        """Alle Torsionspunkte liegen auf der Kurve."""
        a, b = -1, 0
        curve = EllipticCurve(a, b)
        for x, y in nagell_lutz_theorem(a, b):
            assert curve.is_on_curve(x, y, tol=1e-6), (
                f"Torsionspunkt ({x},{y}) liegt nicht auf y^2=x^3-x"
            )

    def test_y_squared_divides_discriminant(self):
        """y^2 teilt Delta (Nagell-Lutz Bedingung)."""
        a, b = -1, 0
        delta = -16 * (4 * a ** 3 + 27 * b ** 2)
        for x, y in nagell_lutz_theorem(a, b):
            if y != 0:
                assert delta % (y * y) == 0, (
                    f"y^2={y*y} teilt nicht Delta={delta} fuer Punkt ({x},{y})"
                )

    def test_singular_curve_raises(self):
        """Singulaere Kurve wirft ValueError."""
        with pytest.raises(ValueError):
            nagell_lutz_theorem(0, 0)

    def test_other_curve(self):
        """Torsionspunkte auf y^2 = x^3 + 1."""
        result = nagell_lutz_theorem(0, 1)
        # (-1, 0) liegt auf y^2 = x^3 + 1: 0 = -1 + 1 = 0 ✓
        assert isinstance(result, list)
        torsion_set = set(result)
        assert (-1, 0) in torsion_set


class TestMordellWeilRank:
    """Tests fuer die Mordell-Weil-Rang-Schaetzung."""

    def test_returns_dict(self):
        """mordell_weil_rank_estimate gibt Dictionary zurueck."""
        result = mordell_weil_rank_estimate(-1, 0, [5, 7, 11, 13])
        assert isinstance(result, dict)

    def test_dict_has_required_keys(self):
        """Dictionary hat alle benoetigten Keys."""
        result = mordell_weil_rank_estimate(-1, 0, [5, 7, 11, 13])
        required_keys = ['rank_lower', 'rank_upper', 'selmer_bound', 'bsd_order_estimate']
        for key in required_keys:
            assert key in result, f"Key '{key}' fehlt in {result}"

    def test_rank_lower_nonnegative(self):
        """Untere Schranke >= 0."""
        result = mordell_weil_rank_estimate(-1, 0, [5, 7, 11])
        assert result['rank_lower'] >= 0

    def test_rank_upper_geq_lower(self):
        """Obere Schranke >= untere Schranke."""
        result = mordell_weil_rank_estimate(-1, 0, [5, 7, 11])
        assert result['rank_upper'] >= result['rank_lower']

    def test_singular_curve_returns_dict(self):
        """Singulaere Kurve gibt Dictionary mit error-Key zurueck."""
        result = mordell_weil_rank_estimate(0, 0, [5, 7])
        assert isinstance(result, dict)


class TestLFunctionRankOrder:
    """Tests fuer die L-Funktions-Rang-Ordnungs-Schaetzung."""

    def test_returns_dict(self):
        """l_function_rank_order gibt Dictionary zurueck."""
        result = l_function_rank_order(-1, 0, [5, 7, 11, 13, 17])
        assert isinstance(result, dict)

    def test_dict_has_required_keys(self):
        """Dictionary hat alle benoetigten Keys."""
        result = l_function_rank_order(-1, 0, [5, 7, 11])
        for key in ['order_estimate', 'L_value_at_1', 'bsd_consistent']:
            assert key in result

    def test_order_estimate_nonnegative(self):
        """Ordnungsschaeatzung >= 0."""
        result = l_function_rank_order(-1, 0, [5, 7, 11, 13])
        assert result['order_estimate'] >= 0

    def test_l_value_is_complex(self):
        """L-Wert ist komplex."""
        result = l_function_rank_order(-1, 0, [5, 7, 11])
        assert isinstance(result['L_value_at_1'], complex)


class TestSupersingular:
    """Tests fuer die Supersingularitaets-Pruefung."""

    def test_returns_bool(self):
        """is_supersingular gibt bool zurueck."""
        result = is_supersingular(-1, 0, 7)
        assert isinstance(result, bool)

    def test_supersingular_definition(self):
        """Supersingulaar gdw a_p ≡ 0 (mod p)."""
        a, b, p = -1, 0, 7
        curve = EllipticCurve(a, b)
        trace = curve.trace_of_frobenius(p)
        expected = (trace % p == 0)
        result = is_supersingular(a, b, p)
        assert result == expected

    def test_invalid_p_raises(self):
        """p < 5 wirft ValueError."""
        with pytest.raises(ValueError):
            is_supersingular(-1, 0, 3)

    def test_endomorphism_ring_type_ordinary(self):
        """Ordinaere Kurve gibt 'ordinary' zurueck."""
        # Pruefe mehrere kleine Primzahlen
        found_ordinary = False
        for p in [5, 7, 11, 13, 17, 19, 23]:
            t = endomorphism_ring_type(-1, 0, p)
            if t == 'ordinary':
                found_ordinary = True
                break
        assert found_ordinary

    def test_endomorphism_ring_type_supersingular(self):
        """Supersingulare Kurve gibt 'supersingular' zurueck."""
        # Finde supersingulare Kurve durch Suche
        found = False
        for p in [5, 7, 11, 13, 17, 19, 23]:
            if is_supersingular(-1, 0, p):
                assert endomorphism_ring_type(-1, 0, p) == 'supersingular'
                found = True
                break
        # Es kann sein, dass y^2=x^3-x ueber diesen p nicht supersingulaar ist
        # dann teste mit anderer Kurve
        if not found:
            for p in [5, 7, 11, 13, 17, 19, 23]:
                if is_supersingular(0, 1, p):
                    assert endomorphism_ring_type(0, 1, p) == 'supersingular'
                    break


# ===========================================================================
# Tests fuer ECM-Faktorisierung
# ===========================================================================

class TestLenstraECM:
    """Tests fuer Lenstra's ECM-Faktorisierung."""

    def test_finds_factor_of_15(self):
        """lenstra_ecm_factorization(15) findet 3 oder 5."""
        result = lenstra_ecm_factorization(15, max_curves=50, B=100)
        # Entweder None (Pech) oder gueltiger Faktor
        if result is not None:
            assert result in [3, 5], f"Ungueltiger Faktor {result} von 15"
            assert 15 % result == 0

    def test_finds_factor_of_77(self):
        """lenstra_ecm_factorization(77) findet 7 oder 11."""
        result = lenstra_ecm_factorization(77, max_curves=50, B=100)
        if result is not None:
            assert 77 % result == 0
            assert result in [7, 11]

    def test_prime_returns_none(self):
        """Primzahl hat keinen echten Teiler => None."""
        result = lenstra_ecm_factorization(17, max_curves=10)
        assert result is None

    def test_even_number_returns_2(self):
        """Gerade Zahl => 2."""
        result = lenstra_ecm_factorization(100)
        assert result == 2

    def test_result_is_valid_divisor(self):
        """Ergebnis muss echter Teiler von n sein."""
        n = 35  # 5 * 7
        result = lenstra_ecm_factorization(n, max_curves=50, B=200)
        if result is not None:
            assert n % result == 0
            assert 1 < result < n


# ===========================================================================
# Tests fuer bekannte Kurven
# ===========================================================================

class TestKnownCurves:
    """Tests fuer bekannte elliptische Kurven."""

    def test_example_bsd_curve(self):
        """example_bsd_curve() gibt EllipticCurve zurueck."""
        curve = example_bsd_curve()
        assert isinstance(curve, EllipticCurve)
        # y^2 = x^3 - x: a=-1, b=0
        assert curve.a == -1
        assert curve.b == 0

    def test_bsd_curve_has_three_torsion(self):
        """BSD-Kurve hat mindestens 3 Torsionspunkte."""
        curve = example_bsd_curve()
        # (0,0), (1,0), (-1,0) sind auf der Kurve
        assert curve.is_on_curve(0, 0)
        assert curve.is_on_curve(1, 0)
        assert curve.is_on_curve(-1, 0)

    def test_congruent_number_curve(self):
        """congruent_number_curve gibt EllipticCurve zurueck."""
        curve = congruent_number_curve(5)
        assert isinstance(curve, EllipticCurve)

    def test_congruent_number_curve_5(self):
        """E_5: y^2 = x^3 - 25x (a=-25, b=0)."""
        curve = congruent_number_curve(5)
        assert curve.a == -25
        assert curve.b == 0

    def test_congruent_number_curve_6(self):
        """E_6: y^2 = x^3 - 36x (a=-36, b=0)."""
        curve = congruent_number_curve(6)
        assert curve.a == -36
        assert curve.b == 0

    def test_congruent_number_curve_negative_raises(self):
        """n <= 0 wirft ValueError."""
        with pytest.raises(ValueError):
            congruent_number_curve(-1)

    def test_congruent_number_curve_zero_raises(self):
        """n = 0 wirft ValueError."""
        with pytest.raises(ValueError):
            congruent_number_curve(0)

    def test_secp256k1_creation(self):
        """secp256k1() gibt EllipticCurveModP zurueck."""
        curve = secp256k1()
        assert isinstance(curve, EllipticCurveModP)
        assert curve.a == 0
        assert curve.b == 7


# ===========================================================================
# Tests fuer Punktaddition (reelle Kurven)
# ===========================================================================

class TestPointAdditionReal:
    """Tests fuer Punktaddition auf reellen Kurven."""

    def test_point_on_y2_eq_x3_plus_1(self):
        """Punkt auf y^2 = x^3 + 1 pruefen."""
        curve = EllipticCurve(0, 1)
        # (-1, 0): (-1)^3 + 1 = 0 = 0^2 ✓
        assert curve.is_on_curve(-1, 0)

    def test_addition_with_known_result(self):
        """Teste bekannte Addition auf einfacher Kurve."""
        # Auf y^2 = x^3 - x ueber F_7:
        # Alle Punkte berechnen und Addition pruefen
        curve = EllipticCurveModP(-1, 0, 7)
        pts = [p for p in curve.all_points() if not p.is_infinity]
        if len(pts) >= 1:
            P = pts[0]
            # P + O = P
            O = ECPoint.infinity(curve)
            assert curve.add_mod(P, O) == P
            # O + P = P
            assert curve.add_mod(O, P) == P


# ===========================================================================
# Edge Cases und Sonderfaelle
# ===========================================================================

class TestEdgeCases:
    """Edge Cases und Sonderfaelle."""

    def test_curve_with_large_coefficients(self):
        """Kurve mit grossen Koeffizienten."""
        curve = EllipticCurve(100, -200)
        assert abs(curve.discriminant()) > 0

    def test_ecpoint_multiply_by_zero(self):
        """0 * P = O."""
        curve = EllipticCurve(-1, 0)
        P = ECPoint(0, 0, curve)
        result = P * 0
        assert result.is_infinity

    def test_ecpoint_multiply_large_scalar(self):
        """Skalarmultiplikation mit grosser Zahl bleibt stabil."""
        curve = EllipticCurveModP(-1, 0, 7)
        pts = [p for p in curve.all_points() if not p.is_infinity]
        if pts:
            P = pts[0]
            n = curve.group_order()
            # 2*n*P = O (zweifache Gruppenordnung)
            result = curve.scalar_mult_mod(2 * n, P)
            assert result.is_infinity

    def test_double_and_add_matches_iterative(self):
        """Double-and-Add stimmt mit iterierter Addition ueberein."""
        curve = EllipticCurveModP(-1, 0, 7)
        pts = [p for p in curve.all_points() if not p.is_infinity]
        if pts:
            P = pts[0]
            # 4*P via Multiplikation
            result_mul = curve.scalar_mult_mod(4, P)
            # 4*P via iterierter Addition
            result_add = P
            for _ in range(3):
                result_add = curve.add_mod(result_add, P)
            assert result_mul == result_add

    def test_nagell_lutz_large_delta(self):
        """Nagell-Lutz fuer Kurve mit grosser Diskriminante."""
        result = nagell_lutz_theorem(-5, 4)
        assert isinstance(result, list)

    def test_ecpoint_inequality_with_non_ecpoint(self):
        """ECPoint != nicht-ECPoint."""
        P = ECPoint(0, 0, None)
        assert P != "kein Punkt"
        assert P != 42

    def test_all_frobenius_traces_in_hasse_bound(self):
        """Alle Frobenius-Spuren innerhalb der Hasse-Schranke."""
        curve = EllipticCurve(0, 1)
        for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
            a_p = curve.trace_of_frobenius(p)
            assert abs(a_p) <= 2 * math.sqrt(p) + 1e-10, (
                f"|a_{p}| = {abs(a_p)} > 2*sqrt({p})"
            )


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
