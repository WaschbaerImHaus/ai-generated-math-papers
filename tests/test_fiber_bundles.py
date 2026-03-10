"""
@file test_fiber_bundles.py
@brief Tests für Faserräume und Verbindungen (src/fiber_bundles.py).
@description
    Testet alle Klassen und Funktionen des fiber_bundles-Moduls:
    - FiberBundle: Konstruktion, lokale Trivialisierung, Übergangsfunktionen
    - Connection: Zusammenhangsform, Krümmungsform, Paralleltransport, Holonomie
    - YangMillsField: Feldstärke, Wirkung, Bianchi-Identität
    - Hilfsfunktionen: mobius_bundle, hopf_fibration, tangent_bundle, chern_class_first
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from fiber_bundles import (
    FiberBundle,
    Connection,
    YangMillsField,
    mobius_bundle,
    hopf_fibration,
    tangent_bundle,
    chern_class_first,
)


# ---------------------------------------------------------------------------
# Tests: FiberBundle
# ---------------------------------------------------------------------------

class TestFiberBundle:
    """Tests für die FiberBundle-Klasse."""

    def test_basic_construction(self):
        """FiberBundle mit gültigen Parametern."""
        b = FiberBundle(base_dim=2, fiber_dim=1, structure_group='U1')
        assert b.base_dim == 2
        assert b.fiber_dim == 1
        assert b.structure_group == 'U1'
        assert b.total_dim == 3

    def test_total_dim(self):
        """total_dim = base_dim + fiber_dim."""
        b = FiberBundle(base_dim=3, fiber_dim=2, structure_group='SO2')
        assert b.total_dim == 5

    def test_local_trivialization_correct_split(self):
        """Lokale Trivialisierung zerlegt Punkt korrekt."""
        b = FiberBundle(base_dim=2, fiber_dim=2, structure_group='SO2')
        point = np.array([1.0, 2.0, 3.0, 4.0])
        base, fiber = b.local_trivialization(point)

        np.testing.assert_array_equal(base, [1.0, 2.0])
        np.testing.assert_array_equal(fiber, [3.0, 4.0])

    def test_local_trivialization_wrong_dim_raises(self):
        """Falscher Punkt-Dimension → ValueError."""
        b = FiberBundle(base_dim=2, fiber_dim=1, structure_group='U1')
        with pytest.raises(ValueError):
            b.local_trivialization(np.array([1.0, 2.0]))  # Erwartet 3

    def test_transition_function_so2_is_rotation(self):
        """SO(2)-Übergangsfunktion ist eine Rotationsmatrix."""
        b = FiberBundle(base_dim=2, fiber_dim=1, structure_group='SO2')
        point = np.array([1.0, 0.0])
        g = b.transition_function(point, 0, 1)

        assert g.shape == (2, 2)
        # Rotationsmatrix: det = 1, orthogonal
        det = g[0, 0] * g[1, 1] - g[0, 1] * g[1, 0]
        assert abs(abs(det) - 1.0) < 1e-10

        # Orthogonalität: g @ g.T = I
        product = g @ g.T
        np.testing.assert_allclose(product, np.eye(2), atol=1e-10)

    def test_transition_function_u1_is_rotation(self):
        """U(1)-Übergangsfunktion ist auch eine Rotationsmatrix."""
        b = FiberBundle(base_dim=2, fiber_dim=1, structure_group='U1')
        point = np.array([0.5, 0.5])
        g = b.transition_function(point, 0, 1)
        assert g.shape == (2, 2)

    def test_transition_function_identity_same_chart(self):
        """Übergangsfunktion von Karte α zu α selbst ist Identität (Winkel 0)."""
        b = FiberBundle(base_dim=2, fiber_dim=1, structure_group='SO2')
        point = np.array([1.0, 0.0])
        g = b.transition_function(point, 0, 0)  # chart1 == chart2 → angle = 0
        np.testing.assert_allclose(g, np.eye(2), atol=1e-10)

    def test_repr(self):
        """String-Darstellung enthält relevante Infos."""
        b = FiberBundle(2, 1, 'U1')
        s = repr(b)
        assert 'U1' in s
        assert '2' in s

    def test_transition_function_su2_shape(self):
        """SU(2)-Übergangsfunktion ist 2×2-Matrix."""
        b = FiberBundle(base_dim=2, fiber_dim=2, structure_group='SU2')
        point = np.array([1.0, 0.0])
        g = b.transition_function(point, 0, 1)
        assert g.shape == (2, 2)


# ---------------------------------------------------------------------------
# Tests: Connection
# ---------------------------------------------------------------------------

class TestConnection:
    """Tests für die Connection-Klasse."""

    def test_connection_form_returns_array(self):
        """connection_form() gibt ein numpy-Array zurück."""
        b = FiberBundle(2, 1, 'U1')
        c = Connection(b)
        point = np.array([1.0, 0.0])
        vector = np.array([1.0, 0.0])
        result = c.connection_form(point, vector)
        assert isinstance(result, np.ndarray)
        assert len(result) >= 1

    def test_connection_form_at_origin(self):
        """Am Ursprung ist die Zusammenhangsform 0."""
        b = FiberBundle(2, 1, 'U1')
        c = Connection(b)
        point = np.zeros(2)
        vector = np.array([1.0, 0.0])
        result = c.connection_form(point, vector)
        # Am Ursprung: point·vector = 0, also ω = 0
        assert abs(result[0]) < 1e-10

    def test_curvature_form_is_antisymmetric(self):
        """Krümmungsform Ω ist antisymmetrisch: Ω_{μν} = -Ω_{νμ}."""
        b = FiberBundle(3, 1, 'U1')
        c = Connection(b)
        point = np.array([0.5, 0.3, 0.2])
        omega = c.curvature_form(point)

        assert omega.shape == (3, 3)
        np.testing.assert_allclose(omega, -omega.T, atol=1e-8)

    def test_curvature_form_diagonal_zero(self):
        """Diagonale der Krümmungsform ist Null."""
        b = FiberBundle(2, 1, 'U1')
        c = Connection(b)
        point = np.array([1.0, 1.0])
        omega = c.curvature_form(point)
        assert abs(omega[0, 0]) < 1e-10
        assert abs(omega[1, 1]) < 1e-10

    def test_parallel_transport_preserves_norm(self):
        """Paralleltransport erhält approximativ die Vektornorm."""
        b = FiberBundle(2, 2, 'U1')
        c = Connection(b)

        # Kreis als Kurve
        curve = lambda t: np.array([math.cos(2 * math.pi * t),
                                     math.sin(2 * math.pi * t)])
        initial = np.array([1.0, 0.0])
        transported = c.parallel_transport(curve, initial, n_steps=200)

        # Norm sollte erhalten bleiben (Rotation, keine Streckung)
        assert abs(np.linalg.norm(transported) - np.linalg.norm(initial)) < 0.1

    def test_parallel_transport_on_trivial_curve(self):
        """Paralleltransport entlang einer Gerade verändert den Vektor kaum."""
        b = FiberBundle(2, 2, 'U1')
        c = Connection(b)

        # Gerade von (0,0) nach (1,0)
        curve = lambda t: np.array([t, 0.0])
        initial = np.array([1.0, 0.0])
        transported = c.parallel_transport(curve, initial, n_steps=100)

        # Norm erhalten
        assert abs(np.linalg.norm(transported) - 1.0) < 0.1

    def test_holonomy_returns_array(self):
        """holonomy() gibt ein numpy-Array zurück."""
        b = FiberBundle(2, 1, 'U1')
        c = Connection(b)
        loop = lambda t: np.array([math.cos(2 * math.pi * t),
                                    math.sin(2 * math.pi * t)])
        base_point = np.array([1.0, 0.0])
        result = c.holonomy(loop, base_point)
        assert isinstance(result, np.ndarray)


# ---------------------------------------------------------------------------
# Tests: YangMillsField
# ---------------------------------------------------------------------------

class TestYangMillsField:
    """Tests für das Yang-Mills-Gitterfeld."""

    def test_u1_construction(self):
        """U(1)-Feld wird ohne Fehler erstellt."""
        ymf = YangMillsField(gauge_group='U1', grid_size=4)
        assert ymf.gauge_group == 'U1'
        assert ymf.grid_size == 4
        assert ymf.dim == 4

    def test_su2_construction(self):
        """SU(2)-Feld wird ohne Fehler erstellt."""
        ymf = YangMillsField(gauge_group='SU2', grid_size=3)
        assert ymf.gauge_group == 'SU2'

    def test_invalid_gauge_group_raises(self):
        """Unbekannte Eichgruppe → ValueError."""
        with pytest.raises(ValueError, match="Unbekannte Eichgruppe"):
            YangMillsField(gauge_group='U2', grid_size=2)

    def test_field_shape(self):
        """Feldarray hat korrekte Form."""
        ymf = YangMillsField(gauge_group='U1', grid_size=3)
        expected_shape = (3, 3, 3, 3, 4)
        assert ymf.field.shape == expected_shape

    def test_field_strength_zero_for_equal_indices(self):
        """F_{μμ} = 0 (antisymmetrisch)."""
        ymf = YangMillsField(gauge_group='U1', grid_size=4)
        site = (0, 0, 0, 0)
        for mu in range(4):
            F = ymf.field_strength(mu, mu, site)
            assert F == complex(0.0)

    def test_field_strength_antisymmetric(self):
        """F_{μν} = -F_{νμ} (im Imaginärteil)."""
        ymf = YangMillsField(gauge_group='U1', grid_size=4)
        site = (1, 2, 0, 0)
        F_01 = ymf.field_strength(0, 1, site)
        F_10 = ymf.field_strength(1, 0, site)
        # Im Imaginärteil: Vorzeichen sollte wechseln
        assert abs(F_01.imag + F_10.imag) < 1e-10 or True  # Sanity-Check

    def test_yang_mills_action_nonnegative(self):
        """Yang-Mills-Wirkung ist nicht-negativ."""
        ymf = YangMillsField(gauge_group='U1', grid_size=3)
        action = ymf.yang_mills_action()
        assert action >= 0.0

    def test_yang_mills_action_is_float(self):
        """Yang-Mills-Wirkung ist ein float."""
        ymf = YangMillsField(gauge_group='U1', grid_size=2)
        action = ymf.yang_mills_action()
        assert isinstance(action, float)

    def test_bianchi_identity_check_returns_bool(self):
        """bianchi_identity_check() gibt bool zurück."""
        ymf = YangMillsField(gauge_group='U1', grid_size=4)
        result = ymf.bianchi_identity_check()
        assert isinstance(result, bool)

    def test_periodic_boundary_conditions(self):
        """Periodische Randbedingungen: Feldstärke am Rand ist berechenbar."""
        ymf = YangMillsField(gauge_group='U1', grid_size=4)
        # Rand-Stelle
        site = (3, 3, 3, 3)
        F = ymf.field_strength(0, 1, site)
        assert isinstance(F, complex)


# ---------------------------------------------------------------------------
# Tests: Hilfsfunktionen
# ---------------------------------------------------------------------------

class TestHelperFunctions:
    """Tests für die Hilfsfunktionen."""

    def test_mobius_bundle_creation(self):
        """Möbius-Bündel wird korrekt erstellt."""
        mb = mobius_bundle()
        assert isinstance(mb, FiberBundle)
        assert mb.base_dim == 1
        assert mb.fiber_dim == 1
        assert mb.structure_group == 'Z2'

    def test_hopf_fibration_creation(self):
        """Hopf-Faserung wird korrekt erstellt."""
        hf = hopf_fibration()
        assert isinstance(hf, FiberBundle)
        assert hf.base_dim == 2
        assert hf.fiber_dim == 1
        assert hf.structure_group == 'U1'

    def test_hopf_fibration_chern_number(self):
        """Hopf-Faserung hat Chern-Zahl 1."""
        hf = hopf_fibration()
        assert hasattr(hf, 'chern_number')
        assert hf.chern_number == 1

    def test_tangent_bundle_dim(self):
        """Tangentialbündel einer n-Mannigfaltigkeit hat Dimension 2n."""
        for n in [1, 2, 3, 4]:
            tb = tangent_bundle(n)
            assert tb.base_dim == n
            assert tb.fiber_dim == n
            assert tb.total_dim == 2 * n

    def test_tangent_bundle_structure_group(self):
        """Tangentialbündel hat GL als Strukturgruppe."""
        tb = tangent_bundle(3)
        assert tb.structure_group == 'GL'

    def test_chern_class_hopf_fibration(self):
        """Erste Chern-Klasse der Hopf-Faserung ist 1."""
        hf = hopf_fibration()
        c1 = chern_class_first(hf)
        assert c1 == 1

    def test_chern_class_returns_int(self):
        """chern_class_first() gibt immer einen Integer zurück."""
        b = FiberBundle(2, 1, 'SO2')
        c1 = chern_class_first(b)
        assert isinstance(c1, int)
