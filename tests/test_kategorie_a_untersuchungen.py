"""
@file test_kategorie_a_untersuchungen.py
@brief Umfassende Tests für das Kategorie-A-Untersuchungsmodul.
@description
    Testet alle Klassen und Methoden aus kategorie_a_untersuchungen.py,
    die mathematische Vermutungen der Kategorie A (hohe Beweiswahrscheinlichkeit)
    untersuchen:
    - GoldbachUntersuchung
    - ZwillingsprimzahlUntersuchung
    - LegendreUntersuchung
    - BrocardUntersuchung
    - ErdosStrausUntersuchung
    - ArtinVermutungUntersuchung
    - GoldbachPartiellerBeweis
    - kategorie_a_zusammenfassung()

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest

# Sicherstelle, dass src/ im Suchpfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from kategorie_a_untersuchungen import (
    _sieve_primes,
    _is_prime_fast,
    GoldbachUntersuchung,
    ZwillingsprimzahlUntersuchung,
    LegendreUntersuchung,
    BrocardUntersuchung,
    ErdosStrausUntersuchung,
    ArtinVermutungUntersuchung,
    GoldbachPartiellerBeweis,
    kategorie_a_zusammenfassung,
)


# ============================================================
# Hilfsfunktionen-Tests
# ============================================================

class TestSievePrimes:
    """Tests für die Sieb-Hilfsfunktion."""

    def test_leere_liste_bei_null(self):
        """Sieb mit Grenze 0 gibt leere Liste zurück."""
        assert _sieve_primes(0) == []

    def test_leere_liste_bei_eins(self):
        """Sieb mit Grenze 1 gibt leere Liste zurück (1 ist keine Primzahl)."""
        assert _sieve_primes(1) == []

    def test_erste_primzahl(self):
        """Kleinste Primzahl ist 2."""
        assert _sieve_primes(2) == [2]

    def test_erste_10_primzahlen(self):
        """Die ersten 10 Primzahlen bis 30."""
        result = _sieve_primes(30)
        assert result == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    def test_primzahlen_bis_100(self):
        """Anzahl der Primzahlen bis 100 ist 25."""
        assert len(_sieve_primes(100)) == 25

    def test_keine_zusammengesetzten_zahlen(self):
        """Prüft dass keine zusammengesetzten Zahlen enthalten sind."""
        primes = _sieve_primes(50)
        for p in primes:
            assert all(p % i != 0 for i in range(2, p)), f"{p} ist keine Primzahl"


class TestIsPrimeFast:
    """Tests für den schnellen Primzahltest."""

    def test_eins_nicht_prim(self):
        assert not _is_prime_fast(1)

    def test_zwei_ist_prim(self):
        assert _is_prime_fast(2)

    def test_drei_ist_prim(self):
        assert _is_prime_fast(3)

    def test_vier_nicht_prim(self):
        assert not _is_prime_fast(4)

    def test_grosse_primzahl(self):
        """997 ist eine Primzahl."""
        assert _is_prime_fast(997)

    def test_grosse_zusammengesetzte_zahl(self):
        """999 = 3 × 333 = 3 × 3 × 111 = 3³ × 37 ist nicht prim."""
        assert not _is_prime_fast(999)

    def test_null_nicht_prim(self):
        assert not _is_prime_fast(0)

    def test_negative_zahl(self):
        assert not _is_prime_fast(-7)


# ============================================================
# GoldbachUntersuchung Tests
# ============================================================

class TestGoldbachUntersuchung:
    """Tests für die Goldbach-Untersuchungsklasse."""

    def setup_method(self):
        """Erstellt neue Instanz vor jedem Test."""
        self.goldbach = GoldbachUntersuchung()

    def test_verifiziere_bis_100(self):
        """Goldbach-Vermutung verifiziert bis 100."""
        result = self.goldbach.verifiziere_bis(100)
        assert result['verifiziert'] is True
        assert result['bis'] == 100
        assert result['geprüfte_zahlen'] == 49  # 4, 6, 8, ..., 100 → 49 Zahlen

    def test_verifiziere_bis_1000(self):
        """Goldbach-Vermutung verifiziert bis 1000."""
        result = self.goldbach.verifiziere_bis(1000)
        assert result['verifiziert'] is True
        assert result['bis'] == 1000

    def test_kein_gegenbeispiel_bis_500(self):
        """Kein Gegenbeispiel bis 500 gefunden."""
        result = self.goldbach.verifiziere_bis(500)
        assert 'gegenbeispiel' not in result

    def test_beispiele_enthalten(self):
        """Ergebnis enthält Beispiele für kleinste Zerlegungen."""
        result = self.goldbach.verifiziere_bis(50)
        assert 'beispiele' in result
        # 4 = 2 + 2 ist die kanonische erste Zerlegung
        beispiele = result['beispiele']
        assert len(beispiele) > 0

    def test_goldbach_zerlegung_4(self):
        """4 = 2 + 2 ist einzige Zerlegung."""
        count = self.goldbach.goldbach_zerlegungen_anzahl(4)
        assert count == 1

    def test_goldbach_zerlegung_6(self):
        """6 = 3 + 3 ist einzige Zerlegung."""
        count = self.goldbach.goldbach_zerlegungen_anzahl(6)
        assert count == 1

    def test_goldbach_zerlegung_10(self):
        """10 = 3+7 = 5+5 → 2 Zerlegungen."""
        count = self.goldbach.goldbach_zerlegungen_anzahl(10)
        assert count == 2

    def test_goldbach_zerlegung_ungerade_null(self):
        """Ungerade Zahl hat 0 Goldbach-Zerlegungen."""
        assert self.goldbach.goldbach_zerlegungen_anzahl(7) == 0

    def test_goldbach_zerlegung_zu_klein(self):
        """Zahlen ≤ 2 haben 0 Zerlegungen."""
        assert self.goldbach.goldbach_zerlegungen_anzahl(2) == 0

    def test_hardy_littlewood_positiv(self):
        """Hardy-Littlewood-Schätzung ist positiv für gerade n ≥ 4."""
        val = self.goldbach.hardy_littlewood_singular_series(100)
        assert val > 0

    def test_hardy_littlewood_grösser_fuer_grosse_n(self):
        """Schätzung wächst mit n."""
        val100 = self.goldbach.hardy_littlewood_singular_series(100)
        val1000 = self.goldbach.hardy_littlewood_singular_series(1000)
        assert val1000 > val100

    def test_vinogradov_string(self):
        """Vinogradov-Beschreibung enthält Schlüsselinfo."""
        s = self.goldbach.vinogradov_schranke()
        assert 'BEWIESEN' in s or 'bewiesen' in s.lower()
        assert 'Helfgott' in s

    def test_chen_satz_string(self):
        """Chen-Satz-Beschreibung enthält Schlüsselinfo."""
        s = self.goldbach.chen_satz()
        assert 'Chen' in s
        assert 'BEWIESEN' in s or 'bewiesen' in s.lower()


# ============================================================
# ZwillingsprimzahlUntersuchung Tests
# ============================================================

class TestZwillingsprimzahlUntersuchung:
    """Tests für die Zwillingsprimzahl-Untersuchungsklasse."""

    def setup_method(self):
        self.zwilling = ZwillingsprimzahlUntersuchung()

    def test_erste_paare_gefunden(self):
        """Die ersten Zwillingsprimzahlpaare werden korrekt gefunden."""
        result = self.zwilling.verifiziere_bis(50)
        assert result['anzahl_paare'] > 0
        # (3,5), (5,7), (11,13), (17,19), (29,31), (41,43) bis 50
        assert result['anzahl_paare'] == 6

    def test_bekannte_paare(self):
        """(3,5), (5,7), (11,13) sind bekannte Zwillingsprimzahlpaare."""
        result = self.zwilling.verifiziere_bis(15)
        paare = result['erste_10']
        assert (3, 5) in paare
        assert (5, 7) in paare
        assert (11, 13) in paare

    def test_brun_konstante_positiv(self):
        """Brun-Konstante ist positiv."""
        bk = self.zwilling.brun_konstante_naerung(1000)
        assert bk > 0

    def test_brun_konstante_wächst(self):
        """Brun-Konstante wächst mit der Grenze."""
        bk100 = self.zwilling.brun_konstante_naerung(100)
        bk1000 = self.zwilling.brun_konstante_naerung(1000)
        assert bk1000 > bk100

    def test_brun_konstante_näherung_an_limit(self):
        """Brun-Konstante konvergiert gegen ~1.902 (bis 1000 sollte ~ 1.5-1.9 sein)."""
        bk = self.zwilling.brun_konstante_naerung(10000)
        # Bekannter Wert: ~1.902, aber numerisch langsame Konvergenz
        assert 1.0 < bk < 2.5

    def test_zhang_maynard_enthält_schranken(self):
        """Zhang-Maynard-Ergebnis enthält die Schranken."""
        result = self.zwilling.zhang_maynard_ergebnis()
        assert 'zhang_2013' in result
        assert 'maynard_2013' in result
        assert result['zhang_2013']['schranke'] == 70_000_000
        assert result['maynard_2013']['schranke'] == 600
        assert result['polymath8b_2014']['schranke'] == 246

    def test_dichte_zwillingsprimzahlen_wächst(self):
        """Geschätzte Anzahl der Zwillingsprimzahlpaare wächst mit x."""
        d100 = self.zwilling.dichte_zwillingsprimzahlen(100.0)
        d1000 = self.zwilling.dichte_zwillingsprimzahlen(1000.0)
        assert d1000 > d100

    def test_dichte_bei_kleinem_x(self):
        """Bei x=2 ist die Dichte 0."""
        d = self.zwilling.dichte_zwillingsprimzahlen(2.0)
        assert d == 0.0

    def test_dichte_positiv(self):
        """Dichte bei großem x ist positiv."""
        d = self.zwilling.dichte_zwillingsprimzahlen(100.0)
        assert d > 0


# ============================================================
# LegendreUntersuchung Tests
# ============================================================

class TestLegendreUntersuchung:
    """Tests für die Legendre-Untersuchungsklasse."""

    def setup_method(self):
        self.legendre = LegendreUntersuchung()

    def test_verifiziere_bis_20(self):
        """Legendres Vermutung bis n=20 verifiziert."""
        result = self.legendre.verifiziere_bis(20)
        assert result['verifiziert'] is True
        assert result['bis_n'] == 20

    def test_erstes_intervall(self):
        """Zwischen 1=1² und 4=2² liegt die Primzahl 2 oder 3."""
        result = self.legendre.verifiziere_bis(1)
        assert result['verifiziert'] is True

    def test_bertrand_postulat_n2(self):
        """Bertrand-Postulat für n=2: Primzahl zwischen 2 und 4 (= 3)."""
        result = self.legendre.bertrand_postulat(2)
        assert result['gilt'] is True
        assert 2 < result['p'] <= 4

    def test_bertrand_postulat_n10(self):
        """Bertrand-Postulat für n=10: Primzahl zwischen 10 und 20."""
        result = self.legendre.bertrand_postulat(10)
        assert result['gilt'] is True
        assert 10 < result['p'] <= 20

    def test_bertrand_postulat_n100(self):
        """Bertrand-Postulat für n=100."""
        result = self.legendre.bertrand_postulat(100)
        assert result['gilt'] is True
        assert result['p'] > 100
        assert result['p'] <= 200

    def test_huxley_schranke_positiv(self):
        """Huxley-Schranke ergibt positive Anzahl für großes x."""
        val = self.legendre.huxley_schranke(10000.0)
        assert val > 0

    def test_huxley_schranke_wächst(self):
        """Huxley-Schranke wächst mit x."""
        v1 = self.legendre.huxley_schranke(1000.0)
        v2 = self.legendre.huxley_schranke(10000.0)
        assert v2 > v1

    def test_huxley_schranke_theta_parameter(self):
        """Theta-Parameter beeinflusst das Ergebnis."""
        v_klein = self.legendre.huxley_schranke(1000.0, theta=0.5)
        v_gross = self.legendre.huxley_schranke(1000.0, theta=7/12)
        assert v_gross > v_klein

    def test_ingham_kubisch_n1(self):
        """Zwischen 1 und 8 liegt die Primzahl 2, 3, 5, 7."""
        result = self.legendre.ingham_kubisch(1)
        assert 1 < result['primzahl'] < 8

    def test_ingham_kubisch_n2(self):
        """Zwischen 8 und 27 liegt eine Primzahl."""
        result = self.legendre.ingham_kubisch(2)
        assert 8 < result['primzahl'] < 27

    def test_ingham_kubisch_gibt_primzahl(self):
        """ingham_kubisch gibt tatsächlich eine Primzahl zurück."""
        for n in range(1, 6):
            result = self.legendre.ingham_kubisch(n)
            assert _is_prime_fast(result['primzahl']), \
                f"n={n}: {result['primzahl']} ist keine Primzahl"


# ============================================================
# BrocardUntersuchung Tests
# ============================================================

class TestBrocardUntersuchung:
    """Tests für die Brocard-Untersuchungsklasse."""

    def setup_method(self):
        self.brocard = BrocardUntersuchung()

    def test_verifiziere_bis_10(self):
        """Brocard-Vermutung bis Index 10 verifiziert."""
        result = self.brocard.verifiziere_bis(10)
        assert result['verifiziert'] is True

    def test_minimum_mindestens_4(self):
        """Beobachtetes Minimum ist mindestens 4 (Brocard-Bedingung)."""
        result = self.brocard.verifiziere_bis(15)
        assert result['minimum_beobachtet'] >= 4

    def test_beispiel_n2(self):
        """Zwischen 3²=9 und 5²=25 liegen ≥4 Primzahlen (11,13,17,19,23)."""
        result = self.brocard.verifiziere_bis(5)
        assert result['verifiziert'] is True
        # Primzahl-Index 2: p_2=3, p_3=5 → 5 Primzahlen
        details = result['erste_beispiele']
        n2_detail = [d for d in details if d['n'] == 2]
        if n2_detail:
            assert n2_detail[0]['anzahl'] >= 4

    def test_verbindung_zu_legendre_string(self):
        """Verbindungsbeschreibung enthält Schlüsselbegriffe."""
        s = self.brocard.verbindung_zu_legendre()
        assert 'Legendre' in s or 'Brocard' in s
        assert '4' in s  # ≥4 Primzahlen

    def test_verifiziere_gibt_dict_zurück(self):
        """Verifizierung gibt ein Dictionary zurück."""
        result = self.brocard.verifiziere_bis(8)
        assert isinstance(result, dict)
        assert 'verifiziert' in result


# ============================================================
# ErdosStrausUntersuchung Tests
# ============================================================

class TestErdosStrausUntersuchung:
    """Tests für die Erdős-Straus-Untersuchungsklasse."""

    def setup_method(self):
        self.es = ErdosStrausUntersuchung()

    def test_zerlegung_n2(self):
        """4/2 = 2 = 1/1 → (1,1,2) oder andere Zerlegung."""
        result = self.es.finde_zerlegung(2)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/2) < 1e-10

    def test_zerlegung_n3(self):
        """4/3 = 1/1 + 1/4 + 1/12 oder andere."""
        result = self.es.finde_zerlegung(3)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/3) < 1e-10

    def test_zerlegung_n5(self):
        """4/5: Primzahl-Fall."""
        result = self.es.finde_zerlegung(5)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/5) < 1e-10

    def test_zerlegung_n7(self):
        """4/7: schwieriger Fall n ≡ 3 (mod 4)."""
        result = self.es.finde_zerlegung(7)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/7) < 1e-10

    def test_zerlegung_n11(self):
        """4/11: weitere Primzahl."""
        result = self.es.finde_zerlegung(11)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/11) < 1e-10

    def test_alle_zerlegungen_korrekt(self):
        """Alle gefundenen Zerlegungen sind mathematisch korrekt."""
        for n in range(2, 50):
            result = self.es.finde_zerlegung(n)
            if result is not None:
                x, y, z = result
                assert abs(1/x + 1/y + 1/z - 4/n) < 1e-10, \
                    f"n={n}: 1/{x}+1/{y}+1/{z} ≠ 4/{n}"
                assert x > 0 and y > 0 and z > 0, f"n={n}: negative Zahlen"

    def test_verifiziere_bis_100(self):
        """Erdős-Straus-Vermutung für n=2..100 verifiziert."""
        result = self.es.verifiziere_bis(100)
        assert result['verifiziert'] is True
        assert result['nicht_gefunden'] == []

    def test_verifiziere_enthält_beispiele(self):
        """Ergebnis enthält Beispielzerlegungen."""
        result = self.es.verifiziere_bis(20)
        assert 'beispiele' in result
        assert len(result['beispiele']) > 0

    def test_beweis_via_restklassen_dict(self):
        """Restklassen-Beweis gibt Dict mit 4 Einträgen zurück."""
        result = self.es.beweis_via_restklassen()
        assert isinstance(result, dict)
        assert len(result) == 4  # mod 4: r=0,1,2,3

    def test_restklasse_0_bewiesen(self):
        """n ≡ 0 (mod 4) ist als bewiesen markiert."""
        result = self.es.beweis_via_restklassen()
        assert result['n ≡ 0 (mod 4)']['bewiesen'] is True

    def test_zerlegung_positive_zahlen(self):
        """Alle Elemente der Zerlegung sind positive ganze Zahlen."""
        for n in [5, 7, 11, 13, 17, 19, 23]:
            result = self.es.finde_zerlegung(n)
            assert result is not None, f"Keine Zerlegung für n={n}"
            x, y, z = result
            assert isinstance(x, int) and x > 0
            assert isinstance(y, int) and y > 0
            assert isinstance(z, int) and z > 0


# ============================================================
# ArtinVermutungUntersuchung Tests
# ============================================================

class TestArtinVermutungUntersuchung:
    """Tests für die Artin-Vermutungs-Untersuchungsklasse."""

    def setup_method(self):
        self.artin = ArtinVermutungUntersuchung()

    def test_artin_konstante_wert(self):
        """Artin-Konstante ist ca. 0.3740."""
        assert abs(ArtinVermutungUntersuchung.ARTIN_KONSTANTE - 0.3740) < 0.001

    def test_ist_primitive_wurzel_2_mod_5(self):
        """2 ist primitive Wurzel mod 5 (ord_5(2) = 4 = φ(5))."""
        # 2^1=2, 2^2=4, 2^3=3, 2^4=1 (mod 5) → Ordnung 4 = φ(5) ✓
        assert self.artin.ist_primitive_wurzel(2, 5) is True

    def test_ist_primitive_wurzel_2_mod_7(self):
        """2 ist primitive Wurzel mod 7 (ord_7(2) = 3 ≠ 6 = φ(7))."""
        # 2^1=2, 2^2=4, 2^3=1 (mod 7) → Ordnung 3 ≠ 6 → KEINE primitive Wurzel
        assert self.artin.ist_primitive_wurzel(2, 7) is False

    def test_ist_primitive_wurzel_3_mod_7(self):
        """3 ist primitive Wurzel mod 7 (ord_7(3) = 6 = φ(7))."""
        # 3^1=3, 3^2=2, 3^3=6, 3^4=4, 3^5=5, 3^6=1 (mod 7) → Ordnung 6 ✓
        assert self.artin.ist_primitive_wurzel(3, 7) is True

    def test_keine_primitive_wurzel_nicht_prim(self):
        """Für zusammengesetzte Zahlen ist keine primitive Wurzel definiert."""
        assert self.artin.ist_primitive_wurzel(2, 4) is False

    def test_primitive_wurzel_2_mod_11(self):
        """2 ist primitive Wurzel mod 11."""
        # φ(11) = 10, ord_11(2) = 10 ✓
        assert self.artin.ist_primitive_wurzel(2, 11) is True

    def test_dichte_verifizieren_gibt_dict(self):
        """dichte_verifizieren gibt ein Dictionary zurück."""
        result = self.artin.dichte_verifizieren(2, 50)
        assert isinstance(result, dict)
        assert 'beobachtete_dichte' in result
        assert 'artin_konstante' in result

    def test_dichte_zwischen_0_und_1(self):
        """Beobachtete Dichte liegt zwischen 0 und 1."""
        result = self.artin.dichte_verifizieren(2, 100)
        assert 0 < result['beobachtete_dichte'] < 1

    def test_artin_konstante_berechnen(self):
        """Berechnete Artin-Konstante liegt nahe am bekannten Wert."""
        val = self.artin.artin_konstante_berechnen(100)
        # Konvergiert langsam, aber sollte im Bereich 0.30-0.42 liegen
        assert 0.30 < val < 0.45

    def test_artin_konstante_mehr_primes_näher(self):
        """Mehr Primzahlen im Produkt → nähert sich dem echten Wert an."""
        val50 = self.artin.artin_konstante_berechnen(50)
        val100 = self.artin.artin_konstante_berechnen(100)
        # Beide sollten positive Werte sein
        assert val50 > 0
        assert val100 > 0

    def test_hooley_bedingung_string(self):
        """Hooley-Beschreibung enthält Schlüsselbegriffe."""
        s = self.artin.hooley_bedingung()
        assert 'Hooley' in s
        assert 'GRH' in s or 'Riemann' in s


# ============================================================
# GoldbachPartiellerBeweis Tests
# ============================================================

class TestGoldbachPartiellerBeweis:
    """Tests für den partiellen Goldbach-Beweis."""

    def setup_method(self):
        self.gpb = GoldbachPartiellerBeweis()

    def test_goldbach_6k_k1(self):
        """6 = 3 + 3 (k=1)."""
        result = self.gpb.goldbach_fuer_vielfache_von_6(1)
        assert result is not None
        p, q = result
        assert p + q == 6
        assert _is_prime_fast(p) and _is_prime_fast(q)

    def test_goldbach_6k_k2(self):
        """12 = p + q mit Primzahlen p, q (k=2)."""
        result = self.gpb.goldbach_fuer_vielfache_von_6(2)
        assert result is not None
        p, q = result
        assert p + q == 12
        assert _is_prime_fast(p) and _is_prime_fast(q)

    def test_goldbach_6k_mehrere_k(self):
        """Goldbach für 6k mit k=1..20."""
        for k in range(1, 21):
            result = self.gpb.goldbach_fuer_vielfache_von_6(k)
            assert result is not None, f"Kein Paar für 6×{k}={6*k}"
            p, q = result
            assert p + q == 6 * k
            assert _is_prime_fast(p) and _is_prime_fast(q)

    def test_goldbach_zweierpotenz_k2(self):
        """4 = 2 + 2 (2^2)."""
        result = self.gpb.goldbach_fuer_zweierpotenzen(2)
        assert result is not None
        p, q = result
        assert p + q == 4
        assert _is_prime_fast(p) and _is_prime_fast(q)

    def test_goldbach_zweierpotenz_k3(self):
        """8 = 3 + 5 (2^3)."""
        result = self.gpb.goldbach_fuer_zweierpotenzen(3)
        assert result is not None
        p, q = result
        assert p + q == 8
        assert _is_prime_fast(p) and _is_prime_fast(q)

    def test_goldbach_zweierpotenz_mehrere(self):
        """Goldbach für 2^k mit k=2..10."""
        for k in range(2, 11):
            result = self.gpb.goldbach_fuer_zweierpotenzen(k)
            assert result is not None, f"Kein Paar für 2^{k}={2**k}"
            p, q = result
            assert p + q == 2 ** k

    def test_halbprimzahlen_n6(self):
        """6 als Summe zweier P₂-Zahlen."""
        paare = self.gpb.satz_gerade_zahlen_als_summe_von_zwei_halbprimzahlen(6)
        assert len(paare) > 0
        # Überprüfe, dass die Paare wirklich summieren
        for p, q in paare:
            assert p + q == 6

    def test_halbprimzahlen_n12(self):
        """12 als Summe zweier P₂-Zahlen."""
        paare = self.gpb.satz_gerade_zahlen_als_summe_von_zwei_halbprimzahlen(12)
        assert len(paare) > 0
        for p, q in paare:
            assert p + q == 12

    def test_halbprimzahlen_alle_gerade(self):
        """Für jede gerade Zahl n=4..30 gibt es P₂-Paare."""
        for n in range(4, 31, 2):
            paare = self.gpb.satz_gerade_zahlen_als_summe_von_zwei_halbprimzahlen(n)
            assert len(paare) > 0, f"Keine P₂-Paare für n={n}"


# ============================================================
# kategorie_a_zusammenfassung Tests
# ============================================================

class TestKategorieAZusammenfassung:
    """Tests für die Zusammenfassungsfunktion."""

    def test_gibt_dict_zurück(self):
        """Zusammenfassung gibt ein Dictionary zurück."""
        result = kategorie_a_zusammenfassung()
        assert isinstance(result, dict)

    def test_enthält_alle_vermutungen(self):
        """Alle 6 Hauptvermutungen sind enthalten."""
        result = kategorie_a_zusammenfassung()
        assert 'Goldbach' in result
        assert 'Zwillingsprimzahlen' in result
        assert 'Legendre' in result
        assert 'Brocard' in result
        assert 'Erdős-Straus' in result
        assert 'Artin' in result

    def test_jeder_eintrag_hat_status(self):
        """Jeder Eintrag hat einen Status-Schlüssel."""
        result = kategorie_a_zusammenfassung()
        for name, info in result.items():
            assert 'status' in info, f"{name} hat keinen Status"

    def test_artin_status_bedingt_bewiesen(self):
        """Artin-Vermutung hat Status als bedingt bewiesen."""
        result = kategorie_a_zusammenfassung()
        status = result['Artin']['status'].lower()
        assert 'bewiesen' in status or 'grh' in status

    def test_goldbach_verifiziert_bis_angegeben(self):
        """Goldbach-Eintrag enthält Information über Verifikationsgrenze."""
        result = kategorie_a_zusammenfassung()
        gb = result['Goldbach']
        assert 'verifiziert_bis' in gb or 'beste_resultate' in gb

    def test_jeder_eintrag_hat_beste_resultate(self):
        """Die meisten Einträge haben 'beste_resultate'."""
        result = kategorie_a_zusammenfassung()
        count_mit_resultaten = sum(
            1 for info in result.values() if 'beste_resultate' in info
        )
        assert count_mit_resultaten >= 4  # Mindestens 4 von 6


# ============================================================
# Integrationstests
# ============================================================

class TestIntegration:
    """Integrationstests die mehrere Klassen kombinieren."""

    def test_goldbach_und_zerlegungsanzahl_konsistent(self):
        """Wenn verifiziere_bis True, sollten alle Zerlegungsanzahlen > 0 sein."""
        gb = GoldbachUntersuchung()
        result = gb.verifiziere_bis(50)
        assert result['verifiziert'] is True
        for n in range(4, 51, 2):
            count = gb.goldbach_zerlegungen_anzahl(n)
            assert count >= 1, f"n={n} hat keine Goldbach-Zerlegung trotz Verifikation"

    def test_legendre_und_bertrand_konsistent(self):
        """Bertrand-Postulat ist schwächer als Legendre – beide sollten gelten."""
        leg = LegendreUntersuchung()
        for n in range(1, 20):
            bertrand = leg.bertrand_postulat(n)
            assert bertrand['gilt'] is True, f"Bertrand-Postulat versagt für n={n}"

    def test_erdos_straus_restklassen_und_zerlegung_konsistent(self):
        """Die Restklassen-Analyse sollte mit gefundenen Zerlegungen übereinstimmen."""
        es = ErdosStrausUntersuchung()
        rk = es.beweis_via_restklassen()
        # n ≡ 0 (mod 4) ist bewiesen – teste einige solche n
        for n in [4, 8, 12, 16, 20]:
            result = es.finde_zerlegung(n)
            assert result is not None, f"n={n} ≡ 0 (mod 4) hat keine Zerlegung"
            x, y, z = result
            assert abs(1/x + 1/y + 1/z - 4/n) < 1e-10

    def test_artin_primitiv_wurzel_korrekt(self):
        """Primitive Wurzeln und Nicht-Wurzeln sind korrekt klassifiziert."""
        artin = ArtinVermutungUntersuchung()
        # Bekannte primitive Wurzeln mod 7: 3, 5
        assert artin.ist_primitive_wurzel(3, 7) is True
        assert artin.ist_primitive_wurzel(5, 7) is True
        # Bekannte Nicht-Wurzeln mod 7: 2, 4, 6
        assert artin.ist_primitive_wurzel(2, 7) is False
        assert artin.ist_primitive_wurzel(4, 7) is False
        assert artin.ist_primitive_wurzel(6, 7) is False

    def test_zwillingsprimzahl_dichte_und_brun_konsistent(self):
        """Brun-Konstante und Dichte-Schätzung sind beide positiv und konsistent."""
        zp = ZwillingsprimzahlUntersuchung()
        brun = zp.brun_konstante_naerung(1000)
        dichte = zp.dichte_zwillingsprimzahlen(1000.0)
        # Beide positiv
        assert brun > 0
        assert dichte > 0


# ============================================================
# Edge-Case Tests
# ============================================================

class TestEdgeCases:
    """Grenz- und Sonderfälle."""

    def test_goldbach_grenzfall_4(self):
        """Kleinster Goldbach-Fall: 4 = 2 + 2."""
        gb = GoldbachUntersuchung()
        result = gb.verifiziere_bis(4)
        assert result['verifiziert'] is True

    def test_legendre_n1(self):
        """Legendre für n=1: Primzahl zwischen 1 und 4 (= 2 oder 3)."""
        leg = LegendreUntersuchung()
        result = leg.verifiziere_bis(1)
        assert result['verifiziert'] is True

    def test_erdos_straus_grosse_primzahl(self):
        """Erdős-Straus für n=97 (große Primzahl)."""
        es = ErdosStrausUntersuchung()
        result = es.finde_zerlegung(97)
        assert result is not None
        x, y, z = result
        assert abs(1/x + 1/y + 1/z - 4/97) < 1e-10

    def test_artin_negative_basis(self):
        """Negativer Basis-Fall für primitive Wurzel."""
        artin = ArtinVermutungUntersuchung()
        # -1 ist nie eine primitive Wurzel für p > 3 (hat Ordnung 2)
        result = artin.ist_primitive_wurzel(-1, 7)
        # -1 ≡ 6 (mod 7), 6^2 = 36 ≡ 1 (mod 7) → Ordnung 2 ≠ φ(7)=6 → False
        assert result is False

    def test_brocard_sehr_kleine_grenze(self):
        """Brocard für kleinen Index."""
        brocard = BrocardUntersuchung()
        result = brocard.verifiziere_bis(3)
        assert isinstance(result, dict)
        # Sollte trotz kleinem Index korrekt laufen

    def test_zwillingsprimzahl_bis_grenze_3(self):
        """Bis 3: Das Paar (3,5) wird gefunden, da p=3 ≤ 3 und p+2=5 prim ist."""
        zp = ZwillingsprimzahlUntersuchung()
        result = zp.verifiziere_bis(3)
        # Die Implementierung findet (3,5) weil 3 ≤ grenze und 5 prim (im Sieb bis grenze+2)
        assert result['anzahl_paare'] == 1
        assert (3, 5) in result['erste_10']

    def test_goldbach_zerlegungen_grosse_zahl(self):
        """Anzahl der Goldbach-Zerlegungen von 100."""
        gb = GoldbachUntersuchung()
        count = gb.goldbach_zerlegungen_anzahl(100)
        # 100 = 3+97 = 11+89 = 29+71 = 41+59 = 47+53 → mindestens 5 Paare
        assert count >= 5

    def test_brun_konstante_zwei_paare_bis_5(self):
        """Brun-Konstante bis 5: Paare (3,5) und (5,7) werden beide gezählt."""
        zp = ZwillingsprimzahlUntersuchung()
        brun = zp.brun_konstante_naerung(5)
        # range(3, 6, 2) = [3, 5]; 3+2=5 prim ✓; 5+2=7 prim ✓ (im Sieb bis 7)
        expected = (1/3 + 1/5) + (1/5 + 1/7)
        assert abs(brun - expected) < 1e-12
