"""
@file test_erdos_brocard.py
@brief Tests für ErdosStrausExt und BrocardExtension (TDD).
@description
    Umfassende Testsammlung für:
    1. ErdosStrausExt (erdos_straus_ext.py)
       - Lösungsalgorithmen (Unit Fractions)
       - Konstruktive Teilbeweise für Primzahlklassen
       - Vollständige Analyse und Ausnahmensuche
    2. BrocardExtension (brocard_extension.py)
       - Modulare Ausschlussanalyse
       - Numerische Suche (exakte Quadratwurzel-Prüfung)
       - Schranken- und p-adische Analyse

    Mathematische Korrektheit wird durch:
    - Bruchverifikation über Kreuzprodukt (exakt, ohne Gleitkomma)
    - Exakte Quadratwurzelprüfung mit math.isqrt
    - Vergleich gegen bekannte Lösungen

    Abdeckung: Happy Path, Edge Cases, Grenzwerte, mathematische Invarianten.

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import sys
import os

# Suchpfad für src-Module setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest

from erdos_straus_ext import (
    ErdosStrausExt,
    _bruch_pruefen,
    _bruch_pruefen_2terme,
    _gcd,
)
from brocard_extension import (
    BrocardExtension,
    _quadratische_reste,
    _ist_quadrat,
    _fakultaet_mod,
)


# ===========================================================
# Hilfsfunktionen für Tests
# ===========================================================

def _bruch_exakt(n: int, x: int, y: int, z: int) -> bool:
    """
    Unabhängige Verifikation: 4/n = 1/x + 1/y + 1/z
    über Kreuzprodukt (identisch zur Implementierung, aber unabhängig kodiert).
    """
    return 4 * x * y * z == n * (y * z + x * z + x * y)


# ===========================================================
# Tests für Hilfsfunktionen (erdos_straus_ext)
# ===========================================================

class TestHilfsfunktionenErdos:
    """Tests für die Hilfsfunktionen in erdos_straus_ext.py."""

    def test_bruch_pruefen_korrekte_zerlegung(self):
        """4/5 = 1/2 + 1/4 + 1/20 ist eine bekannte korrekte Zerlegung."""
        # 1/2 + 1/4 + 1/20 = 10/20 + 5/20 + 1/20 = 16/20 = 4/5
        assert _bruch_pruefen(5, 2, 4, 20) is True

    def test_bruch_pruefen_falsche_zerlegung(self):
        """1/2 + 1/3 + 1/4 ≠ 4/5."""
        assert _bruch_pruefen(5, 2, 3, 4) is False

    def test_bruch_pruefen_negative_werte(self):
        """Negative Nenner müssen False ergeben."""
        assert _bruch_pruefen(5, -1, 2, 3) is False
        assert _bruch_pruefen(5, 2, -3, 4) is False
        assert _bruch_pruefen(5, 2, 3, -4) is False

    def test_bruch_pruefen_null_nenner(self):
        """Nullnenner müssen False ergeben."""
        assert _bruch_pruefen(5, 0, 2, 3) is False

    def test_bruch_pruefen_2terme(self):
        """4/7 = 1/2 + 1/14: 2-Term-Zerlegung testen."""
        # 1/2 + 1/14 = 7/14 + 1/14 = 8/14 = 4/7 ✓
        assert _bruch_pruefen_2terme(7, 2, 14) is True

    def test_bruch_pruefen_2terme_falsch(self):
        """Falsche 2-Term-Zerlegung wird korrekt abgelehnt."""
        assert _bruch_pruefen_2terme(7, 2, 10) is False

    def test_gcd_grundfall(self):
        """ggT(12, 8) = 4."""
        assert _gcd(12, 8) == 4

    def test_gcd_teilerfremd(self):
        """ggT(7, 13) = 1 für teilerfremde Zahlen."""
        assert _gcd(7, 13) == 1

    def test_gcd_gleiche_zahlen(self):
        """ggT(n, n) = n."""
        assert _gcd(15, 15) == 15

    def test_gcd_mit_null(self):
        """ggT(0, n) = n (Python-Standard)."""
        assert _gcd(0, 7) == 7


# ===========================================================
# Tests für ErdosStrausExt.loese_unit_fraction
# ===========================================================

class TestLoesungUnitFraction:
    """Tests für die Zerlegung 4/n = 1/x + 1/y + 1/z."""

    def setup_method(self):
        """Initialisierung vor jedem Test."""
        self.es = ErdosStrausExt()

    def test_n_gleich_2(self):
        """4/2 = 2 = 1/1 + 1/1... Spezialfall: einfache Zerlegung."""
        loes = self.es.loese_unit_fraction(2)
        assert loes is not None
        assert _bruch_exakt(2, *loes), f"4/2 ≠ 1/{loes[0]}+1/{loes[1]}+1/{loes[2]}"

    def test_n_gleich_3(self):
        """4/3: bekannte Zerlegung existiert."""
        loes = self.es.loese_unit_fraction(3)
        assert loes is not None
        assert _bruch_exakt(3, *loes)

    def test_n_gleich_4(self):
        """4/4 = 1: Zerlegung existiert."""
        loes = self.es.loese_unit_fraction(4)
        assert loes is not None
        assert _bruch_exakt(4, *loes)

    def test_n_gleich_5(self):
        """4/5: bekannte Zerlegung 1/2 + 1/4 + 1/20."""
        loes = self.es.loese_unit_fraction(5)
        assert loes is not None
        assert _bruch_exakt(5, *loes)

    def test_n_gleich_7(self):
        """4/7: p ≡ 3 (mod 4), konstruktive Lösung via Formel."""
        loes = self.es.loese_unit_fraction(7)
        assert loes is not None
        assert _bruch_exakt(7, *loes)

    def test_n_gleich_11(self):
        """4/11: p ≡ 3 (mod 4)."""
        loes = self.es.loese_unit_fraction(11)
        assert loes is not None
        assert _bruch_exakt(11, *loes)

    def test_n_gleich_13(self):
        """4/13: p ≡ 1 (mod 4)."""
        loes = self.es.loese_unit_fraction(13)
        assert loes is not None
        assert _bruch_exakt(13, *loes)

    def test_mod4_klassen_alle_abgedeckt(self):
        """Alle 4 Restklassen mod 4 werden korrekt gelöst."""
        testfaelle = {
            0: [4, 8, 12, 20],   # n ≡ 0 mod 4
            1: [5, 9, 13, 17],   # n ≡ 1 mod 4
            2: [6, 10, 14, 18],  # n ≡ 2 mod 4
            3: [7, 11, 15, 19],  # n ≡ 3 mod 4
        }
        for r, ns in testfaelle.items():
            for n in ns:
                assert n % 4 == r
                loes = self.es.loese_unit_fraction(n)
                assert loes is not None, f"Keine Lösung für n={n} (n≡{r} mod 4)"
                assert _bruch_exakt(n, *loes), f"Falsche Lösung für n={n}: {loes}"

    def test_alle_n_bis_100(self):
        """Für alle n von 2 bis 100 muss eine Lösung existieren (Verifikation)."""
        for n in range(2, 101):
            loes = self.es.loese_unit_fraction(n)
            assert loes is not None, f"Keine Lösung gefunden für n={n}"
            assert _bruch_exakt(n, *loes), f"Falsche Lösung für n={n}: {loes}"

    def test_positive_nenner(self):
        """Alle Nenner der Lösung müssen positiv sein."""
        for n in range(2, 51):
            loes = self.es.loese_unit_fraction(n)
            if loes:
                x, y, z = loes
                assert x > 0 and y > 0 and z > 0, (
                    f"Nicht-positive Nenner für n={n}: {loes}"
                )

    def test_ungueltige_eingabe(self):
        """n < 2 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.es.loese_unit_fraction(1)
        with pytest.raises(ValueError):
            self.es.loese_unit_fraction(0)

    def test_primzahlen_bis_50(self):
        """Für alle Primzahlen bis 50 muss eine Lösung existieren."""
        import sympy
        primes = list(sympy.primerange(2, 51))
        for p in primes:
            loes = self.es.loese_unit_fraction(p)
            assert loes is not None, f"Keine Lösung für Primzahl p={p}"
            assert _bruch_exakt(p, *loes), f"Falsche Lösung für p={p}: {loes}"


# ===========================================================
# Tests für konstruktive Teilbeweise (mod 4)
# ===========================================================

class TestBeweisKlasseMod4_3:
    """Tests für den konstruktiven Beweis p ≡ 3 (mod 4)."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_p_gleich_3(self):
        """p=3: Grundfall p≡3(mod 4)."""
        res = self.es.beweis_klasse_mod4_3(3)
        assert res['2term_verifiziert'] is True
        assert res['3term_verifiziert'] is True

    def test_p_gleich_7(self):
        """p=7: k=1, x=2, y=14. 4/7 = 1/2 + 1/14."""
        res = self.es.beweis_klasse_mod4_3(7)
        assert res['k'] == 1
        assert res['x'] == 2
        assert res['y_2term'] == 14
        assert res['2term_verifiziert'] is True
        assert res['3term_verifiziert'] is True

    def test_p_gleich_11(self):
        """p=11: k=2, x=3, y=33."""
        res = self.es.beweis_klasse_mod4_3(11)
        assert res['k'] == 2
        assert res['x'] == 3
        assert res['y_2term'] == 33
        assert res['2term_verifiziert'] is True

    def test_p_gleich_19(self):
        """p=19: k=4, x=5, y=95."""
        res = self.es.beweis_klasse_mod4_3(19)
        assert res['k'] == 4
        assert res['x'] == 5
        assert res['y_2term'] == 95
        assert res['2term_verifiziert'] is True

    def test_alle_primzahlen_mod4_3(self):
        """Alle Primzahlen ≡ 3 (mod 4) bis 200: 2-Term-Beweis gilt."""
        import sympy
        for p in sympy.primerange(3, 201):
            p = int(p)
            if p % 4 != 3:
                continue
            res = self.es.beweis_klasse_mod4_3(p)
            assert res['2term_verifiziert'] is True, (
                f"2-Term-Beweis fehlgeschlagen für p={p}"
            )
            assert res['3term_verifiziert'] is True, (
                f"3-Term-Erweiterung fehlgeschlagen für p={p}"
            )

    def test_nicht_anwendbar_fuer_mod4_1(self):
        """Für p ≡ 1 (mod 4) ist die Methode nicht direkt anwendbar."""
        res = self.es.beweis_klasse_mod4_3(5)  # 5 ≡ 1 mod 4
        assert res['status'] == 'NICHT-ANWENDBAR (p ≡ 1 mod 4)'

    def test_formelkorrektheit_manuell(self):
        """Manueller Beweis: x=(p+1)/4 für p≡3(mod 4), Bruchverifikation."""
        p = 23  # 23 ≡ 3 (mod 4), k=5, x=6
        k = (p - 3) // 4
        x = k + 1
        y = x * p
        # 2-Term: 4/p = 1/x + 1/(x*p)?
        # 1/x + 1/(x*p) = (p+1)/(x*p) = (p+1)/((k+1)*p) = 4(k+1)/((k+1)*(4k+3)) = 4/p ✓
        assert _bruch_pruefen_2terme(p, x, y) is True


class TestBeweisKlasseMod4_1:
    """Tests für den konstruktiven Beweis p ≡ 1 (mod 4)."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_p_gleich_5(self):
        """p=5: k=1, x=2."""
        res = self.es.beweis_klasse_mod4_1(5)
        assert res['k'] == 1
        assert res['x'] == 2
        assert res['loesung'] is not None
        if res['loesung']:
            assert _bruch_exakt(5, *res['loesung'])

    def test_p_gleich_13(self):
        """p=13: k=3, x=4."""
        res = self.es.beweis_klasse_mod4_1(13)
        assert res['k'] == 3
        assert res['x'] == 4
        assert res['loesung'] is not None
        if res['loesung']:
            assert _bruch_exakt(13, *res['loesung'])

    def test_p_gleich_17(self):
        """p=17: k=4, x=5."""
        res = self.es.beweis_klasse_mod4_1(17)
        assert res['k'] == 4
        assert res['x'] == 5
        assert res['loesung'] is not None

    def test_alle_primzahlen_mod4_1(self):
        """Alle Primzahlen ≡ 1 (mod 4) bis 100 werden gelöst."""
        import sympy
        for p in sympy.primerange(5, 101):
            p = int(p)
            if p % 4 != 1:
                continue
            res = self.es.beweis_klasse_mod4_1(p)
            assert res['loesung'] is not None, f"Keine Lösung für p={p}"
            assert res['loesung_verifiziert'] is True, (
                f"Falsche Lösung für p={p}: {res['loesung']}"
            )


# ===========================================================
# Tests für beweis_klasse_mod3 und beweis_klasse_mod_12
# ===========================================================

class TestBeweisKlasseMod3:
    """Tests für die mod-3-Klassen-Analyse."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_p_mod3_1(self):
        """p=7: 7 ≡ 1 (mod 3)."""
        res = self.es.beweis_klasse_mod3(7, 1)
        assert res['p_mod_3'] == 1
        assert res['klasse_passend'] is True
        assert res['loesung'] is not None

    def test_p_mod3_2(self):
        """p=5: 5 ≡ 2 (mod 3)."""
        res = self.es.beweis_klasse_mod3(5, 2)
        assert res['p_mod_3'] == 2
        assert res['klasse_passend'] is True

    def test_p_mod3_0(self):
        """p=3: 3 ≡ 0 (mod 3)."""
        res = self.es.beweis_klasse_mod3(3, 0)
        assert res['p_mod_3'] == 0
        assert res['klasse_passend'] is True

    def test_klasse_nicht_passend(self):
        """Wenn p nicht zur angefragten Klasse gehört, wird das korrekt angezeigt."""
        res = self.es.beweis_klasse_mod3(7, 2)  # 7 ≡ 1, nicht 2
        assert res['klasse_passend'] is False

    def test_ungueltige_klasse(self):
        """r ∉ {0,1,2} löst Fehlermeldung aus."""
        res = self.es.beweis_klasse_mod3(7, 5)
        assert 'fehler' in res


class TestBeweisKlasseMod12:
    """Tests für die mod-12-Klassen-Analyse."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_p_gleich_7_mod12_7(self):
        """p=7: 7 ≡ 7 (mod 12), also ≡ 3 (mod 4) und ≡ 1 (mod 3)."""
        res = self.es.beweis_klasse_mod_12(7, 7)
        assert res['p_mod_12'] == 7
        assert res['p_mod_4'] == 3
        assert res['p_mod_3'] == 1
        assert res['klasse_passend'] is True

    def test_p_gleich_11_mod12_11(self):
        """p=11: 11 ≡ 11 (mod 12), also ≡ 3 (mod 4) und ≡ 2 (mod 3)."""
        res = self.es.beweis_klasse_mod_12(11, 11)
        assert res['p_mod_12'] == 11
        assert res['p_mod_4'] == 3
        assert res['p_mod_3'] == 2

    def test_loesung_verifiziert(self):
        """Für alle getesteten Primzahlen ist die Lösung korrekt."""
        import sympy
        for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
            r = p % 12
            res = self.es.beweis_klasse_mod_12(p, r)
            if res['loesung']:
                assert res['loesung_verifiziert'] is True, (
                    f"Falsche Lösung für p={p}: {res['loesung']}"
                )


# ===========================================================
# Tests für vollstaendige_analyse und suche_ausnahmen
# ===========================================================

class TestVollstaendigeAnalyse:
    """Tests für die vollständige Analyse und Ausnahmensuche."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_vollstaendig_bis_50(self):
        """Für n von 2 bis 50 müssen alle Fälle gelöst sein."""
        res = self.es.vollstaendige_analyse(50)
        assert res['vollstaendig'] is True
        assert res['nicht_gelöst'] == []
        assert res['gelöst_total'] == 49  # n=2 bis 50

    def test_ausnahmen_bis_100_leer(self):
        """Ausnahmensuche bis 100: keine Ausnahmen erwartet."""
        ausnahmen = self.es.suche_ausnahmen(100)
        assert ausnahmen == [], f"Unerwartete Ausnahmen: {ausnahmen}"

    def test_beispiele_korrekt(self):
        """Die ersten 20 Beispiele in vollstaendige_analyse sind korrekt."""
        res = self.es.vollstaendige_analyse(30)
        for bsp in res['beispiele']:
            n, loes = bsp['n'], bsp['loesung']
            assert _bruch_exakt(n, *loes), f"Beispiel n={n} falsch: {loes}"
            assert bsp['verifiziert'] is True

    def test_statistik_methoden(self):
        """Methodenstatistik enthält nicht-negative Zahlen."""
        res = self.es.vollstaendige_analyse(40)
        assert res['methode_restklasse'] >= 0
        assert res['methode_teilersuche'] >= 0
        assert res['methode_bruteforce'] >= 0
        total = (res['methode_restklasse']
                 + res['methode_teilersuche']
                 + res['methode_bruteforce'])
        assert total == res['gelöst_total']


# ===========================================================
# Tests für parametrische_loesungen
# ===========================================================

class TestParametrischeLoesungen:
    """Tests für die Tabellarische Lösungsübersicht."""

    def setup_method(self):
        self.es = ErdosStrausExt()

    def test_anzahl_eintraege(self):
        """Es müssen genau 4 Einträge (für mod 0, 1, 2, 3) vorhanden sein."""
        formeln = self.es.parametrische_loesungen()
        assert len(formeln) == 4

    def test_alle_klassen_vorhanden(self):
        """Alle 4 Restklassen mod 4 sind vorhanden."""
        formeln = self.es.parametrische_loesungen()
        klassen = [f['klasse'] for f in formeln]
        assert 'n ≡ 0 (mod 4)' in klassen
        assert 'n ≡ 1 (mod 4)' in klassen
        assert 'n ≡ 2 (mod 4)' in klassen
        assert 'n ≡ 3 (mod 4)' in klassen

    def test_beispiele_korrekt(self):
        """Alle Beispiele in der Formelübersicht sind mathematisch korrekt."""
        formeln = self.es.parametrische_loesungen()
        for f in formeln:
            bsp = f['beispiel']
            n, x, y, z = bsp['n'], bsp['x'], bsp['y'], bsp['z']
            assert _bruch_exakt(n, x, y, z), (
                f"Formelbeispiel falsch für Klasse '{f['klasse']}': n={n}, {(x,y,z)}"
            )

    def test_status_felder(self):
        """Jeder Eintrag hat ein Status-Feld."""
        formeln = self.es.parametrische_loesungen()
        for f in formeln:
            assert 'status' in f
            assert len(f['status']) > 0


# ===========================================================
# Tests für Hilfsfunktionen (brocard_extension)
# ===========================================================

class TestHilfsfunktionenBrocard:
    """Tests für die Hilfsfunktionen in brocard_extension.py."""

    def test_quadratische_reste_mod_5(self):
        """QR(5) = {0, 1, 4}."""
        qr = _quadratische_reste(5)
        assert qr == {0, 1, 4}

    def test_quadratische_reste_mod_8(self):
        """QR(8) = {0, 1, 4}."""
        qr = _quadratische_reste(8)
        assert qr == {0, 1, 4}

    def test_quadratische_reste_mod_7(self):
        """QR(7) = {0, 1, 2, 4}."""
        qr = _quadratische_reste(7)
        assert qr == {0, 1, 2, 4}

    def test_ist_quadrat_perfekte_zahlen(self):
        """Perfekte Quadrate werden korrekt erkannt."""
        assert _ist_quadrat(0) is True
        assert _ist_quadrat(1) is True
        assert _ist_quadrat(4) is True
        assert _ist_quadrat(9) is True
        assert _ist_quadrat(16) is True
        assert _ist_quadrat(25) is True
        assert _ist_quadrat(100) is True

    def test_ist_quadrat_keine_quadrate(self):
        """Nicht-Quadratzahlen werden korrekt abgelehnt."""
        assert _ist_quadrat(2) is False
        assert _ist_quadrat(3) is False
        assert _ist_quadrat(5) is False
        assert _ist_quadrat(7) is False
        assert _ist_quadrat(10) is False
        assert _ist_quadrat(99) is False

    def test_ist_quadrat_negative_zahl(self):
        """Negative Zahlen sind keine Quadrate."""
        assert _ist_quadrat(-1) is False
        assert _ist_quadrat(-4) is False

    def test_fakultaet_mod_einfach(self):
        """5! = 120, 120 mod 7 = 1."""
        assert _fakultaet_mod(5, 7) == 120 % 7

    def test_fakultaet_mod_grosse_n(self):
        """Für n ≥ p (Primzahl): n! ≡ 0 (mod p)."""
        assert _fakultaet_mod(7, 5) == 0   # 7! enthält Faktor 5
        assert _fakultaet_mod(10, 7) == 0  # 10! enthält Faktor 7

    def test_fakultaet_mod_1(self):
        """Alles mod 1 ist 0."""
        assert _fakultaet_mod(5, 1) == 0
        assert _fakultaet_mod(100, 1) == 0


# ===========================================================
# Tests für BrocardExtension.numerische_suche
# ===========================================================

class TestNumerischeSuche:
    """Tests für die numerische Verifikation der Brocard-Gleichung."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_bekannte_loesungen_gefunden(self):
        """Alle drei bekannten Lösungen (4,5,7) werden gefunden."""
        res = self.ext.numerische_suche(10)
        gefundene_n = {s['n'] for s in res['loesungen']}
        assert 4 in gefundene_n, "Lösung n=4 nicht gefunden"
        assert 5 in gefundene_n, "Lösung n=5 nicht gefunden"
        assert 7 in gefundene_n, "Lösung n=7 nicht gefunden"

    def test_korrekte_m_werte(self):
        """Die m-Werte der gefundenen Lösungen stimmen."""
        res = self.ext.numerische_suche(10)
        loes_map = {s['n']: s['m'] for s in res['loesungen']}
        assert loes_map.get(4) == 5, f"Erwartet m=5 für n=4, erhalten {loes_map.get(4)}"
        assert loes_map.get(5) == 11, f"Erwartet m=11 für n=5, erhalten {loes_map.get(5)}"
        assert loes_map.get(7) == 71, f"Erwartet m=71 für n=7, erhalten {loes_map.get(7)}"

    def test_verifikation_n_fak_plus1(self):
        """Für jede gefundene Lösung gilt n!+1 = m² exakt."""
        res = self.ext.numerische_suche(10)
        for s in res['loesungen']:
            n, m = s['n'], s['m']
            fak = math.factorial(n)
            assert fak + 1 == m * m, f"n!+1 ≠ m² für n={n}, m={m}"

    def test_keine_unerwarteten_loesungen_bis_100(self):
        """Zwischen n=8 und n=100 darf keine neue Lösung auftauchen."""
        res = self.ext.numerische_suche(100)
        unerwartete = [
            s for s in res['loesungen']
            if s['n'] not in BrocardExtension.BEKANNTE_LOESUNGEN
        ]
        assert unerwartete == [], f"Unerwartete Lösung: {unerwartete}"

    def test_vermutungsstatus_bestaetigt(self):
        """Vermutungsstatus ist 'BESTÄTIGT' für n bis 100."""
        res = self.ext.numerische_suche(100)
        assert 'BESTÄTIGT' in res['vermutung_status']

    def test_n_gleich_4_5_7_verifikation_direkt(self):
        """Direkte Verifikation der drei bekannten Lösungen."""
        assert math.factorial(4) + 1 == 5 * 5
        assert math.factorial(5) + 1 == 11 * 11
        assert math.factorial(7) + 1 == 71 * 71

    def test_laufzeit_messung_vorhanden(self):
        """Laufzeit wird korrekt gemessen und zurückgegeben."""
        res = self.ext.numerische_suche(50)
        assert 'laufzeit_sekunden' in res
        assert res['laufzeit_sekunden'] >= 0


# ===========================================================
# Tests für BrocardExtension.modular_ausschluss
# ===========================================================

class TestModularAusschluss:
    """Tests für die modulare Ausschlussanalyse."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_bekannte_loesungen_nicht_ausgeschlossen(self):
        """Bekannte Lösungen n=4,5,7 dürfen NICHT ausgeschlossen werden."""
        moduli = [5, 7, 11, 13, 8, 9, 25]
        for n in [4, 5, 7]:
            res = self.ext.modular_ausschluss(n, moduli)
            assert res['ausgeschlossen'] is False, (
                f"Bekannte Lösung n={n} fälschlich ausgeschlossen durch "
                f"m={res['ausschluss_modulus']}"
            )

    def test_ergebnis_struktur(self):
        """Das Ergebnis hat alle erforderlichen Felder."""
        res = self.ext.modular_ausschluss(10, [5, 7, 11])
        assert 'n' in res
        assert 'moduli' in res
        assert 'ausgeschlossen' in res
        assert 5 in res['moduli']
        assert 7 in res['moduli']
        assert 11 in res['moduli']

    def test_ausschluss_korrektheit(self):
        """Wenn ausgeschlossen=True, muss rest ∉ QR(modulus) sein."""
        res = self.ext.modular_ausschluss(10, list(range(2, 30)))
        if res['ausgeschlossen']:
            m = res['ausschluss_modulus']
            r = res['ausschluss_rest']
            qr = _quadratische_reste(m)
            assert r not in qr, (
                f"Widerspruch: {r} sollte kein Quadrat mod {m} sein, aber {r} ∈ QR"
            )

    def test_leere_moduli_liste(self):
        """Leere Moduli-Liste: kein Ausschluss möglich."""
        res = self.ext.modular_ausschluss(10, [])
        assert res['ausgeschlossen'] is False

    def test_modulus_1_ignoriert(self):
        """Modulus 1 wird ignoriert (nicht sinnvoll)."""
        res = self.ext.modular_ausschluss(10, [1, 5])
        assert 1 not in res['moduli']


# ===========================================================
# Tests für BrocardExtension.finde_ausschlusskandidaten
# ===========================================================

class TestFindeAusschlusskandidaten:
    """Tests für die Suche nach Ausschluss-Moduli."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_bekannte_loesungen_nicht_ausgeschlossen(self):
        """Bekannte Lösungen n=4,5,7 dürfen durch keinen Modulus ausgeschlossen werden."""
        for n in [4, 5, 7]:
            res = self.ext.finde_ausschlusskandidaten(n)
            assert res['ausgeschlossen'] is False, (
                f"Bekannte Lösung n={n} fälschlich ausgeschlossen"
            )

    def test_ergebnis_struktur(self):
        """Ergebnis hat korrekte Struktur."""
        res = self.ext.finde_ausschlusskandidaten(10)
        assert 'n' in res
        assert 'ausschluss_moduli' in res
        assert 'kein_ausschluss' in res
        assert 'ausgeschlossen' in res

    def test_konsistenz_ausgeschlossen_flag(self):
        """ausgeschlossen=True iff ausschluss_moduli nicht leer."""
        for n in range(8, 25):
            res = self.ext.finde_ausschlusskandidaten(n)
            assert res['ausgeschlossen'] == (len(res['ausschluss_moduli']) > 0)


# ===========================================================
# Tests für BrocardExtension.schranken_analyse
# ===========================================================

class TestSchrankenAnalyse:
    """Tests für die Schrankenanalyse."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_bekannte_loesungen_in_beispielen(self):
        """n=4, 5, 7 müssen als Lösungen markiert sein."""
        res = self.ext.schranken_analyse()
        loesungen_n = {b['n'] for b in res['beispiele'] if b['ist_loesung']}
        assert 4 in loesungen_n
        assert 5 in loesungen_n
        assert 7 in loesungen_n

    def test_m_werte_korrekt(self):
        """m-Werte für n=4,5,7 stimmen mit bekannten Lösungen überein."""
        res = self.ext.schranken_analyse()
        for b in res['beispiele']:
            n = b['n']
            if n in BrocardExtension.BEKANNTE_LOESUNGEN:
                m_erwartet = BrocardExtension.BEKANNTE_LOESUNGEN[n]
                assert b['m_exakte_wurzel'] == m_erwartet, (
                    f"Falsches m für n={n}: {b['m_exakte_wurzel']} ≠ {m_erwartet}"
                )

    def test_n_fak_stellen_wachsend(self):
        """n! hat für größere n mehr Stellen (monoton)."""
        res = self.ext.schranken_analyse()
        stellen = [b['n_fak_stellen'] for b in res['beispiele']]
        for i in range(1, len(stellen)):
            assert stellen[i] >= stellen[i - 1], (
                f"n! Stellen nicht monoton an Position {i}"
            )

    def test_formeln_vorhanden(self):
        """Alle Formelfelder sind vorhanden."""
        res = self.ext.schranken_analyse()
        assert 'formel_untere_schranke' in res
        assert 'formel_obere_schranke' in res
        assert 'wachstumsrate' in res


# ===========================================================
# Tests für BrocardExtension.p_adische_analyse
# ===========================================================

class TestPAdischeAnalyse:
    """Tests für die p-adische Analyse."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_bekannte_loesungen_korrekt_markiert(self):
        """In allen Beispielen sind n=4,5,7 korrekt als Lösungen markiert."""
        for p_test in [2, 3, 5, 7]:
            res = self.ext.p_adische_analyse(p_test)
            for e in res['beispiele']:
                n = e['n']
                loesung_erwartet = n in BrocardExtension.BEKANNTE_LOESUNGEN
                assert e['ist_loesung'] == loesung_erwartet, (
                    f"p={p_test}, n={n}: ist_loesung={e['ist_loesung']} "
                    f"aber erwartet={loesung_erwartet}"
                )

    def test_keine_bekannte_loesung_ausgeschlossen(self):
        """Bekannte Lösungen dürfen durch p-adische Bewertung nicht ausgeschlossen werden."""
        for p_test in [2, 3, 5]:
            res = self.ext.p_adische_analyse(p_test)
            for e in res['beispiele']:
                if e['ist_loesung']:
                    assert e['ausschluss'] is False, (
                        f"Bekannte Lösung n={e['n']} durch p={p_test} "
                        f"fälschlich ausgeschlossen (v_p={e['v_p(n!+1)']})"
                    )

    def test_vp_gerade_konsistenz(self):
        """v_p_gerade ist korrekt aus v_p(n!+1) abgeleitet."""
        res = self.ext.p_adische_analyse(3)
        for e in res['beispiele']:
            assert e['v_p_gerade'] == (e['v_p(n!+1)'] % 2 == 0)

    def test_ungueltige_eingabe(self):
        """Nicht-Primzahl erzeugt Fehlermeldung."""
        res = self.ext.p_adische_analyse(4)
        assert 'fehler' in res

    def test_ergebnis_struktur(self):
        """Ergebnis hat alle Pflichtfelder."""
        res = self.ext.p_adische_analyse(5)
        assert 'p' in res
        assert 'beispiele' in res
        assert 'legendre_formel' in res
        assert 'fazit' in res

    def test_legendre_formel_korrekt(self):
        """Legendres Formel: v_5(10!) = ⌊10/5⌋ + ⌊10/25⌋ = 2."""
        res = self.ext.p_adische_analyse(5)
        e10 = next((e for e in res['beispiele'] if e['n'] == 10), None)
        assert e10 is not None
        assert e10['v_p(n!)'] == 2, f"v_5(10!) = {e10['v_p(n!)']} ≠ 2"


# ===========================================================
# Tests für analysiere_restklassen
# ===========================================================

class TestAnalyseRestklassen:
    """Tests für die Restklassenanalyse (welche Moduli sind nützlich)."""

    def setup_method(self):
        self.ext = BrocardExtension()

    def test_grundstruktur(self):
        """Ergebnis hat korrekte Felder."""
        res = self.ext.analysiere_restklassen(15)
        assert 'max_mod' in res
        assert 'moduli_analyse' in res
        assert 'nuetzliche_moduli' in res
        assert res['max_mod'] == 15

    def test_alle_moduli_analysiert(self):
        """Alle Moduli von 2 bis max_mod sind analysiert."""
        res = self.ext.analysiere_restklassen(10)
        for m in range(2, 11):
            assert m in res['moduli_analyse'], f"Modulus {m} fehlt"

    def test_qr_korrekt(self):
        """Quadratische Reste für m=5 sind korrekt."""
        res = self.ext.analysiere_restklassen(5)
        qr = set(res['moduli_analyse'][5]['quadratische_reste'])
        assert qr == {0, 1, 4}


# ===========================================================
# Integrationstests (end-to-end)
# ===========================================================

class TestIntegration:
    """Integrationstests, die beide Module zusammen testen."""

    def test_erdos_straus_und_brute_force_konsistent(self):
        """Restklassenformeln und Brute-Force liefern für gleiche n korrekte Ergebnisse."""
        es = ErdosStrausExt()
        for n in [7, 11, 13, 17, 19, 23]:
            loes_formel = es._strategie_restklasse(n)
            loes_brute = es._strategie_brute_force(n, limit=200)
            if loes_formel:
                assert _bruch_exakt(n, *loes_formel), f"Formel falsch für n={n}"
            if loes_brute:
                assert _bruch_exakt(n, *loes_brute), f"Brute-Force falsch für n={n}"

    def test_brocard_bekannte_loesungen_vollstaendig(self):
        """Alle bekannten Lösungen werden numerisch verifiziert."""
        ext = BrocardExtension()
        for n, m in BrocardExtension.BEKANNTE_LOESUNGEN.items():
            fak = math.factorial(n)
            assert fak + 1 == m * m, f"Bekannte Lösung ({n},{m}) falsch"

    def test_erdos_straus_primzahlklassen_konsistenz(self):
        """beweis_klasse_mod4_3 und loese_unit_fraction liefern beide gültige Lösungen."""
        es = ErdosStrausExt()
        import sympy
        for p in sympy.primerange(3, 50):
            p = int(p)
            if p % 4 != 3:
                continue
            res_beweis = es.beweis_klasse_mod4_3(p)
            loes_allg = es.loese_unit_fraction(p)

            # Beweis liefert korrekte Lösung
            assert res_beweis['3term_verifiziert'] is True

            # Allgemeine Methode liefert ebenfalls korrekte Lösung
            assert loes_allg is not None
            assert _bruch_exakt(p, *loes_allg)
