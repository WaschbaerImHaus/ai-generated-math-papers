"""
@file test_giuga_4prim.py
@brief pytest-Tests für giuga_4prim.py — Giuga 4-Prim-Beweisversuch.
@description
    Testet alle Kernfunktionen des 4-Prim-Giuga-Moduls:

    1. Erkennung bekannter Giuga-Zahlen (30, 858, 1722, 66198) als Giuga-Zahlen
       (nur schwache Bedingung), aber NICHT als Giuga-Pseudoprimes.
    2. Schrankenberechnungen für die Fälle p=2, p=3 und p≥5.
    3. Numerische Suche bis 10^6 ohne Pseudoprime-Fund.
    4. Korrekte Untere-Schranken-Analyse für p=2.
    5. Einzelne Giuga-Bedingungsprüfungen (schwach/stark, quadratfrei, prim).

    Edge-Cases:
    - n < 4 (zu klein)
    - Primzahlen (kein Pseudoprime möglich)
    - Nicht-quadratfreie Zahlen (p² | n → kein Giuga-Kandidat)
    - Korrekte Abgrenzung Giuga-Zahl vs. Giuga-Pseudoprime

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
"""

import pytest
import sympy
from sympy import factorint, isprime

# Modul importieren — liegt unter src/
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from giuga_4prim import (
    berechne_giuga_bedingung,
    numerische_suche_4prim,
    schranke_fuer_4prim_fall_p_gleich_2,
    schranke_fuer_4prim_fall_p_gleich_3,
    schranke_fuer_4prim_allgemeiner_fall,
    Giuga4PrimBeweis,
)


# =============================================================================
# Tests: berechne_giuga_bedingung
# =============================================================================


class TestBerechnGiugaBedingung:
    """
    @brief Testet die Funktion berechne_giuga_bedingung(n).
    @lastModified 2026-03-12
    """

    # --- Bekannte Giuga-Zahlen (nur schwache Bedingung) ---

    def test_giuga_zahl_30(self):
        """
        @brief 30 = 2·3·5 ist eine Giuga-Zahl (schwach), kein Giuga-Pseudoprime.
        @description
            Schwache Bedingung für 30 = 2·3·5:
            - 2 | (15-1) = 14 ✓
            - 3 | (10-1) = 9 ✓
            - 5 | (6-1) = 5 ✓
            Starke Bedingung:
            - (2-1)=1 | 14 ✓
            - (3-1)=2 | 9 = 9/2 nicht ganzzahlig ✗ → kein Pseudoprime
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(30)
        assert ergebnis['ist_giuga_zahl'] is True, \
            "30 muss als Giuga-Zahl (schwach) erkannt werden"
        assert ergebnis['ist_giuga_pseudoprime'] is False, \
            "30 darf NICHT als Giuga-Pseudoprime erkannt werden"
        assert ergebnis['primfaktoren'] == [2, 3, 5]

    def test_giuga_zahl_858(self):
        """
        @brief 858 = 2·3·11·13 ist eine Giuga-Zahl (schwach), kein Pseudoprime.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(858)
        assert ergebnis['ist_giuga_zahl'] is True, \
            "858 muss als Giuga-Zahl erkannt werden"
        assert ergebnis['ist_giuga_pseudoprime'] is False, \
            "858 darf NICHT als Pseudoprime erkannt werden"
        assert set(ergebnis['primfaktoren']) == {2, 3, 11, 13}

    def test_giuga_zahl_1722(self):
        """
        @brief 1722 = 2·3·7·41 ist eine Giuga-Zahl, kein Pseudoprime.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(1722)
        assert ergebnis['ist_giuga_zahl'] is True, \
            "1722 muss als Giuga-Zahl erkannt werden"
        assert ergebnis['ist_giuga_pseudoprime'] is False, \
            "1722 darf NICHT als Pseudoprime erkannt werden"

    def test_giuga_zahl_66198(self):
        """
        @brief 66198 = 2·3·11·17·59 ist eine Giuga-Zahl, kein Pseudoprime.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(66198)
        assert ergebnis['ist_giuga_zahl'] is True, \
            "66198 muss als Giuga-Zahl erkannt werden"
        assert ergebnis['ist_giuga_pseudoprime'] is False, \
            "66198 darf NICHT als Pseudoprime erkannt werden"

    # --- Primzahlen: kein Pseudoprime ---

    def test_prim_2(self):
        """
        @brief n=2 ist keine Giuga-Zahl (zu klein oder prim, je nach Code-Pfad).
        @description
            n=2 ist prim und < 4. Die Funktion lehnt beide Bedingungen ab.
            Der genaue Grund ('n < 4' oder 'prim') hängt von der Implementierung ab —
            entscheidend ist nur, dass kein Giuga-Ergebnis zurückgegeben wird.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(2)
        assert ergebnis['ist_giuga_zahl'] is False
        assert ergebnis['ist_giuga_pseudoprime'] is False
        # Prüfe, dass ein 'grund' vorhanden und nicht leer ist
        assert 'grund' in ergebnis
        assert len(ergebnis['grund']) > 0

    def test_prim_97(self):
        """@brief Primzahl 97 ist kein Pseudoprime. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(97)
        assert ergebnis['ist_giuga_zahl'] is False
        assert ergebnis['ist_giuga_pseudoprime'] is False

    # --- Zu kleine Zahlen ---

    def test_n_kleiner_4(self):
        """@brief n < 4 wird als nicht relevant zurückgegeben. @lastModified 2026-03-12"""
        for n in [1, 2, 3]:
            ergebnis = berechne_giuga_bedingung(n)
            assert ergebnis['ist_giuga_zahl'] is False
            assert ergebnis['ist_giuga_pseudoprime'] is False

    # --- Nicht-quadratfreie Zahlen ---

    def test_nicht_quadratfrei_4(self):
        """
        @brief 4 = 2² ist nicht quadratfrei → kein Giuga-Kandidat.
        @description
            Satz 1 aus beweisversuche.py: Alle Giuga-Pseudoprimes sind quadratfrei.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(4)
        assert ergebnis['ist_giuga_zahl'] is False
        assert ergebnis['ist_giuga_pseudoprime'] is False

    def test_nicht_quadratfrei_12(self):
        """@brief 12 = 2²·3 ist nicht quadratfrei. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(12)
        assert ergebnis['ist_giuga_zahl'] is False

    def test_nicht_quadratfrei_18(self):
        """@brief 18 = 2·3² ist nicht quadratfrei. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(18)
        assert ergebnis['ist_giuga_zahl'] is False

    # --- Normale zusammengesetzte Zahlen (nicht Giuga) ---

    def test_n_15_kein_giuga(self):
        """
        @brief 15 = 3·5 ist keine Giuga-Zahl.
        @description
            Schwache Bedingung: 3 | (5-1)=4? Nein → keine Giuga-Zahl.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(15)
        assert ergebnis['ist_giuga_zahl'] is False

    def test_n_6_kein_giuga(self):
        """
        @brief 6 = 2·3 ist keine Giuga-Zahl.
        @description
            Schwache Bedingung: 2 | (3-1)=2 ✓, aber 3 | (2-1)=1? Nein.
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(6)
        assert ergebnis['ist_giuga_zahl'] is False

    # --- Rückgabestruktur ---

    def test_rückgabe_enthält_primfaktoren_für_30(self):
        """@brief Rückgabe enthält 'primfaktoren' für zusammengesetztes quadratfreies n. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(30)
        assert 'primfaktoren' in ergebnis
        assert isinstance(ergebnis['primfaktoren'], list)
        assert len(ergebnis['primfaktoren']) > 0

    def test_rückgabe_enthält_details(self):
        """@brief Rückgabe enthält schwach- und stark-Details. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(30)
        assert 'schwach_details' in ergebnis
        assert 'stark_details' in ergebnis
        assert len(ergebnis['schwach_details']) == 3  # 3 Primfaktoren von 30


# =============================================================================
# Tests: numerische_suche_4prim
# =============================================================================


class TestNumerischeSuche4Prim:
    """
    @brief Testet die Funktion numerische_suche_4prim(grenze).
    @lastModified 2026-03-12
    """

    def test_kein_4prim_pseudoprime_bis_1_000_000(self):
        """
        @brief Kein 4-Prim-Giuga-Pseudoprime bis 10^6.
        @description
            Haupttest: Bestätigt numerisch, dass bis zur Schranke von 1.000.000
            kein 4-Prim-Giuga-Pseudoprime existiert.
        @lastModified 2026-03-12
        """
        ergebnis = numerische_suche_4prim(1_000_000)
        assert ergebnis['kein_pseudoprime_gefunden'] is True, \
            f"Unerwarteter Fund: {ergebnis['giuga_pseudoprimes_4prim']}"
        assert ergebnis['giuga_pseudoprimes_4prim'] == []

    def test_suchgrenze_wird_respektiert(self):
        """
        @brief Alle gefundenen n liegen ≤ grenze.
        @lastModified 2026-03-12
        """
        grenze = 100_000
        ergebnis = numerische_suche_4prim(grenze)
        for item in ergebnis['giuga_zahlen_4prim']:
            assert item['n'] <= grenze, f"n={item['n']} überschreitet grenze={grenze}"

    def test_giuga_zahlen_haben_4_faktoren(self):
        """
        @brief Alle gefundenen Giuga-Zahlen haben genau 4 Primfaktoren.
        @lastModified 2026-03-12
        """
        ergebnis = numerische_suche_4prim(100_000)
        for item in ergebnis['giuga_zahlen_4prim']:
            assert len(item['faktoren']) == 4, \
                f"n={item['n']} hat {len(item['faktoren'])} Faktoren, erwartet 4"

    def test_bekannte_giuga_zahlen_858_1722_nicht_in_4prim(self):
        """
        @brief 858 und 1722 sind keine 4-Prim-Zahlen (sie haben 4 Primfaktoren).
        @description
            858 = 2·3·11·13 (4 Faktoren), 1722 = 2·3·7·41 (4 Faktoren).
            Diese sind Giuga-Zahlen (schwach), aber kein Pseudoprime.
            Die Suche sollte sie als Giuga-Zahlen (schwach) aber NICHT als
            Pseudoprimes finden.
        @lastModified 2026-03-12
        """
        ergebnis = numerische_suche_4prim(2000)
        # Prüfe, dass sie in giuga_zahlen sind (schwach), nicht in pseudoprimes
        gefundene_n = [item['n'] for item in ergebnis['giuga_zahlen_4prim']]
        assert 858 in gefundene_n, \
            "858 sollte als 4-Prim-Giuga-Zahl (schwach) gefunden werden"
        assert 1722 in gefundene_n, \
            "1722 sollte als 4-Prim-Giuga-Zahl (schwach) gefunden werden"
        # Beide dürfen NICHT in den Pseudoprimes sein
        pseudoprime_n = [item['n'] for item in ergebnis['giuga_pseudoprimes_4prim']]
        assert 858 not in pseudoprime_n
        assert 1722 not in pseudoprime_n

    def test_rückgabe_struktur(self):
        """@brief Rückgabe-Dictionary hat alle erwarteten Schlüssel. @lastModified 2026-03-12"""
        ergebnis = numerische_suche_4prim(1000)
        assert 'grenze' in ergebnis
        assert 'geprueft' in ergebnis
        assert 'giuga_zahlen_4prim' in ergebnis
        assert 'giuga_pseudoprimes_4prim' in ergebnis
        assert 'kein_pseudoprime_gefunden' in ergebnis
        assert 'fazit' in ergebnis

    def test_geprueft_positiv(self):
        """@brief Mindestens eine Zahl wurde geprüft. @lastModified 2026-03-12"""
        ergebnis = numerische_suche_4prim(10_000)
        assert ergebnis['geprueft'] > 0


# =============================================================================
# Tests: schranke_fuer_4prim_fall_p_gleich_2
# =============================================================================


class TestSchrankeFall2:
    """
    @brief Testet die Schrankenanalyse für Fall p=2.
    @lastModified 2026-03-12
    """

    def test_kein_pseudoprime_im_suchbereich(self):
        """
        @brief Für q_max=50 wird kein 4-Prim-Pseudoprime der Form 2·q·r·s gefunden.
        @lastModified 2026-03-12
        """
        ergebnis = schranke_fuer_4prim_fall_p_gleich_2(q_max=50)
        assert ergebnis['kein_pseudoprime'] is True
        assert ergebnis['starke_kandidaten_pseudoprimes'] == []

    def test_rückgabe_enthält_fall_bezeichnung(self):
        """@brief Rückgabe enthält 'fall' Beschreibung. @lastModified 2026-03-12"""
        ergebnis = schranke_fuer_4prim_fall_p_gleich_2(q_max=20)
        assert 'fall' in ergebnis
        assert '2' in ergebnis['fall']

    def test_untere_schranke_s_kleiner_als_sqrt_2qr(self):
        """
        @brief Schranken-Beschreibung enthält das korrekte Argument s·(s-1) | (2qr-1).
        @description
            Aus der Analyse: s(s-1) | (2qr-1) → s ≤ √(2qr).
            Dies ist das zentrale Schranken-Argument für Fall p=2.
        @lastModified 2026-03-12
        """
        ergebnis = schranke_fuer_4prim_fall_p_gleich_2(q_max=20)
        assert 'schranken_beschreibung' in ergebnis
        # Das Argument muss in der Beschreibung erwähnt sein
        arg = ergebnis['schranken_beschreibung']['aus_III_und_III_strich']
        assert 's' in arg and 'q' in arg and 'r' in arg

    def test_schwache_kandidaten_erfüllen_schwache_bedingung(self):
        """
        @brief Alle gefundenen schwachen Kandidaten erfüllen die schwache Giuga-Bedingung.
        @lastModified 2026-03-12
        """
        ergebnis = schranke_fuer_4prim_fall_p_gleich_2(q_max=50)
        for kand in ergebnis['schwache_kandidaten']:
            n = kand['n']
            faktoren = [2, kand['q'], kand['r'], kand['s']]
            for p in faktoren:
                assert (n // p - 1) % p == 0, \
                    f"Schwache Bedingung verletzt: p={p}, n={n}"

    def test_q_max_parameter_wirkt(self):
        """@brief Mit größerem q_max werden mehr Kandidaten gesucht. @lastModified 2026-03-12"""
        ergebnis_klein = schranke_fuer_4prim_fall_p_gleich_2(q_max=20)
        ergebnis_gross = schranke_fuer_4prim_fall_p_gleich_2(q_max=60)
        # Größerer Suchraum → mindestens so viele schwache Kandidaten
        assert len(ergebnis_gross['schwache_kandidaten']) >= \
               len(ergebnis_klein['schwache_kandidaten'])


# =============================================================================
# Tests: Giuga4PrimBeweis Klasse
# =============================================================================


class TestGiuga4PrimBeweisKlasse:
    """
    @brief Testet die Klasse Giuga4PrimBeweis.
    @lastModified 2026-03-12
    """

    @pytest.fixture
    def beweis(self):
        """@brief Erstellt eine Instanz von Giuga4PrimBeweis. @lastModified 2026-03-12"""
        return Giuga4PrimBeweis()

    # --- fall1_p_gleich_2_analyse ---

    def test_fall1_status_enthält_offen(self, beweis):
        """@brief Fall 1 meldet Status 'offen' (kein vollständiger Beweis). @lastModified 2026-03-12"""
        ergebnis = beweis.fall1_p_gleich_2_analyse()
        assert 'offen' in ergebnis['status'].lower() or \
               'OFFEN' in ergebnis['status'] or \
               'kein vollständiger' in ergebnis['status'].lower()

    def test_fall1_kein_pseudoprime_numerisch(self, beweis):
        """@brief Numerischer Teil von Fall 1 findet keinen Pseudoprime. @lastModified 2026-03-12"""
        ergebnis = beweis.fall1_p_gleich_2_analyse()
        assert ergebnis['numerisch']['bestätigt'] is True
        assert ergebnis['numerisch']['pseudoprimes_gefunden'] == 0

    def test_fall1_schlüsselbedingungen_vorhanden(self, beweis):
        """@brief Fall 1 enthält Schlüsselbedingungen. @lastModified 2026-03-12"""
        ergebnis = beweis.fall1_p_gleich_2_analyse()
        assert 'schlüsselbedingungen' in ergebnis
        assert 'aus_s' in ergebnis['schlüsselbedingungen']

    # --- fall2_p_gleich_3_analyse ---

    def test_fall2_kein_pseudoprime_numerisch(self, beweis):
        """@brief Fall 2 findet keinen Pseudoprime numerisch. @lastModified 2026-03-12"""
        ergebnis = beweis.fall2_p_gleich_3_analyse()
        assert ergebnis['numerisch']['bestätigt'] is True

    def test_fall2_enthält_kongruenz_analyse(self, beweis):
        """@brief Fall 2 enthält Kongruenz-Analyse (mod 3). @lastModified 2026-03-12"""
        ergebnis = beweis.fall2_p_gleich_3_analyse()
        assert 'kongruenz_analyse' in ergebnis

    # --- fall3_p_ungerade_analyse ---

    def test_fall3_kein_pseudoprime_numerisch(self, beweis):
        """@brief Fall 3 (p≥5) findet keinen Pseudoprime. @lastModified 2026-03-12"""
        ergebnis = beweis.fall3_p_ungerade_analyse()
        assert ergebnis['numerisch']['bestätigt'] is True

    def test_fall3_enthält_erklärung_warum_3prim_versagt(self, beweis):
        """
        @brief Fall 3 erklärt warum der 3-Prim-Beweisansatz nicht greift.
        @lastModified 2026-03-12
        """
        ergebnis = beweis.fall3_p_ungerade_analyse()
        assert 'warum_3prim_beweis_versagt' in ergebnis
        erkl = ergebnis['warum_3prim_beweis_versagt']
        assert len(erkl) > 50  # substantielle Erklärung

    # --- schranken_analyse ---

    def test_schranken_enthält_externe_literatur(self, beweis):
        """@brief Schrankenanalyse referenziert externe Literatur (Borwein, Bednarek). @lastModified 2026-03-12"""
        ergebnis = beweis.schranken_analyse()
        assert 'externe_schranken_literatur' in ergebnis
        literatur = ergebnis['externe_schranken_literatur']
        assert any('Borwein' in k for k in literatur.keys())
        assert any('Bednarek' in k for k in literatur.keys())

    def test_schranken_enthält_eigene_resultate(self, beweis):
        """@brief Schrankenanalyse nennt eigene Resultate aus beweisversuche.py. @lastModified 2026-03-12"""
        ergebnis = beweis.schranken_analyse()
        assert 'eigene_resultate' in ergebnis
        eigene = ergebnis['eigene_resultate']
        assert len(eigene) >= 3  # Mindestens Sätze 1, 2, Korollar

    def test_schranken_beispiel_n_1365(self, beweis):
        """
        @brief Schrankenanalyse enthält ein Beispiel (3·5·7·s-Kandidat).
        @description
            3·5·7 - 1 = 104. Die Primteiler > 7 von 104 sind: 13.
            Also s=13 ist der einzige Primzahlteiler > 7 von 104.
            n = 3·5·7·13 = 1365 wäre ein Kandidat (aber kein Pseudoprime).
        @lastModified 2026-03-12
        """
        ergebnis = beweis.schranken_analyse()
        assert 'beispiel_schranke' in ergebnis
        bsp = ergebnis['beispiel_schranke']
        assert bsp['pqr_minus_1'] == 104  # 3*5*7 - 1 = 104
        assert 13 in bsp['s_kandidaten']  # 13 ist Primteiler > 7 von 104

    # --- numerische_verifikation ---

    def test_numerische_verifikation_bis_500k(self, beweis):
        """@brief Numerische Verifikation bis 500.000 findet keinen Pseudoprime. @lastModified 2026-03-12"""
        ergebnis = beweis.numerische_verifikation(500_000)
        assert ergebnis['kein_pseudoprime_gefunden'] is True

    # --- korollar_kein_4prim_pseudoprime ---

    def test_korollar_status_ist_offen(self, beweis):
        """@brief Korollar meldet den Status als offen (noch kein vollständiger Beweis). @lastModified 2026-03-12"""
        korollar = beweis.korollar_kein_4prim_pseudoprime()
        assert 'OFFEN' in korollar['status'] or 'offen' in korollar['status'].lower()

    def test_korollar_nennt_vorherige_sätze(self, beweis):
        """@brief Korollar referenziert die bewiesenen Sätze 2–4 aus beweisversuche.py. @lastModified 2026-03-12"""
        korollar = beweis.korollar_kein_4prim_pseudoprime()
        assert 'vorherige_sätze' in korollar
        saetze = korollar['vorherige_sätze']
        assert 'satz_2' in saetze
        assert 'satz_3' in saetze
        assert 'satz_4' in saetze

    def test_korollar_enthält_empfehlung(self, beweis):
        """@brief Korollar enthält eine Empfehlung für weitere Arbeit. @lastModified 2026-03-12"""
        korollar = beweis.korollar_kein_4prim_pseudoprime()
        assert 'empfehlung' in korollar
        assert len(korollar['empfehlung']) > 30

    def test_korollar_erwähnt_literaturschranken(self, beweis):
        """@brief Korollar nennt externe Schranken (Borwein, Bednarek). @lastModified 2026-03-12"""
        korollar = beweis.korollar_kein_4prim_pseudoprime()
        assert 'literatur_schranken' in korollar


# =============================================================================
# Tests: Edge-Cases und Integration
# =============================================================================


class TestEdgeCases:
    """
    @brief Edge-Case-Tests und Integrationstests.
    @lastModified 2026-03-12
    """

    def test_giuga_bedingung_für_alle_bekannten_giuga_zahlen(self):
        """
        @brief Alle vier bekannten Giuga-Zahlen werden korrekt als solche erkannt.
        @description
            Die bekannten Giuga-Zahlen {30, 858, 1722, 66198} müssen:
            - ist_giuga_zahl = True  (schwache Bedingung)
            - ist_giuga_pseudoprime = False (nicht beide Bedingungen)
        @lastModified 2026-03-12
        """
        bekannte_giuga_zahlen = [30, 858, 1722, 66198]
        for n in bekannte_giuga_zahlen:
            ergebnis = berechne_giuga_bedingung(n)
            assert ergebnis['ist_giuga_zahl'] is True, \
                f"n={n} sollte Giuga-Zahl sein"
            assert ergebnis['ist_giuga_pseudoprime'] is False, \
                f"n={n} sollte KEIN Giuga-Pseudoprime sein"

    def test_858_hat_4_primfaktoren(self):
        """@brief 858 = 2·3·11·13 hat genau 4 Primfaktoren. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(858)
        assert len(ergebnis['primfaktoren']) == 4
        assert sorted(ergebnis['primfaktoren']) == [2, 3, 11, 13]

    def test_1722_hat_4_primfaktoren(self):
        """@brief 1722 = 2·3·7·41 hat genau 4 Primfaktoren. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(1722)
        assert len(ergebnis['primfaktoren']) == 4
        assert sorted(ergebnis['primfaktoren']) == [2, 3, 7, 41]

    def test_30_hat_3_primfaktoren(self):
        """@brief 30 = 2·3·5 hat genau 3 Primfaktoren. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(30)
        assert len(ergebnis['primfaktoren']) == 3

    def test_numerische_suche_bis_210_nur_bekannte_kandidaten(self):
        """
        @brief Bis n=210 (= 2·3·5·7, kleinstes 4-Prim-Produkt) gibt es nur wenige Kandidaten.
        @lastModified 2026-03-12
        """
        ergebnis = numerische_suche_4prim(210)
        # 210 = 2·3·5·7 könnte ein Kandidat sein (schwach prüfen)
        # Es soll kein Pseudoprime vorhanden sein
        assert ergebnis['kein_pseudoprime_gefunden'] is True

    def test_210_schwache_bedingung_direkt(self):
        """
        @brief 210 = 2·3·5·7 direkte Prüfung: kein Giuga-Pseudoprime.
        @description
            Prüfe manuell für n=210, p=2: 2|(105-1)=104? 104/2=52 ✓
            p=3: 3|(70-1)=69? 69/3=23 ✓
            p=5: 5|(42-1)=41? 41/5=8.2 ✗ → keine Giuga-Zahl
        @lastModified 2026-03-12
        """
        ergebnis = berechne_giuga_bedingung(210)
        # 210/5 - 1 = 42 - 1 = 41, und 41 % 5 = 1 ≠ 0 → schwache Bedingung verletzt
        assert ergebnis['ist_giuga_zahl'] is False
        assert ergebnis['ist_giuga_pseudoprime'] is False

    def test_schranke_fall_2_gibt_korrekte_fall_bezeichnung(self):
        """@brief Fall-2-Funktion gibt korrekte Bezeichnung zurück. @lastModified 2026-03-12"""
        ergebnis = schranke_fuer_4prim_fall_p_gleich_3(q_max=20)
        assert 'fall' in ergebnis
        assert '3' in ergebnis['fall']

    def test_allgemeiner_fall_p_max(self):
        """@brief Allgemeiner Fall schließt korrekt bis p_max ab. @lastModified 2026-03-12"""
        ergebnis = schranke_fuer_4prim_allgemeiner_fall(p_max=11)
        assert ergebnis['p_max'] == 11
        assert ergebnis['kein_pseudoprime'] is True

    def test_giuga_bedingung_kein_crash_für_große_n(self):
        """@brief berechne_giuga_bedingung() stürzt nicht bei großem n ab. @lastModified 2026-03-12"""
        # n = 2·3·5·7·11·13 = 30030 (6 Primfaktoren)
        ergebnis = berechne_giuga_bedingung(30030)
        assert 'ist_giuga_zahl' in ergebnis
        assert 'ist_giuga_pseudoprime' in ergebnis

    def test_giuga_bedingung_für_n_1(self):
        """@brief n=1 ist kein Giuga-Kandidat. @lastModified 2026-03-12"""
        ergebnis = berechne_giuga_bedingung(1)
        assert ergebnis['ist_giuga_zahl'] is False
        assert ergebnis['ist_giuga_pseudoprime'] is False
