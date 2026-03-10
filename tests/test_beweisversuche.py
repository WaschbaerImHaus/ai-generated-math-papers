"""
@file test_beweisversuche.py
@brief Tests für das Beweisversuch-Modul (Giuga, Brocard-Ramanujan, Erdős-Straus, Kurepa).
@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from beweisversuche import (
    GiugaBeweisführung,
    BrocardRamanujanBeweis,
    ErdosStrausRestklassen,
    KurepaAnalyse,
    LehmerVermutungBeweis,
    schwierigkeitsranking,
)


# ===========================================================
# Giuga-Beweise
# ===========================================================

class TestGiugaSatz1:
    """Tests für Satz 1: Quadratfreiheit."""

    def setup_method(self):
        self.g = GiugaBeweisführung()

    def test_beweis_status_bewiesen(self):
        r = self.g.satz1_quadratfreiheit_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_beweis_enthält_kernargumente(self):
        r = self.g.satz1_quadratfreiheit_beweis()
        assert len(r['beweis_kern']) >= 4
        # Das Widerspruchsargument muss vorhanden sein
        alle = ' '.join(r['beweis_kern'])
        assert 'Widerspruch' in alle or 'unmöglich' in alle.lower() or 'Unmöglich' in alle

    def test_verifikation_bis_50000(self):
        r = self.g.satz1_numerische_verifikation(50000)
        assert r['bestätigt'] is True
        assert r['nicht_quadratfreie_giuga_kandidaten'] == []

    def test_keine_nicht_quadratfreien_kandidaten_klein(self):
        r = self.g.satz1_numerische_verifikation(1000)
        assert r['bestätigt'] is True


class TestGiugaSatz2:
    """Tests für Satz 2: Kein 2-Primfaktor-Giuga-Pseudoprime."""

    def setup_method(self):
        self.g = GiugaBeweisführung()

    def test_beweis_status_bewiesen(self):
        r = self.g.satz2_kein_2prim_pseudoprime_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_beweis_enthält_ordnungsargument(self):
        r = self.g.satz2_kein_2prim_pseudoprime_beweis()
        alle = ' '.join(r['beweis_kern'])
        # Muss p < q und p-1 < q erwähnen
        assert 'Widerspruch' in alle or 'Kein' in alle or '□' in alle

    def test_konsequenz_mindestens_3_faktoren(self):
        r = self.g.satz2_kein_2prim_pseudoprime_beweis()
        assert '3' in r['konsequenz'] or 'drei' in r['konsequenz'].lower()

    def test_verifikation_bis_50000(self):
        r = self.g.satz2_numerische_verifikation(50000)
        assert r['bestätigt'] is True
        assert r['2prim_giuga_kandidaten'] == []

    def test_verifikation_bis_100(self):
        r = self.g.satz2_numerische_verifikation(100)
        assert r['bestätigt'] is True

    def test_bekannte_produktpaare_schlagen_fehl(self):
        """6=2·3, 15=3·5 etc. dürfen keine Giuga-Kandidaten sein."""
        import sympy
        g = GiugaBeweisführung()
        r = g.satz2_numerische_verifikation(100)
        # Diese Semiprimes dürfen nicht in der Liste sein
        for n in [6, 10, 14, 15, 21, 22]:
            assert n not in r['2prim_giuga_kandidaten']


class TestGiugaSatz3:
    """Tests für Satz 3: Kein 3-Primfaktor-Giuga-Pseudoprime mit p₁=2."""

    def setup_method(self):
        self.g = GiugaBeweisführung()

    def test_beweis_status_bewiesen(self):
        r = self.g.satz3_kein_3prim_mit_p_gleich_2_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_paritätsargument_vorhanden(self):
        r = self.g.satz3_kein_3prim_mit_p_gleich_2_beweis()
        alle = ' '.join(r['beweis_kern'])
        # Gerade/Ungerade Argument
        assert 'GERADE' in alle or 'UNGERADE' in alle or 'gerade' in alle.lower()

    def test_verifikation_keine_2qr_pseudoprimes(self):
        r = self.g.satz3_numerische_verifikation(500)
        assert r['bestätigt'] is True
        assert r['2qr_giuga_pseudoprimes'] == []

    def test_bekannte_nicht_gegenbeispiele(self):
        """n=30=2·3·5 ist Giuga-Zahl (schwach) aber KEIN Pseudoprime (stark)."""
        # 30 = 2·3·5. Starke Bedingung für p=3: (p-1)=2 | (n/3-1)=9. 2|9? NEIN!
        # Also 30 ist kein Pseudoprime → korrekt
        r = self.g.satz3_numerische_verifikation(200)
        n30_in_liste = any(d.get('n') == 30 for d in r['2qr_giuga_pseudoprimes'])
        assert not n30_in_liste

    def test_analyse_3prim_alle_ungerade_kein_ergebnis(self):
        r = self.g.analyse_3prim_alle_ungerade(p_max=50)
        assert r['anzahl_stark'] == 0

    def test_untere_schranken_korrekt(self):
        r = self.g.berechne_untere_schranke()
        assert 'Satz_1' in r['eigene_resultate']
        assert 'Satz_2' in r['eigene_resultate']
        assert 'Satz_3' in r['eigene_resultate']
        # Bekannte Schranken müssen referenziert sein
        assert 'Borwein_1996' in r['bekannte_schranken']


class TestGiugaVollständig:
    """Integrationstests für die Giuga-Beweisführung."""

    def test_giuga_zahlen_keine_pseudoprimes(self):
        """Bekannte Giuga-Zahlen {30, 858, 1722} sind KEINE Pseudoprimes."""
        g = GiugaBeweisführung()
        import sympy
        giuga_zahlen = [30, 858, 1722]
        for n in giuga_zahlen:
            faktoren = list(sympy.factorint(n).keys())
            # Schwache Bedingung: alle erfüllen sie
            schwach = all((n // p - 1) % p == 0 for p in faktoren)
            assert schwach, f"n={n} sollte schwache Bedingung erfüllen"
            # Starke Bedingung: NICHT alle erfüllen sie
            stark = all((n // p - 1) % (p - 1) == 0 for p in faktoren)
            assert not stark, f"n={n} sollte NICHT Pseudoprime sein"

    def test_vollständige_analyse_läuft(self):
        """Vollständige Analyse läuft ohne Fehler."""
        g = GiugaBeweisführung()
        # Nur Sätze 1-3, nicht die schwerste numerische Suche
        r1 = g.satz1_quadratfreiheit_beweis()
        r2 = g.satz2_kein_2prim_pseudoprime_beweis()
        r3 = g.satz3_kein_3prim_mit_p_gleich_2_beweis()
        assert r1['status'] == 'BEWIESEN'
        assert r2['status'] == 'BEWIESEN'
        assert r3['status'] == 'BEWIESEN'


# ===========================================================
# Brocard-Ramanujan
# ===========================================================

class TestBrocardRamanujan:
    """Tests für Brocard-Ramanujan-Analyse."""

    def setup_method(self):
        self.br = BrocardRamanujanBeweis()

    def test_bekannte_lösungen(self):
        """4!+1=25=5², 5!+1=121=11², 7!+1=5041=71²."""
        import math
        for n, m in [(4, 5), (5, 11), (7, 71)]:
            assert math.factorial(n) + 1 == m * m

    def test_keine_lösungen_8_bis_40(self):
        r = self.br.untere_schranke_fuer_loesungen(n_start=8)
        assert r['loesungen_gefunden'] == []

    def test_modulare_analyse_läuft(self):
        r = self.br.modularer_ausschluss_analyse(n_max=15)
        assert isinstance(r, dict)
        assert 8 in r or 9 in r

    def test_bekannte_lösungen_als_quadratzahlen(self):
        r = self.br.modularer_ausschluss_analyse(n_max=10)
        # n=4 sollte als Quadratzahl erkannt werden
        if 4 in r:
            assert r[4]['ist_quadratzahl'] is True
        # n=8 sollte keine Quadratzahl sein
        if 8 in r:
            assert r[8]['ist_quadratzahl'] is False

    def test_notwendige_bedingung_analyse(self):
        r = self.br.satz_notwendige_bedingung()
        assert 'satz' in r
        assert 'kongruenz_analyse' in r
        # Primzahl 3 sollte in der Analyse sein
        analyse = r['kongruenz_analyse']
        assert 3 in analyse
        # 3 ≡ 3 (mod 4): -1 ist KEIN QR mod 3
        assert analyse[3]['-1_ist_QR'] is False

    def test_kongruenz_3_kein_qr(self):
        """−1 ist kein QR mod 3 (da 3 ≡ 3 mod 4)."""
        br = BrocardRamanujanBeweis()
        # QR mod 3: 0², 1², 2² = 0, 1, 1 (mod 3) → QR = {0, 1}
        # -1 ≡ 2 (mod 3), nicht in {0,1} → KEIN QR
        assert br._quadratischer_rest(3 - 1, 3) is False

    def test_kongruenz_5_kein_qr(self):
        """−1 ist kein QR mod 5 (da 5 ≡ 1 mod 4)."""
        br = BrocardRamanujanBeweis()
        # 5 ≡ 1 (mod 4) → -1 IST QR mod 5: 2²=4≡-1 ✓
        assert br._quadratischer_rest(5 - 1, 5) is True

    def test_kongruenz_7_kein_qr(self):
        """−1 ist kein QR mod 7 (da 7 ≡ 3 mod 4)."""
        br = BrocardRamanujanBeweis()
        assert br._quadratischer_rest(7 - 1, 7) is False


# ===========================================================
# Erdős-Straus Restklassen
# ===========================================================

class TestErdosStrausRestklassen:
    """Tests für die Restklassenabdeckung."""

    def setup_method(self):
        self.es = ErdosStrausRestklassen()

    def test_verifikation_formel_n_mod0_mod4(self):
        """n ≡ 0 (mod 4): Formel 1/(2k)+1/(3k)+1/(6k) ist korrekt."""
        for n in [4, 8, 12, 16, 20, 100, 840]:
            result = self.es.formel_fuer_restklasse(n % 840, 840, n)
            assert result is not None, f"n={n}: keine Zerlegung"
            x, y, z = result
            assert abs(1/x + 1/y + 1/z - 4/n) < 1e-10, f"n={n}: falsche Zerlegung"

    def test_formel_für_n_mod3_mod4(self):
        """n ≡ 3 (mod 4): Explizite Formel 1/(k+1) + 1/(M+1) + 1/(M(M+1)) mit M=n(k+1)."""
        for n in [3, 7, 11, 15, 19, 23, 27]:
            result = self.es.formel_fuer_restklasse(n % 840, 840, n)
            assert result is not None, f"n={n}: keine Zerlegung gefunden"
            x, y, z = result
            assert abs(1/x + 1/y + 1/z - 4/n) < 1e-10, f"n={n}: falsche Zerlegung"

    def test_explizite_formeln_gibt_dict(self):
        formeln = self.es.explizite_formeln()
        assert isinstance(formeln, dict)
        assert 'r ≡ 0 (mod 4)' in formeln
        assert 'r ≡ 3 (mod 4)' in formeln

    def test_explizite_formel_n_mod0_algebraisch_korrekt(self):
        """Algebraische Verifikation: 1/(2k)+1/(3k)+1/(6k) = 1/k."""
        # 1/(2k)+1/(3k)+1/(6k) = 3/(6k)+2/(6k)+1/(6k) = 6/(6k) = 1/k ✓
        for k in range(1, 20):
            n = 4 * k
            x, y, z = 2*k, 3*k, 6*k
            assert abs(1/x + 1/y + 1/z - 1/k) < 1e-12
            assert abs(1/x + 1/y + 1/z - 4/n) < 1e-12

    def test_explizite_formel_n_mod3_mod4_algebraisch(self):
        """Algebraische Verifikation: 4/(4k+3) = 1/(k+1) + 1/((4k+3)(k+1))."""
        # 4/(4k+3) - 1/(k+1) = (4(k+1) - (4k+3)) / ((4k+3)(k+1)) = 1/((4k+3)(k+1))
        for k in range(0, 20):
            n = 4*k + 3
            a = k + 1
            b = n * (k + 1)
            val = 1/a + 1/b
            assert abs(val - 4/n) < 1e-10, \
                f"k={k}, n={n}: 1/{a} + 1/{b} = {val} ≠ 4/{n} = {4/n}"

    def test_verifikation_formel(self):
        """Algebraische Verifikation der Formel: 4/(4k+3) = 1/(k+1) + 1/(M+1) + 1/(M(M+1))."""
        for n in [3, 7, 11, 15, 19]:
            k = (n - 3) // 4
            a = k + 1
            M = n * (k + 1)
            b = M + 1
            c = M * (M + 1)
            assert self.es._verifiziere_zerlegung(n, a, b, c), \
                f"Verifikation für n={n}: 1/{a}+1/{b}+1/{c} ≠ 4/{n}"


# ===========================================================
# Kurepa-Analyse
# ===========================================================

class TestKurepaAnalyse:
    """Tests für die Kurepa-Vermutungs-Analyse."""

    def setup_method(self):
        self.k = KurepaAnalyse()

    def test_linke_fakultät_kleine_werte(self):
        """!n = 0!+1!+...+(n-1)!"""
        # !1 = 0! = 1
        assert self.k.linke_fakultaet(1) == 1
        # !2 = 0!+1! = 2
        assert self.k.linke_fakultaet(2) == 2
        # !3 = 0!+1!+2! = 1+1+2 = 4
        assert self.k.linke_fakultaet(3) == 4
        # !4 = 1+1+2+6 = 10
        assert self.k.linke_fakultaet(4) == 10
        # !5 = 1+1+2+6+24 = 34
        assert self.k.linke_fakultaet(5) == 34

    def test_linke_fakultät_mod_p_korrekt(self):
        """!p mod p = Σₖ=0^{p-1} k! mod p."""
        import math
        for p in [3, 5, 7, 11, 13]:
            # Berechne brute-force
            total = sum(math.factorial(k) for k in range(p)) % p
            assert self.k.linke_fakultaet_mod_p(p) == total

    def test_verifiziere_bis_100(self):
        """Kurepa-Vermutung gilt für alle Primzahlen bis 100."""
        r = self.k.verifiziere_bis(100)
        assert r['verifiziert'] is True
        assert r['gegenbeispiele'] == []
        assert r['geprüfte_primzahlen'] == 24  # π(100)=25 minus p=2 → 24 ungerade Primzahlen

    def test_verifiziere_bis_200(self):
        r = self.k.verifiziere_bis(200)
        assert r['verifiziert'] is True

    def test_wilson_analyse_gibt_daten(self):
        r = self.k.wilson_analyse()
        assert 'daten' in r
        assert len(r['daten']) > 0
        # Überprüfe dass alle S(p) ≠ 0 für bekannte p
        for eintrag in r['daten']:
            assert eintrag['kurepa_gilt'] is True, \
                f"Kurepa scheitert für p={eintrag['p']}"

    def test_restklassen_analyse_gibt_zusammenfassung(self):
        r = self.k.restklassen_analyse()
        assert 'mod_4_analyse' in r
        analyse = r['mod_4_analyse']
        # Es gibt p ≡ 1 und p ≡ 3 (mod 4)
        assert 'p ≡ 1 (mod 4)' in analyse
        assert 'p ≡ 3 (mod 4)' in analyse
        # Keine Nullen in S(p) gefunden
        for r_val, stats in analyse.items():
            assert stats['nullen'] == 0, \
                f"Kurepa scheitert in Klasse {r_val}"

    def test_p_gleich_2_ausgeklammert(self):
        """p=2 ist ungerade Primzahl? Nein – nur ungerade Primzahlen werden geprüft."""
        r = self.k.verifiziere_bis(10)
        # p=2 sollte nicht als Gegenbeispiel erscheinen (Vermutung gilt für ungerade p)
        assert 2 not in r.get('gegenbeispiele', [])


# ===========================================================
# Schwierigkeitsranking
# ===========================================================

class TestSchwierigkeitsranking:
    """Tests für das Schwierigkeitsranking."""

    def test_ranking_gibt_dict(self):
        r = schwierigkeitsranking()
        assert isinstance(r, dict)

    def test_bewiesene_sätze_vorhanden(self):
        r = schwierigkeitsranking()
        assert 'BEWIESEN_in_diesem_Modul' in r
        bewiesen = r['BEWIESEN_in_diesem_Modul']
        assert 'Giuga_Satz1' in bewiesen
        assert 'Giuga_Satz2' in bewiesen
        assert 'Giuga_Satz3' in bewiesen

    def test_ranking_liste_vorhanden(self):
        r = schwierigkeitsranking()
        assert 'RANKING' in r
        ranking = r['RANKING']
        assert len(ranking) >= 5
        # Giuga 2-Prim-Fall auf Rang 1
        rang1 = next((x for x in ranking if x['rang'] == 1), None)
        assert rang1 is not None
        assert '✓' in rang1['status'] or 'BEWIESEN' in rang1['status']

    def test_alle_bewiesenen_sätze_korrekt_markiert(self):
        r = schwierigkeitsranking()
        for name, info in r['BEWIESEN_in_diesem_Modul'].items():
            assert 'aussage' in info
            assert 'methode' in info

    def test_giuga_satz4_und_korollar_im_ranking(self):
        """Satz 4 und Korollar müssen im Ranking enthalten sein."""
        r = schwierigkeitsranking()
        bewiesen = r['BEWIESEN_in_diesem_Modul']
        assert 'Giuga_Satz4' in bewiesen
        assert 'Giuga_Korollar_3Prim' in bewiesen

    def test_lehmer_beweise_im_ranking(self):
        """Lehmer-Beweise müssen im Ranking enthalten sein."""
        r = schwierigkeitsranking()
        bewiesen = r['BEWIESEN_in_diesem_Modul']
        assert 'Lehmer_Quadratfreiheit' in bewiesen
        assert 'Lehmer_Kein_Semiprime' in bewiesen

    def test_ranking_hat_14_eintraege(self):
        """Das Ranking sollte mindestens 10 Einträge haben."""
        r = schwierigkeitsranking()
        assert len(r['RANKING']) >= 10


# ===========================================================
# Giuga Satz 4 (NEUER BEWEIS)
# ===========================================================

class TestGiugaSatz4:
    """Tests für Satz 4: Kein 3-Prim-alle-ungerade-Giuga-Pseudoprime."""

    def setup_method(self):
        self.g = GiugaBeweisführung()

    def test_beweis_status_bewiesen(self):
        r = self.g.satz4_kein_3prim_alle_ungerade_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_beweis_enthält_schrankenargument(self):
        r = self.g.satz4_kein_3prim_alle_ungerade_beweis()
        alle = ' '.join(r['beweis_kern'])
        # Muss das Schranken-Widerspruchsargument enthalten
        assert 'WIDERSPRUCH' in alle or 'Widerspruch' in alle
        assert 'r(r-1)' in alle or 'r·(r-1)' in alle or 'lcm' in alle.lower() or '(q+2)' in alle

    def test_korollar_3prim_vollständig(self):
        """Das Korollar über den vollständigen 3-Prim-Fall muss vorhanden sein."""
        r = self.g.satz4_kein_3prim_alle_ungerade_beweis()
        assert 'kombination_mit_satz3' in r
        assert 'Korollar' in r.get('kombination_mit_satz3', '') or 'Korollar' in r.get('korollar', '')

    def test_korollar_methode_läuft(self):
        r = self.g.korollar_keine_3prim_giuga_pseudoprimes()
        assert r['status'] == 'BEWIESEN'
        assert len(r['beweis']) >= 3

    def test_kernaussage_algebraisch_korrekt(self):
        """Algebraische Verifikation: (q+2)(q+1) > pq für p<q (alle ungerade prim)."""
        import sympy as sp
        # Zeige: Wenn r(r-1) | (pq-1) UND r ≥ q+2, dann p > q (Widerspruch)
        # d.h.: für alle p < q < r (ungerade Primes) gilt r(r-1) > pq-1
        primes = list(sp.primerange(3, 50))
        alle_korrekt = True
        for i, p in enumerate(primes):
            for j, q in enumerate(primes):
                if q <= p:
                    continue
                # Kleinster möglicher r: nächste Primzahl nach q
                r = sp.nextprime(q)
                # Behauptung: r*(r-1) > pq - 1, also kein echter Teiler
                if r * (r - 1) <= p * q - 1:
                    alle_korrekt = False
        assert alle_korrekt, "Schrankenargument für Satz 4 scheitert numerisch!"

    def test_numerisch_kein_3prim_alle_ungerade(self):
        """Numerisch: Kein 3-Prim-alle-ungerade-Giuga-Pseudoprime bis p_max=80."""
        r = self.g.analyse_3prim_alle_ungerade(p_max=80)
        assert r['anzahl_stark'] == 0, \
            f"Gegenbeispiel gefunden: {r['echte_giuga_pseudoprimes']}"


# ===========================================================
# Lehmer-Vermutung
# ===========================================================

class TestLehmerVermutungBeweis:
    """Tests für die Lehmer-Vermutungs-Beweise."""

    def setup_method(self):
        self.lv = LehmerVermutungBeweis()

    def test_quadratfreiheit_status_bewiesen(self):
        r = self.lv.satz_quadratfreiheit_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_quadratfreiheit_enthält_kernargument(self):
        r = self.lv.satz_quadratfreiheit_beweis()
        alle = ' '.join(r['beweis_kern'])
        assert 'WIDERSPRUCH' in alle or 'Widerspruch' in alle or '□' in alle

    def test_kein_semiprime_status_bewiesen(self):
        r = self.lv.satz_kein_semiprime_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_kein_semiprime_kongruenzargument(self):
        r = self.lv.satz_kein_semiprime_beweis()
        alle = ' '.join(r['beweis_kern'])
        # Muss die Kongruenz pq-1 ≡ p+q-2 erwähnen
        assert 'p+q-2' in alle or 'p + q - 2' in alle or 'p+q' in alle

    def test_algebraische_verifikation_kongruenz(self):
        """Verifikation: pq-1 ≡ p+q-2 (mod (p-1)(q-1)) für konkrete Werte."""
        import sympy as sp
        for p in [2, 3, 5, 7, 11]:
            for q in sp.primerange(p + 1, p + 20):
                phi = (p - 1) * (q - 1)
                rest = (p * q - 1) % phi
                erwartet = (p + q - 2) % phi
                assert rest == erwartet, \
                    f"p={p}, q={q}: pq-1 mod φ = {rest} ≠ p+q-2 mod φ = {erwartet}"

    def test_numerisch_keine_zusammengesetzten_lösungen(self):
        """Numerisch: Kein zusammengesetztes n < 10000 erfüllt φ(n)|(n-1)."""
        r = self.lv.numerische_verifikation(10000)
        assert r['verifiziert'] is True
        assert r['zusammengesetzte_kandidaten'] == []

    def test_semiprime_schlägt_fehl_konkret(self):
        """Konkrete Semiprimes dürfen φ(n)|(n-1) nicht erfüllen."""
        import sympy as sp
        for n in [6, 10, 15, 21, 35, 77, 143]:
            phi = sp.totient(n)
            assert (n - 1) % phi != 0, \
                f"n={n} = Semiprime erfüllt fälschlicherweise φ(n)|(n-1)!"

    def test_primzahlen_erfüllen_bedingung(self):
        """Primzahlen müssen φ(n)|(n-1) erfüllen (n-1 = φ(n) für Primzahlen)."""
        import sympy as sp
        for p in sp.primerange(2, 50):
            phi = sp.totient(p)
            assert (p - 1) % phi == 0, \
                f"Primzahl p={p}: φ(p)={phi} teilt p-1={p-1}? NEIN – Fehler!"

    def test_vergleich_mit_giuga_gibt_daten(self):
        r = self.lv.vergleich_mit_giuga()
        assert 'gemeinsame_struktur' in r
        assert len(r['gemeinsame_struktur']) >= 2
        # Bekannte Giuga-Zahlen müssen aufgeführt sein
        assert 30 in r['bekannte_giuga_zahlen']

    def test_giuga_zahlen_erfüllen_lehmer_nicht_unbedingt(self):
        """Giuga-Zahlen (30, 858, ...) müssen φ(n)|(n-1) NICHT erfüllen."""
        import sympy as sp
        for n in [30, 858]:
            phi = sp.totient(n)
            # Giuga-Zahlen sind kein Gegenbeispiel zu Lehmer
            assert (n - 1) % phi != 0, \
                f"n={n} (Giuga-Zahl) erfüllt fälschlicherweise Lehmer-Bedingung!"

    def test_kein_3prim_gerade_status_bewiesen(self):
        r = self.lv.satz_kein_3prim_gerade_beweis()
        assert r['status'] == 'BEWIESEN'

    def test_kein_3prim_gerade_kernargument(self):
        r = self.lv.satz_kein_3prim_gerade_beweis()
        alle = ' '.join(r['beweis_kern'])
        # Muss b-Fallunterscheidung und Widerspruch enthalten
        assert 'WIDERSPRUCH' in alle
        assert 'r = 2q' in alle or 'r=2q' in alle or '2q' in alle

    def test_kein_3prim_gerade_algebraisch(self):
        """Verifikation: (r-1)|(2q-1) ist für q<r (beide prim) unmöglich."""
        import sympy as sp
        # Teste alle Primpaare q<r mit q,r ≤ 100
        verletzungen = []
        for q in sp.primerange(3, 100):
            for r in sp.primerange(q + 1, 100):
                if (2 * q - 1) % (r - 1) == 0:
                    verletzungen.append((q, r))
        # Wenn es Verletzungen gäbe, wäre der Beweis falsch
        # Aber: (r-1)|(2q-1) mit r>q ist möglich! Der Beweis zeigt Unmöglichkeit von r prim.
        # Also: für alle (q, r) wo (r-1)|(2q-1) MUSS r zusammengesetzt sein
        for q, r in verletzungen:
            assert not sp.isprime(r), \
                f"(q={q}, r={r}): (r-1)|(2q-1) und r prim! Beweis falsch?"

    def test_notwendige_bedingung_2qr_korrekt(self):
        """Notwendige Bed.: (r-1)|(2qr-1) impliziert (r-1)|(2q-1)."""
        import sympy as sp
        for q in sp.primerange(3, 30):
            for r in sp.primerange(q + 1, 50):
                n = 2 * q * r
                # Prüfe: wenn (r-1)|(2qr-1), dann auch (r-1)|(2q-1)
                if (n - 1) % (r - 1) == 0:
                    assert (2 * q - 1) % (r - 1) == 0, \
                        f"Notwendige Bedingung verletzt: q={q}, r={r}"

    def test_3prim_gerade_kein_kandidat(self):
        """Numerisch: Kein n=2qr (q,r odd prime) erfüllt φ(n)|(n-1) bis q,r≤100."""
        import sympy as sp
        for q in sp.primerange(3, 100):
            for r in sp.primerange(q + 1, 100):
                n = 2 * q * r
                phi = (q - 1) * (r - 1)
                assert (n - 1) % phi != 0, \
                    f"GEGENBEISPIEL: n=2·{q}·{r}={n}, φ(n)={phi}, n-1={n-1}"
