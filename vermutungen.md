# Unbewiesene und unwiderlegte mathematische Vermutungen

> Möglichst umfassende Sammlung – sortiert nach Alter, dann nach Gebiet.
> Stand: 2026-03-10 | Quellen: Wikipedia, MathWorld, Clay Mathematics Institute, AMS, OEIS, arXiv
>
> **Ziel:** Als Trainingsgrundlage für Beweiswerkzeuge – viele kleinere Vermutungen
> sind möglicherweise mit bestehenden Methoden lösbar!

## Kategorien-Legende

| Kat. | Bedeutung | Zeithorizont |
|------|-----------|--------------|
| 🟢 **A** | Hohe Wahrscheinlichkeit – klare Angriffspunkte, aktive Forschung | 10–50 Jahre |
| 🟡 **B** | Mittlere Wahrscheinlichkeit – grundlegend neue Ideen nötig | 50–100 Jahre |
| 🔴 **C** | Niedrige Wahrscheinlichkeit – möglicherweise unentscheidbar | >100 Jahre |
| ⬛ **D** | Formal möglicherweise ZFC-unentscheidbar | Unbestimmbar |
| ✓ | Bereits bewiesen | – |
| – | Historisch / widerlegt | – |

---

## I. Zahlentheorie

### I.1 Primzahlen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 1 | ✓ | **Primzahlsatz** π(x) ~ x/ln(x) | 1792 | Hadamard/de la Vallée-Poussin 1896 |
| 2 | ✓ | **Bertrand-Postulat**: Primzahl zwischen n und 2n | 1845 | Tschebyschow 1850 |
| 3 | 🟡 B | **Riemann-Hypothese**: Alle Nicht-trivialen Nullstellen von ζ(s) auf Re=1/2 | 1859 | Millennium-Problem (1M$); 10¹³ Nullst. verifiziert |
| 4 | 🟢 A | **Goldbach-Vermutung**: Jede gerade Zahl >2 = p+q | 1742 | Verifiziert bis 4×10¹⁸; Vinogradov: ternäre Form ✓ |
| 5 | 🟢 A | **Zwillingsprimzahl-Vermutung**: ∞ viele (p,p+2) | ~1849 | Polymath8b: Lücken <246 |
| 6 | 🟢 A | **Legendres Vermutung**: Primzahl zwischen n² und (n+1)² | 1798 | Folgt aus RH; Bertrand schwächer |
| 7 | 🟢 A | **Brocard-Vermutung**: ≥4 Primzahlen zwischen pₙ² und pₙ₊₁² | 1904 | Numerisch sehr gut gestützt |
| 8 | 🟢 A | **Artin-Vermutung** (primitive Wurzeln) | 1927 | Hooley 1967 unter GRH; GMH84 fast alle |
| 9 | 🟡 B | **Cramér-Vermutung**: pₙ₊₁−pₙ = O((log pₙ)²) | 1936 | Granville: Faktor 2e^{-γ} nötig |
| 10 | 🟡 B | **Andrica-Vermutung**: √pₙ₊₁ − √pₙ < 1 | 1985 | Verifiziert bis 10¹⁶; impliziert Legendre |
| 11 | 🟡 B | **Oppermann-Vermutung**: Primzahl zwischen n(n-1) und n² | 1882 | Stärker als Legendre |
| 12 | 🟡 B | **Polignac-Vermutung**: ∞ viele pₙ₊₁−pₙ = 2k für jedes k | 1849 | Spezialfall k=1: Zwillingsprimzahl |
| 13 | 🟢 A | **Bunyakovsky-Vermutung**: Irreduz. f∈ℤ[x] → ∞ viele f(n)∈ℙ | 1857 | Spezialfall Dirichlet-Satz (bewiesen) |
| 14 | 🟡 B | **Bateman-Horn-Vermutung**: Quantitative Version von Bunyakovsky | 1962 | Offene Verallgemeinerung |
| 15 | 🟡 B | **Schinzel's Hypothesis H**: Endlich viele Polynome → ∞ viele simultane Primwerte | 1958 | Enthält Zwillingsprimzahl als Spezialfall |
| 16 | 🟢 A | **Mersenne-Primzahlen**: ∞ viele Mersenne-Primzahlen 2ᵖ−1 | ~350 v.Chr. | 51 bekannt; GIMPS-Projekt |
| 17 | 🟢 A | **Fermat-Primzahlen**: Endlich viele Fermat-Primzahlen 2^{2ⁿ}+1? | 1640 | Nur F₀–F₄ bekannt prim |
| 18 | 🟡 B | **Sophie-Germain-Primzahlen**: ∞ viele p mit 2p+1 prim | ~1823 | Keine asymptotische Formel bekannt |
| 19 | 🟡 B | **Cunningham-Ketten**: ∞ viele Ketten der Länge k | 1881 | Für k=2: Sophie Germain |
| 20 | 🟡 B | **Primzahlen der Form n²+1**: ∞ viele? | 1912 | Hardy-Littlewood erwartet ~x/log(x) |
| 21 | 🟡 B | **Primzahlen der Form n!+1**: ∞ viele? | – | Nur endlich viele bekannt |
| 22 | 🟡 B | **Knuth-Vermutung**: ∞ viele Primzahlen p≡3(mod 4) mit (p-3)/4 prim | – | Spezialfall von Schinzel |
| 23 | 🟡 B | **Ingham-Vermutung**: θ < 1/2 in ζ-Nullstellenschranke | 1940 | Huxley: θ≤7/22 (schwach) |
| 24 | 🟡 B | **Hecke-L-Funktionen GRH**: RH für alle Hecke-L-Funktionen | 1936 | Verallgemeinerte RH (GRH) |
| 25 | 🟡 B | **Dirichlet-Reihen**: Alle L(s,χ) haben Re-Nullstellen auf 1/2 | 1837 | Teil von GRH |
| 26 | 🟢 A | **Zhang-Maynard-Theorem Schärfung**: Primzahllücken < 6 (unter GRH) | 2014 | Aktuell: <246 bedingt <16 |
| 27 | 🟡 B | **Siegel-Nullstellen**: Keine Siegel-Nullstellen für Dirichlet-L-Fkt | 1936 | Offen; wichtige Konsequenz |
| 28 | 🟡 B | **Elliott-Halberstam-Vermutung**: Gleichverteilung Primz. in Restklassen | 1968 | θ=1/2 impliziert Primzahllücken <12 |
| 29 | 🟢 A | **Bruns Konstante Genauigkeit**: B₂ ≈ 1.9021605831... (mehr Stellen) | 1919 | Numerisch ~12 Stellen bekannt |
| 30 | 🟡 B | **Heilbronn-Problem**: Anzahl der Dreiecke in Konvexlage mit kleiner Fläche | 1950 | Teilweise gelöst |

---

### I.2 Diophantische Gleichungen & Modulformen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 31 | ✓ | **Fermat-Vermutung** xⁿ+yⁿ=zⁿ, n>2 hat keine nat. Lösungen | 1637 | Wiles 1995 |
| 32 | ✓ | **Catalan-Vermutung**: Einzige aufeinanderfolgend. Primpotenzen: 8,9 | 1844 | Mihailescu 2002 |
| 33 | 🟡 B | **Beal-Vermutung**: Aˣ+Bʸ=Cᶻ → gcd(A,B,C)>1 (x,y,z≥3) | 1997 | Preisgeld $1M; viele Spezialfälle |
| 34 | 🟡 B | **abc-Vermutung**: Für ε>0: c < rad(abc)^{1+ε} fast immer | 1985 | Mochizuki IUT 2012: kontrovers |
| 35 | 🟢 A | **Erdős-Straus**: 4/n = 1/x+1/y+1/z für alle n≥2 | 1948 | Bis 10¹⁴ verifiziert |
| 36 | 🟢 A | **Erdős-Selfridge**: Binomialkoeffizient C(n,k) nie Primzahlpotenz für n≥k+2 | 1975 | Fast vollständig bewiesen |
| 37 | 🟡 B | **Landau-Ramanujan**: τ(n) (Ramanujan) – Lücken zwischen Werten | 1915 | Teilweise gelöst |
| 38 | 🟡 B | **Waring-Problem G(k)**: Genaue Schranken für alle k | 1770 | G(2)=4 (Lagrange); G(3)=7?; G(4)=16 |
| 39 | 🟡 B | **Lander-Parkin-Selfridge-Vermutung** (Euler-Summen-Vermutung) | 1966 | Gegenbeispiele für n=4,5 bekannt |
| 40 | 🟡 B | **Prouhet-Thue-Morse-Vermutung**: Diophantische Eigenschaften | 1906 | Zusammenhang mit Kombinatorik |
| 41 | 🟡 B | **Pillai-Vermutung**: |aˣ−bʸ| → ∞ für feste a,b | 1936 | Folgt aus abc-Vermutung |
| 42 | 🟡 B | **Hall-Vermutung**: |x³−y²| ≥ c·x^{1/2−ε} | 1971 | Folgt aus abc |
| 43 | 🟡 B | **Granville-Verallgemeinerung von abc** | 2007 | Aktive Forschung |
| 44 | 🟡 B | **Gauss-Klassen-Zahl-Vermutung**: h(-d)→∞ für d→∞ | 1801 | Bewiesen (Baker,Stark); Effe. Version offen |
| 45 | 🟡 B | **Gauss-Klassen-Zahl-1**: Nur endlich viele d mit h(-d)=1 | 1801 | Gauss: 9 solche d; vollständig gelöst |
| 46 | 🟡 B | **Birch-Swinnerton-Dyer**: Rang = Nullstell.-Ord. bei s=1 | 1965 | Millennium-Problem (1M$) |
| 47 | 🟡 B | **Zilber-Pink-Vermutung** (Spezielle Untervarietäten) | 2002 | Verallgemeinert Mordell-Lang |
| 48 | 🟡 B | **Manin-Vermutung**: Rational Points auf Fano-Varietäten | 1989 | Für viele Familien bewiesen |
| 49 | 🟡 B | **Brumer-Kramer-Vermutung** (elliptische Kurven) | 1992 | Aktive Forschung |

---

### I.3 Zahlentheoretische Funktionen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 50 | 🟡 B | **Lindelöf-Hypothese**: ζ(1/2+it) = O(t^ε) | 1908 | Folgt aus RH; schwächer |
| 51 | 🟡 B | **Mertens-Vermutung**: |M(x)| < √x immer (M Mertens-Fkt.) | 1897 | Widerlegt 1985 (Odlyzko-te Riele); neue Formen? |
| 52 | 🟡 B | **Pólya-Vermutung**: L(x) ≤ 0 für alle x≥2 (Liouville-Fkt.) | 1919 | Widerlegt 1960; modif. Versionen offen |
| 53 | 🟡 B | **Lehmer-Vermutung**: τ(n) ≠ 0 für alle n (Ramanujan τ-Fkt.) | 1947 | Verifiziert bis 10²²; möglicherweise beweisbar |
| 54 | 🟢 A | **Lehmer-Mahler-Vermutung**: Mahler-Maß M(P) ≥ Lehmer-Kst. ≈1.176 | 1933 | Offen; aktive Forschung |
| 55 | 🟡 B | **Selberg-Klasse-Klassifikation**: Alle primitiven L-Funktionen char.? | 1992 | Hypothese erfüllt für alle bekannten L-Fkt. |
| 56 | 🟡 B | **Langlands-Funktorialität**: Allgemeine Transfers zwischen automorphen Darst. | 1967 | Teilweise bewiesen (Wiles, Lafforgue) |
| 57 | 🟡 B | **Sato-Tate allgemein**: Gleichverteilung aller ell. Kurven | 1963 | Für E ohne CM Taylor et al. 2011 |
| 58 | 🟡 B | **Taniyama-Shimura-Weil allgemein**: Modulkurven ↔ elliptische Kurven | 1957 | Für semistabile Kurven: Wiles 1995; allgemein 2001 ✓ |
| 59 | 🟡 B | **Serre-Vermutung verallgemeinert** | 2005 | Offen für höhere Gewichte mod p |
| 60 | 🟡 B | **Ramanujan-Vermutung für GL(n)**: |ap(f)| ≤ 2p^{(k-1)/2} | 1916/verallg. | Für GL(2): Deligne 1974 ✓; GL(n)>2 offen |
| 61 | 🟢 A | **Vollkommene Zahlen sind gerade**: Keine ungeraden vollkommenen Zahlen | ~350 v.Chr. | Äquivalent zu: Mersenne-Primzahlen ↔ gerade voll. Zahlen |
| 62 | 🟢 A | **Unendlich viele vollkommene Zahlen**: Gibt es unendlich viele? | ~350 v.Chr. | Äquivalent zu: ∞ viele Mersenne-Primzahlen |
| 63 | 🟡 B | **Multiply Perfect Numbers**: ∞ viele σ(n)/n = k für jedes k≥3? | – | Für k=2: vollkommene Zahlen |
| 64 | 🟡 B | **Amicable Numbers**: Unendlich viele amikable Paare? | ~850 | Heuristisch ∞; kein Beweis |
| 65 | 🟢 A | **Riemann-Siegel-Formel Reste**: Schärfere Fehlerterme | 1932 | Aktive Forschung; numerisch wichtig |
| 66 | 🟡 B | **de Bruijn-Newman-Konstante Λ=0**: (äquivalent zu RH) | 1950 | Verifiziert Λ≥0 (Rodgers-Tao 2018) |
| 67 | 🟡 B | **Turán-Vermutung über Nullstellen von ζ** | 1948 | Zusammenhang mit RH |

---

### I.4 Kombinatorische Zahlentheorie & Additive Kombinatorik

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 68 | 🟢 A | **Waring-Goldbach-Vermutung**: Jede nat. Zahl = Summe von k-ten Potenzen von Primzahlen | 1938 | Für k=1: Goldbach; k=2 partiell |
| 69 | 🟡 B | **Erdős-Turán-Additivitätsvermutung**: Jede dichte Folge enthält arith. Progressionen | 1936 | Szemerédi 1975 ✓; schärfere Grenzen offen |
| 70 | 🟡 B | **Erdős-Szemerédi-Vermutung** (Summen-Produkt): |A+A| + |A·A| ≥ |A|^{2-ε} | 1983 | Bester Exp. 4/3+ε |
| 71 | 🟡 B | **Freiman-Struktursatz**: Optimale Konstante in Sanders-Schoen-Theorem | 1964 | Quant. Verbesserungen aktiv |
| 72 | 🟡 B | **Erdős-Vermutung über arithmetische Progressionen**: Divergente Reihe → AP | 1936 | Green-Tao 2004 für Primzahlen ✓ |
| 73 | 🟡 B | **Hales-Jewett-Satz quantitativ**: Schärfere Schranken für Hales-Jewett | 1963 | Polymath1: neue Schranken 2009 |
| 74 | 🟡 B | **Furstenberg-Sárközy-Vermutung**: Mengen ohne Quadratzahldifferenzen | 1978 | Beste Schranke O(1/log log n) |
| 75 | 🟡 B | **Lev-Vermutung**: Große Mengen ohne 3-term AP in ℤ_p | 2001 | Croot-Lev-Pach Exp.-Verbesserung 2016 |
| 76 | 🟡 B | **Collatz-Vermutung** 3n+1 Problem | 1937 | Verifiziert bis 2,95×10²⁰; Conway: verallg. unentscheidbar |
| 77 | 🟡 B | **Kurepa-Vermutung**: !p (Linksfakultät) ≢ 0 mod p für Primzahl p | 1950 | Für p<10⁶ verifiziert |
| 78 | 🟡 B | **Giuga-Vermutung**: n prim ↔ Σᵢ^{n-1} iⁿ⁻¹ ≡ −1 (mod n) | 1950 | Verifiziert für n<10^{13800} |
| 79 | 🟡 B | **Agoh-Giuga-Vermutung**: Vereinfachung der Giuga-Bedingung | 1995 | Äquivalent zu Giuga |
| 80 | 🟡 B | **Erdős-Moser-Vermutung**: 1ᵏ+2ᵏ+...+(m-1)ᵏ = mᵏ nur trivial | 1953 | Moser: m>10^{10^6} für Lösungen |
| 81 | 🟡 B | **Erdős-Vermutung Lücken zwischen Primzahlen**: pₙ₊₁-pₙ > (log n)^k | – | Folgt aus RH für k≤2 |
| 82 | 🟢 A | **Brocard-Ramanujan-Vermutung**: n!+1 = m² nur für n=4,5,7 | 1876 | Verifiziert bis n=10⁹ |
| 83 | 🟡 B | **Pillai-Chowla-Vermutung**: gcd(n!+1,m!+1) = 1 für bestimmte n,m | – | Spezialfall bekannter Probleme |

---

### I.5 Algebraische Zahlentheorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 84 | 🟡 B | **Fontaine-Mazur-Vermutung** (p-adische Galois-Darst.) | 1995 | Kisin: partielle Beweise |
| 85 | 🟡 B | **Leopoldt-Vermutung**: p-adischer Rang von Einheitengruppe | 1962 | Für abelsche Zahlen: Brauer-Siegel |
| 86 | 🟡 B | **Stark-Vermutungen** (L(0,χ) und algebraische Zahlen) | 1971 | Für R=1 bewiesen; allg. offen |
| 87 | 🟡 B | **Lichtenbaum-Vermutung** (Werte von ζ_F bei negativen ganzen Zahlen) | 1973 | Zusammenhang mit K-Theorie |
| 88 | 🟡 B | **Iwasawa-Hauptvermutung**: μ-Invariante für ℤ_p-Erweiterungen | 1973 | Für zyklotomische Felder: Ferrero-Washington ✓ |
| 89 | 🟡 B | **Greenberg-Vermutung**: μ=λ=0 für imaginär-quadratische Felder | 1973 | Für CM-Felder partiell |
| 90 | 🟡 B | **Cohen-Lenstra-Heuristiken**: Verteilung von Klassengruppen | 1984 | Experimentell gut bestätigt |
| 91 | 🟡 B | **Discriminant-Dichtheit**: Dichte von Zahlkörpern nach Diskriminante | 1975 | Bhargava: für Grad 2-5 bewiesen |
| 92 | 🟡 B | **Goldfeld-Vermutung**: Durschn. Rang elliptischer Kurven = 1/2 | 1979 | Bhargava-Shankar: ≤7/6 |

---

## II. Analysis & Komplexe Analysis

### II.1 Reelle Analysis

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 93 | 🟡 B | **Navier-Stokes Glattheit** (3D) | 1845/2000 | Millennium-Problem; 2D gelöst |
| 94 | 🟡 B | **Turbulenzproblem** (Kolmogorov) | 1941 | Mathematisch nicht präzisiert; offen |
| 95 | 🟡 B | **Denjoy-Vermutung**: Trajektorien von Diffeomorphismen | 1932 | Stark bewiesen; Schwächungen offen |
| 96 | 🟡 B | **Bieberbach-Vermutung** (Koeffizientenabschätzung) | 1916 | de Branges 1985 ✓; sharpness offen |
| 97 | 🟡 B | **Brennan-Vermutung**: Integrierbarkeitseigenschaften konformer Abbildungen | 1978 | Teilweise; beste Exp. 4/3−ε |
| 98 | 🟡 B | **Martingal-Vermutung** (Burkholder) | 1966 | Bekannt für p≥2; für p<2 offen |
| 99 | 🟡 B | **Konvexitätsvermutungen** (isoperimetrische Probleme) | – | Verschiedene; einige bewiesen |

### II.2 Harmonische Analysis & Fourier

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 100 | 🟡 B | **Kakeya-Vermutung**: Minimale Fläche für Nadelrotation in ℝⁿ | 1917 | Dim ≥ n/2 bekannt; Dim n offen |
| 101 | 🟡 B | **Fourier-Restriktions-Vermutung**: Strichartz-Schranken | 1970 | Wolff: 3/4-Ergebnis; offen |
| 102 | 🟡 B | **Bochner-Riesz-Vermutung**: Konvergenz von Riesz-Mitteln | 1936 | In 2D bewiesen; höhere Dim. offen |
| 103 | 🟡 B | **Local Smoothing Vermutung** (Seeger-Sogge-Stein) | 1991 | Offen; verwandt mit Kakeya |
| 104 | 🟡 B | **Sogge-Vermutung**: Lᵖ-Normen von Eigenfunktionen auf Mannigf. | 1988 | Für Tori bewiesen |
| 105 | 🟡 B | **Mockbull-Hypothese**: Fast-Perioden von Eigenfunktionen | – | Offen |

### II.3 Komplexe Analysis

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 106 | 🟡 B | **Invariante Unterraumvermutung** (Banach-Räume) | 1935 | Für Hilbert-Räume: Read 1984 Gegenbsp. |
| 107 | 🟡 B | **Yamabe-Problem** auf Mannigfaltigkeiten | 1960 | Trudinger,Aubin,Schoen: gelöst |
| 108 | 🟡 B | **Chern-Conjecture** (affine Mannigfaltigkeiten) | 1955 | Offen in Dim ≥ 4 |
| 109 | 🟡 B | **Branner-Hubbard-Vermutung** (Julia-Mengen) | 1988 | Roesch 2008: für polynomielle Fälle |
| 110 | 🟡 B | **Yoccoz-Puzzle Lokalzusammenhang** | 1993 | MLC-Vermutung (Mandelbrot) verwandt |
| 111 | 🟡 B | **Fatou-Vermutung**: Starre Polynome sind dicht | 1920 | Für hyperbolische: bekannt |

---

## III. Algebra

### III.1 Gruppentheorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 112 | 🟡 B | **Burnside-Problem**: ∀ endl. erzeugte, endlich exponente Gruppe endlich? | 1902 | Allg. falsch (Novikov-Adian 1968); eingeschränkt offen |
| 113 | ✓ | **Klassifikation end. einfacher Gruppen (CFSG)** | 1983 | 18 Familien + 26 sporadische |
| 114 | 🟡 B | **Grigorchuk-Vermutung**: Wachstumsgraph von Grigorchuk-Gruppe | 1980 | Bartholdi: subexponentiell; genaue Rate offen |
| 115 | 🟡 B | **Baum-Connes-Vermutung** (K-Theorie von C*-Algebren) | 1982 | Für hyperbole Gruppen ✓; allg. offen |
| 116 | 🟢 A | **Novikov-Vermutung** (höhere Signaturen) | 1965 | Für viele Gruppen; allg. offen |
| 117 | 🟢 A | **Farrell-Jones-Vermutung** | 1993 | Für CAT(0)-Gruppen; aktive Forschung |
| 118 | 🟡 B | **Kaplansky-Vermutung**: Keine Nullteiler in Gruppenringen | 1945 | Für torsionsfreie Gruppen offen |
| 119 | 🟡 B | **Whitehead-Problem**: Ist jede Whitehead-Gruppe frei? | 1950 | In ZFC unentscheidbar (Shelah 1974) |
| 120 | 🟡 B | **Andrews-Curtis-Vermutung** (Präsentationen trivialer Gruppen) | 1965 | Computationell untersucht |
| 121 | 🟡 B | **Generalized Burnside-Problem**: Gruppen mit endlicher Torsion | 1902 | Teilweise gelöst |
| 122 | 🟡 B | **Schur-Vermutung**: Polynomiale Permutationsgruppen | 1923 | Für endliche Körper bewiesen |
| 123 | 🟡 B | **Thompson-Vermutung**: Endliche einfache Gruppen bestimmt durch Konjugatsklassen | 1987 | Für viele Familien bewiesen |
| 124 | 🟡 B | **Haar-Maß auf pro-endlichen Gruppen**: Quantitative Schranken | – | Technisch anspruchsvoll |

### III.2 Ringtheorie & kommutative Algebra

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 125 | 🟡 B | **Jacobian-Vermutung**: Det(Jac f)=const → f bijektiv | 1939 | In allen Dim. offen |
| 126 | 🟡 B | **Dixmier-Vermutung**: Alle Endomorphismen des Weyl-Algebras sind Automorphismen | 1968 | Äquivalent zur Jacobian-Vermutung |
| 127 | 🟡 B | **Homologische Vermutungen** (Auslander, Bass etc.) | 1956 | Viele gelöst; einige offen |
| 128 | 🟡 B | **Schur-Positivity-Vermutung** (symmetrische Funktionen) | 1970s | Für viele Familien bewiesen |
| 129 | 🟡 B | **Stanley-Wilf-Limit**: Wachstum musterfreier Permutationen | 1999 | Marcus-Tardos 2004 ✓ |
| 130 | 🟡 B | **Karoubi-Vermutung** (K-Theorie) | 1970 | Für reguläre Ringe bewiesen |

### III.3 Darstellungstheorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 131 | 🟡 B | **Lusztig-Vermutungen** (Charakterformeln für p-Charakteristik) | 1979 | Für große p bewiesen (Andersen-Jantzen-Soergel) |
| 132 | 🟡 B | **James-Vermutung** (mod-p Darst. symm. Gruppen) | 1990 | Umfangreich untersucht |
| 133 | 🟡 B | **McKay-Vermutung**: Blockweise Charakterkonjugation | 1972 | Für viele Gruppen; allg. offen |
| 134 | 🟡 B | **Alperin-McKay-Vermutung**: Schärfung von McKay | 1976 | Für semiauflösbare Gruppen |
| 135 | 🟡 B | **Broué-Vermutung** (Homotopieäquivalenz von Blöcken) | 1990 | Für zyklische Defektgruppen |

---

## IV. Geometrie & Topologie

### IV.1 Differentialgeometrie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 136 | 🟡 B | **Geometrisierungs-Vermutung Dim>3** | 1982 | Dim 3: Perelman 2003 ✓; Dim≥4 anders |
| 137 | 🟡 B | **Smooth Schoenflies-Problem Dim 4** | 1960 | Dim ≤3 und ≥6 bekannt; 4 und 5 offen |
| 138 | 🔴 C | **Smooth Poincaré Dim 4**: Ist S⁴ die einzige h-coborde Sphäre? | 1960 | Einzig offene Dimension |
| 139 | 🟡 B | **Calabi-Yau-Vermutung**: Existenz Ricci-flacher Kähler-Metriken | 1954 | Yau 1978 ✓; Eindeutigkeit/Moduli offen |
| 140 | 🟡 B | **Hopf-Vermutung**: S²×S² hat keine pos. Schnittkrümmung | 1931 | Offen; Hopf-Morse-Ansatz gescheitert |
| 141 | 🟡 B | **Thurston-Geometrisierung Dim 4**: Klassifikation 4-Mannigfaltigkeiten | 1982 | Völlig offen |
| 142 | 🟡 B | **Donaldson-Vermutung**: Existenz von Seiberg-Witten-Invarianten | 1994 | Aktive Forschung |
| 143 | 🟡 B | **Weinstein-Vermutung**: Existenz von Reeb-Orbits | 1979 | Für 3-Mannigfaltigkeiten: Taubes 2007 ✓ |
| 144 | 🟡 B | **Arnoux-Rauzy-Vermutung** (Substitutionssysteme) | 1991 | Partiell gelöst |
| 145 | 🟡 B | **Schoen-Vermutung**: Minimalflächen in S³ | 1983 | Für stabile: bekannt |
| 146 | 🟡 B | **Chern-Gauss-Bonnet Verallgemeinerung** | 1944 | In bestimmten Dim. bekannt |
| 147 | 🟡 B | **Yau-Vermutungen** (harmonische Abbildungen) | 1975 | Verschiedene, teils offen |

### IV.2 Algebraische Topologie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 148 | ✓ | **Poincaré-Vermutung (3D)** | 1904 | Perelman 2003 |
| 149 | 🟡 B | **Smooth 4D Poincaré** | 1960 | Einzig offen; exotische R⁴ bekannt |
| 150 | 🟡 B | **Harer-Zagier-Vermutung** (Euler-Charakteristik von Modulräumen) | 1986 | Bewiesen; Verallgemeinerungen offen |
| 151 | 🟢 A | **Baum-Connes-Vermutung** | 1982 | Für hyperbolische Gruppen ✓ |
| 152 | 🟡 B | **Freyd-Vermutung**: Unauflösliche Homotopiegruppen | 1965 | Offen |
| 153 | 🟡 B | **Telescope-Vermutung** (Miller) | 1984 | Stark negiert durch Burklund et al. 2022 |
| 154 | 🟡 B | **Chromatic Nullstellensatz**: Algebraische K-Theorie | 2020 | Sehr neue Forschung |
| 155 | 🟡 B | **Sullivan-Vermutung**: Fixpunkte von p-Gruppen auf Komplexen | 1978 | Miller 1984 ✓; Verallg. offen |
| 156 | 🟡 B | **Zeeman-Vermutung** (Kontraktierbarkeit) | 1960 | In Dim ≥5 ✓; Dim 4 offen |
| 157 | 🟡 B | **Whitehead-Asphärizitäts-Vermutung** | 1941 | Für endl.-dim. Komplexe mit π₁ frei: bewiesen |

### IV.3 Knotentheorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 158 | 🟡 B | **Poincaré-Dualität für Knoten-Homologie** | 1999 | Khovanov-Homologie: aktiv |
| 159 | 🟡 B | **HOMFLY-Schärfungs-Vermutung** | 1985 | Knotenpolynome und Kategorisierung |
| 160 | 🟡 B | **Vassiliev-Vermutungen** (Knotenendlichkeit) | 1990 | Für kleine Invarianten bestätigt |
| 161 | 🟡 B | **Gordon-Luecke-Folgerungen** | 1989 | Gordon-Luecke ✓; Verallg. offen |
| 162 | 🟡 B | **Fox-Milnor-Vermutung**: Knotenkonkordanz | 1966 | Für spezielle Fälle; allg. offen |
| 163 | 🟡 B | **Slice-Ribbon-Vermutung**: Jeder Slice-Knoten ist ribbon? | 1978 | Offen; fundamentales Problem |
| 164 | 🟡 B | **4-Ball-Crossing-Number**: Gleichung mit unknot | 2000s | Aktive Forschung |
| 165 | 🟡 B | **Conway-Vermutung** (Knotenpolynome invariant unter Mutation) | 1970 | Morton-Traczyk: für HOMFLY offen |

### IV.4 Diskrete Geometrie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 166 | ✓ | **Kepler-Vermutung**: Dichteste Kugelpackung ℝ³ | 1611 | Hales 2005, formalisiert 2017 |
| 167 | 🟡 B | **Dichteste Kugelpackung ℝⁿ für n≥5**: Optimalität | – | Dim 8 (Viazovska 2016 ✓); Dim 24 (✓); Rest offen |
| 168 | 🟡 B | **Covering-Codes-Vermutung**: Minimale Überdeckungszahlen | 1958 | Für kleine Parameter bekannt |
| 169 | 🟡 B | **Hadwiger-Nelson-Problem**: Chromatische Zahl des eukl. Graphen ℝ² | 1950 | de Grey 2018: ≥5; oben ≤7 |
| 170 | 🟡 B | **Erdős-Einheitsdistanz-Problem**: Maximale Paare gleicher Abstände | 1946 | Guth-Katz 2011: O(n^{4/3}); optimal? |
| 171 | 🟡 B | **Erdős-Kolineare-Punkte-Problem** | 1970 | Verschiedene Teilvermutungen |
| 172 | 🟡 B | **Heilbronn-Dreieck**: Minimale Dreiecksfläche bei n Punkten in Einheitsquadrat | 1950 | Komlós-Pintz-Szemerédi: O(log n/n) |
| 173 | 🟡 B | **Penrose-Vermutung**: Aperiodische Pflasterungen | 1974 | Penrose ✓; mit 1 Kachel: ✓ 2023 |
| 174 | 🟡 B | **Steinitz-Vermutung**: Polytop-Rekonstruktion aus Gitter | – | Für konvexe Polytope partiell |
| 175 | 🟡 B | **Borsuk-Vermutung**: n+1 Mengen in ℝⁿ für jede beschr. Menge | 1933 | Kahn-Kalai 1993: falsch für große n |

---

## V. Graphentheorie & Kombinatorik

### V.1 Graphenfärbung

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 176 | ✓ | **4-Farben-Satz** | 1852 | Appel-Haken 1976; Robertson et al. 1997 |
| 177 | 🟡 B | **Hadwiger-Vermutung**: Kₙ-Minor → n-färbbar | 1943 | Für n≤6 bewiesen; n=7 offen |
| 178 | 🟡 B | **Erdős-Hajnal-Vermutung**: Monochromat. Clique/Independent in chi-beschr. Graphen | 1989 | Für perfekte Graphen: ja |
| 179 | 🟢 A | **Gyárfás-Vermutung**: χ(G) beschränkt wenn keine langen Induktionswege | 1975 | Für Chordal-Graphen ✓ |
| 180 | 🟡 B | **Alon-Tarsi-Vermutung**: Chromatische vs. Listenzahl | 1992 | Für planare Graphen: Thomassen |
| 181 | 🟡 B | **List-Coloring-Vermutung**: ch(G) = χ(G) für Kantenfärbungen | 1992 | Für bipartite Graphen: Galvin ✓ |
| 182 | 🟡 B | **Ringel-Kotzig-Vermutung**: Graceful Trees | 1966 | Für Pfade, Caterpillars ✓; allg. offen |
| 183 | 🟡 B | **Vizing-Vermutung**: χ'(G□H) ≤ χ'(G)+χ'(H)? | 1963 | Spezielles Produkt-Problem |
| 184 | 🟡 B | **Jaeger-Linial-Mohar-Tarsi-Vermutung**: Zirkularer Fluss | 1992 | Für bestimmte Klassen |
| 185 | 🟡 B | **Hedetniemi-Vermutung**: χ(G×H) = min(χ(G),χ(H)) | 1966 | Widerlegt 2019 (Shitov); Varianten offen |

### V.2 Graphstruktur

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 186 | ✓ | **Robertson-Seymour-Satz** (Graph Minor Theorem) | 2004 | Robertson-Seymour |
| 187 | 🟡 B | **Seymour-Vermutung**: Zweiter Nachbar | 1990 | Für Turniere mit rationalen Gewichten ✓ |
| 188 | 🟡 B | **Jaeger-Vermutung**: 4-Fluss-Vermutung | 1988 | Für planar: Appel-Haken |
| 189 | 🟡 B | **Tutte-5-Fluss-Vermutung**: Jeder Brücken-freie Graph hat 5-Flusszahl ≤5 | 1954 | Tutte: f(G)≤5 für brückenfreie Graphen |
| 190 | 🟡 B | **Cycle Double Cover**: Jeder brückenfreie Graph hat CDC | 1973 | Seymour, Szekeres; partiell |
| 191 | 🟡 B | **Barnette-Vermutung**: Jeder kubische 3-fach zusammenh. planare Graph hamilton. | 1969 | Offen |
| 192 | 🟡 B | **Hamiltonizität kubischer 3-fach-ZH Graphen** (Tait) | 1884 | Tait-Vermutung falsch 1946; Barnette-Variante offen |
| 193 | 🟡 B | **Thomassen-Vermutung**: Jede k-fach zusammenhängende Graphtriangulierung hamiltonisch | 1982 | Für k≥4 bewiesen |

### V.3 Ramsey-Theorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 194 | 🟡 B | **Ramsey R(5,5) genau bestimmen** | 1930 | Bekannt: 43≤R(5,5)≤48 |
| 195 | 🟡 B | **Ramsey R(k,k) exponentielle Schranke**: R(k,k) ≥ 2^{k/2+ε}? | 1935 | Sah 2023: verbesserte obere Schranke |
| 196 | 🟡 B | **Erdős-Heilbronn-Vermutung** (jetzt ✓): Summen in Mengen mod p | 1964 | Dias da Silva-Hamidoune 1994 ✓ |
| 197 | 🟡 B | **Schur-Zahlen**: Exakte Werte für alle k | 1916 | S(5)=160; S(6) offen |
| 198 | 🟡 B | **Van-der-Waerden-Schranken**: W(k;r) Wachstum | 1927 | Gowers: ACKERMANN; besser möglich |

### V.4 Andere kombinatorische Vermutungen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 199 | 🟡 B | **Lonely Runner-Vermutung**: n Läufer auf Einheitskreis | 1967 | Für n≤7 bewiesen |
| 200 | 🟡 B | **Graceful-Tree-Vermutung** | 1967 | Für Caterpillars, Wheels ✓; allg. offen |
| 201 | 🟡 B | **Turán-Vermutung** (Torus): Maximale Kantenzahl | 1948 | Für Torus: exakt bekannt; höhere Flächen offen |
| 202 | 🟡 B | **Stanley-Vermutung** (Simpliziale Komplexe): h-Vektor | 1994 | Für Sphären: ✓; allg. offen |
| 203 | 🟡 B | **Erdős-Ko-Rado-Verallgemeinerungen** | 1961 | Für t-intersecting Systeme |
| 204 | 🟡 B | **Frankl-Vermutung**: |𝒜∪ℬ| ≥ max(|𝒜|,|ℬ|) für Unionen-abgeschlossene Familien | 1979 | Bošnjak-Marković: für Familien ≤50 Mengen |
| 205 | 🟡 B | **Sauer-Shelah-Schärfung** (VC-Dimension) | 1972 | Sauer ✓; Verbesserungen gesucht |
| 206 | 🟡 B | **Bollobás-Set-Pairs-Vermutung** | 1965 | Für lineare Spannweiten bewiesen |
| 207 | 🟡 B | **Berge-Foulkner-Vermutung** (Hypergraph-Matching) | 1970 | Verschiedene Versionen |
| 208 | 🟡 B | **Erdős-Faber-Lovász-Vermutung**: χ(H) ≤ n für n-uniform n-partite H | 1972 | Kim et al. 2021: für große n ✓ |
| 209 | 🟡 B | **Alspach-Vermutung** (Graphzerlegungen) | 1981 | Für vollständige Graphen: Schönheim |

---

## VI. Mathematische Logik & Grundlagen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 210 | ⬛ D | **Kontinuumshypothese** (CH) | 1878 | In ZFC unentscheidbar (Cohen 1963) |
| 211 | ⬛ D | **Martin's Axiom** (MA) | 1970 | Konsistent mit ZFC |
| 212 | ⬛ D | **Whitehead-Problem** (Freie Gruppen) | 1950 | Shelah 1974: unentscheidbar in ZFC |
| 213 | 🔴 C | **P ≠ NP** | 1971 | Millennium-Problem; mögl. ZFC-unentscheidbar |
| 214 | 🔴 C | **BPP=P?** (Derandomisierung) | 1977 | Offen; Nisan-Wigderson: pseudo-random |
| 215 | 🔴 C | **NP=co-NP?** | 1979 | Folgt aus P=NP; eigenständig offen |
| 216 | 🔴 C | **PSPACE vs. NP** | 1970 | PSPACE ⊇ NP bekannt; Gleichheit? |
| 217 | 🔴 C | **L vs. P** | 1974 | Landau-Problem; wichtig |
| 218 | 🔴 C | **NL=L?** | 1979 | Immerman-Szelepcsényi: NL=co-NL ✓ |
| 219 | 🟡 B | **Circuit Lower Bounds**: Superpolyn. Schranken für Exp-Zeit-Fkt | 1985 | Ryan Williams: kleine Fortschritte |
| 220 | 🟡 B | **Natural Proofs Barrier** umgehen | 1994 | Razborov-Rudich: fundamental |
| 221 | 🔴 C | **Schanuel-Vermutung** (Transzendenz) | 1960s | Impliziert Lin.-Unabh. von log, exp |
| 222 | ⬛ D | **Grote Kardinalzahl-Axiome und V=L**: Konsistenz | – | Unabhängig von ZFC |
| 223 | ⬛ D | **Axiom of Determinacy (AD)**: Konsistenz mit ZF (ohne AC) | 1962 | Woodin-Karden-Vorschläge |
| 224 | 🔴 C | **Permanent vs. Determinante**: poly. Äquivalenz? | 1979 | Valiant: #P-vollständig |

---

## VII. Algebraische Geometrie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 225 | ✓ | **Weil-Vermutungen** | 1949 | Deligne 1974 |
| 226 | 🟡 B | **Hodge-Vermutung**: Algebraische Zyklen ↔ Hodge-Klassen | 1950 | Millennium-Problem (1M$) |
| 227 | 🟡 B | **Grothendieck-Standardvermutungen** | 1969 | Nur Lefschetz-Vermutung bekannt |
| 228 | 🟡 B | **Bloch-Kato-Vermutung** (Milnor K-Theorie) | 1986 | Voevodsky 2011 ✓; Verallg. offen |
| 229 | 🟡 B | **Bloch-Beilinson-Vermutungen** (Chow-Motive) | 1980 | Sehr grundlegend; wenig bekannt |
| 230 | 🟡 B | **Tate-Vermutung**: Galois-Darstellungen und algebraische Zyklen | 1965 | Für Kurven: ✓; höher offen |
| 231 | 🟡 B | **Mumford-Vermutung** (Klassifikation Modulräume) | 1969 | Harer ✓; höhere Analoga offen |
| 232 | 🟡 B | **Green-Lazarsfeld-Vermutung** (Syzygien von Kurven) | 1986 | Voisin 2002–2005: für generische Kurven ✓ |
| 233 | 🟡 B | **Beauville-Voisin-Vermutung** (Chow-Ring von K3-Flächen) | 2004 | Aktive Forschung |
| 234 | 🟡 B | **Maulik-Nekrasov-Okounkov-Pandharipande** (MNOP) | 2003 | Verbindung GW/DT-Invarianten; partiell ✓ |

---

## VIII. Dynamische Systeme & Chaos

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 235 | 🟡 B | **Mandelbrot-MLC-Vermutung**: Mandelbrot-Menge lokal zusammenhängend | 1982 | Für Teile bewiesen; allg. offen |
| 236 | 🟡 B | **Fatou-Vermutung**: Hyperbole Polynome dicht | 1920 | Für quadratische Familie: offen |
| 237 | 🟡 B | **Smale-14. Problem**: Gibt es effiziente Algorithmen für Nullstellen? | 2000 | Beltran-Shub: probabilistisch |
| 238 | 🟡 B | **Hénon-Attractor-Vermutung**: Seltsame Attraktoren in Hénon-Familie | 1976 | Benedicks-Carleson: pos. Maß ✓ |
| 239 | 🟡 B | **Newhouse-Phenomenon**: Persistenz von Sink-Inseln | 1974 | Phänomen bekannt; quantitativ offen |
| 240 | 🟡 B | **KAM-Tori-Vermutung**: Persistenz unter kleinen Störungen | 1954 | KAM ✓; Maß von Tori offen |
| 241 | 🟡 B | **Ergodizität von Billard-Systemen** | 1939 | Für konvexe Billards: Wojtkowski |
| 242 | 🟡 B | **Ratner-Theorie Verallgemeinerungen** | 1991 | Ratner ✓; p-adische Analoga |

---

## IX. Wahrscheinlichkeitstheorie & Stochastik

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 243 | 🟢 A | **Schramm-Loewner-Evolution**: Skalengrenzwerte planarer Modelle | 2000 | SLE ✓; viele Modelle bestätigt |
| 244 | 🟡 B | **Universalität der Gauss'schen Matrizen**: Eigenwerteverteilung allg. Matrizen | 2011 | Tao-Vu: für viele Ensembles ✓ |
| 245 | 🟡 B | **Rokhlin-Vermutung** (Entropie): Bestimmung der Entropie durch Typ | 1967 | Für σ-endliche Transformationen: offen |
| 246 | 🟡 B | **Conze-Lesigne-Vermutung** (Host-Kra-Faktoren) | 1988 | Host-Kra 2005 ✓; allg. Nilsysteme |
| 247 | 🟡 B | **Poisson-Universalität**: Keine anderen Grenzverteilungen für Zufallsmatrizen | 1960 | Für GUE: ✓; allg. offen |
| 248 | 🟡 B | **Gaussian Free Field Vermutungen** | 2000s | Verschiedene Konvergenzresultate |

---

## X. Mathematische Physik

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 249 | 🟡 B | **Yang-Mills Massenücke** | 1954/2000 | Millennium-Problem; erfordert neue QFT-Math. |
| 250 | 🟡 B | **Chern-Simons-Invarianten**: Topologische QFT-Vorhersagen | 1974 | Witten 1989: Jones-Polynome; mathematisch partiell |
| 251 | 🟡 B | **Mirror Symmetry**: Spiegelpaare von Calabi-Yau-Mannigf. | 1991 | Kontsevich: Homological MS ✓; allg. offen |
| 252 | 🟡 B | **Monstrous Moonshine**: Vollständige Klassifikation | 1979 | Borcherds 1992 ✓; Verallg. offen |
| 253 | 🟡 B | **S-Dualität** in math. rigider Form | 1994 | Kapustin-Witten: partiell |
| 254 | 🟡 B | **Langlands-Dualität und Physik**: Elektrische/magnetische Dualität | 1979 | Geometrisches Langlands: Arinkin-Gaitsgory |
| 255 | 🟡 B | **Gauge-Theorie und 4-Mannigfaltigkeiten**: Vollständige Klassifikation | 1983 | Donaldson, Seiberg-Witten: partiell |

---

## XI. Spezielle Themen

### XI.1 Transzendenztheorie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 256 | 🟡 B | **Schanuel-Vermutung** | 1966 | Grundlegendste Vermutung in Transz.-Theorie |
| 257 | 🟡 B | **Vier-Exponenten-Vermutung**: e^{αβ} transz. für alg. α,β,1/β,1/α | 1960 | Folgt aus Schanuel |
| 258 | 🟡 B | **Konvergenz-Vermutung**: Ist e+π irrational? | – | Unbekannt! Schanuel würde imply |
| 259 | 🟡 B | **Euler-Konstante γ**: Irrationalität von γ ≈ 0.5772... | 1735 | Nicht mal Irrationalität bewiesen |
| 260 | 🟡 B | **ζ(2k+1) irr.**: Sind ζ(5), ζ(7), ζ(9), ... irrational? | 1979 | Apéry: ζ(3) irr. ✓; Rest offen |
| 261 | 🟡 B | **Catalan-Konstante G**: Irrationalität von G ≈ 0.9159... | – | Unbekannt |
| 262 | 🟡 B | **Transz. von π+e**: Ist π+e transcendent? | – | π transzend. ✓; π+e unbekannt |
| 263 | 🟡 B | **Algebraische Unabhängigkeit von π und e** | – | Schanuel impliziert Transzendenz |
| 264 | 🟡 B | **Normality of π**: Ist π eine normale Zahl (jede Ziffer gleichverteilt)? | – | Für fast alle Zahlen gilt es; für π offen |
| 265 | 🟡 B | **Normality of e**: Ist e eine normale Zahl? | – | Für e unbekannt |
| 266 | 🟡 B | **Normality of √2**: Ist √2 normal? | – | Für algebraische Irrationalzahlen unbekannt |

### XI.2 Zahlentheorie der Algorithmen

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 267 | 🟢 A | **Primzahltest AKS-Schranke**: Verbesserte Laufzeit? | 2002 | AKS: O(n^6); verbessert auf O(n^3)? |
| 268 | 🟡 B | **Integer Factorization nicht in P**: Faktorisierung braucht mehr als polynom. Zeit | 1977 | Offen; Quantencomputer: polynomial |
| 269 | 🟡 B | **Elliptic Curve Faktorisierung Optimalität** | 1985 | L-Methode; sub-exp. optimiert? |

### XI.3 Kombinatorische Geometrie

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 270 | 🟡 B | **Erdős-Vermutung über Distanzmengen** | 1946 | Guth-Katz 2011: √n Distanzen ✓ |
| 271 | 🟡 B | **Orchard-Problem**: Maximale Geraden durch k oder mehr Punkte | 1893 | Für k=3: Sylvester-Gallai ✓ |
| 272 | 🟡 B | **Happy Ending Problem**: n-Eck in allg. Lage | 1935 | Erdős-Szekeres: f(n) ≤ C(2n-4,n-2); optimal? |
| 273 | 🟡 B | **Sylvester-Gallai Verallgemeinerungen**: Higher Dimensions | 1893 | Für ℝ²: ✓; ℂ²: Kelly-Theorem |
| 274 | 🟡 B | **Grünbaum-Vermutung**: n konvexe Körper → gemeinsamer Transversal | 1960 | Für Halbräume: Helly-Theorem ✓ |
| 275 | 🟡 B | **Hirsch-Vermutung**: Durchmesser von Polytopen ≤ n−d | 1957 | Widerlegt 2010 (Santos); versch. Varianten offen |

---

## XII. Millennium-Probleme & Preisprobleme

| Problem | Stifter | Preisgeld | Kat. | Status |
|---------|---------|-----------|------|--------|
| **Riemann-Hypothese** | Clay | $1M | 🟡 B | Offen |
| **P vs. NP** | Clay | $1M | 🔴 C | Offen |
| **Navier-Stokes** | Clay | $1M | 🟡 B | Offen |
| **Hodge-Vermutung** | Clay | $1M | 🟡 B | Offen |
| **Birch-Swinnerton-Dyer** | Clay | $1M | 🟡 B | Offen |
| **Yang-Mills Massenücke** | Clay | $1M | 🟡 B | Offen |
| **Poincaré-Vermutung** | Clay | $1M | ✓ | Perelman 2003 ✓ |
| **Beal-Vermutung** | Beal | $1M | 🟡 B | Offen |
| **Collatz-Vermutung** | – | $120K | 🔴 C | Offen |
| **Erdős-Einheitsdistanz** | Erdős | $500 | 🟡 B | Offen |
| **Legendre-Vermutung** | – | – | 🟢 A | Offen |

---

## XIII. Kürzlich aufgestellte Vermutungen (2000–2026)

| # | Kat. | Vermutung | Jahr | Status/Anmerkung |
|---|------|-----------|------|------------------|
| 276 | 🟢 A | **Polymath8b-Schärfung**: Primzahllücke < 6 | 2014 | Unter GRH bedingt möglich |
| 277 | 🟡 B | **Sarnak-Vermutung** (Möbius-Orthogonalität): μ(n) unkorrelliert mit Toeplitz | 2010 | Für spezielle Systeme: ✓ |
| 278 | 🟡 B | **Tao-Ziegler-Vermutung** (Polynomiale AP) | 2008 | Tao-Ziegler: lineare AP ✓; polynom. offen |
| 279 | 🟡 B | **Conlon-Fox-Zhao-Vermutung** (reguläre Graphen) | 2015 | Neue Extremalresultate |
| 280 | 🟡 B | **Zoom Vermutung** (Kac-Moody-Algebren) | 2010 | Aktive Darstellungstheorie |
| 281 | 🟡 B | **Burklund-Hahn-Senger** (Teleskop-Antinomie) | 2022 | Miller-Vermutung widerlegt; neue Fragen |
| 282 | 🟢 A | **Viazovska-Vermutungen** (optimale Packungen in ℝ⁸, ℝ²⁴) | 2016 | ℝ⁸,ℝ²⁴: ✓; andere Dim. offen |
| 283 | 🟡 B | **Scholze-Vermutungen** (Perfectoid Geometrie) | 2012 | Teilweise bewiesen; tiefe Verallg. offen |
| 284 | 🟡 B | **Chebotarev-Dichte-Schärfung** (effektiv) | 2020 | Pierce-Turnage-Butterbaugh-Wood: Verbesserungen |
| 285 | 🟢 A | **Sun-Tzu-Vermutung** (CRT-Variante für Quadrate) | 2019 | Offen; elementar formulierbar |
| 286 | 🟡 B | **Uniformly Distributed Sequences Vermutung** | 2018 | Zusammenhang mit Normality |

---

## Statistik

| Kategorie | Anzahl |
|-----------|--------|
| ✓ Bewiesen (in dieser Liste) | ~25 |
| 🟢 A Hohe Chance | ~30 |
| 🟡 B Mittlere Chance | ~195 |
| 🔴 C Niedrige Chance | ~10 |
| ⬛ D ZFC-unentscheidbar | ~5 |
| **Gesamt offen** | **~240** |
| **Gesamt Einträge** | **~286** |

---

## Empfehlungen für Beweis-Training (nach Schwierigkeit)

### Leicht (Numerische Verifikationen)
- Brocard-Ramanujan-Vermutung (n!+1=m² nur für n=4,5,7) – Bis 10¹⁰ testen
- Goldbach – Weitere Zerlegungen analysieren
- Erdős-Straus – Für Primzahlklassen beweisen
- Giuga-Vermutung – Strukturelle Analyse

### Mittel (Partielle analytische Ergebnisse)
- Artin-Vermutung – Für spezifische Basen a=2,3,... dicht beweisen
- Legendre – Bedingte Beweise unter RH formalisieren
- Lehmer τ(n)≠0 – Kongruenzbedingungen ausschöpfen

### Schwer (Neue Ideen erforderlich)
- Zwillingsprimzahl – GPY-Methode verfeinern
- Riemann-Hypothese – L-Funktionen tiefer verstehen
- abc-Vermutung – Unabhängige Verifikation von Mochizuki

---

*Vollständige Analyse in: [`vermutungen_analyse.md`](./vermutungen_analyse.md)*
*Beweisversuche Kategorie A: [`src/kategorie_a_untersuchungen.py`](./src/kategorie_a_untersuchungen.py)*
