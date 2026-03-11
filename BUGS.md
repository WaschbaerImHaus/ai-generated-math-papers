# BUGS.md - specialist-maths

## Offene Bugs
_(keine bekannten Bugs)_

## Behobene Bugs

### BUG-B4-P18-EN-001 (Build 92): paper18 EN — Vereinfachung Remark 4.3
- **Datei**: papers/batch4/paper18_vinogradov_three_primes.tex
- **Problem**: In Remark 4.3 stand `$1 + \mu(2)^3 c_2(n)/\phi(2)^3 = 1 - c_2(n)$`
  — das `\mu`-Vorzeichen war implizit, was zu Verwirrung führen konnte
- **Behoben**: Auf `$1 - c_2(n)/(2-1)^3 = 1 - c_2(n)$` vereinfacht (explizites Vorzeichen)

### BUG-B4-P18-DE-001 (Build 92): paper18 DE — Gleiche Korrektur Bemerkung 4.3
- **Datei**: papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: Analog zu BUG-B4-P18-EN-001
- **Behoben**: Analog korrigiert

### BUG-B4-P18-DE-002 (Build 91): paper18 DE — Fallunterscheidung Satz 5.1
- **Datei**: papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: Beweis enthielt ursprünglich nur den Fall n≤N
- **Behoben**: Alle drei Fälle (n≤N, N<n≤2N, 2N<n≤3N) ergänzt

### BUG-B4-P19-EN-001 (Build 91): paper19 EN — Querverweise auf Paper 17
- **Datei**: papers/batch4/paper19_waring_problem.tex
- **Problem**: Falsche Referenzen Definition 5.1 und Lemma 5.3
- **Behoben**: Auf Definition 5.2 und Lemma 5.5 korrigiert

### BUG-B4-P19-DE-001 (Build 91): paper19 DE — Fehlendes Korollar Weyl-Schranke G(k)
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: Korollar G(k) ≤ (k-2)2^{k-1}+5 fehlte
- **Behoben**: Korollar mit Beweis in Abschnitt 4 ergänzt

### BUG-B4-P19-DE-002 (Build 91): paper19 DE — Wooley-Korollar ohne Beweis
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: G(k) ≤ k(log k + 4.20032) hatte keinen Beweis
- **Behoben**: Kurzen Beweis ergänzt

### BUG-B4-P19-DE-003 (Build 91): paper19 DE — Fehlender formaler Satz J_{s,k}
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: Nur Bemerkungen statt formalem Satz
- **Behoben**: Formaler Satz für J_{s,k}-Schranke aus EN adaptiert

### BUG-B4-P20-EN-001 (Build 91): paper20 EN — Satzfragment in Remark "Odd n"
- **Datei**: papers/batch4/paper20_goldbach_singular_series.tex
- **Problem**: "odd only for..." Satzfragment
- **Behoben**: Fragment entfernt, Remark vollständig formuliert

### BUG-B4-P20-EN-002 (Build 91): paper20 EN — Falsche Tabellenwerte S(n)
- **Datei**: papers/batch4/paper20_goldbach_singular_series.tex
- **Problem**: Tabellenwerte enthielten extra Faktor (p-1)/(p-2) für p∤n
- **Behoben**: Korrekte Werte n=4: 2C_2; n=6: 4C_2; n=8: 2C_2; n=10: 8C_2/3; n=12: 4C_2; n=30: 16C_2/3

### BUG-B4-P20-DE-001 (Build 91): paper20 DE — Identische Tabellenfehler
- **Datei**: papers/batch4/paper20_goldbach_singulaere_reihe_de.tex
- **Problem**: Analog zu BUG-B4-P20-EN-002
- **Behoben**: Korrekte Werte wie EN eingetragen

### BUG-010 (Build 91): paper2 DE — Redundanz in Abschnitt 3
- **Datei**: papers/reviewed/batch1/paper2_lehmer_drei_primfaktoren_de.tex
- **Problem**: Zwei vollständige Beweise ohne Kennzeichnung
- **Behoben**: `\medskip\noindent\textit{Alternativer Beweis:}` vor zweiten Beweis eingefügt

### BUG-011-Papers (Build 91): paper7 DE — Schwacher Zwischenschritt
- **Datei**: papers/reviewed/batch1/paper7_lehmer_kein_semiprim_de.tex
- **Problem**: `≥ 1·(p-1) - 1 = p-2 ≥ 1` war zu schwach
- **Behoben**: Auf `≥ 1 \cdot p - 1 = p-1 \ge 2 > 0` verschärft

### BUG-011 (Build 62): test_config.py – Falscher Autor-Wert
- **Datei**: tests/test_config.py
- **Problem**: TestMetadata::test_author_value prüfte auf 'Kurt Ingwer' statt 'Michael Fuhrmann'
- **Behoben**: Test korrigiert auf korrekten Autor 'Michael Fuhrmann'
