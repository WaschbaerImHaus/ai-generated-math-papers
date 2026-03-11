# Gutachten: Paper 6 — Jede Lösung von Lehmers Totient-Problem ist quadratfrei

**Papers:** `paper6_lehmer_squarefree.tex` (EN) · `paper6_lehmer_quadratfrei_de.tex` (DE)
**Gutachter:** Claude Code (automatisiertes Peer-Review)
**Datum:** 2026-03-11
**Gesamturteil:** ✅ Annahme empfohlen (beide Versionen)

---

## 1. Zusammenfassung

Beide Papers beweisen, dass jede hypothetische **Lehmer-Zahl** (zusammengesetztes $n$ mit $\varphi(n) \mid (n-1)$) notwendigerweise **quadratfrei** ist. Der Beweis nutzt zwei Divisibilitätsketten: $p^2 \mid n \Rightarrow p \mid \varphi(n) \Rightarrow p \mid (n-1)$, kombiniert mit $p \mid n$, um $p \mid 1$ zu erzwingen — Widerspruch.

---

## 2. Mathematische Korrektheit

### Hauptsatz (Theorem 2.1 / Satz 2.1)
**Korrekt.** Detailprüfung:

**Schritt 1:** $p^2 \mid n \Rightarrow p \mid \varphi(n)$.

Begründung: Für $p^a \| n$ (genau $p^a$ teilt $n$) mit $a \ge 2$ gilt:
$$\varphi(n) = \varphi(p^a) \cdot \varphi(n/p^a) = p^{a-1}(p-1) \cdot \varphi(n/p^a)$$
Für $a \ge 2$: $p^{a-1} \ge p$, also $p \mid \varphi(n)$. Für $a = 2$: $p(p-1) \mid \varphi(n)$. ✓

**Schritt 2:** $p \mid \varphi(n)$ und $\varphi(n) \mid (n-1)$ $\Rightarrow$ $p \mid (n-1)$ (Transitivität). ✓

**Schritt 3:** $p \mid n$ (da $p^2 \mid n$) und $p \mid (n-1)$ $\Rightarrow$ $p \mid n - (n-1) = 1$. ✓

**Schritt 4:** $p \ge 2 \Rightarrow p \nmid 1$. Widerspruch. ✓

Die Klammerbemerkung im Beweis ("Recall: if $p^a \| n$ with $a \ge 2$...") ist korrekt und erhöht die Selbständigkeit des Papers.

### Korollare
**Korollar 1:** Form $n = p_1 \cdots p_k$ mit $k \ge 2$ (nicht $k = 1$, da Lehmer-Zahlen zusammengesetzt sind). ✓

**Korollar 2:** Totient-Formel $\varphi(n) = (p_1-1)\cdots(p_k-1)$ für quadratfreies $n$. ✓
Die explizite Lehmer-Bedingung $(p_1-1)\cdots(p_k-1) \mid p_1 \cdots p_k - 1$ ist korrekt. ✓

### Bemerkung zu $M \cdot \varphi(n) = n-1$
Die Bemerkung, dass der Beweis wörtlich für die Verallgemeinerung gilt, ist korrekt:
$p \mid \varphi(n) \Rightarrow p \mid M \cdot \varphi(n) = n-1 \Rightarrow p \mid 1$. ✓

### Vergleich mit Giuga (Abschnitt 4)
Der Parallelismus wird korrekt herausgearbeitet:
- Giuga: $p^2 \mid n \Rightarrow p \mid \frac{n}{p} \Rightarrow \frac{n}{p}-1 \equiv -1 \pmod{p} \Rightarrow p \mid -1$
- Lehmer: $p^2 \mid n \Rightarrow p \mid \varphi(n) \Rightarrow p \mid (n-1) \Rightarrow p \mid 1$

Beide Argumente nutzen die Unverträglichkeit $p \mid n$ und $p \mid n \pm 1$. Sehr gut beobachtet. ✓

---

## 3. Literatur und Zitate

| Referenz | Bewertung |
|---|---|
| Lehmer 1932 | ✅ Standardreferenz, korrekt |
| Cohen & Hagis 1980 (≥13 Primfaktoren) | ✅ Korrekt |
| Hagis 1988 ($M \cdot \varphi(n) = n-1$) | ✅ Korrekt |
| Banks & Luca 2009 (≥15 Primfaktoren) | ✅ Korrekt |
| Ingwer 2026b, 2026lehmer3 | Interne Referenzen, korrekt |

Alle Literaturangaben sind etablierte, ververifizierbare Quellen. Keine Inkonsistenzen.

---

## 4. Qualität und Darstellung

- **Einleitung:** Klar und präzise, Definition von Lehmer-Zahlen explizit gegeben
- **Beweis:** Selbständig und gut kommentiert; die Klammerbemerkung zu $\varphi(p^a)$ ist hilfreich
- **Abschnitt 3 (Consequences):** Korollar 2 macht die Lehmer-Bedingung für quadratfreie Zahlen explizit — wertvoll für folgende Papers
- **Abschnitt 4 (Comparison):** Exzellent — zeigt die strukturelle Analogie Giuga/Lehmer, motiviert das Gesamtprojekt

---

## 5. Vergleich EN ↔ DE

| Aspekt | EN | DE |
|---|---|---|
| Mathematik | Identisch ✓ | Identisch ✓ |
| Struktur | Identisch | Identisch |
| Bezeichnungen | φ = "totient" | φ = "Totient-Funktion" |
| Keine Fehler | ✓ | ✓ |

---

## 6. Befunde

**Keine Fehler gefunden.** Beide Versionen sind mathematisch korrekt, vollständig und gut geschrieben.

---

## 7. Empfehlung

**Annahme ohne Auflagen** (beide Sprachversionen).
Das Paper ist das Lehmer-Analogon zu Paper 4 und liefert die nötige Grundlage für Paper 7.
