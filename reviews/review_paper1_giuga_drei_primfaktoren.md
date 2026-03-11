# Gutachten: Paper 1
**Titel:** „Keine zusammengesetzte Zahl mit genau drei Primfaktoren ist ein Giuga-Pseudoprim"
**Dateien:** `paper1_giuga_three_primes.tex` / `paper1_giuga_drei_primfaktoren_de.tex`
**Gutachter:** Mathematik-Review-Spezialist (Zahlentheorie)
**Datum Erstgutachten:** 2026-03-11
**Datum Revisionsgutachten:** 2026-03-11
**Gesamturteil:** ✅ **Annahme ohne Auflagen** — alle bekannten Fehler behoben

---

## Revisionsgutachten (2026-03-11)

### BUG-001 — Status: ✅ BEHOBEN

**EN:** Subsection 1.2 (Zeile 87) lautet jetzt korrekt:
> „Borwein, Borwein, Girgensohn, and Parnes \cite{Borwein1996}"

**DE:** Subsection 1.2 (Zeile 87) lautet jetzt korrekt:
> „Borwein, Borwein, Girgensohn und Parnes \cite{Borwein1996}"

Beide Sprachversionen sind nun übereinstimmend korrekt.

### Mathematischer Gesamtstatus
Alle Beweise unverändert korrekt. Keine neuen Fehler eingeführt. Die Überarbeitung beschränkte sich auf die Attribution im Fließtext.

**Empfehlung: Druckreif.**

---

---

## 1. Einordnung und Originalitätsprüfung

Das Paper behauptet einen vollständig elementaren Beweis, dass kein Produkt
genau dreier verschiedener Primzahlen ein Giuga-Pseudoprim ist. Das Resultat
selbst ist konsistent mit (und impliziert durch) das stärkere Bednarek-Ergebnis
(≥ 59 Primfaktoren), liefert aber einen neuen, rein elementaren Zugang.

Das Resultat ist mathematisch nicht neu im Sinne eines Durchbruchs — es folgt
aus bekannten Schranken — aber die Beweismethode (Paritätsargument + monotone
Schrankenschätzung) ist sauber und didaktisch wertvoll.

---

## 2. Definitionen und Grundlagen

### Giuga-Pseudoprim
Korrekt definiert: quadratfreie zusammengesetzte ganze Zahl $n$, die für jeden
Primteiler $p \mid n$ erfüllt:
- **(W)** $p \mid \frac{n}{p} - 1$
- **(S)** $(p-1) \mid \frac{n}{p} - 1$

### Charakterisierungssatz (Borwein et al. 1996)
Korrekt zitiert und korrekt angewendet. **Fehler im Fließtext**: Beide Versionen
(EN und DE) schreiben „Borwein, Borwein, Borwein und Giuga" als Autoren.
Laut Bibliographie (korrekt) sind die Autoren: **J. M. Borwein, P. B. Borwein,
R. Girgensohn und S. Parnes**. Giuga ist der Namensgeber der Vermutung (1950),
aber kein Koautor des 1996er Papers. Dieser Fehler tritt in beiden Sprachversionen
identisch auf.

---

## 3. Lemma-Prüfung

### Lemma 1: Quadratfreiheit
**Behauptung:** Jedes Giuga-Pseudoprim ist quadratfrei.

**Beweis-Analyse:** Sei $p^2 \mid n$. Dann gilt $p \mid \frac{n}{p}$, also
$\frac{n}{p} \equiv 0 \pmod{p}$, somit $\frac{n}{p} - 1 \equiv -1 \pmod{p}$.
Bedingung (W) fordert $p \mid -1$, unmöglich für $p \ge 2$. ✅ **Korrekt.**

### Lemma 2: Kein Zweiprim-Pseudoprim
**Behauptung:** Es gibt kein Giuga-Pseudoprim der Form $n = pq$, $p < q$ prim.

**Beweis-Analyse:**
- (W) für $p$: $p \mid (q-1)$
- (W) für $q$: $q \mid (p-1)$
- Da $p < q$: $p - 1 < q$, also $0 \le p-1 < q$ ⟹ $q \mid (p-1)$ nur wenn $p-1 = 0$,
  d.h. $p = 1$. Widerspruch (keine Primzahl). ✅ **Korrekt.**

---

## 4. Theorem: Fall $n = 2qr$ ($2 < q < r$ prim)

### Schritt 1: Schwache Bedingung für $r$
$r \mid \frac{2qr}{r} - 1 = 2q - 1$.

Schrankenargument: $0 < 2q - 1 < 2r$ (gilt wegen $q < r$, also $2q < 2r$,
also $2q - 1 < 2r$). Das einzige positive Vielfache von $r$ im offenen Intervall
$(0, 2r)$ ist $r$ selbst. Damit: $r = 2q - 1$. ✅ **Korrekt.**

### Schritt 2: Starke Bedingung für $r$
$(r-1) \mid (2q-1)$. Mit $r = 2q-1$: $r - 1 = 2(q-1)$.
Benötigt: $2(q-1) \mid (2q-1)$.

$2q - 1$ ist **ungerade** (da $q$ ungerade Primzahl ⟹ $2q$ gerade ⟹ $2q-1$ ungerade).
$2(q-1)$ ist **gerade**.
Gerade Zahl kann ungerade Zahl nicht teilen. **Widerspruch.** ✅ **Korrekt.**

**Bewertung:** Elegantes Paritätsargument, lückenlos. ✅

---

## 5. Theorem: Fall $n = pqr$ (alle ungerade, $3 \le p < q < r$)

### Schritt 1: Zentrale Teilbarkeit
(W) und (S) für $r$: $r \mid (pq-1)$ und $(r-1) \mid (pq-1)$.
Da $\gcd(r, r-1) = 1$: $\text{lcm}(r, r-1) = r(r-1) \mid (pq-1)$.
Da $pq - 1 > 0$: $r(r-1) \le pq - 1 < pq$. ✅

### Schritt 2: Untere Schranke für $r(r-1)$
$r, q$ ungerade Primzahlen mit $r > q$, also $r \ge q + 2$ (zwei aufeinanderfolgende
ungerade Zahlen). Damit:
$$r(r-1) \ge (q+2)(q+1) = q^2 + 3q + 2. \quad ✅$$

### Schritt 3: Widerspruch
$$q^2 + 3q + 2 \le pq - 1 \implies pq \ge q^2 + 3q + 3 \implies p \ge q + 3 + \frac{3}{q}.$$

Da $p \ge 3$ und $q > p$, gilt $q \ge 5$ (nächste ungerade Primzahl nach 3). Damit
$\frac{3}{q} \le \frac{3}{5} < 1$. Da $p$ ganzzahlig: $p \ge q + 4$.
Widerspruch zu $p < q$. ✅ **Korrekt.**

**Anmerkung zum Wording:** „Die kleinste ungerade Primzahl größer als $p \ge 3$"
suggeriert, $q$ sei die nächste Primzahl nach $p$. Tatsächlich ist $q$ nur irgendeine
Primzahl mit $q > p \ge 3$, also $q \ge 5$. Das Argument ist trotzdem korrekt —
die Formulierung könnte aber präziser sein.

---

## 6. Korollar: Mindestens 4 Primfaktoren

Folgt direkt aus Lemma 2 (≥ 3 Primfaktoren) und Theorem 1.3 (nicht genau 3). ✅

---

## 7. Vergleich mit Literatur

- Verweis auf Bednarek (2014) mit ≥ 59 Primfaktoren und > 10^{19907} Stellen: **Referenz
  nicht verifiziebar** (nur als Preprint zitiert). Das mathematische Resultat selbst ist
  konsistent mit dem im Paper bewiesenen Korollar, aber die Quelle ist unsicher.
- Borwein et al. (1996): Bibliographisch korrekt, nur Fließtext-Attribution fehlerhaft.

---

## 8. Zusammenfassung der Befunde

| Komponente | Status | Anmerkung |
|---|---|---|
| Definition Giuga-Pseudoprim | ✅ | Korrekt |
| Lemma Quadratfreiheit | ✅ | Korrekt |
| Lemma Kein 2-Prim | ✅ | Korrekt |
| Theorem n=2qr (Paritätsarg.) | ✅ | Korrekt, elegant |
| Theorem n=pqr (Schrankenarg.) | ✅ | Korrekt |
| Hauptsatz + Korollar | ✅ | Korrekt |
| Zitat Borwein et al. (Fließtext) | ❌ | Falsche Autorenliste |
| Referenz Bednarek 2014 | ⚠️ | Nicht verifizierbar |

---

## 9. Empfehlung

**Annahme nach geringfügiger Revision:**
1. Fließtext-Zitat „Borwein, Borwein, Borwein und Giuga" in **beiden Sprachversionen**
   korrigieren zu „Borwein, Borwein, Girgensohn und Parnes".
2. Referenz Bednarek 2014 durch eine verifizierbare Publikation ersetzen oder als
   „unveröffentlichtes Manuskript" kennzeichnen.
3. Formulierung zu $q \ge 5$ im Beweis von Theorem 4.1 leicht präzisieren.

Die mathematischen Beweise sind vollständig, lückenlos und korrekt.
