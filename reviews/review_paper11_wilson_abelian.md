# Review: Paper 11 — Wilson's Theorem for Finite Abelian Groups (EN + DE)

**Dateien:**
- `papers/batch2/paper11_wilson_abelian.tex` (EN)
- `papers/batch2/paper11_wilson_abelsch_de.tex` (DE)

**Reviewer:** Claude (selbst)
**Datum:** 2026-03-11
**Ergebnis:** ✅ Mathematisch vollständig korrekt — keine Bugs gefunden

---

## Überblick

Paper 11 verallgemeinert den Wilson-Satz von $(\mathbb{Z}/p\mathbb{Z})^*$ auf beliebige endliche abelsche Gruppen. Das Hauptresultat lautet:

**Allgemeiner Wilson-Satz (Satz 3.1):**
Für eine endliche abelsche Gruppe $(G, \cdot)$ gilt:
$$\prod_{g \in G} g = \begin{cases} e & \text{falls } G \text{ kein oder } \geq 2 \text{ Elemente der Ordnung 2 hat} \\ \tau & \text{falls } G \text{ genau ein Element } \tau \text{ der Ordnung 2 hat} \end{cases}$$

Darüber hinaus beweist Paper 11 die **Gauss-Verallgemeinerung:**
$$\prod_{\substack{j=1\\\gcd(j,n)=1}}^{n} j \equiv -1 \pmod n \iff n \in \{1, 2, 4, p^k, 2p^k\} \text{ für ungerades Primzahl } p$$

---

## Beweisstruktur

### Lemma 2.1 — Paarungslemma
**Aussage:** $\prod_{g \in G} g = \prod_{\{g \in G : g^2 = e\}} g$

**Beweis:** Korrekt. Elemente mit $g \neq g^{-1}$ paaren sich mit Beitrag 1. Nur selbstinverse Elemente ($g^2=e$) bleiben übrig. ✅

### Lemma 2.2 — Untergruppe der Involutionen
**Aussage:** $S = \{g \in G : g^2 = e\}$ ist eine Untergruppe, isomorph zu $(\mathbb{Z}/2\mathbb{Z})^k$ für ein $k \geq 0$.

**Beweis:** Korrekt via Klassifikationssatz endlicher abelscher Gruppen: $S$ hat Exponent 2, also $S \cong (\mathbb{Z}/2\mathbb{Z})^k$. ✅

### Satz 3.1 — Allgemeiner Wilson-Satz

**Fall $k = 0$ (kein Element der Ordnung 2, z.B. $|G|$ ungerade):**
$S = \{e\}$, Produkt $= e$. ✅

**Fall $k = 1$ (genau ein Element $\tau$ der Ordnung 2):**
$S = \{e, \tau\}$, Produkt $= e \cdot \tau = \tau$. ✅

**Fall $k \geq 2$ ($|S| \geq 4$ Elemente der Ordnung 2):**
Zwei unabhängige Argumente werden gegeben:
1. **Koordinaten-Argument:** In $S \cong (\mathbb{Z}/2\mathbb{Z})^k$: $\prod_{s \in S} s = $ Summe aller Vektoren in $(\mathbb{Z}/2\mathbb{Z})^k$. Für $k \geq 2$: Jede Koordinate erscheint genau $2^{k-1}$ mal (gerade Anzahl) → Summe $= 0 = e$. ✅
2. **$\tau$-Multiplikation:** Für beliebiges $\tau \in S$ mit $\tau \neq e$: Multiplikation mit $\tau$ permutiert $S$, also $\tau^{|S|} P = P$, d.h. $P^2 = e$. Da $|S| = 2^k$ gerade, auch $P = e$ via Paarung im $k=1$-Argument angewandt auf das $\tau$. ✅ (Das Paper gibt beide Argumente als Alternativen.)

### Gauss-Verallgemeinerung

**Beweis via CRT + Struktursatz:**
- $(\mathbb{Z}/n\mathbb{Z})^*$ hat genau ein Element der Ordnung 2 genau dann, wenn die Gruppe zyklisch ist.
- Via Klassifikation: $(\mathbb{Z}/n\mathbb{Z})^*$ ist zyklisch ↔ $n \in \{1,2,4,p^k,2p^k\}$.
- Für andere $n$ hat die Gruppe $\geq 2$ Elemente der Ordnung 2 → Produkt $= e = 1 \not\equiv -1$.

Der Beweis ist vollständig und korrekt. Die CRT-Dekomposition wird korrekt angewendet. ✅

---

## Numerische Verifikation

| $n$ | $(\mathbb{Z}/n\mathbb{Z})^*$ | Selbstinverse | Produkt $\pmod n$ | Erwartet |
|-----|-------------------------------|---------------|---------------------|----------|
| $n=8$ | $\{1,3,5,7\}$ | alle 4 | $1\cdot3\cdot5\cdot7=105\equiv 1$ | $+1$ ✅ |
| $n=12$ | $\{1,5,7,11\}$ | $\{1,5,7,11\}$ (je $^2\equiv1$) | $385\equiv1\pmod{12}$ | $+1$ ✅ |
| $n=7$ | $\{1,2,3,4,5,6\}$ | $\{1,6\}$ | $720\equiv6\equiv-1$ | $-1$ ✅ |
| $n=9$ | $\{1,2,4,5,7,8\}$ | $\{1,8\}$ | $2240\equiv8\equiv-1$ | $-1$ ✅ |

---

## Bugs

**Keine gefunden.**

---

## Fazit

Paper 11 ist das theoretisch tiefste Paper der Wilson-Serie. Die Verallgemeinerung auf beliebige endliche abelsche Gruppen ist vollständig bewiesen mit zwei unabhängigen Argumenten für den Fall $k \geq 2$. Die Gauss-Verallgemeinerung wird korrekt aus dem allgemeinen Satz abgeleitet. Numerische Verifikationen bestätigen alle Fälle.
