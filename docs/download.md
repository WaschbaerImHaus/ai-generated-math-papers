---
layout: default
title: Download
---

# Download

## Aktuelle Papers — Build 209

Alle 97 mathematischen Papers (je EN + DE) als kompilierte PDFs.

---

<div style="margin: 2em 0; text-align: center;">
  <a href="https://github.com/WaschbaerImHaus/ai-generated-math-papers/releases/latest"
     style="display:inline-block; padding: 0.8em 2em; background:#2d9f4e; color:#fff;
            border-radius:6px; font-size:1.1em; text-decoration:none; font-weight:bold;">
    ⬇ Neuestes Release herunterladen
  </a>
  &nbsp;
  <a href="https://github.com/WaschbaerImHaus/ai-generated-math-papers/releases"
     style="display:inline-block; padding: 0.8em 2em; background:#555; color:#fff;
            border-radius:6px; font-size:1.1em; text-decoration:none;">
    Alle Releases
  </a>
</div>

Das Release enthält eine ZIP-Datei mit allen kompilierten PDFs der 97 Papers (195 Dateien, Batches 1–25).

---

## Gesamtes Repository klonen

```bash
git clone https://github.com/WaschbaerImHaus/ai-generated-math-papers.git
```

Nur das `papers/`-Verzeichnis (LaTeX-Quellen, sparse checkout):

```bash
git clone --no-checkout https://github.com/WaschbaerImHaus/ai-generated-math-papers.git
cd ai-generated-math-papers
git sparse-checkout init --cone
git sparse-checkout set papers
git checkout main
```

---

## PDFs selbst kompilieren

Voraussetzung: `texlive-full` oder `pdflatex`

```bash
git clone https://github.com/WaschbaerImHaus/ai-generated-math-papers.git
cd ai-generated-math-papers
bash build_pdfs.sh
# PDFs landen in papers-pdf/
```

---

## Inhalt (97 Papers, Batches 1–25)

| Batch | Papers | Thema |
|-------|--------|-------|
| 1–4 | 1–20 | Giuga, Lehmer, Wilson, Siebmethoden, Goldbach, Kreismethode |
| 5–6 | 21–28 | Riemann-Hypothese, Elliptische Kurven, BSD |
| 7–8 | 29–36 | Collatz (Tao), Modulformen, abc-Vermutung, Navier-Stokes |
| 9–10 | 37–39 | Algebraische Zahlentheorie, Iwasawa, Langlands |
| 11–18 | 40–71 | Gruppentheorie, Topologie, Spez. Funktionen, Yang-Mills, Hodge, P vs NP |
| 19–25 | 72–97 | Gruppe B: 26 weitere offene Vermutungen |

---

## Lizenz

**MIT-Lizenz** — frei nutzbar, modifizierbar und weitergabe-berechtigt unter Nennung des Autors **Michael Fuhrmann**.

[→ LICENSE](https://github.com/WaschbaerImHaus/ai-generated-math-papers/blob/main/LICENSE)

---

*Autor: Michael Fuhrmann | Generiert von [Claude Code](https://claude.ai/claude-code)*
