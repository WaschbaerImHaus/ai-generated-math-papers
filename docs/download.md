---
layout: default
title: Download
---

# Download

## Aktuelle Version — Build 210

Alle 97 mathematischen Papers (je EN + DE) als LaTeX-Quellen und kompilierte PDFs.

---

## Papers-ZIP (aktuelles Release)

<div style="margin: 2em 0; text-align: center;">
  <a href="https://github.com/WaschbaerImHaus/ai-generated-math-papers/releases/latest"
     style="display:inline-block; padding: 0.8em 2em; background:#2d9f4e; color:#fff;
            border-radius:6px; font-size:1.1em; text-decoration:none; font-weight:bold;">
    ⬇ Neuestes Release herunterladen
  </a>
</div>

Das Release enthält:
- **`papers-latex.zip`** — Alle 97 LaTeX-Quelldateien (EN + DE, Batches 1–25)
- **`papers-pdf.zip`** — Alle kompilierten PDFs (195 Dateien)

> Releases werden bei jedem neuen Build automatisch über GitHub Actions erstellt.

---

## Gesamtes Repository

Das vollständige Repository (inkl. Python-Module, Tests, Research-Dateien):

```
https://github.com/WaschbaerImHaus/ai-generated-math-papers/archive/refs/heads/main.zip
```

<a href="https://github.com/WaschbaerImHaus/ai-generated-math-papers/archive/refs/heads/main.zip"
   style="display:inline-block; padding:0.5em 1.5em; background:#159?; color:#fff;
          border:1px solid #ccc; border-radius:4px; text-decoration:none;">
  ⬇ Repository als ZIP
</a>

---

## Nur Papers-Verzeichnis (via Git)

```bash
# Komplettes Repo klonen
git clone https://github.com/WaschbaerImHaus/ai-generated-math-papers.git

# Nur papers/ auschecken (sparse checkout):
git clone --no-checkout https://github.com/WaschbaerImHaus/ai-generated-math-papers.git
cd ai-generated-math-papers
git sparse-checkout init --cone
git sparse-checkout set papers
git checkout main
```

---

## PDFs selbst kompilieren

Voraussetzung: `texlive-full` oder `pdflatex` installiert.

```bash
git clone https://github.com/WaschbaerImHaus/ai-generated-math-papers.git
cd ai-generated-math-papers
bash build_pdfs.sh
# PDFs landen in papers-pdf/
```

---

## Inhalt der Papers (97 Papers, Batches 1–25)

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

Alle Inhalte stehen unter der **MIT-Lizenz** — kostenlos nutzbar, modifizierbar und weitergabe-berechtigt unter Nennung des Autors **Michael Fuhrmann**.

[→ LICENSE anzeigen](https://github.com/WaschbaerImHaus/ai-generated-math-papers/blob/main/LICENSE)

---

*Autor: Michael Fuhrmann | Generiert von [Claude Code](https://claude.ai/claude-code)*
