# Van der Veen crevasse mechanics as a damage law

This repository rewrites the van der Veen surface and basal crevasse criteria
as a prognostic law for the fractured thickness fraction

$$
D=D_s+D_b=\frac{d_s+d_b}{H}.
$$

The literal LEFM form is

$$
\frac{\mathrm{D}D}{\mathrm{D}t}
=\frac{v_s}{H}\mathcal F\!\left(\frac{K_{\mathrm I}^{\mathrm s}}{K_{\mathrm{Ic}}}\right)
+\frac{v_b}{H}\mathcal F\!\left(\frac{K_{\mathrm I}^{\mathrm b}}{K_{\mathrm{Ic}}}\right),
\qquad D<1.
$$

The Freund-style kinetic function is

$$
\mathcal F(q)=
\begin{cases}
0, & q<1,\\
1-q^{-2}, & q\geq 1.
\end{cases}
$$

The stress-intensity factors and fracture threshold follow van der Veen. The
kinetic closure follows the Freund formulation used by Lipovsky (2018). The
note also gives a reduced algebraic zero-stress version and observational
predictions that could test the assumed kinetics.

## Files

- damage_law.tex — model assumptions, derivation, predictions, and limitations
- damage_law.pdf — compiled technical note
- damage_law_demo.ipynb — minimal NumPy implementation with illustrative figures
- references.bib — bibliography
- AGENTS.md — repository research instructions
- brad-lipovsky-academic-style-guide.md — academic prose guide
- Makefile — reproducible PDF build

## Build

Run make. The build requires latexmk, pdflatex, and bibtex.

Run the notebook from top to bottom with Jupyter. It requires only NumPy and
Matplotlib.
