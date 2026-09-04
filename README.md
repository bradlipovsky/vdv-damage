# Van der Veen crevasse mechanics as a damage law

This repository rewrites the van der Veen surface and basal crevasse criteria
as a prognostic law for the fractured thickness fraction

$$
D=D_s+D_b=\frac{d_s+d_b}{H}.
$$

The literal LEFM form is

\[
\frac{\mathrm{D}D}{\mathrm{D}t}
=\frac{v_s}{H}\left\langle\frac{K_s}{K_{\mathrm{Ic}}}-1\right\rangle_+
+\frac{v_b}{H}\left\langle\frac{K_b}{K_{\mathrm{Ic}}}-1\right\rangle_+,
\qquad D<1.
\]

The stress-intensity factors and fracture threshold follow van der Veen. The
linear relation between fracture excess and crack-tip speed is a new,
explicitly stated constitutive assumption. The note also gives a reduced
algebraic zero-stress version and observational predictions that could test
the assumed kinetics.

## Files

- damage_law.tex — model assumptions, derivation, predictions, and limitations
- damage_law.pdf — compiled technical note
- references.bib — bibliography
- AGENTS.md — repository research instructions
- brad-lipovsky-academic-style-guide.md — academic prose guide
- Makefile — reproducible PDF build

## Build

Run make. The build requires latexmk, pdflatex, and bibtex.
