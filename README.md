# Van der Veen hydrofracture as a damage law

This repository rewrites surface and basal crevasse criteria as an explicit,
irreversible evolution law for the fractured thickness fraction

\[
D = D_s + D_b = \frac{d_s+d_b}{H}.
\]

The technical note distinguishes the literal van der Veen LEFM criterion from
the cheaper zero-crack-tip-stress closure commonly used in large-scale
crevasse-depth models. Both use the same kinetic wrapper, whose total-damage
source is

\[
\frac{\mathrm{D}D}{\mathrm{D}t}
=(1-D)^p\left[
\frac{\langle\Psi_s\rangle_+^{m_s}}{\tau_s}
+\frac{\langle\Psi_b\rangle_+^{m_b}}{\tau_b}
\right].
\]

## Files

- `damage_law.tex` — complete derivation, assumptions, model variants, and implementation notes
- `damage_law.pdf` — compiled technical note
- `references.bib` — bibliography
- `Makefile` — reproducible PDF build

## Build

```sh
make
```

This requires `latexmk`, `pdflatex`, and `bibtex`.

## Key modeling choice

The timescales \(\tau_s,\tau_b\), rate exponents \(m_s,m_b\), and saturation
exponent \(p\) are new constitutive parameters. The original van der Veen
criterion specifies when a crack propagates or arrests; it does not specify a
propagation rate.
