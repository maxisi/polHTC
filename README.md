polHTC
======

_PolHTC_ is a package for the analysis of nontensorial continuous gravitational waves (CWs) from known galactic pulsars.
The targeted signals are nonstandard in that the are allowed to carry polarization content beyond what is allowed in Einstein's general relativity (GR), including vector and tensor modes in addition to the tensor modes predicted by GR.

The name _polHTC_ refers to the fact the code provides utilities to set up distributed workflows in high-throughput computing (HTC) clusters using [HTCondor](https://htcondor.org).

This code was used in the publication _Detecting Beyond-Einstein Polarizations of Continuous Gravitational Waves_ by M. Isi et al (2015) [[arXiv:1502.00333](https://arxiv.org/abs/1502.00333)].

The code takes as input LIGO-Virgo data preheterodyned by the code described in Pitkin et al [[arXiv:1705.08978](https://arxiv.org/abs/1705.08978)].

### Attribution

If you find this code useful, please cite the following paper:

```
@article{Isi:2015cva,
    author = "Isi, Maximiliano and Weinstein, Alan J. and Mead, Carver and Pitkin, Matthew",
    title = "{Detecting Beyond-Einstein Polarizations of Continuous Gravitational Waves}",
    eprint = "1502.00333",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "LIGO-P1400169",
    doi = "10.1103/PhysRevD.91.082002",
    journal = "Phys. Rev. D",
    volume = "91",
    number = "8",
    pages = "082002",
    year = "2015"
}
```
