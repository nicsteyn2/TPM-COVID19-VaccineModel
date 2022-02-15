# A COVID-19 Vaccination Model for Aotearoa New Zealand

Two implementations of a model for COVID-19 spread in a partially vaccinated population. Default parameters reflect New Zealand's population, vaccination programme, and non-pharmaceutical interventions used.

All code is implemented in MATLAB.

The first implementation is a deterministic SEIR-style model (see `runSEIR.m`). This is more useful for considering population-level spread, the kind of which might be expected when borders are re-opened. The second implementation is a stochastic branching process model (see `runBranching.m`). This is more useful for considering small outbreaks, seeded by a single re-incursion.

Default parameters are included in `getPar.m`. Simple examples are provided in `example.m`. Code to reproduce the results from the paper is found in `resultsSEIR.m` and `resultsBranching.m`.

For more details:
> [Steyn. N., Plank. M.J., Binny. R.N., Hendy. S.C., Lustig. A., Ridings. K. A COVID-19 Vaccination Model for Aotearoa New Zealand. Scientific Reports. (2022).](https://www.tepunahamatatini.ac.nz/2021/06/30/a-covid-19-vaccination-model-for-aotearoa-new-zealand/)
