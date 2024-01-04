# CureAuxSP
This is an R package used to fit cure models with auxiliary survival Information.
The underlying methods are based on the paper titled "Efficient auxiliary information synthesis for cure rate model", which has been published in Journal of the Royal Statistial Society Series C with DOI: 10.1093/jrsssc/qlad106.

More detailed description is still in editing. Thank you for you attention!

### Available models

We main focus on cure models, including:
- Mixture cure models
- Mixture cure models with external survival rates as auxiliary information

### Introduction to main functions
- **SMC.AuxSP.AIS**: fit the semiparametric mixture cure model with homogeneous auxiliary subgroup survival rates using AIS
- **SMC.AuxSP.PAIS**: fit the semiparametric mixture cure model using with heterogeneous auxiliary subgroup survival rates using PAIS
