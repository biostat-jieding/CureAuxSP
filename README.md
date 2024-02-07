# CureAuxSP
This is an R package used to fit cure models with auxiliary survival Information.
The underlying methods are based on the paper titled "Efficient auxiliary information synthesis for cure rate model", which has been published in Journal of the Royal Statistial Society Series C with DOI: 10.1093/jrsssc/qlad106.

### Available models

We main focus on cure models, including:
- Both PH mixture cure model and AFT mixture cure model (Semiparametric) with and without auxiliary information

### Introduction to main functions
- **SMC.AuxSP**: fit the semiparametric mixture cure model with auxiliary subgroup survival probabilities (with various options of arguments)
- **print.SMC.AuxSP**: print outputted results from SMC.AuxSP().
