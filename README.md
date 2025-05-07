# Radigen

**Radigen** is a Python-based kernel for simulating radical oxidation mechanisms in complex organic mixtures, including **edible oils**, **biofuels**, and **polymer materials**. It supports the generation, structuring, and numerical resolution of large-scale reaction networks arising from reactive functions and radical-mediated pathways.



## 🔬 Scientific Foundation

This project is grounded in two recent peer-reviewed publications:

1. Touffet M., Smith P., Vitrac O.  
   *A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures*,  
   **Food Research International**, 173(1), 2023, 113289.  
   [https://doi.org/10.1016/j.foodres.2023.113289](https://doi.org/10.1016/j.foodres.2023.113289)

2. Touffet M., Vitrac O.  
   *Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms*,  
   **Food Chemistry**, 481, 2025, 143952.  
   [https://doi.org/10.1016/j.foodchem.2025.143952](https://doi.org/10.1016/j.foodchem.2025.143952)

   

## 🧠 Core Module: `radigen3.oxidation.py`

This module defines the **combinatorial reaction kernel** and related infrastructure for:

- Reactive species definition with radical and functional group metadata.
- Automatic generation of chemical products using known transformation rules.
- Recursive construction of reaction networks (mono- and bimolecular).
- Canonical fingerprinting of reactions for parameter lookup.
- Association of known reaction rates (`k0`, `Ea`) from literature.
- Internal registry for `reactionRateDB` with confidence intervals and structured querying.
- Support for homolytic and heterolytic decompositions.
- 

## 🧪 Quick Start

```python
from radigen3.oxidation import mixture

oil = mixture()
oil.add("L1H", concentration=3000)       # fatty ester (use P1H for polylers, R1H, R2H, R3H for a generic substance where x in RxH indicates the allicity of the labile proton)
oil.add("L1OOH", concentration=100)      # hydroperoxide
oil.add("O2", concentration=10)          # dissolved oxygen

oil.addProducts()                        # infer and add all products
oil.addReactions()                       # build A+B -> C+D reaction scheme
oil.populateReactionRates()             # populate Ea and k0 from registry
````

`oil.reactions` yields

```
 *R0: L1H + HO• -> L1• + H2O
 *R1: L1H + L1O• -> L1• + L1OH
 *R2: L1H + L1OO• -> L1• + L1OOH
 *R3: L1OOH -> L1O• + HO•
 *R4: L1OOH + L1OOH -> L1O• + L1OO•
 *R5: L1• + O2 -> L1OO•
 *R6: L1• + L1• -> L1-polymer + L1-polymer
 *R7: L1• + L1OO• -> L1-polymer + L1-polymer
 *R8: L1O• -> L1=O
 *R9: L1OO• + L1OO• -> L1-polymer + L1-polymer
```

where `*` indicates that the reaction has been parameterized



## 📘 Features in Progress

* [x] Combinatorial generation of products and reaction networks
* [x] Canonical reaction fingerprinting for tabulated rates
* [x] Detection of monomolecular and bimolecular topologies
* [ ] Diffusion-limited rate correction for fast radical reactions
* [ ] Time-resolved kinetic solver with temperature coupling
* [ ] Dynamic simulation under variable O₂ conditions

## 📂 Project Structure

```
radigen/
├── radigen3/
│   ├── oxidation.py           # main module for species, reactions, rates
│   ├── ...
├── examples/
│   └── simulate_oil.py        # sample kinetic setup for edible oil oxidation (coming soon)
├── tests/
│   └── test_oxidation.py      # unit tests for mechanism generation (coming soon)
├── README.md
```

## 📜 License

MIT License – see [`LICENSE`](./LICENSE) file.

## 🧠 Credits

Developed by the **Generative Simulation Initiative**, led by Olivier Vitrac
Contact: *[olivier.vitrac@gmail.com](mailto:olivier.vitrac@gmail.com)*

---

Last updated:  $2025-05-07$