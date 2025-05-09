# Radigen

**Radigen** is a Python-based kernel for simulating radical oxidation mechanisms in complex organic mixtures, including **edible oils**, **biofuels**, and **polymer materials**.

It supports the **generative construction**, **parameterization**, and **numerical resolution** of complex reaction networks involving reactive functions, radicals, and temperature/transport effects. This kernel is also the basis for **LLM integration** in the **Generative Simulation Initiative**.

------

## 🔬 Scientific Foundation

Radigen builds on two mechanistic modeling frameworks published in peer-reviewed journals:

1. **Touffet M., Smith P., Vitrac O.**
    *A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures*,
    *Food Research International*, 173(1), 2023, 113289.
    https://doi.org/10.1016/j.foodres.2023.113289
2. **Touffet M., Vitrac O.**
    *Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms*,
    *Food Chemistry*, 481, 2025, 143952.
    https://doi.org/10.1016/j.foodchem.2025.143952

------

## 🧠 Core Capabilities

- ✅ **Declarative species and reaction construction**
- ✅ **Combinatorial generation** of mono- and bimolecular reactions with canonical ordering
- ✅ **Cross-reaction inheritance** and parameter inference from self-reactions
- ✅ **Rate constant assignment** (`k₀`, `Eₐ`) from a structured database (`reactionRateDB`)
- ✅ **Temperature and viscosity effects** via Arrhenius + Smoluchowski models
- ✅ **Equilibrium modeling** for cage vs. free hydroperoxide decomposition
- ✅ **O₂ transport coupling** via Henry's law and first-order exchange kinetics
- ✅ **Sparse/dense stoichiometric matrix generation**
- ✅ **Time-resolved ODE integration** (based on `scipy.integrate.solve_ivp`)
- ✅ **Plotting and DataFrame export** for concentration profiles
- ✅ **Construction of lumped species** based on chemical structures for rapid interpretation and reporting

------

## 🤖 LLM Integration

Radigen is designed as a kernel for **LLM-driven chemistry**. Its architecture supports:

- Modular inspection and manipulation of species and reactions
- Natural language mapping to chemical entities (`C1H`, `L1OOH`, etc.)
- Symbolic access to mechanisms, rate equations, and observables
- Prompt-based generation of mixtures, mechanisms, and simulations

This kernel supports ongoing developments in generative scientific simulation and **LLM fine-tuning on physical–chemical representations**.

------

## 🧪 Quick Start

```python
from radigen3.oxidation import mixture

oil = mixture()
oil.add("L1H", concentration=3000)
oil.add("L1OOH", concentration=100)
oil.add("O2", concentration=10)

oil.addProducts()
oil.addReactions()
oil.populateReactionRates()

oilmodel = mixtureKinetics(oil)
oilmodel.solve(3600*24*10, T=60)  # 10 days at 60°C
oilmodel.plot()
df = oilmodel.results_as_dataframe(["L1H", "L1OOH", "L1O•"])
print(df.head())
```

------

```python
oil = mixture()
oil.add("L1H",concentration=3000)
oil.add("L2H",concentration=1000)
oil.add("L3H",concentration=500)
oil.add("L1OOH",concentration=100)
oil.add("O2",concentration=10)
oil.addProducts()
oil.addReactions()
oil.populateReactionRates()

oilmodel = mixtureKinetics(oil) # kinetic model
oilmodel.solve(10*24*3600,60)
oilmodel.plot()
df = oilmodel.results_as_dataframe(["L1H","L2H","L3H","L1OOH","L2OOH","L3OOH"])
print(df)
```

<details>
	<summary>Click here to see results</summary>

⏱ 59 reactions involving 27 species were generated, including:

```
 *R0: L1H + HO• → L1• + H2O,
 *R1: L1H + L1O• → L1• + L1OH,
 *R2: L1H + L1OO• → L1• + L1OOH,
 *R3: L1H + L2OO• → L1• + L2OOH,
 *R4: L1H + L3OO• → L1• + L3OOH,
 *R5: L1H + L2O• → L1• + L2OH,
 *R6: L1H + L3O• → L1• + L3OH,
 *R7: L2H + HO• → L2• + H2O,
 *R8: L2H + L1O• → L2• + L1OH,
 *R9: L2H + L1OO• → L2• + L1OOH,
 *R10: L2H + L2OO• → L2• + L2OOH,
 *R11: L2H + L3OO• → L2• + L3OOH,
 *R12: L2H + L2O• → L2• + L2OH,
 *R13: L2H + L3O• → L2• + L3OH,
 *R14: L3H + HO• → L3• + H2O,
 *R15: L3H + L1O• → L3• + L1OH,
 *R16: L3H + L1OO• → L3• + L1OOH,
 *R17: L3H + L2OO• → L3• + L2OOH,
 *R18: L3H + L3OO• → L3• + L3OOH,
 *R19: L3H + L2O• → L3• + L2OH,
 *R20: L3H + L3O• → L3• + L3OH,
 *R21: L1OOH → L1O• + HO•,
 *R22: L1OOH + L1OOH → L1O• + L1OO•,
 *R23: L1OOH + L2OOH → L1O• + L2OO•,
 *R24: L1OOH + L2OOH → L2O• + L1OO•,
 *R25: L1OOH + L3OOH → L1O• + L3OO•,
 *R26: L1OOH + L3OOH → L3O• + L1OO•,
 *R27: L1• + O2 → L1OO•,
 *R28: L1• + L1• → L1-polymer + L1-polymer,
 *R29: L1• + L2• → L1-polymer + L2-polymer,
 *R30: L1• + L3• → L1-polymer + L3-polymer,
 *R31: L1• + L1OO• → L1-polymer + L1-polymer,
 *R32: L1• + L2OO• → L1-polymer + L2-polymer,
 *R33: L1• + L3OO• → L1-polymer + L3-polymer,
 *R34: L2• + O2 → L2OO•,
 *R35: L2• + L2• → L2-polymer + L2-polymer,
 *R36: L2• + L3• → L2-polymer + L3-polymer,
 *R37: L2• + L1OO• → L2-polymer + L1-polymer,
 *R38: L2• + L2OO• → L2-polymer + L2-polymer,
 *R39: L2• + L3OO• → L2-polymer + L3-polymer,
 *R40: L3• + O2 → L3OO•,
 *R41: L3• + L3• → L3-polymer + L3-polymer,
 *R42: L3• + L1OO• → L3-polymer + L1-polymer,
 *R43: L3• + L2OO• → L3-polymer + L2-polymer,
 *R44: L3• + L3OO• → L3-polymer + L3-polymer,
 *R45: L1O• → L1=O,
 *R46: L1OO• + L1OO• → L1-polymer + L1-polymer,
 *R47: L1OO• + L2OO• → L1-polymer + L2-polymer,
 *R48: L1OO• + L3OO• → L1-polymer + L3-polymer,
 *R49: L2OO• + L2OO• → L2-polymer + L2-polymer,
 *R50: L2OO• + L3OO• → L2-polymer + L3-polymer,
 *R51: L3OO• + L3OO• → L3-polymer + L3-polymer,
 *R52: L2OOH → L2O• + HO•,
 *R53: L2OOH + L2OOH → L2O• + L2OO•,
 *R54: L2OOH + L3OOH → L2O• + L3OO•,
 *R55: L2OOH + L3OOH → L3O• + L2OO•,
 *R56: L3OOH → L3O• + HO•,
 *R57: L3OOH + L3OOH → L3O• + L3OO•,
 *R58: L2O• → L2=O,
 *R59: L3O• → L3=O
```
*A star indicates that the reaction has been correctly parameterized*.
</details>

------

## 📘 Feature Progress

| Feature                                                      | Status             |
| ------------------------------------------------------------ | ------------------ |
| Reaction fingerprinting + registry lookup                    | ✅ done             |
| Monomolecular & bimolecular pathway logic                    | ✅ done             |
| Cross-reaction inheritance                                   | ✅ done             |
| Transport coupling (O₂ interface kinetics)                   | ✅ done             |
| Cage ↔ free hydroperoxide equilibrium                        | ✅ done             |
| Time-resolved simulation (solve_ivp)                         | ✅ done             |
| Diffusion-limited reaction modeling                          | ✅ done             |
| Data export (DataFrame, plotting)                            | ✅ done             |
| Creation of arbitrary “lumped” species from initial conditions and simulated results | ✅ done             |
| GUI or LLM prompt wrappers                                   | 🚧 work in progress |

------

## 📂 Project Structure

```
radigen/
├── radigen3/
│   ├── oxidation.py         # main kernel
│   └── ...
├── examples/
│   └── simulate_oil.py      # full oxidation scenario
├── tests/
│   └── test_oxidation.py    # unit tests
├── README.md
```

------


## 🧩 Overview of Core Classes
| Class             | Purpose                                                                 | Key Features                                                                                    |
| ----------------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `species`         | Represents a chemical species, radical or non-radical                   | Tracks identity, concentration, radical status, reactive functions, origin and products         |
| `reaction`        | Encodes a single chemical reaction with stoichiometry and type          | Handles bimolecular/monomolecular logic, canonical fingerprinting, and directionality           |
| `reactionRateDB`  | Stores and queries rate constants (`k₀`, `Eₐ`) for known reaction types | Supports fingerprint lookups, priority ranks, redox heuristics, and confidence tagging          |
| `mixture`         | Container for species and reactions in a chemical system                | Builds product and reaction networks automatically from a few initial species                   |
| `mixtureKinetics` | Numerical kinetic solver based on a `mixture` instance                  | Integrates ODEs, applies temperature & diffusion corrections, and returns time-resolved results |
| `lumped`          | Groups multiple species under a shared name or tag (e.g., LOOH)         | Supports species aliases, output simplification, and aggregate quantification                   |

------

## 📦 Lumped Species

Radigen supports the **automatic grouping of chemically related species** into *lumped species*, allowing users to visualize and extract **aggregate concentrations**. These are especially useful for:

- Hydroperoxides (e.g., `L1OOH`, `L2OOH`, `L3OOH`)
- Aldehydes, ketones, alcohols
- Radicals centered on carbon or oxygen
- Any class of polar compounds (with at least one `O` atom)

<details>
	<summary>Click here to expand</summary>

This feature relies on two mechanisms:

### 1 | **Auto-Grouping from `mixture`**

```python
# Get a lumped object containing all hydroperoxides
rooh = oil.lumped_hydroperoxides()

# Similarly for aldehydes, radicals, etc.
radC = oil.lumped_radicals_on_C()
polar = oil.lumped_polar_compounds()
```

Each method returns a `lumped` object with a `.name` and `.species` list. These objects can be passed to `mixtureKinetics`.

### 2 | **Registration in Kinetics Model**

You can register lumped species to be **plotted or exported as if they were real species**:

```python
oilmodel = mixtureKinetics(oil)
oilmodel.register_lumped("LOOH", rooh)           # name can be arbitrary
oilmodel.register_lumped("radicals_C", radC)

oilmodel.solve(10*24*3600, T=60)
oilmodel.plot()                                  # includes lumped if species=None
df = oilmodel.results_as_dataframe()             # same: lumped automatically appended
```

If you only want to include specific lumped species:

```python
oilmodel.plot(["L1H", "LOOH", "radicals_C"])
df = oilmodel.results_as_dataframe(["L1H", "LOOH", "radicals_C"])
```

### 🔧 Custom Grouping by Function or Pattern

Advanced users can create lumped species on the fly:

```python
# Group all with reactive function == 'COOH'
l1 = oil.get_lumped_by_function("COOH")

# Group all species whose name matches a regex
peroxyls = oil.get_lumped_by_pattern(r".*OO•$")
```

These groups can be registered and used identically to individual species.

</details>

------

## 📜 License

MIT License – see [`LICENSE`](https://chatgpt.com/c/LICENSE) file.

------

## 🧠 Credits

Developed by the **Generative Simulation Initiative**
 Lead: Olivier Vitrac
 Contact: *[olivier.vitrac@gmail.com](mailto:olivier.vitrac@gmail.com)*

------

*Last updated: 2025-05-09*

