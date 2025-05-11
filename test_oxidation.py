#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RADIGEN – test_oxidation.py

Basic unit test for the RADIGEN oxidation engine.

This script performs a simple, minimal validation of core RADIGEN components:
    1. Simulates oxidation of a mono-unsaturated species (C1H) with hydroperoxides (C1OOH),
    2. Under two conditions:
        a. Anoxic system (no oxygen)
        b. Oxygenated system (with O₂ transfer, kO₂ > 0)
    3. Automatically generates product species and the reaction network.
    4. Solves the oxidation kinetics over 100 days at 25°C.
    5. Plots concentration profiles (log-log).
    6. Dumps the underlying reaction scheme.

Classes used:
    - `mixture`: container for species, physical conditions, and reactions
    - `mixtureKinetics`: solver for concentration evolution
    - `reactionRateDB`: fingerprint-based reaction rate parameter database

Note:
    The species prefix `C1` (e.g., `C1H`, `C1OOH`) is mapped to standard FAME patterns `L1`
    (e.g., `L1H`, `L1OOH`) to reuse known reaction rate constants.

@author: olivier.vitrac@gmail.com
@date: 2025-05-10
"""

# %% Setup (ensure radigen3 is in the current folder)
from radigen3.oxidation import mixture, mixtureKinetics, reactionRateDB

# Register class fingerprint substitution: C1 → L1 (use L1 kinetic constants)
reactionRateDB.register_fingerprint_substitution(r"C([123])", r"L\1")

# %% Oxidation in anoxia (no oxygen available)
mix1 = mixture(name="oxidation in anoxia", kO2=0)
mix1.add("C1H", concentration=3000)      # base FAME
mix1.add("C1OOH", concentration=300)     # initial hydroperoxide
mix1.addProducts()
mix1.addReactions()
mix1.populateReactionRates()

# Kinetic model (25°C for 100 days)
mix1kin = mixtureKinetics(mix1)
mix1kin.solve(tspan=(100, 'days'), T=25)

# Plot log-log scale
mix1kin.plot(xscale="log", yscale="log")

# Print reaction scheme used
print("Considered reaction Scheme", "-"*30, mix1.reactionScheme, sep="\n")

# %% Oxidation in the presence of oxygen (O₂ transport enabled)
mix2 = mix1.copy()
mix2.name = "oxidation with oxygenation"
mix2.kO2 = 1e-5                     # realistic oxygenation rate
mix2.add("O2", concentration=10)    # introduce oxygen
mix2.addProducts()
mix2.addReactions()
mix2.populateReactionRates()

# Kinetic model
mix2kin = mixtureKinetics(mix2)
mix2kin.solve(tspan=(100, 'days'), T=25) # simulate 100 days @25°C

# Plot log-log
mix2kin.plot(xscale="log", yscale="log")

# End of test
"""
Tested system:
  • 11 species, 10 reactions, 7 lumped groups.

Example reaction scheme:
    *R0: C1H + C1O•     → C1• + C1OH
    *R1: C1H + HO•      → C1• + H2O
    *R2: C1H + C1OO•    → C1• + C1OOH
    *R3: C1OOH          → C1O• + HO•
    *R4: C1OOH + C1OOH  → C1O• + C1OO•
    *R5: C1• + C1•      → C1-polymer + C1-polymer
    *R6: C1• + C1OO•    → C1-polymer + C1-polymer
    *R7: C1O•           → C1=O
    *R8: C1OO• + C1OO•  → C1-polymer + C1-polymer
"""
