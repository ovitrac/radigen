#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RADIGEN – simulate_oil.py

Oxidation modeling of fatty acid methyl esters (FAMEs) mixtures using the RADIGEN simulation kernel.

This script provides four didactic scenarios:
    1. **Constant-temperature oxidation** (1h at 140°C)
    2. **Extended low-temperature storage** (100h at 80°C after scenario 1)
    3. **Dynamic oxidation** using a temperature/oxygenation cycle (180→40°C)
    4. **Partial oil renewal**: simulation of abused oil partially replaced by fresh oil

🔬 **Chemical context**:
    • The species L1H, L2H, and L3H represent labile protons in mono-, di-, and triallylic positions:
        - L1H → mono-allylic (oleic acid equivalent, C18:1)
        - L2H → di-allylic (linoleic acid equivalent, C18:2)
        - L3H → tri-allylic (linolenic acid equivalent, C18:3)

    • We discard glycerol structure: triacylglycerols are modeled as equivalent methyl esters (FAMEs)

📘 **Core classes used**:
    - `mixture`: defines chemical species, concentrations, oxygen transport parameters, and reactions
    - `mixtureKinetics`: integrates oxidation kinetics over time with thermal and transport coupling
    - `TKO2cycle`: constructs dynamic oxidation profiles with time-dependent T(t) and kO₂(t)

📦 **Figures**:
    Figures are generated and saved to `"images/"` using:
        fig.print("filename", "images")

💡 **Example LLM prompts**:
    - "Model linoleic acid oxidation during heating and cooling"
    - "Simulate a FAME mixture under cycling conditions with varying oxygenation"
    - "Compare oxidized oil to a mixture renewed with fresh oil"

@author: olivier.vitrac@gmail.com
@date: 2025-05-11
"""

# %% Setup
from radigen3.oxidation import mixture, mixtureKinetics, TKO2cycle

# %% Define pristine oil mixture (unoxidized)
# Species roots:
#   L1: oleic (monoallylic, monounsaturated) fatty-ester
#   L2: linoleic (diallylic, disaturated) fatty-ester
#   L3: linolenic (triallylic, trisaturated) fatty-ester
#
# List of suffixes:
#    H = labile proton (in alpha position of a double bond)
#    • (alkxyl), O• (alkoxyl), OO• (hydroperoxil), HO• (peroxyl), •OO• (triplet dioxyden) radicals
#    OH (alcohol), HO (aldehyde), =O (ketone),-polymer (cross-linked polymer) stable products
#
# L1H is a pristine oleic ester, L1OOH is its hydroperoxide etc.

pristineoil = mixture(name="virgin oil", description="L1H+L2H+L3H", kO2=100)  # high O₂ bubbling
pristineoil.add("L1H", concentration=2500) # concentration in mol/m³
pristineoil.add("L2H", concentration=2500)
pristineoil.add("L3H", concentration=2500)
pristineoil.add("L1OOH", concentration=10)  # seed ROOH to initiate (here only L1OOH)
pristineoil.add("O2", concentration=5)      # initial oxygen content

# Autogenerate full chemical reaction network
pristineoil.addProducts()
pristineoil.addReactions()
pristineoil.populateReactionRates()

# Save reaction scheme
print("Considered reaction Scheme","-"*30,pristineoil.reactionScheme,sep="\n")
# %% Scenario 1 – Oxidation for 1h @ 140°C
oil1 = pristineoil.copy()
oil1.name = "Oxidized oil (1h @ 140°C)"
model1 = mixtureKinetics(oil1)
result1 = model1.solve(tspan=(1, "hour"), T=140)

# %% Plot all species and selected species
figscenario1 = model1.plot()
figscenario1.print("scenario1_all_species","images")
model1.plot(["L1H", "L2H", "L3H", "L1OOH", "L2OOH", "L3OOH", "O2"])

# %% Scenario 2 – Continue oxidation at 80°C for 100h
# Reuse oil1 composition to simulate long-term oxidation at 80°C
model1.set_C()  # update oil1 concentrations
# Continue oxidation
model2 = mixtureKinetics(oil1)  # reuse oil1
model2.name = "(100h @ 80°C)"
result2 = model2.solve(tspan=(100, "hours"), T=80)
model2.set_C()

# Merge and plot full trajectory
merged_model12 = model1 + model2
figscenario2 = merged_model12.plot()
figscenario2.print("scenario2_chained_oil","images")

# %% Scenario 3 – Dynamic oxidation cycle
# 1h @180°C (fast oxidation), ramp to 40°C in 4h, hold at 40°C for 2 days
cycle = TKO2cycle.constant(T=180, kO2=1.0, duration=(1, "hour")) +\
        TKO2cycle.ramp(T=(180, 40), kO2=(1.0, 0.01), duration=(4, "hours")) +\
        TKO2cycle.constant(T=40, kO2=1e-4, duration=(2, "days"))

oil3 = pristineoil.copy()
oil3.name = "Dynamic oxidation: 180→40°C"
model3 = mixtureKinetics(oil3)
result3 = model3.solve(tspan=(7, "days"), cycle=cycle)
figscenario3 = model3.plot(["L1H", "L2H", "L3H", "L1OOH", "O2"])
figscenario3.print("scenario3_cycles_oil","images")

# %% Scenario 4 – Partial oil renewal after abuse
# Step 1: abuse pristine oil 3h @160°C
oil_abused = pristineoil.copy()
oil_abused.kO2 = 1e-3 # [m/s] oxygen dissolution is limiting oxidation
oil_abused.name = "Oil abused for 3h @160°C"
model_abused = mixtureKinetics(oil_abused)
model_abused.solve(tspan=(3, "hours"), T=160)
model_abused.set_C()

# Step 2: mix with fresh oil (50/50 by volume), keeping total V unchanged
volume_original = pristineoil.V
oil_fresh = pristineoil.copy()
oil_fresh.V = volume_original / 2
oil_abused.V = volume_original / 2
oil_renewed = oil_abused + oil_fresh
oil_renewed.V = volume_original  # restore original volume
oil_renewed.A = pristineoil.A    # restore original interface area

# Step 3: simulate storage @40°C for 3 days
model_renewed = mixtureKinetics(oil_renewed)
model_renewed.name = f"({model_renewed.name}) for 3d @ 40°C"
model_renewed.solve(tspan=(3, "days"), T=40)
model_renewed.plot(["L1H", "L1OOH", "O2"])

merged_model_abused = model_abused + model_renewed
figscenario4 = merged_model_abused.plot(["L1H", "L2H", "L3H", "L1OOH", "L2OOH", "L3OOH", "O2"])
figscenario4.print("scenario4_abused_oil","images")

# end of examples

"""
    27 species and 60 reactions + 7 lumped species

    Details of the reaction scheme (generated and parameterized automatically):
        *R0: L1H + L1O• -> L1• + L1OH
        *R1: L1H + HO• -> L1• + H2O
        *R2: L1H + L1OO• -> L1• + L1OOH
        *R3: L1H + L2OO• -> L1• + L2OOH
        *R4: L1H + L3OO• -> L1• + L3OOH
        *R5: L1H + L2O• -> L1• + L2OH
        *R6: L1H + L3O• -> L1• + L3OH
        *R7: L2H + L1O• -> L2• + L1OH
        *R8: L2H + HO• -> L2• + H2O
        *R9: L2H + L1OO• -> L2• + L1OOH
        *R10: L2H + L2OO• -> L2• + L2OOH
        *R11: L2H + L3OO• -> L2• + L3OOH
        *R12: L2H + L2O• -> L2• + L2OH
        *R13: L2H + L3O• -> L2• + L3OH
        *R14: L3H + L1O• -> L3• + L1OH
        *R15: L3H + HO• -> L3• + H2O
        *R16: L3H + L1OO• -> L3• + L1OOH
        *R17: L3H + L2OO• -> L3• + L2OOH
        *R18: L3H + L3OO• -> L3• + L3OOH
        *R19: L3H + L2O• -> L3• + L2OH
        *R20: L3H + L3O• -> L3• + L3OH
        *R21: L1OOH -> L1O• + HO•
        *R22: L1OOH + L1OOH -> L1O• + L1OO•
        *R23: L1OOH + L2OOH -> L1O• + L2OO•
        *R24: L1OOH + L2OOH -> L2O• + L1OO•
        *R25: L1OOH + L3OOH -> L1O• + L3OO•
        *R26: L1OOH + L3OOH -> L3O• + L1OO•
        *R27: L1• + O2 -> L1OO•
        *R28: L1• + L1• -> L1-polymer + L1-polymer
        *R29: L1• + L2• -> L1-polymer + L2-polymer
        *R30: L1• + L3• -> L1-polymer + L3-polymer
        *R31: L1• + L1OO• -> L1-polymer + L1-polymer
        *R32: L1• + L2OO• -> L1-polymer + L2-polymer
        *R33: L1• + L3OO• -> L1-polymer + L3-polymer
        *R34: L2• + O2 -> L2OO•
        *R35: L2• + L2• -> L2-polymer + L2-polymer
        *R36: L2• + L3• -> L2-polymer + L3-polymer
        *R37: L2• + L1OO• -> L2-polymer + L1-polymer
        *R38: L2• + L2OO• -> L2-polymer + L2-polymer
        *R39: L2• + L3OO• -> L2-polymer + L3-polymer
        *R40: L3• + O2 -> L3OO•
        *R41: L3• + L3• -> L3-polymer + L3-polymer
        *R42: L3• + L1OO• -> L3-polymer + L1-polymer
        *R43: L3• + L2OO• -> L3-polymer + L2-polymer
        *R44: L3• + L3OO• -> L3-polymer + L3-polymer
        *R45: L1O• -> L1=O
        *R46: L1OO• + L1OO• -> L1-polymer + L1-polymer
        *R47: L1OO• + L2OO• -> L1-polymer + L2-polymer
        *R48: L1OO• + L3OO• -> L1-polymer + L3-polymer
        *R49: L2OO• + L2OO• -> L2-polymer + L2-polymer
        *R50: L2OO• + L3OO• -> L2-polymer + L3-polymer
        *R51: L3OO• + L3OO• -> L3-polymer + L3-polymer
        *R52: L2OOH -> L2O• + HO•
        *R53: L2OOH + L2OOH -> L2O• + L2OO•
        *R54: L2OOH + L3OOH -> L2O• + L3OO•
        *R55: L2OOH + L3OOH -> L3O• + L2OO•
        *R56: L3OOH -> L3O• + HO•
        *R57: L3OOH + L3OOH -> L3O• + L3OO•
        *R58: L2O• -> L2=O
        *R59: L3O• -> L3=O
"""
