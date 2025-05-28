#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module: radigen3.oxidation.py

Core module for radical oxidation modeling within the Generative Simulation Initiative.

This kernel supports the generative construction, parameterization, and simulation
of oxidation reaction networks in chemically diverse mixtures, such as edible oils,
FAMEs, fuels, and polymers. The oxidation model is built from chemical species annotated
with reactive functions and product maps. Radical propagation, decomposition, and
termination pathways are handled automatically, including temperature effects and
diffusion-limited corrections.

Features:
    - Auto-generation of all relevant reaction species and products from labile H donors.
    - Rule-based reaction network construction with automatic cross-reaction inference.
    - Temperature-dependent Arrhenius and Smoluchowski kinetics.
    - Support for hydroperoxide decomposition equilibrium (free vs. cage).
    - Stoichiometric matrix generation and fast vectorized ODE solving.
    - Oxygen dissolution modeled as a transport-limited source term.
    - Support for hundreds of reactions and species from text prompts.
    - Support arbitrary lumped species.

LLM compatibility:
    This module is designed to be a kernel for Large Language Models (LLMs) that can
    reason about and manipulate chemical mixtures and reaction systems from natural
    language. The full system is composable, traceable, and extensible.

Key References:
[2] Touffet M., Vitrac O. (2025), "Temperature-dependent kinetics of unsaturated fatty acid methyl esters:
    Modeling autoxidation mechanisms", Food Chemistry, 481, 143952. https://doi.org/10.1016/j.foodchem.2025.143952
[1] Touffet M., Smith P., Vitrac O. (2023), "A comprehensive two-scale model for predicting the oxidizability
    of fatty acid methyl ester mixtures", Food Research International, 173(1), 113289. https://doi.org/10.1016/j.foodres.2023.113289


### 🧩 Overview of Core Classes
+---------------------+---------------------------------------------------------------+
| Class               | Purpose                                                       |
+=====================+===============================================================+
| `species`           | Represents chemical species (radicals, peroxides, stable)     |
+---------------------+---------------------------------------------------------------+
| `reaction`          | Canonical reaction with fingerprint and type inference        |
+---------------------+---------------------------------------------------------------+
| `reactionRateDB`    | Registry of rate constants (k₀, Eₐ) searchable by fingerprint  |
+---------------------+---------------------------------------------------------------+
| `mixture`           | Encodes physical/chemical system: species, V, A, T, kO₂, pO₂  |
+---------------------+---------------------------------------------------------------+
| `mixtureKinetics`   | Numerical integrator for the chemical system (`solve_ivp`)    |
+---------------------+---------------------------------------------------------------+
| `lumped`            | Group of species for aggregate observables                    |
+---------------------+---------------------------------------------------------------+
| `TKO2cycle`         | Dynamic temperature and kO₂ definition for advanced scenarios |
+---------------------+---------------------------------------------------------------+


Author: Olivier Vitrac — olivier.vitrac@gmail.com
Revision: 2025-05-13
"""

# Revision history
# 2025-05-11 production release (v0.51)
# 2025-05-13 add reaction.ispolymerization (True if C-polymer, D-polymer)
#            Stoichiometry increment (+0.5 instead of 1) (v0.512)

# %% Indentication
__project__ = "Radigen"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2025"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@gmail.com"
__version__ = "0.512"


# %% Dependencies
import os, re, math, copy
import pandas as pd
import numpy as np
from types import SimpleNamespace
from datetime import datetime
from textwrap import indent
from contextlib import contextmanager
from collections import defaultdict
from scipy.sparse import dok_matrix, csr_matrix
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

__all__ = ['H2O', 'PrintableFigure', 'TKO2cycle', 'custom_plt_figure', 'custom_plt_subplots', 'diAllylicC', 'diAllylicCH', 'diAllylicCHO', 'diAllylicCO', 'diAllylicCOH', 'diAllylicCOO', 'diAllylicCOOH', 'diAllylicCeqO', 'is_valid_figure', 'kB1', 'kB2', 'kB3', 'ka1', 'ka2', 'ka3', 'kb1', 'kb2', 'kb3', 'kbb1', 'kbb2', 'kbb3', 'kd1', 'kd2', 'kd3', 'kdbi1', 'kdbi2', 'kdbi3', 'kp1', 'kp2', 'kp3', 'kt1', 'kt2', 'kt3', 'ktt1', 'ktt2', 'ktt3', 'kttt1', 'kttt2', 'kttt3', 'lumped', 'mixture', 'mixtureKinetics', 'monoAllylicC', 'monoAllylicCH', 'monoAllylicCHO', 'monoAllylicCO', 'monoAllylicCOH', 'monoAllylicCOO', 'monoAllylicCOOH', 'monoAllylicCeqO', 'oxygen', 'peroxyl', 'print_figure', 'print_pdf', 'print_png', 'reaction', 'reactionRateDB', 'species', 'terminationPolymers', 'triAllylicC', 'triAllylicCH', 'triAllylicCHO', 'triAllylicCO', 'triAllylicCOH', 'triAllylicCOO', 'triAllylicCOOH', 'triAllylicCeqO']

# %% Figure Functions
_fig_metadata_atrr_ = "__filename__"
def is_valid_figure(fig):
    """
    Checks if `fig` is a valid and open Matplotlib figure.

    Parameters:
    - fig: object to check

    Returns:
    - bool: True if `fig` is a valid, open Matplotlib figure.
    """
    return isinstance(fig, Figure) and hasattr(fig, 'canvas') and fig.canvas is not None

def _generate_figname(fig, extension):
    """
    Generate a clean filename based on metadata or current date/time.

    Parameters:
    - fig: Matplotlib figure object.
    - extension: File extension ('.pdf' or '.png').

    Returns:
    - str: Cleaned filename with correct extension.
    """
    # Try to retrieve the hidden filename metadata
    if hasattr(fig, _fig_metadata_atrr_):
        filename = getattr(fig, _fig_metadata_atrr_)
    else:
        # Default: Use date-time format if metadata is missing
        filename = "fig" + datetime.now().strftime("%Y%m%d_%H%M%S")
    # Clean filename (replace spaces, trim, remove special characters)
    filename = filename.strip().replace(" ", "_")
    # Ensure correct file extension
    if not filename.lower().endswith(extension):
        filename += extension
    return filename

def print_pdf(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a given figure as a PDF.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, PDF filename (auto-generated if empty).
    - destinationfolder: str, folder to save the file.
    - overwrite: bool, overwrite existing file.
    - dpi: int, resolution (default=300).
    """
    if not is_valid_figure(fig):
        print("no valid figure")
        return
    # Generate filename if not provided
    if not filename:
        filename = _generate_figname(fig, ".pdf")
    # add extension if missing
    if not filename.endswith(".pdf"):
        filename = filename+".pdf"
    # Ensure full path
    filename = os.path.join(destinationfolder, filename)
    # Prevent overwriting unless specified
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use overwrite=True to replace it.")
        return
    # Save figure as PDF
    fig.savefig(filename, format="pdf", dpi=dpi, bbox_inches="tight")
    print(f"Saved PDF: {filename}")

def print_png(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
    """
    Save a given figure as a PNG.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, PNG filename (auto-generated if empty).
    - destinationfolder: str, folder to save the file.
    - overwrite: bool, overwrite existing file.
    - dpi: int, resolution (default=300).
    """
    if not is_valid_figure(fig):
        print("no valid figure")
        return
    # Generate filename if not provided
    if not filename:
        filename = _generate_figname(fig, ".png")
    # add extension if missing
    if not filename.endswith(".png"):
        filename = filename+".png"
    # Ensure full path
    filename = os.path.join(destinationfolder, filename)
    # Prevent overwriting unless specified
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use overwrite=True to replace it.")
        return
    # Save figure as PNG
    fig.savefig(filename, format="png", dpi=dpi, bbox_inches="tight")
    print(f"Saved PNG: {filename}")

def print_figure(fig, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi={"png":150,"pdf":300}):
    """
    Save the figure in both PDF and PNG formats.

    Parameters:
    - fig: Matplotlib figure object to be saved.
    - filename: str, base filename (auto-generated if empty).
    - destinationfolder: str, folder to save the files.
    - overwrite: bool, overwrite existing files.
    - dpi: int, resolution (default=300).
    """
    if is_valid_figure(fig):
        print_pdf(fig, filename, destinationfolder, overwrite, dpi["pdf"])
        print_png(fig, filename, destinationfolder, overwrite, dpi["png"])
    else:
        print("no valid figure")

# Define PrintableFigure class
class PrintableFigure(Figure):
    """Custom Figure class with show and print methods."""

    def show(self, display_mode=None):
        """
        Show figure based on the environment:
        - In Jupyter: Uses display(fig)
        - In scripts or GUI mode: Uses plt.show()
        """
        try:
            get_ipython  # Check if running inside Jupyter
            if display_mode is None or display_mode == "auto":
                display(self)  # Preferred for Jupyter
            else:
                super().show()  # Use Matplotlib’s default show
        except NameError:
            super().show()  # Use Matplotlib’s default show in scripts

    def print(self, filename="", destinationfolder=os.getcwd(), overwrite=True, dpi={"png":150,"pdf":300}):
        print_figure(self, filename, destinationfolder, overwrite, dpi)

    def print_png(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_png(self, filename, destinationfolder, overwrite, dpi)

    def print_pdf(self, filename="", destinationfolder=os.getcwd(), overwrite=False, dpi=300):
        print_pdf(self, filename, destinationfolder, overwrite, dpi)

# ✅ Store original references
original_plt_figure = plt.figure
original_plt_subplots = plt.subplots

# ✅ Override `plt.figure()` to use PrintableFigure
def custom_plt_figure(*args, **kwargs):
    """Ensure all figures use PrintableFigure."""
    kwargs.setdefault("FigureClass", PrintableFigure)
    return original_plt_figure(*args, **kwargs)

# ✅ Override `plt.subplots()` to return PrintableFigure
def custom_plt_subplots(*args, **kwargs):
    """Ensure plt.subplots() returns a PrintableFigure."""
    fig, ax = original_plt_subplots(*args, **kwargs)  # Create normal figure
    fig.__class__ = PrintableFigure  # Explicitly convert to PrintableFigure
    return fig, ax

# ✅ Apply overrides
plt.figure = custom_plt_figure
plt.subplots = custom_plt_subplots
plt.FigureClass = PrintableFigure  # Ensure all figures default to PrintableFigure globally
plt.rcParams['figure.figsize'] = (8, 6)  # Optional: Default size

# %% Reaction Class
class reaction:
    """
    Represents a chemical reaction of the form A (+ B) → C (+ D).

    This class handles the logical structure and parameterization of elementary reactions in
    radical oxidation mechanisms. It supports Arrhenius and Smoluchowski rate models, diffusion-limited
    corrections, and fingerprinting for unique identification and cross-term inference.

    It distinguishes between monomolecular and bimolecular reactions, between self- and cross-reactions,
    and supports equilibrium treatment for specific classes (e.g. hydroperoxides).

    Attributes:
        A, B (species): Reactants; B can be None (monomolecular reaction).
        C, D (species): Products; D can be None.
        index (int): Reaction index within the mixture.
        k0 (float): Pre-exponential rate constant at T0.
        Ea (float): Activation energy [J/mol].
        T0 (float): Reference temperature [°C].
        kunits, Eaunits, Tunits (str): Units for k0, Ea, T0.
        diffusionLimited (bool): Whether reaction is diffusion-limited.
        ratio_rgh (float): Ratio of gyration to hydrodynamic radii.
        eta0, eta_Ea, eta_T0, etaunits: Viscosity model parameters.
        model (str): String describing the rate model origin.
        phi (float): Multiplicative factor for cross-reactions.
        inherit (tuple): Self-reactions from which a cross-reaction inherits its rate.

    Properties:
        ismonomolecular (bool): True if only A is defined.
        isbimolecular (bool): True if both A and B are defined.
        isselfroot (bool): True if A and B share the same root.
        iscrossroot (bool): True if A and B differ by root.
        isdecompositionROOH (bool): True if hydroperoxide decomposition.
        iscagedecompositionROOH (bool): True if bimolecular hydroperoxide cage decomposition.

    Methods:
        k(T): Return rate constant (includes cross-term inference).
        kArrhenius(T): Arrhenius-only rate constant.
        kSmoluchowski(T): Smoluchowski diffusion rate.
        viscosity(T): Temperature-dependent viscosity.
        fingerprint: Human-readable identifier (e.g. "A + B -> C + D").
        crossfingerprints: Tuple of fingerprints for each self-reaction of a cross-term.
        __str__, __repr__: Human-readable and detailed representations.

    Static Methods:
        hash_reagents(A, B): Hash based on reactants.
        hash_reagents_products(A, B, C, D): Hash including products.
        order_reagents(A, B, C, D): Canonical ordering of inputs.

    Notes:
        - Reactions are ordered to ensure uniqueness (e.g., oxidant = B if ranks differ).
        - Rate constants are inferred from registered `reactionRateDB` entries.
        - Cross-reactions use a geometric mean of inherited rate constants.

    See also:
        - species: for chemical participant definitions.
        - mixture: for reaction network generation and storage.
        - reactionRateDB: for assigning kinetic constants.
    """

    def __init__(self, *, A, C, B=None, D=None,
                 index=None, k0=1.0, Ea=45.0e3, T0=25.0,
                 Tunits = "°C", kunits = "m³/mol·s", Eaunits = "J/mol",
                 diffusionLimited = False, ratio_rgh = 0.77,
                 eta0 = 1.9e-3, eta_Ea = 17.0e3, eta_T0 = 80, etaunits = "Pa·s",
                 model="default", Tdefault=40, phi=2.0):
        """
        Initialize a reaction A (+ B) -> C (+ D)

        Parameters:
            A (species): Reactant A
            B (species or None): Reactant B (None for monomolecular)
            C (species): Product C
            D (species or None): Product D (optional)
            index (int or None): Reaction index
            k0 (float): Rate constant at T0 (units depend on reaction order)
            Ea (float): Activation energy (J/mol)
            T0 (float): Reference temperature (°C)
            model (str): Keyword for the reaction rate model
        """
        A,B,C,D = reaction.order_reagents(A, B, C, D)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        # index and inherited reactions
        self.index = index
        self.inherit = None # Returns both self-reactions matching a cross-reaction
        # Arrhenius terms
        self.k0 = k0
        self.kunits = kunits
        self.Ea = Ea
        self.Eaunits = Eaunits
        self.T0 = T0
        self.Tunits = Tunits
        # Smoluchowski terms
        self.diffusionLimited = diffusionLimited # False
        self.ratio_rgh = ratio_rgh # 0.77 # ratio between gyration (0.52 nm ± 0.08 nm) and hydrodynamic radii
        # Viscosity parameters
        self.eta0 = eta0      # 1.9e-3 # Pa·s (C18:1)
        self.eta_Ea = eta_Ea  # 17.0e3 # J/mol
        self.eta_T0 = eta_T0  # 80 # °C
        # other attributes
        self.model = model
        self.assigned = False
        self.Tdefault = Tdefault #40 #°C (default temperature if not defined)
        # cross terms
        self.phi = phi # cross-reaction rate mutiplicative constant respec­tively to their geometric mean

    @property
    def isselfroot(self):
        """Returns True if the reaction is bimolecular and involves same roots (e.g., L1OOH+L1OOH)"""
        return ((self.A.root == self.B.root) or self.B.root=="") if self.B is not None else True

    @property
    def iscrossroot(self):
        """Returns True if the reaction is monomolecular or involves different roots (e.g., L1OOH+L2OOH)"""
        return ((self.A.root != self.B.root) and self.B.root!="") if self.B is not None else False

    @property
    def ismonomolecular(self):
        """Returns True if the reaction is monomolecular"""
        return self.B is None

    @property
    def isbimolecular(self):
        """Returns True if the reaction is bimolecular"""
        return self.B is not None

    @property
    def isdecompositionROOH(self):
        """Returns True is the reaction involves the decomposition of hydroperoxides"""
        if self.ismonomolecular:
            return self.A.reactiveFunction == "COOH"
        else:
            return self.A.reactiveFunction == "COOH" and self.B.reactiveFunction == "COOH"

    @property
    def iscagedecompositionROOH(self):
        """Returns True if the reaction involves a bimolecular cage decomposition of hydroperoxides"""
        return self.isbimolecular and self.isdecompositionROOH

    @property
    def ispolymerization(self):
        """Returns True if the reaction involves a polymerization"""
        if self.ismonomolecular:
            return False
        else:
            return self.C.reactiveFunction == "-polymer" and self.D.reactiveFunction == "-polymer"

    def kArrhenius(self,T=None):
        """Returns the Arrhenian reaction rate"""
        T = self.Tdefault if T is None else T
        R = 8.314  # J/mol·K
        T_kelvin = 273.15 + T
        T0_kelvin = 273.15 + self.T0
        exponent = -self.Ea / R * (1/T_kelvin - 1/T0_kelvin)
        return self.k0 * math.exp(exponent)

    def kSmoluchowski(self,T=None):
        """Returns the Smoluchowski reaction rate"""
        T = self.Tdefault if T is None else T
        R = 8.314  # J/mol·K
        T_kelvin = 273.15 + T
        return 16.0/6.0 * self.ratio_rgh * R * T_kelvin / self.viscosity(T)

    def viscosity(self,T=None):
        """Returns the dynamic viscosity"""
        T = self.Tdefault if T is None else T
        R = 8.314  # J/mol·K
        T_kelvin = 273.15 + T
        T0_kelvin = 273.15 + self.eta_T0
        exponent = self.eta_Ea / R * (1/T_kelvin - 1/T0_kelvin)
        return self.eta0 * math.exp(exponent)

    def k(self,T=None):
        """
        Returns a function k(T) giving the Arrhenius rate constant at temperature T (°C).
        """
        T = self.Tdefault if T is None else T
        if self.inherit:
            return self.phi * math.sqrt( self.inherit[0].k(T) * self.inherit[1].k(T) )
        else:
            if self.diffusionLimited:
                return 1 / (1/self.kSmoluchowski(T) + 1/self.kArrhenius(T))
            else:
                return self.kArrhenius(T)

    @property
    def fingerprint(self):
        """Returns the human-readable fingerprint of a reaction"""
        parts = [self.A.shortname]
        if self.B:
            parts.append("+")
            parts.append(self.B.shortname)
        parts.append("->")
        parts.append(self.C.shortname)
        if self.D:
            parts.append("+")
            parts.append(self.D.shortname)
        return " ".join(parts)

    @property
    def crossfingerprints(self):
        """Returns the human-readable fingerprints for cross-reactions"""
        if not self.iscrossroot:
            return self.fingerprint
        A_rootA = self.A.root + self.A.suffix
        B_rootA = self.A.root + self.B.suffix
        C_rootA = self.A.root + self.C.suffix
        D_rootA = self.A.root + self.D.suffix if self.D else None
        A_rootB = self.B.root + self.A.suffix
        B_rootB = self.B.root + self.B.suffix
        C_rootB = self.B.root + self.C.suffix
        D_rootB = self.B.root + self.D.suffix if self.D else None
        parts_rootA = [A_rootA,"+",B_rootA,"->",C_rootA]
        parts_rootB = [A_rootB,"+",B_rootB,"->",C_rootB]
        if self.D:
            parts_rootA.append("+")
            parts_rootA.append(D_rootA)
            parts_rootB.append("+")
            parts_rootB.append(D_rootB)
        return (" ".join(parts_rootA)," ".join(parts_rootB))

    def __str__(self):
        """Shows reactions as strings, add * if the reaction is documented"""
        prefix = "*" if self.assigned else " "
        return prefix + f"{'R'+str(self.index) if self.index is not None else '?'}: {self.fingerprint}"


    def __repr__(self):
        """Show details of results"""
        k0_units = "s⁻¹" if self.ismonomolecular else "m³·mol⁻¹·s⁻¹"
        Ea_units = "J/mol"
        T0_units = "°C"
        props = {
            'A': self.A.name,
            'B': self.B.name if self.B else None,
            'C': self.C.name,
            'D': self.D.name if self.D else None,
            'index': self.index,
            'model': self.model,
            f'kA(T0={self.T0:.1f} [{T0_units}])': f"{self.kArrhenius(self.T0):.3g} [{k0_units}]",
            'Ea': f"{self.Ea:.3g} [{Ea_units}]"
            }
        if self.isbimolecular and self.diffusionLimited:
            props[f"kSL(T0={self.T0:.1f} [{T0_units}])"] = f"{self.kSmoluchowski(self.T0):.3g} [{k0_units}]"
        width = max(len(k) for k in props)
        print("\n".join(f"{k.rjust(width)}: {v}" for k, v in props.items()))
        return str(self)

    def __hash__(self):
        """
        Returns a hash based on the reactants A and B (unordered for bimolecular).
        """
        return self.hash_reagents(self.A, self.B)

    @staticmethod
    def hash_reagents(A, B=None):
        """
        Returns a hash based on reagent identities.

        Args:
            A (species): First reactant.
            B (species or None): Second reactant (optional).

        Returns:
            int: Unique hash for the reaction based on its reactants.
        """
        A,B,_,_ = reaction.order_reagents(A, B)
        if B is None:
            return hash(A.name)
        else:
            return hash(frozenset([A.name, B.name]))

    @staticmethod
    def hash_reagents_products(A, B=None, C=None, D=None):
        """
        Returns a hash based on reagent and product identities.
        """
        A,B,C,D = reaction.order_reagents(A, B, C, D)
        if B is None:
            return hash(A.name)
        elif D is None:
            return hash(frozenset([A.name, B.name, C.name]))
        else:
            return hash(frozenset([A.name, B.name, C.name, D.name]))

    @staticmethod
    def order_reagents(A, B, C=None, D=None):
        """
        Orders A + B -> C + D for reaction consistency (B is the oxidizer)
        Note that the products C and D must match A and B

        Rules:
        - If A is None and B is not, swap.
        - If both are defined and A._redoxRank > B._redoxRank, swap.

        Returns:
            tuple: (A, B) in canonical order
        """
        if A is None and B is not None:
            return B, None, D, C

        if B is not None:
            A_rank = getattr(A, "_redoxRank", float("inf"))
            B_rank = getattr(B, "_redoxRank", float("inf"))
            if A_rank and B_rank and A_rank > B_rank:
                return B, A, D, C

        return A, B, C, D


# %% Base Mixture and MixtureKinetics classes

# --------------------------------------------------------
# Base mixture class
# --------------------------------------------------------
class mixture:
    """
    Container class for chemical species and reactions involved in oxidation mechanisms.

    The `mixture` class represents a chemical system composed of individual `species` objects,
    along with their physical properties, reaction relationships, and dynamic behaviors.
    It supports prompt-driven population of reactive species, combinatorial generation of reactions,
    and parameterization from a curated reaction rate database.

    This class also encodes physical parameters relevant to kinetic modeling, including
    density, molar mass, volume, area, oxygen pressure, and diffusion properties.

    It is the entry point for model creation in the oxidation kernel, and is tightly coupled with
    the `mixtureKinetics` class, which provides numerical integration of the reaction network.

    Attributes:
        name (str): Name of the mixture.
        description (str): Optional description.
        A (float): Interfacial area for mass transfer [m²].
        V (float): Volume of the reacting domain [m³].
        pO2 (float): Partial pressure of oxygen [Pa].
        M (float): Molar mass of the mixture [kg/mol].
        rho0 (float): Reference density [kg/m³] at `rho_T0` °C.
        rho_T0 (float): Reference temperature for density [°C].
        rho_beta (float): Thermal expansion coefficient [1/K].
        kO2 (float): Oxygen mass transfer coefficient [m/s].

    Internals:
        _substances (list): List of all species in the mixture.
        _reactions (dict): Indexed dictionary of all reactions.
        _reaction_hashes (set): Set of reaction fingerprints for uniqueness.
        _reaction_counter (int): Tracks reaction indices for auto-numbering.

    Properties:
        listReactiveFunctions (list): lists unique reactive functions held by species in the mixture.
        listProductsBySubstance (dict): lists the products it can form for each substance
        listProducts (list): lists the unique names of all product species that can be formed
        listFormattedProductsBySubstance (dict): lists the *named* products it can form for each each substance
        listFormattedProducts (list): lists the unique *named* products that can be formed in the mixture
        species (list): lists species names
        substances (list(species)): lists species objects
        reactions (list(reactions)): lists reactions objects
        reactionScheme (str): Reaction scheme as a string

    Key Methods:
        add(clsname, name=..., **kwargs) → species
            Add a species to the mixture by class or alias name.
        addProducts() → int
            Recursively add all valid product species inferred from reactive pairs.
        addReactions() → int
            Combinatorially create all valid reactions among species in the mixture.
        populateReactionRates() → int
            Assign rate constants using registered reaction rate data.

        buildStoichiometryMatrix(sparse=False) → (S, species_list, reactions_list)
            Build the full stoichiometric matrix (dense or sparse).
        stoichiometryDataFrame(sparse=False) → pd.DataFrame
            Return a human-readable DataFrame version of the stoichiometry.

        getReactionByFingerprint(fingerprint: str) → reaction
            Look up a reaction using its string representation.
        getSpeciesByReactiveFunction(func: str) → list[species]
            Return species that match a given reactive function (e.g., "COOH").

        get_lumped_by_function(function: str) → lumped
            Return a lumped object for all species with the given reactive function.
        get_lumped_by_pattern(pattern: str, attr="name") → lumped
            Return a lumped object for species matching a regex pattern on a given attribute.

        lumped_hydroperoxides() → lumped
            Return a lumped species for all hydroperoxides (reactiveFunction = "COOH").
        lumped_alcohols() → lumped
            Return a lumped species for all alcohols (reactiveFunction = "COH").
        lumped_aldehydes() → lumped
            Return a lumped species for all aldehydes (reactiveFunction = "CHO").
        lumped_ketones() → lumped
            Return a lumped species for all ketones (reactiveFunction = "C=O").
        lumped_radicals_on_C() → lumped
            Return a lumped species for all radicals centered on carbon atoms.
        lumped_radicals_on_O() → lumped
            Return a lumped species for all radicals centered on oxygen atoms.
        lumped_polar_compounds() → lumped
            Return a lumped species for polar compounds (any compound with at least one oxygen).
        lumped_polymers() → lumped
            Return a lumped species for all termination polymers (reactiveFunction = "-polymer").

        copy() → mixture
            Return a copy of a mixture.

        mix1 + mix2 → mixture
            Operator '+' combine mixtures and dilute them accordingly


    Usage Example:
        oil = mixture()
        oil.add("L1H", concentration=3000)
        oil.add("L1OOH", concentration=100)
        oil.add("O2", concentration=10)
        oil.addProducts()
        oil.addReactions()
        oil.populateReactionRates()
        S, species_list, rxn_list = oil.buildStoichiometryMatrix()
        df = oil.stoichiometryDataFrame()

    Notes:
        - Products are named by combining the reagent root and the suffix of the product class.
        - Cross-reactions are inferred automatically but do not get their own parameter sets.
        - This class does not perform any kinetic simulation — see `mixtureKinetics` for that.
        - All parameters are expressed in SI units, unless otherwise noted.

    See also:
        - mixtureKinetics: For kinetic integration and modeling.
        - species: For base and derived chemical species.
        - reaction: For handling A (+B) → C (+D) logic.
        - reactionRateDB: To define and register kinetic data.

    Part of the Generative Simulation Initiative — olivier.vitrac@gmail.com
    """

    def __init__(self, name="mixture", description="",
                 A = 0.45, V = 3e-3, pO2 = 21e3, M = 0.28246,
                 rho0 = 883, rho_T0 = 25, rho_beta = 6.3e-4, kO2 = 1e-5,
                 Aunits = "m²", Vunits = "m³", pO2units = "Pa", Munits = "kg·m⁻¹",
                 rhounits = "kg·m⁻³", Tunits = "°C", betaunits = "K⁻¹", kO2units="m·s⁻¹"):
        """
        Attributes:
            name (str): The name of the mixture.
            description (str): Descriptive text.
            substances (list): List of species instances, each with an assigned index.

        Mixture container:
            Physical properties: default values (units)
                 A: 0.3*0.15 (m²), surface Area
                 V: 3e-3 (m³), volume
               pO2: 21e3 (Pa), O2 partial pressure in equilibrium with the mixture
              rho0: 883 (kg/m³)
            rho_T0: 25 (°C) reference temperature for density
          rho_beta: 6.3e-4 (K⁻¹), thermal expansion coefficient
        """

        self.name = name
        self.description = description
        # physical properties
        self.A = A
        self.V = V
        self.pO2 = pO2
        self.rho0 = rho0
        self.rho_T0 = rho_T0
        self.rho_beta = rho_beta
        self.Aunits = Aunits
        self.Vunits = Vunits
        self.pO2units = pO2units
        self.rhounits = rhounits
        self.Tunits = Tunits
        self.betaunits = betaunits
        self.M = M
        self.Munits = Munits
        self.kO2 = kO2
        self.kO2units = kO2units
        # substances
        self._substances = []
        self._indexmap = {}
        self._namemap = {}
        # reactions
        self._reactions = {}
        self._reaction_hashes = set()
        self._reaction_counter = 0

    def add(self, clsname, *, name=None, **kwargs):
        """
        Add a new species to the mixture using a class name or alias.

        Args:
            clsname (str): Alias or registered species class name.
            name (str, optional): Name to assign to the species instance.
                Defaults to the class's `defaultName`.
            **kwargs: Additional keyword arguments for the species constructor,
                      including e.g. concentration, root, shortname.

        Returns:
            species: The created and registered species object.
        """
        kwargs["_alias_used"] = clsname
        sp = species.create(clsname, **kwargs)

        # Override name if provided, otherwise use class default
        if name:
            sp.name = name
        elif getattr(sp, "name", "") == "":
            sp.name = getattr(sp, "defaultName", clsname)

        sp.index = len(self._substances)
        self._substances.append(sp)
        self._indexmap[sp.index] = sp
        self._namemap[sp.name] = sp
        setattr(self, sp.name, sp)
        self.refresh_indices()
        return sp

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._substances[key]
        elif isinstance(key, slice):
            return self._make_submixture(self._substances[key])
        elif isinstance(key, list):
            return self._make_submixture([self._substances[k] for k in key])
        elif isinstance(key, str):
            return self._namemap[key]
        raise KeyError(f"Invalid key: {key}")

    def __setitem__(self, key, value):
        if value == []:
            self._delete(key)
            return
        self._replace(key, value)

    def __setattr__(self, key, value):
        core_attrs = {"name", "description",
            "A","V","pO2","M","kO2","Aunits","Vunits","pO2units","Munits", "kO2units",
            "rho0","rho_T0","rho_beta","rhounits","Tunits","betaunits",
            "_substances", "_indexmap", "_namemap"}
        if key in core_attrs:
            super().__setattr__(key, value)
        elif isinstance(value, list) and value == []:
            self._delete(key)
        elif key in self._namemap:
            self._replace(key, value)
        else:
            super().__setattr__(key, value)

    def _delete(self, key):
        """delete a species in the mixture"""
        if isinstance(key, int):
            obj = self._substances[key]
        elif isinstance(key, str):
            obj = self._namemap.get(key, None)
        elif isinstance(key, slice):
            for k in range(*key.indices(len(self))):
                self._delete(k)
            return
        elif isinstance(key, list):
            for k in key:
                self._delete(k)
            return
        else:
            raise KeyError(f"Cannot delete key: {key}")

        if obj:
            self._substances.remove(obj)
            self._indexmap.pop(obj.index, None)
            self._namemap.pop(obj.name, None)
            if hasattr(self, obj.name):
                delattr(self, obj.name)
        self.refresh_indices()

    def _replace(self, key, newobj):
        """Replace a species in the mixture"""
        if not isinstance(newobj, species):
            raise TypeError("Replacement must be a species instance.")

        old = None
        idx = None

        if isinstance(key, int) and key < len(self._substances):
            old = self._substances[key]
            idx = key
        elif isinstance(key, str):
            old = self._namemap.get(key)
            idx = old.index if old else None

        if idx is not None and idx < len(self._substances):
            self._substances[idx] = newobj
        else:
            idx = len(self._substances)
            self._substances.append(newobj)

        newobj.index = idx
        self._indexmap[idx] = newobj
        self._namemap[newobj.name] = newobj
        super().__setattr__(newobj.name, newobj)
        self.refresh_indices()

    def _make_submixture(self, sublist):
        """Make a submixture from a mixture"""
        m = mixture(name=self.name + "_sub", description="submixture")
        for sp in sublist:
            m.add(sp.name, concentration=sp.concentration)
        return m

    def __iter__(self):
        return iter(self._substances)

    def __len__(self):
        return len(self._substances)

    def __repr__(self):
        width = max(len(sp.name) for sp in self._substances) if self._substances else 10
        lines = []
        for sp in self._substances:
            conc = f"{sp.concentration} [mol/m³]"
            react = f"{sp.reactiveFunction} + {sp.reactWith if len(sp.reactWith)>1 else sp.reactWith[0]} -> {sp.product if len(sp.product)>1 else sp.product[0]}" if sp.isreactive else "non-reactive"
            lines.append(f"[{sp.index}] {sp.name.rjust(width)}: {conc}")
            lines.append(f"{'R'.rjust(width+len(str(sp.index))+3)}: {react}")
        print("\n".join(lines))
        return(str(self))

    def __str__(self):
        return f"<mixture: {self.name}> including {len(self)} substances"

    @property
    def listReactiveFunctions(self):
        """
        Returns a sorted list of unique reactive functions held by species in the mixture.

        Returns:
            list[str]: List of unique `reactiveFunction` values.
        """
        funcs = {
            sp.reactiveFunction
            for sp in self
            if hasattr(sp, "reactiveFunction") and sp.reactiveFunction is not None
        }
        return sorted(funcs)

    @property
    def listProductsBySubstance(self):
        """
        For each species in the mixture, lists the products it can form
        (filtered by the reactive functions currently present in the mixture).

        Returns:
            dict[str, dict[str, str]]:
                {species_name: {reactWith: product}}
        """
        functions_present = self.listReactiveFunctions
        result = {}

        for sp in self:
            filtered = {
                func: prod
                for func, prod in sp.listProducts.items()
                if func == "self" or func in functions_present
            }
            if filtered:
                result[sp.name] = filtered

        return result

    @property
    def listProducts(self):
        """
        Lists the unique names of all product species that can be formed
        by the species in the mixture under current reactive conditions.

        Returns:
            list[str]: Sorted list of unique product class names.
        """
        functions_present = self.listReactiveFunctions
        product_set = set()

        for sp in self:
            for func, prod in sp.listProducts.items():
                if func != "self" and func not in functions_present:
                    continue
                if prod is None:
                    continue
                if isinstance(prod, (tuple, list, set, frozenset)):
                    product_set.update(prod)
                else:
                    product_set.add(prod)

        return sorted(product_set)

    @property
    def listFormattedProductsBySubstance(self):
        """
        For each species in the mixture, lists the *named* products it can form
        (filtered by the reactive functions currently present in the mixture),
        using root + class suffix logic.

        Returns:
            dict[str, dict[str, str]]:
                {species_name: {reactWith: formatted_product_name}}
        """
        functions_present = self.listReactiveFunctions
        result = {}

        for sp in self:
            filtered = {
                func: prod
                for func, prod in sp.formattedProducts.items()
                if func == "self" or func in functions_present
            }
            if filtered:
                result[sp.name] = filtered

        return result

    @property
    def listFormattedProducts(self):
        """
        Lists the unique *named* products that can be formed in the mixture
        using root + suffix logic (e.g., L1O• instead of monoAllylicCO).

        Returns:
            list[str]: Sorted list of formatted product names.
        """
        functions_present = self.listReactiveFunctions
        product_set = {
            prod
            for sp in self
            for func, prod in sp.formattedProducts.items()
            if (func == "self" or func in functions_present) and prod is not None
        }
        return sorted(product_set)

    def addProducts(self):
        """
        Recursively adds all potential product species (with zero concentration)
        derived from reactions under current reactive conditions.

        Each product species is created with:
            name = reagent.root + target_class.suffix

        Returns:
            int: Number of new species added to the mixture.
        """
        added = 0

        while True:
            known_species = {sp.name for sp in self}
            scheduled_names = set()
            new_to_add = []
            for reagent in self:
                for _, product_entry in reagent.listProducts.items():
                    if product_entry is None:
                        continue
                    # Determine iterable set of product names
                    if isinstance(product_entry, set):
                        product_names = list(product_entry)
                    elif isinstance(product_entry, tuple):
                        product_names = list(product_entry)
                    else:
                        product_names = [product_entry]

                    # Assign root to first in set only
                    for i, pname in enumerate(product_names):
                        product_cls = species._registry.get(pname)
                        if not product_cls:
                            continue
                        suffix = getattr(product_cls, "suffix", "")
                        hasroot = getattr(product_cls, "defaultRoot", "") != ""
                        root = reagent.root if not isinstance(product_entry, set) or hasroot else ""
                        full_name = root + suffix
                        if full_name not in known_species and full_name not in scheduled_names:
                            new_to_add.append((full_name, pname, root))
                            scheduled_names.add(full_name)
            if not new_to_add:
                break
            for new_name, product_cls_name, root in new_to_add:
                try:
                    self.add(product_cls_name, name=new_name, shortname=new_name, root=root, concentration=0.0)
                    added += 1
                except Exception as e:
                    print(f"⚠️ Could not add {product_cls_name} as {new_name}: {e}")
        self.refresh_indices()
        return added

    @property
    def species(self):
        """Return the list of substances in the mixture"""
        return [s.name for s in self._substances]

    def refresh_indices(self):
        """
        Refresh indices of compounds in the mixture.
        Ensures self._indexmap is contiguous and updates .index attributes.
        """
        # Ensure we are working with a list of compounds
        if not hasattr(self, '_substances'):
            raise AttributeError("Mixture must have a 'compounds' attribute (list of species)")

        self._indexmap = {}  # reset index map

        for i, compound in enumerate(self._substances):
            self._indexmap[compound.name] = i  # or compound.ID if preferred
            compound.index = i  # update the index attribute in-place

    def addReaction(self, A, B, C, D=None):
        A2, B2, C2, D2 = reaction.order_reagents(A, B, C, D)
        if B2 and D2 and A2.reactiveFunction==B2.reactiveFunction: # A and B are symmetric
            C2, D2, _, _ = reaction.order_reagents(C2, D2) # order C and D
        rhash = reaction.hash_reagents_products(A2,B2,C2,D2) # #rhash = reaction.hash_reagents(A, B)
        if rhash in self._reaction_hashes:
            return None
        r = reaction(A=A2, B=B2, C=C2, D=D2, index=self._reaction_counter)
        self._reactions[self._reaction_counter] = r
        self._reaction_hashes.add(rhash)
        self._reaction_counter += 1
        return r

    def addReactions(self, allow_self_reaction=True):
        """
        Generates all valid reactions based on current species.

        This method uses the rule:
        - C = A.listProducts[B.reactiveFunction]
        - D = B.listProducts[A.reactiveFunction]

        Products are looked up in the current mixture via their inferred names:
        root + suffix (from species registry).

        Returns:
            int: number of reactions added
        """
        cnt = 0

        for A in self:
            # Monomolecular reactions
            Cc = A.listProducts.get("self")
            if Cc:
                if isinstance(Cc, set):
                    Cc = list(Cc)
                    for i in range(len(Cc)):
                        for j in range(i + 1, len(Cc)):
                            clsC = species._registry.get(Cc[i])
                            clsD = species._registry.get(Cc[j])
                            if not clsC or not clsD:
                                continue
                            hasrootC = getattr(clsC, "defaultRoot", "") != ""
                            hasrootD = getattr(clsD, "defaultRoot", "") != ""
                            nameC = (A.root if hasrootC else "") + getattr(clsC, "suffix", "")
                            nameD = (A.root if hasrootD else "") + getattr(clsD, "suffix", "")
                            C = self._namemap.get(nameC)
                            D = self._namemap.get(nameD)
                            C,D,_,_ = reaction.order_reagents(C, D) # order is required ()
                            if C and D and self.addReaction(A, None, C, D):
                                cnt += 1
                else:
                    Cclasses = list(Cc) if isinstance(Cc, tuple) else [Cc]
                    for Cclass in Cclasses:
                        clsC = species._registry.get(Cclass)
                        if not clsC:
                            continue
                        nameC = A.root + getattr(clsC, "suffix", "")
                        C = self._namemap.get(nameC)
                        if C and self.addReaction(A, None, C):
                            cnt += 1

            # Bimolecular
            for B in self:
                if A is B and not allow_self_reaction:
                    continue

                fA = getattr(A, "reactiveFunction", None)
                fB = getattr(B, "reactiveFunction", None)
                if not fA and not fB:
                    continue

                Aprod = A.listProducts.get(fB)
                Bprod = B.listProducts.get(fA)

                Aps = list(Aprod) if isinstance(Aprod, tuple) else ([Aprod] if Aprod else [])
                Bps = list(Bprod) if isinstance(Bprod, tuple) else ([Bprod] if Bprod else [])

                if not Aps or not Bps:
                    continue

                # if the reactive groups are similar, the product lists is reversed
                if len(Aps) == len(Bps) and len(Aps) > 1 and A.reactiveFunction==B.reactiveFunction: #type(A) == type(B):
                    combos = list(zip(Aps, Bps[::-1]))
                else:
                    maxlen = max(len(Aps), len(Bps))
                    combos = list(zip(
                        Aps * ((maxlen + len(Aps) - 1) // len(Aps)),
                        Bps * ((maxlen + len(Bps) - 1) // len(Bps))
                    ))[:maxlen]

                for pc, pd in combos:
                    clsC = species._registry.get(pc)
                    clsD = species._registry.get(pd) if pd else None
                    if not clsC:
                        continue
                    nameC = A.root + getattr(clsC, "suffix", "")
                    nameD = B.root + getattr(clsD, "suffix", "") if clsD else None
                    C = self._namemap.get(nameC)
                    D = self._namemap.get(nameD) if nameD else None
                    if C and self.addReaction(A, B, C, D):
                        cnt += 1

        return cnt

    def populateReactionRates(self):
        """
        Updates reaction rate constants (k0, Ea, T0) from the reactionRateDB
        based on reaction fingerprints.

        If a matching fingerprint is found, the reaction's parameters are updated.
        If not, the default values are retained.

        Returns:
            int: Number of reactions updated.
        """
        updated = 0
        for rxn in self._reactions.values():
            if rxn.isselfroot: # same self-reaction
                fp = rxn.fingerprint
                entry = reactionRateDB.get_latest(fp)
                if entry:
                    rxn.k0 = entry.k0
                    rxn.Ea = entry.Ea
                    rxn.T0 = entry.T0
                    rxn.kunits = entry.kunits
                    rxn.Eaunits = entry.Eaunits
                    rxn.Tunits = entry.Tunits
                    rxn.diffusionLimited = entry.diffusionLimited
                    rxn.ratio_rgh = entry.ratio_rgh
                    rxn.eta0 = entry.eta0
                    rxn.eta_T0 = entry.eta_T0
                    rxn.eta_Ea = entry.eta_Ea
                    rxn.etaunits = entry.etaunits
                    rxn.phi = entry.phi # constant from cross-reaction terms
                    rxn.model = "values from "+entry.source
                    rxn.assigned = True
                    updated += 1
                else:
                    print("Missing reaction rate constant for:",rxn)
            else: # cross-reaction, we store the references to the self-reactions
                fps = rxn.crossfingerprints
                # last reaction rates data
                entry_rootA = reactionRateDB.get_latest(fps[0])
                entry_rootB = reactionRateDB.get_latest(fps[1])
                # matching self-reactions
                rxn_rootA = self.getReactionByFingerprint(fps[0])
                rxn_rootB = self.getReactionByFingerprint(fps[1])
                if entry_rootA and entry_rootB:
                    # minimum refresh of constants (in case they have not been refreshed)
                    rxn_rootA.k0 = entry_rootA.k0
                    rxn_rootB.k0 = entry_rootB.k0
                    rxn_rootA.T0 = entry_rootA.T0
                    rxn_rootB.T0 = entry_rootB.T0
                    rxn_rootA.Ea = entry_rootA.Ea
                    rxn_rootB.Ea = entry_rootB.Ea
                    # equivalent cross-reaction (required by method k())
                    rxn.inherit = (rxn_rootA,rxn_rootB) # inherit self-reactions
                    # equivalent parameters
                    rxn.phi = math.sqrt(entry_rootA.phi*entry_rootB.phi)
                    rxn.Ea = (entry_rootA.Ea * entry_rootB.Ea)/2
                    rxn.T0 = entry_rootA.T0
                    rxn.k0 = rxn.phi * math.sqrt(rxn_rootA.kArrhenius(rxn.T0) * rxn_rootB.kArrhenius(rxn.T0))
                    rxn.diffusionLimited = entry_rootA.diffusionLimited or entry_rootB.diffusionLimited
                    rxn.ratio_rgh = math.sqrt(entry_rootA.ratio_rgh * entry_rootB.ratio_rgh)
                    rxn.Tunits = entry_rootA.Tunits
                    rxn.Eaunits = entry_rootA.Eaunits
                    rxn.etaunits = entry_rootA.etaunits
                    rxn.eta0 = math.sqrt(entry_rootA.eta0 * entry_rootB.eta0)
                    rxn.model = "cross reactions"
                    rxn.assigned = True
                    updated += 1
        return updated

    @property
    def substances(self):
        """Retuns list of substances as a list"""
        return self._substances

    @property
    def reactions(self):
        """Returns defined reactions in the mixture as a list"""
        return list(self._reactions.values())

    @property
    def reactionScheme(self):
        """Returns reaction scheme as a string"""
        return "\n".join([str(r) for r in self.reactions])

    def getReactionByFingerprint(self, fingerprint):
        """
        Retrieves the reaction object matching a given fingerprint string.

        Args:
            fingerprint (str): A string representation of the reaction, e.g., "L1OOH -> L1HO• + L1O•"

        Returns:
            reaction or None: The matching reaction object, or None if not found.
        """
        for r in self._reactions.values():
            if r.fingerprint == fingerprint:
                return r
        return None

    def getSpeciesByReactiveFunction(self, reactiveFunction):
        """
        Retrieves a list of species matching a given reactiveFunction.

        Args:
            reactiveFunction (str): A string representation of the reactive function, e.g., "L1OOH"

        Returns:
            list of species: An empty list is returned if the reactive function is not found.
        """
        return [sp for sp in self._substances if sp.reactiveFunction==reactiveFunction]

    def buildStoichiometryMatrix(self, sparse=False):
        """
        Constructs the stoichiometric matrix S of the mixture's reaction network.

        Rows: species (ordered by .index)
        Columns: reactions (ordered by .index)

        Reactants A,B → -1 (increment, bimolecular → -2 with A=B)
        Products C,D → +1 (increment, except if polymerization → +0.5)

        Args:
            sparse (bool): If True, returns a scipy.sparse CSR matrix.

        Returns:
            S: Stoichiometry matrix (dense ndarray or sparse CSR)
            species_list (list): Ordered list of species
            reactions_list (list): Ordered list of reactions
        """
        species_list = sorted(self._substances, key=lambda sp: sp.index)
        reactions_list = sorted(self._reactions.values(), key=lambda r: r.index)

        n_species = len(species_list)
        n_reactions = len(reactions_list)

        if sparse:
            S = dok_matrix((n_species, n_reactions), dtype=np.float64)
        else:
            S = np.zeros((n_species, n_reactions), dtype=np.float64)

        for j, rxn in enumerate(reactions_list):
            for sp in [rxn.A, rxn.B]:
                if sp is not None:
                    S[sp.index, j] -= 1
            ispolymerization = rxn.ispolymerization
            for sp in [rxn.C, rxn.D]:
                if sp is not None:
                    S[sp.index, j] += 0.5 if ispolymerization else 1

        if sparse:
            S = S.tocsr()

        return S, species_list, reactions_list


    def stoichiometryDataFrame(self, sparse=False):
        """
        Converts the stoichiometry matrix into a human-readable DataFrame.

        Args:
            sparse (bool): If True, builds sparse stoichiometry matrix.

        Returns:
            pd.DataFrame: Species × Reactions stoichiometry.
        """
        S, species_list, reactions_list = self.buildStoichiometryMatrix(sparse=sparse)

        if sparse:
            S = S.toarray()

        row_labels = [sp.name for sp in species_list]
        col_labels = [rxn.fingerprint for rxn in reactions_list]

        return pd.DataFrame(S, index=row_labels, columns=col_labels)

    def get_lumped_by_function(self, function,color=None,linestyle='-'):
        """Returns a lumped object for all species with a given reactive function."""
        matches = [sp for sp in self if getattr(sp, "reactiveFunction", None) == function]
        return lumped(matches,color=color,linestyle=linestyle) if matches else None

    def get_lumped_by_pattern(self, pattern, attr="name",color=None,linestyle="-"):
        """Returns a lumped object for species whose attribute (default 'name') matches a regex."""
        regex = re.compile(pattern)
        matches = [sp for sp in self if regex.match(getattr(sp, attr, ""))]
        return lumped(matches,color=color,linestyle=linestyle) if matches else None

    def lumped_polymers(self,color="DarkSlateGray",linestyle="-"):
        """Returns a lumped species containing all polymer termination products."""
        return self.get_lumped_by_function("-polymer",color=color,linestyle=linestyle)

    def lumped_hydroperoxides(self,color="Magenta",linestyle="-"):
        """Returns a lumped species containing all hydroperoxides (COOH)."""
        return self.get_lumped_by_function("COOH",color=color,linestyle=linestyle)

    def lumped_alcohols(self,color="Indigo",linestyle="-"):
        """Returns a lumped species containing all aldehydes (COH)."""
        return self.get_lumped_by_function("COH",color=color)

    def lumped_aldehydes(self,color="DarkViolet",linestyle="-"):
        """Returns a lumped species containing all aldehydes (CHO)."""
        return self.get_lumped_by_function("CHO",color=color,linestyle=linestyle)

    def lumped_ketones(self,color="DarkSlateBlue",linestyle="-"):
        """Returns a lumped species containing all aldehydes (CHO)."""
        return self.get_lumped_by_function("C=O",color=color,linestyle=linestyle)

    def lumped_radicals_on_C(self,color="DarkGreen",linestyle=":"):
        """Returns a lumped species with radicals centered on carbon atoms."""
        matches = [sp for sp in self if sp.isradical_onC]
        return lumped(matches,color=color,linestyle=linestyle) if matches else None

    def lumped_radicals_on_O(self,color="Olive",linestyle="-."):
        """Returns a lumped species with radicals centered on oxygen atoms."""
        matches = [sp for sp in self if sp.isradical_onO]
        return lumped(matches,color=color,linestyle=linestyle) if matches else None

    def lumped_polar_compounds(self):
        """Returns a lumped species including compounds with at least one oxygen (based on reactiveFunction/suffix)."""
        matches = [
            sp for sp in self
            if 'O' in getattr(sp, "reactiveFunction", "") or 'O' in getattr(sp, "suffix", "")
        ]
        return lumped(matches) if matches else None

    def copy(self):
        """
        Returns a deep copy (clone) of the current mixture instance.

        This method ensures that all species, reactions, and internal structures
        are duplicated without shared references. Modifying the clone will not
        affect the original mixture.

        Returns:
            mixture: A fully independent clone of the original mixture.
        """
        return copy.deepcopy(self)

    def __add__(self, other):
        """
        Combine two mixture instances via physical mixing.

        This operator enables the addition of two compatible `mixture` objects. The result is a
        new mixture with:

            - The combined volume (V = V₁ + V₂)
            - Volume-weighted concentrations of each species
            - Physical properties inherited from the left operand (self), except volume
            - Same species names and classes (must match exactly)

        Constraints:
            - Species names in both mixtures must match exactly.
            - Reactions are not recomputed; only species and concentrations are merged.

        Returns:
            mixture: A new `mixture` instance representing the combination.

        Raises:
            TypeError: If `other` is not a `mixture`.
            ValueError: If the two mixtures have mismatched species.

        Example:
            oil1 = mixture()
            oil1.add("L1H", concentration=1000)
            oil2 = mixture()
            oil2.add("L1H", concentration=500)
            mixed = oil1 + oil2  # L1H final concentration is volume-weighted average
        """
        if not isinstance(other, mixture):
            raise TypeError("Can only mix with another mixture.")
        # Check that substance names match
        names1 = {sp.name for sp in self._substances}
        names2 = {sp.name for sp in other._substances}
        if names1 != names2:
            raise ValueError("Cannot mix incompatible mixtures (different species).")
        # Addition preserves volume
        total_V = self.V + other.V
        # left operand (self) has a higher precedence
        mix = self.copy()
        mix.name = f"{self.name}+{other.name}"
        mix.volume = total_V
        for sp in mix._substances:
            other_sp = next(s for s in other._substances if s.name == sp.name)
            new_conc = (sp.concentration * self.V + other_sp.concentration * other.V) / total_V
            sp.concentration = new_conc
        return mix


# --------------------------------------------------------
# Mixture representation for kinetic modeling
# --------------------------------------------------------
class mixtureKinetics:
    """
    Solver class for time-resolved oxidation kinetics in reactive mixtures.

    This class constructs and integrates the differential system for species concentrations:

        dC/dt = S · R(C, T) + source_terms(C, T)

    where:
        - C: vector of concentrations [mol/m³]
        - S: stoichiometric matrix
        - R: reaction rate vector
        - T: temperature (can be time-dependent)
        - source_terms includes mass transport (e.g., O2₂ influx)

    Key Features:
        - Handles monomolecular and bimolecular mechanisms
        - Automatically applies Arrhenius/Smoluchowski hybrid kinetics
        - Models ROOH ↔ RO• + HO• with cage/free equilibrium
        - Tracks oxygen flux using Henry's law and film model
        - Solves stiff ODE systems using scipy.integrate.solve_ivp
        - Supports arbitrary T(t) and kO₂(t) using a TKO2cycle object
        - Can merge multiple simulations via the '+' operator

    Lumped Species Support:
        Lumped groups of species can be registered via:
            model.register_lumped("LOOH", mix.lumped_hydroperoxides())
        Registered lumped species can be plotted, tracked, and exported as if they were named species.

    Plotting:
        - Colors and line styles are automatically detected from `_color` and `_linestyle` attributes.
        - Lumped species are rendered with thicker lines and a '†' symbol in legends.
        - Annotations show the curve maximum for each species if enabled.

    Dynamic Simulation Support:
        By passing a `TKO2cycle` object to `solve(cycle=...)`, the solver will use time-varying temperature and oxygen transfer coefficients.

    Methods:
        solve(tspan, T, C0, ..., cycle=None)
            Integrate the ODE system over time.
        plot(species=None, ...)
            Plot time-concentration profiles, including lumped groups.
        results_as_dataframe(species=None)
            Export results as a DataFrame.
        get_interpolator()
            Returns interpolator for concentrations at arbitrary time points.
        get_C0(), set_C(C=None)
            Retrieve or update species concentrations.
        register_lumped(name, lumped_obj)
            Register a lumped species group.
        get_lumped_concentration(name)
            Return total concentration over a lumped group.
        copy()
            Return a deep copy of the entire model and simulation state.
        __add__(other)
            Concatenate two compatible simulations for continuous kinetics.

    Static Methods:
        parse_time_unit(tspan, unit)
            Convert time input (e.g., (10, "days")) to seconds with canonical unit tracking.

    Attributes:
        name (str): Title of the model (e.g., "oxidation of oil")
        mixture (mixture): Associated mixture instance
        species_list (list): List of species in simulation
        reactions_list (list): List of reactions used
        n_species (int): Number of species
        n_reactions (int): Number of reactions
        _results (SimpleNamespace): Output from solve()
        _lumped_registry (dict[str, lumped]): Named groups for tracking output

    Example:
        oil = mixture()
        oil.add("L1H", concentration=3000)
        oil.add("L1OOH", concentration=100)
        oil.add("O2", concentration=10)
        oil.addProducts(); oil.addReactions(); oil.populateReactionRates()

        model = mixtureKinetics(oil)
        model.solve((5, 'days'), T=60)
        model.plot(["L1H", "hydroperoxides"])

        df = model.results_as_dataframe()
        print(df.head())
    """

    def __init__(self, mixture, name=None):
        """mixtureKinetics constructor"""
        self.name = mixture.name if name is None else name
        self.mixture = mixture
        self.species_list = sorted(mixture._substances, key=lambda sp: sp.index)
        self.reactions_list = sorted(mixture._reactions.values(), key=lambda r: r.index)
        self.n_species = len(self.species_list)
        self.n_reactions = len(self.reactions_list)
        # Precompute reactions indices involving hydroperoxide decomposition
        self._freeROOH_decomposition_indices = []
        self._cageROOH_decomposition_indices = []
        for j,rxn in enumerate(self.reactions_list):
            if rxn.isdecompositionROOH:
                if rxn.iscagedecompositionROOH:
                    self._cageROOH_decomposition_indices.append(j)
                else:
                    self._freeROOH_decomposition_indices.append(j)
        # Precompute hydroperoxides indices subjected to cage mechanism
        self._ROOH_indices = []
        self._ROOH_indices_byname = {}
        for i,sp in enumerate(self.species_list):
            if sp.reactiveFunction == "COOH":
                self._ROOH_indices.append(i)
                self._ROOH_indices_byname[sp.name] = i
        # Precompute source-regulated species indices (e.g. O2)
        self._source_indices = {}
        for i, sp in enumerate(self.species_list):
            if sp.name == "O2":
                self._source_indices["O2"] = i
        # Container for integration results
        self._solution = None
        self._results = None
        # Registry for lumped species
        self._lumped_registry = {}

    def register_lumped(self, name, lumped_obj):
        """Register a lumped group under a name for plotting/dataframe output."""
        if lumped_obj is None:
            print(f'WARNING: The lumped species "{name}" cannot be set, it is empty')
            return
        if not isinstance(lumped_obj, lumped):
            raise TypeError("Expected a lumped object.")
        self._lumped_registry[name] = lumped_obj

    def get_lumped_concentration(self, name):
        """Return the concentration (cumulated) of a registered lumped species."""
        if self._results is None:
            raise RuntimeError("No solution available. Run `solve()` first.")
        if name not in self._lumped_registry:
            raise ValueError(f"No lumped group named '{name}' registered.")
        indices = [self._results.species.index(sp.name) for sp in self._lumped_registry[name]]
        return self._results.t, self._results.y[indices].sum(axis=0)

    def get_C0(self):
        """Returns the initial concentration vector C0 from species in the mixture."""
        return np.array([sp.concentration for sp in self.species_list])

    def set_C(self, C=None):
        """
        Updates the concentrations of species in the mixture.

        If C is None, uses the last simulated concentration vector from the final time point.
        """
        if C is None:
            if self._results is None:
                raise ValueError("No simulation results available. Run `solve()` first.")
            C = self._results.y[:, -1]  # Last time point
        for sp, conc in zip(self.species_list, C):
            sp.concentration = conc

    def get_R(self, T, C):
        """
        Computes the vector of reaction rates R at temperature T and concentrations C.

        Args:
            T (float): Temperature in °C.
            C (array): Current concentrations of all species.

        Returns:
            R (np.ndarray): Reaction rate vector.
        """
        R = np.zeros(self.n_reactions, dtype=np.float64)
        for j, rxn in enumerate(self.reactions_list):
            k = rxn.k(T)
            A_idx = rxn.A.index
            if rxn.isbimolecular:
                B_idx = rxn.B.index
                if j in self._cageROOH_decomposition_indices: # cage mechanism subjected to eq
                    rate = k * self.cageROOH(T,C,C[A_idx]) * self.cageROOH(T,C,C[B_idx])
                else:
                    rate = k * C[A_idx] * C[B_idx] # other bimolecular reactions
            else:
                if j in self._freeROOH_decomposition_indices: # freeROOH subjected to eq
                    rate = k * self.freeROOH(T,C,C[A_idx])
                else:
                    rate = k * C[A_idx] # other monomolecular reactions
            R[j] = rate
        return R

    def totalROOH(self,C):
        """Returns the total concentration in hydroperoxides"""
        return sum(C[self._ROOH_indices])

    def KcagefreeROOH(self,T,C):
        """Returns the equilibrium constant between cage and free hydroperoxides at T"""
        TCK = 80.0 + 273.15 # K
        TK = T + 273.15     # K
        DeltaH = 60.0e3     # kJ/mol
        R = 8.314           # J/mol·K
        KTC = 2/self.totalROOH(C) # in mixture we consider an equilibrium between all species
        exponent = DeltaH/R * (1/TK - 1/TCK)
        return KTC * np.exp(exponent)

    def freeROOH(self,T,C,CROOH):
        """Returns the free conentration of hydroperoxides of type COOH at T"""
        K = self.KcagefreeROOH(T,C)
        return CROOH + (1 - np.sqrt(1+4*K*CROOH))/(2*K)

    def cageROOH(self,T,C,CROOH):
        """Returns the free conentration of hydroperoxides of type COOH at T"""
        K = self.KcagefreeROOH(T,C)
        return (np.sqrt(1+4*K*CROOH)-1)/(2*K)

    def oxygen_solubility(self,T):
        """
        Returns oxygen solubility in mol/(m³·Pa) using Arai et al. (1989) constants.
        Valid for FAMEs over a wide T range.

        Args:
            T (float): Temperature in °C.

        Returns:
            float: Solubility S_O2 in mol/(m³·Pa)
        """
        # Arai parameters
        B1, B2, B3, B4 = 7080, 56.603, -0.11064, -309.62
        TK = 273.15 + T  # Kelvin
        exponent = B1 / TK + B2 * np.log(TK) + B3 * TK + B4
        return 1/(1e6 * self.molarVolume(T) * math.exp(exponent))

    def rho(self,T):
        """
        Temperature-dependent density model for mixture.
        Args:
            T (float): Temperature in °C.
        Returns:
            float: density in kg/m³
        """
        exponent = -self.mixture.rho_beta * (T - self.mixture.rho_T0)
        return self.mixture.rho0 * math.exp(exponent)


    def molarVolume(self,T):
        """
        Estimate molar volume (m³/mol) of the mixture
        Args:
            T (float): Temperature in °C.

        Returns:
            float: Molar volume in m³/mol
        """
        return (self.mixture.M / self.rho(T)) # m³/mol

    def sourceO2(self,T,concO2):
        """
        Oxygen transport model
        """
        concO2eq = self.oxygen_solubility(T)*self.mixture.pO2
        keffective = self.mixture.kO2 * self.mixture.A / self.mixture.V
        return keffective * (concO2eq - concO2)


    def get_dCdt(self, T, C):
        """
        Computes the time derivative of species concentrations:

            dC/dt = S · R(T, C) + source_vector

        Includes source term for O2 transport.
        """
        S, _, _ = self.mixture.buildStoichiometryMatrix(sparse=False)
        R = self.get_R(T, C)
        source = np.zeros(self.n_species, dtype=np.float64)
        if "O2" in self._source_indices:
            i = self._source_indices["O2"]
            source[i] = self.sourceO2(T, C[i])
        return S @ R + source

    @contextmanager
    def temporary_kO2(self, new_kO2):
        """Temporarily override current kO2 value"""
        original_kO2 = self.mixture.kO2
        self.mixture.kO2 = new_kO2
        try:
            yield
        finally:
            self.mixture.kO2 = original_kO2

    def get_dCdt_dynamic(self, Tfunc, kO2func):
        """
        Computes the time derivative of species concentrations under dynamic T(t) and  kO2(t) conditions

            dC/dt = S · R(T(t), C) + source_vector(t)
        """
        def rhs(t, C):
            T = Tfunc(t)
            kO2 = kO2func(t)
            with self.temporary_kO2(kO2):
                S, _, _ = self.mixture.buildStoichiometryMatrix(sparse=False)
                R = self.get_R(T, C)
                source = np.zeros(self.n_species)
                if "O2" in self._source_indices:
                    i = self._source_indices["O2"]
                    source[i] = self.sourceO2(T, C[i])
                return S @ R + source
        return rhs

    def solve(self, tspan, T=None, C0=None, species=None, cycle=None, method="BDF", rtol=1e-6, atol=1e-7, max_order=5, **kwargs):
        """
        Integrate the system of ODEs dC/dt = S·R(C,T) with optional oxygen source and dynamic conditions.

        Supports both constant and time-dependent temperature and oxygenation profiles. When a `TKO2cycle`
        object is passed to `cycle`, the simulation will use the temperature and kO₂ values defined by the cycle
        for each time point.

        Args:
            tspan (float or tuple): Final time `tf` or `(t0, tf)` in seconds.
                                    Can also be passed as (value, unit) tuples (e.g., (10, "days")).
            T (float, optional): Constant temperature [°C] if no cycle is specified.
            C0 (np.ndarray, optional): Initial concentration vector. If None, uses current species concentrations.
            species (list[str], optional): Species to include in the result. Defaults to all species.
            cycle (TKO2cycle, optional): Time-dependent temperature and oxygenation cycle object.
            method (str): Integration method (e.g., 'BDF', 'LSODA').
                          'BDF' recommended if the concentrations in radicals are initially 0
            rtol (float): Relative tolerance.
            atol (float): Absolute tolerance.
            max_order (int): Maximum solver order (for implicit methods).
            **kwargs: Additional arguments passed to `scipy.integrate.solve_ivp`.

        Returns:
            types.SimpleNamespace:
                - t (np.ndarray): Time points [unit-scaled].
                - y (np.ndarray): Concentrations at each time step.
                - tunit (str): Original time unit (e.g., "days").
                - tscale (float): Scale factor applied to time for output.
                - concunits (str): Concentration units ("mol/m³").
                - species (list[str]): Names of tracked species.
                - success (bool): Whether the integration succeeded.
                - message (str): Solver status message.
                - solver (str): Solver method used.
                - cycle (TKO2cycle or None): Applied cycle, if any.
                - isCycle (bool): True if dynamic conditions were used.
        """
        # Allow single time value
        tunit = "s"
        if isinstance(tspan, (int, float)):
            tfinal,tunit,tscale = mixtureKinetics.parse_time_unit(tspan,tunit)
            tinit = 0.0
        elif isinstance(tspan,tuple):
            if len(tspan)!=2:
                raise ValueError('tspan must be (tinitial,tfinal) or (tfinal,"time unit")')
            if isinstance(tspan[1],str):
                tfinal,tunit,tscale = mixtureKinetics.parse_time_unit(tspan,tunit)
                tinit = 0.0
            else:
                tinit,_,_ = mixtureKinetics.parse_time_unit(tspan[1],tunit)
                tfinal,tunit,tscale = mixtureKinetics.parse_time_unit(tspan[1],tunit)
        tspan = (float(tinit), float(tfinal))

        # Default temperature and initial conditions (no cycle)
        if T is None:
            T = self.mixture.T if hasattr(self,'T') else None

        if C0 is None:
            C0 = self.get_C0()

        # RHSmodel Cycle conditions (if any)
        if cycle is not None:
            if not isinstance(cycle,TKO2cycle):
                raise ValueError(f"cyle must be a TKO2cycle not a {type(cycle).__name__}")
            isCycle = True
            # Cycle RHS Model
            rhs = self.get_dCdt_dynamic(cycle.T, cycle.kO2)
        else:
            isCycle = False
            # Default RHS Model
            def rhs(t, C): return self.get_dCdt(T, C)

        # Options (use BDF instead of LSODA if radicals are 0 at t=0)
        default_opts = dict(method=method, rtol=rtol, atol=atol,
                            max_order=max_order,max_step=np.inf)
        default_opts.update(kwargs)

        # integration
        sol = solve_ivp(rhs, tspan, C0, **default_opts)

        # Store solution
        self._solution = sol

        # Select species
        tracked_species = species or [sp.name for sp in self.species_list]
        tracked_indices = [sp.index for sp in self.species_list if sp.name in tracked_species]
        tracked_conc = sol.y[tracked_indices, :]

        # Save and return result
        result = SimpleNamespace(
            t=sol.t,
            tunit = tunit,
            tscale = tscale,
            y=tracked_conc,
            concunits = "mol/m³",
            species=tracked_species,
            species_names = tracked_species,
            success=sol.success,
            message=sol.message,
            solver=method,
            cycle=cycle,
            isCycle=isCycle
        )
        self._results = result  # optional container

        # Register default lumped species
        self.register_lumped("hydroperoxides", self.mixture.lumped_hydroperoxides())
        self.register_lumped("aldehydes", self.mixture.lumped_aldehydes())
        self.register_lumped("ketones", self.mixture.lumped_ketones())
        self.register_lumped("alcohols", self.mixture.lumped_alcohols())
        self.register_lumped("polar", self.mixture.lumped_polar_compounds())
        self.register_lumped("radicals_on_C", self.mixture.lumped_radicals_on_C())
        self.register_lumped("radicals_on_O", self.mixture.lumped_radicals_on_C())
        self.register_lumped("Term-polymers", self.mixture.lumped_polymers())

        return result

    def get_interpolator(self):
        """
        Returns an interpolator object to estimate concentrations at arbitrary times.

        Returns:
            interp (callable): Function interp(t) that returns concentrations at time t (in seconds).
        """
        if self._results is None:
            raise ValueError("No simulation results available. Run `solve()` first.")
        return interp1d(self._results.t, self._results.y, kind='linear', axis=-1, bounds_error=False, fill_value="extrapolate")


    def plot(self, species=None, ax=None, figsize=None, ncol=None,
             legend_loc="center left", bbox_to_anchor=(1.0, 0.5),
             annotate_max=True, xscale="linear", yscale="linear", **kwargs):
        """
        Plot the concentration profiles of selected species over time.

        Species are rendered with custom colors if their `_color` attribute is set.
        Lumped species (`_lumped = True`) are plotted with thicker lines and shown with a † symbol in the legend.
        The maximum value of each curve can be labeled for clarity.

        Args:
            species (list[str], optional): List of species or lumped names to plot.
                                           If None, all species and registered lumped groups are shown.
            ax (matplotlib.axes.Axes, optional): Axes to draw on. If None, creates a new figure.
            figsize (tuple, optional): Size of the figure.
            ncol (int, optional): Number of columns in the legend.
            legend_loc (str): Location of the legend.
            bbox_to_anchor (tuple): Anchor position for the legend.
            annotate_max (bool): Whether to annotate the curve maximum (default: True).
            xscale (str): Scale for x-axis ('linear' or 'log'). Default: 'linear'.
            yscale (str): Scale for y-axis ('linear' or 'log'). Default: 'linear'.

        Raises:
            RuntimeError: If no solution is available.
            ValueError: If a species or lumped group is unrecognized.
        """
        if not hasattr(self, "_results") or self._results is None:
            raise RuntimeError("No solution available. Run `solve()` first.")

        sol = self._results
        t = sol.t * sol.tscale
        y = sol.y

        if species is None:
            selected_species = sol.species + list(self._lumped_registry.keys())
        else:
            selected_species = species

        resolved_names = []
        y_sel_list = []
        line_styles = []
        colors = []
        labels = []

        for name in selected_species:
            if name in self._lumped_registry:
                lump_obj = self._lumped_registry[name]
                spnames = [sp.name for sp in lump_obj.species]
                try:
                    lump_indices = [sol.species.index(n) for n in spnames]
                    y_sel_list.append(y[lump_indices].sum(axis=0))
                except ValueError as e:
                    raise ValueError(f"One or more species in lumped group '{name}' not in solution.") from e
                resolved_names.append(name)
                colors.append(getattr(lump_obj, "_color", None))
                linestyle = getattr(lump_obj, "_linestyle", "-")
                line_styles.append({'linewidth': 3.5, 'linestyle': linestyle})
                labels.append(f"{name}")
            elif name in sol.species:
                sp_idx = sol.species.index(name)
                y_sel_list.append(y[sp_idx])
                sp = next(sp for sp in self.species_list if sp.name == name)
                resolved_names.append(name)
                colors.append(getattr(sp, "_color", None))
                linestyle = getattr(sp, "_linestyle", "-")
                line_styles.append({'linewidth': 2.0, 'linestyle': linestyle})
                labels.append(name)
            else:
                raise ValueError(f"Unknown species or lumped group: '{name}'")

        y_sel = np.vstack(y_sel_list)
        n_species = len(resolved_names)

        if figsize is None:
            base_width = 8
            figsize = (base_width + 0.5 * n_species, 16)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        # Plot curves and annotate
        for i, name in enumerate(resolved_names):
            plot_kwargs = {
                "label": labels[i],
                "color": colors[i],
                **line_styles[i],
                **kwargs
            }
            line, = ax.plot(t, y_sel[i], **plot_kwargs)

        if annotate_max:
            x1 = ax.get_xlim()[1]
            x0 = ax.get_xlim()[0]
            y1 = ax.get_ylim()[1]
            y0 = ax.get_ylim()[0]

            yspan = y1 - y0
            for i, name in enumerate(resolved_names):
                y_i = y_sel[i]
                t_i = t
                ymax = y_i.max()
                ymin = y_i.min()
                i_max = y_i.argmax()
                if ((ymax-ymin) > (0.02 * yspan)) or yscale=="log":
                    x = max(x0,t_i[i_max])
                    y = max(y0,y_i[i_max])
                    label = labels[i]

                    if i_max == 0:
                        ax.annotate(label, xy=(x, y), xytext=(x + 0.02 * t[-1], y),
                                    textcoords='data', arrowprops=dict(arrowstyle='->', lw=0.5),
                                    fontsize="medium", fontweight="bold", color = colors[i], ha="left", va="center")
                    elif i_max == len(t_i) - 1:
                        ax.annotate(label, xy=(x, y), xytext=(x - 0.02 * t[-1], y),
                                    textcoords='data', arrowprops=dict(arrowstyle='->', lw=0.5),
                                    fontsize="medium", fontweight="bold", color = colors[i], ha="right", va="center")
                    else:
                        ax.text(x, y, label, fontsize="medium", fontweight="bold", color = colors[i], ha="center", va="bottom")

        ax.set_xlabel(f"Time [{sol.tunit}]", fontsize="large")
        ax.set_ylabel(f"Concentration [{sol.concunits}]", fontsize="large")
        ax.set_title(self.name, fontsize="large", fontweight="bold")
        ax.grid(True)

        if ncol is None:
            ncol = math.ceil(n_species / 16)

        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(
            handles,
            labels,
            loc=legend_loc,
            bbox_to_anchor=bbox_to_anchor,
            ncol=ncol,
            fontsize="medium",
            frameon=False
        )

        for text in legend.get_texts():
            if text.get_text() in self._lumped_registry:
                text.set_text(text.get_text() + " †")
                text.set_fontweight("bold")

        plt.tight_layout(rect=[0, 0, 0.85, 1.0])
        plt.show()

        return fig


    def results_as_dataframe(self, species=None):
        """
        Return the time evolution of species as a pandas DataFrame.

        Args:
            species (list[str], optional): If provided, limits the output to these species.

        Returns:
            pd.DataFrame: DataFrame with time and concentrations.
        """
        if not hasattr(self, "_results") or self._results is None:
            raise RuntimeError("No solution available. Run `solve()` first.")

        sol = self._results
        species_all = sol.species
        selected_species = species or species_all
        resolved_names = []
        resolved_data = {}
        for name in selected_species:
            if name in self._lumped_registry:
                spnames = [sp.name for sp in self._lumped_registry[name].species]
                try:
                    lump_indices = [species_all.index(n) for n in spnames]
                    resolved_data[name] = sol.y[lump_indices].sum(axis=0)
                except ValueError as e:
                    raise ValueError(f"One or more species in lumped group '{name}' not in solution.") from e
            elif name in species_all:
                resolved_data[name] = sol.y[species_all.index(name)]
            else:
                raise ValueError(f"Species not in solution or lumped group unknown: '{name}'")
            resolved_names.append(name)

        df = pd.DataFrame(resolved_data)
        df.insert(0, "time [s]", sol.t)

        return df

    def copy(self):
        """
        Returns a deep copy (clone) of the current mixture instance.

        This method ensures that all species, reactions, and internal structures
        are duplicated without shared references. Modifying the clone will not
        affect the original mixture.

        Returns:
            mixture: A fully independent clone of the original mixture.
        """
        return copy.deepcopy(self)

    @staticmethod
    def parse_time_unit(tspan, unit='s'):
        """
        Convert a time span with units to seconds, return also the canonical unit and inverse scale.

        Args:
            tspan (float or tuple): time duration (e.g., 3600 or (10.0, "d"))
            unit (str, optional): time unit if not using tuple form.

        Returns:
            time_in_seconds (float): Time span converted to seconds.
            canonical_unit (str): Canonical form ('s', 'min', 'hours', 'days', 'weeks', 'months', 'years').
            inverse_scale (float): Scaling factor to convert seconds back to canonical unit.

        Raises:
            ValueError: If the unit is unknown or invalid.

        Canonical units and definitions:
            - 's': seconds
            - 'min': 60 seconds
            - 'hours': 3600 seconds
            - 'days': 86400 seconds
            - 'weeks': 7 × 86400
            - 'months': 30 × 86400
            - 'years': 365 × 86400
        """
        unit_map = {
        's': ('s', 1),'sec': ('s', 1),'second': ('s', 1),'seconds': ('s', 1),
        'min': ('min', 60),'mn': ('min', 60),'minutes': ('min', 60),
        'h': ('hours', 3600), 'hr': ('hours', 3600), 'hour': ('hours', 3600), 'hours': ('hours', 3600),
        'd': ('days', 86400), 'day': ('days', 86400),'days': ('days', 86400),
        'w': ('weeks', 7 * 86400),'week': ('weeks', 7 * 86400),'weeks': ('weeks', 7 * 86400),
        'm': ('months', 30.44 * 86400), 'mo': ('months', 30.44 * 86400),'month': ('months', 30.44 * 86400),          'months': ('months', 30.44 * 86400),
        'y': ('years', 365 * 86400),'yr': ('years', 365 * 86400),'year': ('years', 365 * 86400),'years': ('years', 365 * 86400),}
        # Accept tuple format: (value, unit)
        if isinstance(tspan, tuple):
            if len(tspan) != 2:
                raise ValueError("Tuple input for tspan must be (value, unit).")
            tspan, unit = tspan
        if unit is None:
            raise ValueError("A time unit must be provided, either in a tuple or as a separate argument.")
        unit = unit.lower().strip()
        if unit not in unit_map:
            raise ValueError(f"Unknown or unsupported time unit '{unit}'.")
        canonical_unit, seconds_per_unit = unit_map[unit]
        time_in_seconds = float(tspan) * seconds_per_unit
        inverse_scale = 1.0 / seconds_per_unit
        return time_in_seconds, canonical_unit, inverse_scale


    def __repr__(self):
        """repr() method"""
        header = "[mixtureKinetics]\n"
        species_info = f"• {self.n_species} species\n"
        reactions_info = f"• {self.n_reactions} reactions\n"
        sim_info = "• Simulation: "
        if self._results:
            sim_info += f"✓ completed (t = 0 to {self._results.t[-1]:.1f} {self._results.tunit})"
        else:
            sim_info += "× not run yet"
        tracked = self._results.species if self._results else [sp.name for sp in self.species_list]
        tracked_str = ", ".join(tracked[:6]) + ("..." if len(tracked) > 6 else "")
        lumped_keys = list(self._lumped_registry.keys())
        lumped_str = ", ".join(lumped_keys[:6]) + ("..." if len(lumped_keys) > 6 else "")
        lumped_info = f"• Lumped groups: {lumped_str or 'none'}"
        print(
            header
            + indent(species_info + reactions_info + sim_info + "\n", "  ")
            + indent(f"Tracked species: {tracked_str}\n", "  ")
            + indent(lumped_info, "  ")
        )
        return str(self)

    def __str__(self):
        """str() method"""
        n_species = self.n_species
        n_reactions = self.n_reactions
        sim_status = "available" if self._results else "not yet run"
        return (
            f"<mixtureKinetics with {n_species} species and {n_reactions} reactions "
            f"— simulation {sim_status}>"
        )

    def __add__(self, other):
        """Concatenate two mixtureKinetics simulations if compatible."""
        if not isinstance(other, mixtureKinetics):
            return NotImplemented

        # Ensure simulations are solved
        if self._results is None or other._results is None:
            raise ValueError("Both mixtureKinetics instances must have completed simulations.")

        # Check that species and reactions match
        if self.n_species != other.n_species:
            raise ValueError("Species count mismatch.")
        if self.n_reactions != other.n_reactions:
            raise ValueError("Reactions count mismatch.")

        # Check species names in same order
        if [sp.name for sp in self.species_list] != [sp.name for sp in other.species_list]:
            raise ValueError("Species do not match in order or identity.")
        if [r.fingerprint for r in self.reactions_list] != [r.fingerprint for r in other.reactions_list]:
            raise ValueError("Reactions do not match in order or identity.")

        # Extract results
        r1, r2 = self._results, other._results

        # Determine time shift (0 if already continuous)
        dt = r1.t[-1]
        t2_shifted = r2.t + (0.0 if np.isclose(r2.t[0], dt, atol=1e-9) else dt)

        # Concatenate time and concentrations
        t_combined = np.concatenate((r1.t, t2_shifted))
        y_combined = np.concatenate((r1.y, r2.y), axis=1)

        from types import SimpleNamespace
        result = SimpleNamespace(
            t=t_combined,
            tunit=r1.tunit,
            tscale=r1.tscale,
            y=y_combined,
            concunits="mol/m³",
            species=r1.species,
            species_names=r1.species_names,
            success=r1.success and r2.success,
            message=r1.message + " | " + r2.message,
            solver=r1.solver
        )

        # Create a new instance with copied mixture and set the results
        combined = mixtureKinetics(self.mixture.copy())  # deep copy
        combined._results = result

        # Combine lumped registry
        for k in set(self._lumped_registry.keys()).union(other._lumped_registry.keys()):
            reg1 = self._lumped_registry.get(k)
            reg2 = other._lumped_registry.get(k)
            if reg1 and reg2:
                reg1names = [r.name for r in reg1.species]
                reg2names = [r.name for r in reg2.species]
                if reg1names != reg2names:
                    raise ValueError(f"Lumped group '{k}' mismatch between instances.")
            combined._lumped_registry[k] = reg1 or reg2

        # merged names
        combined.name = self.name + "→" + other.name

        return combined


# -----------------------------------------------------------
# class prototype to define heating and oxygentation cycles
# -----------------------------------------------------------
class TKO2cycle:
    """
    Class to define temperature (T) and oxygenation (kO2) cycles over time.

    Supports:
        - Constant values
        - Linear ramps between two values
        - Multiple steps combined using +

    Each cycle defines:
        - T(t): temperature in Celsius
        - kO2(t): oxygen mass transfer coefficient [m/s]

    Example:
        cycle = TKO2cycle.constant(60, 1e-5, duration=2*3600) + \
                TKO2cycle.ramp(60, 80, 1e-5, 2e-5, duration=1*3600)
    """

    def __init__(self, times, T_values, kO2_values):
        """TKO2cycle constructor (do not use it directly, call constant/ramp instead)"""
        if isinstance(times,list):
            for i,t in enumerate(times):
                times[i],_,_ = mixtureKinetics.parse_time_unit(t,'s')
        elif not isinstance(times,np.ndarray):
            times,_,_ = mixtureKinetics.parse_time_unit(times,'s')
        self.times = np.asarray(times)
        self.T_values = np.asarray(T_values)
        self.kO2_values = np.asarray(kO2_values)

    @classmethod
    def constant(cls, T, kO2, duration):
        """generate an holding step for T and kO2"""
        times = [0, duration]
        T_values = [T, T]
        kO2_values = [kO2, kO2]
        return cls(times, T_values, kO2_values)

    @classmethod
    def ramp(cls, T, kO2, duration):
        """generate a ramp from pair values of T and kO2"""
        times = [0, duration]
        if not isinstance(T,(list,tuple)) or len(T)!=2:
            raise ValueError(f"T must be [Tstart, Tstop] or (Tstart,Tstop) not a {type(T).__name__}")
        if not isinstance(kO2,(list,tuple)) or len(kO2)!=2:
            raise ValueError(f"kO2 must be [kO2start, kO2stop] or (kO2start,kO2stop) not a {type(kO2).__name__}")
        return cls(times, T, kO2)

    def __add__(self, other):
        """add a step to an existing cycle"""
        if not isinstance(other, TKO2cycle):
            return NotImplemented
        t_shift = self.times[-1]
        new_times = np.concatenate([self.times, other.times[1:] + t_shift])
        new_T = np.concatenate([self.T_values, other.T_values[1:]])
        new_kO2 = np.concatenate([self.kO2_values, other.kO2_values[1:]])
        return TKO2cycle(new_times, new_T, new_kO2)

    def __call__(self, t):
        """return interpolated T and kO2 values at t"""
        T = np.interp(t, self.times, self.T_values)
        kO2 = np.interp(t, self.times, self.kO2_values)
        return T, kO2

    @property
    def T(self):
        """T(t) interpolant at t as an anonymous function"""
        return lambda t: self(t)[0]

    @property
    def kO2(self):
        """kO2(t) interpolation at t as an anonymous function"""
        return lambda t: self(t)[1]

    def __repr__(self):
        segments = []
        for i in range(len(self.times) - 1):
            seg = f"[{self.times[i]:.0f}s → {self.times[i+1]:.0f}s]: T={self.T_values[i]}→{self.T_values[i+1]}°C, " \
                  f"kO2={self.kO2_values[i]:.1e}→{self.kO2_values[i+1]:.1e} m/s"
            segments.append(seg)
        return "<TKO2cycle>\n" + "\n".join(segments)

    def __str__(self):
        return f"<TKO2cycle with {len(self.times) - 1} segments, duration {self.times[-1]} s>"

# %% Base species and lumped classes
class species:
    """
    Base class for chemical species involved in radical oxidation.

    This class supports structural and kinetic representation of organic species,
    including radicals, hydroperoxides, oxygen, alcohols, aldehydes, and stable products.
    Derived classes should specify suffix, aliases, reactive functions, and redox ranks.

    Attributes:
        name (str): Long name.
        shortname (str): Abbreviation or label.
        root (str): Identifier used to group related species.
        concentration (float): Current concentration in mol/m³.
        index (int): Index in the mixture (set automatically).

    Class Attributes:
        _registry (dict): Mapping from names and aliases to subclasses.

    Core Properties:
        isradical (bool): True if species has nonzero free valence.
        ishydroperoxide (bool): True if reactiveFunction is "COOH".
        isreactive (bool): True if species has a defined reaction partner/product.
        islumped (bool): True if species is a lumped group of species.
        listProducts (dict): Maps reactive functions → product class names.
        formattedProducts (dict): Maps reactive functions → formatted product names (root + suffix).

    Behavior:
        __or__(other): Lump two species (returns a `lumped` object).
        __add__(other): Placeholder for defining reaction A + B.
        __iter__(): Iterator over main attributes.
        __repr__(), __str__(): Print-friendly formatting.

    Class Methods:
        register(cls): Decorator to register subclasses (with aliases).
        create(name, **kwargs): Factory method to instantiate from alias or class name.

    Notes:
        - `reactiveFunction`, `reactWith`, and `product` must be defined in subclasses.
        - `reactWith` may be a list of functional groups or None (monomolecular).
        - Products can be a string, tuple (ordered), or set (unordered).
        - Naming conventions follow: name = root + suffix (e.g., L1OO• from root L1 + suffix OO•).

    Example:
        @species.register
        class monoAllylicCH(species):
            suffix = "H"
            reactiveFunction = "CH"
            reactWith = ["COO•", "CO•"]
            product = ["monoAllylicC"]
            defaultRoot = "L1"
            defaultName = "L1H"

    See also:
        - mixture: For managing collections of species.
        - reaction: For defining A + B → C + D.
        - lumped: For combining species into pseudo-compounds.
    """

    _registry = {}

    def __init__(self, *, root=None, name='', shortname='', concentration=0.0, index=None, _alias_used=None):
        """
        Initialize a species object with structural, chemical, and contextual attributes.

        This constructor supports flexible initialization via direct arguments or alias-based inference,
        particularly when creating species from shorthand names (e.g., "L1OOH"). Class-level defaults
        (`defaultName`, `defaultRoot`, `suffix`) are used as fallback when arguments are not provided.

        Parameters:
            root (str, optional): Structural identifier prefix (e.g., "L1"). Defaults to `defaultRoot`.
            name (str, optional): Full name of the species. If empty, derived from `_alias_used` or `defaultName`.
            shortname (str, optional): Abbreviated label for display. Defaults to `name` or `defaultName`.
            concentration (float, optional): Initial concentration in mol/m³. Default is 0.0.
            index (int, optional): Numerical index in the mixture or solver system.
            _alias_used (str, internal):
                Optional alias (e.g. "L1OOH") used to parse `root` and assign labels.
                If set, and no explicit name is given, `name`, `shortname`, and `root`
                will be inferred using the class's `suffix`.

        Attributes:
            name (str): Fully qualified species name.
            shortname (str): Shorthand label for output and display.
            root (str): Structural root used for tracking transformations.
            concentration (float): Initial concentration.
            index (int or None): Solver index or None if unassigned.
            _color (str or RGB tuple or None): Color used in plots. If None, automatic coloring is applied.
            _linestyle (str, default="-"): Linestyle used in plots.
            _lumped (bool): Flag indicating whether this species is part of a lumped group (default: False).

        Notes:
            - If `_alias_used` matches a class suffix (e.g., "L1OOH" for suffix "OOH"), the `root` is inferred.
            - Species can be later grouped via `lumped` objects using the `_lumped` flag.
            - The `_color` attribute can be set at the class level (e.g., in `oxygen`) or per instance.

        Example:
            sp = monoAllylicCH(name="L1H", concentration=2.5)
            print(sp.name, sp.root, sp.concentration)
        """

        # Alias-based derivation
        if _alias_used:
            if name == '': name = _alias_used
            if shortname == '': shortname = _alias_used
            if root == '' or root is None:
                if hasattr(self, "suffix") and _alias_used.endswith(self.suffix):
                    root = _alias_used[: -len(self.suffix)]
                else:
                    match = re.match(r'^([A-Z]+\d)', _alias_used)
                    if match: root = match.group(1)
        # Set fallback defaults from class attributes
        self.name = name if name else getattr(self, "defaultName", "")
        self.shortname = shortname if shortname else getattr(self, "defaultName", "")
        self.root = root if root else getattr(self, "defaultRoot", "")
        # other Parameters
        self.concentration = concentration
        self.index = index
        self._color = getattr(self.__class__, "_color", None)
        self._linestyle = getattr(self.__class__, "_linestyle", '-')
        self._lumped = True

    def __iter__(self):
        yield from {
            'name': self.name,
            'shortname': self.shortname,
            'root': self.root,
            'concentration': self.concentration,
            'index': self.index
        }.items()

    def __repr__(self):
        cls_attrs = {
            k: getattr(self, k) for k in dir(self)
            if not k.startswith('_')
            and not callable(getattr(self, k))
            and not isinstance(getattr(self, k), property)
        }
        inst_attrs = self.__dict__
        all_attrs = {**inst_attrs, **cls_attrs}
        width = max(len(k) for k in all_attrs)
        lines = [f"{k.rjust(width)}: {v}" for k, v in all_attrs.items()]
        print("\n".join(lines))
        return str(self)

    def __str__(self):
        return f"<[{self.shortname}]={self.concentration} mol/m³> of class {self.__class__.__name__}, idx={self.index}"

    def __or__(self, other):
        return lumped([self, other])

    def __add__(self, other):
        return (self, other)  # placeholder for reaction object

    @classmethod
    def register(cls, derived_cls):
        """
        Class decorator to register a derived species class.

        Adds both the class name and any class-defined aliases to the registry.

        Args:
            derived_cls (type): The class to register.

        Returns:
            type: The same class, unmodified.
        """
        cls._registry[derived_cls.__name__] = derived_cls
        for alias in getattr(derived_cls, "classAliases", []):
            cls._registry[alias] = derived_cls
        return derived_cls

    @staticmethod
    def create(name, **kwargs):
        """
        Factory method to create a species by name or alias.

        Args:
            name (str): Name or alias of the species.
            **kwargs: Passed to the constructor of the target class.

        Returns:
            species: An instance of the requested class.

        Raises:
            ValueError: If no class matches the name or alias.
        """
        cls = species._registry.get(name)
        kwargs["_alias_used"] = name
        if cls is None:
            raise ValueError(f"Unknown species or alias: {name}")
        return cls(**kwargs)

    @property
    def isradical(self):
        """True if species has free valence (i.e., is a radical)."""
        if not hasattr(self, "freevalence") or self.freevalence is None:
            return None
        return self.freevalence > 0 if hasattr(self, "freevalence") else False

    @property
    def isradical_onC(self):
        """True if radical center is on carbon."""
        if not hasattr(self, "freevalenceHolder") or self.freevalenceHolder is None:
            return None
        return self.freevalenceHolder == "C" if self.freevalence else False

    @property
    def isradical_onO(self):
        """True if radical center is on oxygen."""
        if not hasattr(self, "freevalenceHolder") or self.freevalenceHolder is None:
            return None
        return self.freevalenceHolder == "O" if self.freevalence else False

    @property
    def ishydroperoxide(self):
        """True if reactiveFunction is 'COOH', interpreted as a hydroperoxide."""
        if not hasattr(self,"reactiveFunction") or self.reactiveFunction is None:
            return None
        return self.reactiveFunction == "COOH"

    @property
    def isreactive(self):
        """True if the species leads to at least one product."""
        if not hasattr(self,"reactiveFunction") or self.reactiveFunction is None:
            return None
        return hasattr(self, "product") and self.product is not None

    @property
    def ismonomolecular(self):
        """
        Returns True if the species reacts without a partner (reactWith is None).
        If reactWith is a list, returns a list of bools for each reactant type.
        Returns False if not reactive.
        """
        if not self.isreactive:
            return False
        rw = getattr(self, "reactWith", None)
        if rw is None:
            return True
        if isinstance(rw, list):
            return [r is None for r in rw]
        return rw is None

    @property
    def isbimolecular(self):
        """
        Returns True if the species reacts with another reactant (reactWith is not None).
        If reactWith is a list, returns a list of bools for each reactant type.
        Returns False if not reactive.
        """
        if not self.isreactive:
            return False
        rw = getattr(self, "reactWith", None)
        if rw is None:
            return False
        if isinstance(rw, list):
            return [r is not None for r in rw]
        return rw is not None

    @property
    def islumped(self):
        """True if this species is a lumped group of multiple species."""
        return isinstance(self, lumped)

    @property
    def listProducts(self):
        """
        Returns a dictionary mapping the type of reaction to the corresponding product.

        Keys:
            'self' → product of monomolecular decomposition (i.e., if `None` is in reactWith)
            <functional_group> → product when reacting with that function (from reactWith)

        Returns:
            dict[str, str or None]: Maps reaction type to product name.
        """
        if not self.isreactive: return {}
        rw = self.reactWith if isinstance(self.reactWith, list) else [self.reactWith]
        prod = self.product
        # normalize prod to list same length as rw
        if isinstance(prod, tuple):
            prod_list = list(prod)
        elif isinstance(prod, list):
            prod_list = prod[:]  # may be shorter or longer
        else:
            prod_list = [prod]
        if len(prod_list)==1 and len(rw)>1:
            prod_list *= len(rw)
        if len(prod_list)!=len(rw):
            prod_list = (prod_list + [None]*len(rw))[:len(rw)]
        result = {}
        for r,p in zip(rw, prod_list):
            key = 'self' if r is None else r
            result[key] = p
        return result

    @property
    def formattedProducts(self):
        """
        Returns a dictionary mapping reaction types to formatted product names
        using root + suffix (e.g., 'L1OO•' instead of 'monoAllylicCOO').

        Keys:
            'self' → product of monomolecular decomposition
            <functional_group> → product from bimolecular reaction

        Returns:
            dict[str, str or list[str]]: Maps reaction type to formatted name(s)
        """
        if not self.isreactive:
            return {}

        rw = self.reactWith
        prod = self.product

        rw = [] if rw is None else rw
        if not isinstance(rw, list):
            rw = [rw]

        # Normalize product list to match reactWith length
        if isinstance(prod, (str, tuple, set)):
            prod = [prod] * len(rw)
        elif isinstance(prod, list):
            if len(prod) == 1 and len(rw) > 1:
                prod = prod * len(rw)
            elif len(prod) != len(rw):
                prod = [None] * len(rw)

        result = {}
        root = self.root or ""

        for r, p in zip(rw, prod):
            key = "self" if r is None else r

            # Handle set of unordered products
            if isinstance(p, set):
                formatted = []
                for i, name in enumerate(p):
                    cls = species._registry.get(name)
                    suffix = getattr(cls, "suffix", "") if cls else ""
                    has_root = getattr(cls, "defaultRoot", "") != ""
                    formatted_name = (root + suffix) if has_root or i == 0 else suffix
                    formatted.append(formatted_name)
                result[key] = formatted
            elif isinstance(p, (tuple, list)):
                formatted = []
                for name in p:
                    cls = species._registry.get(name)
                    suffix = getattr(cls, "suffix", "") if cls else ""
                    formatted.append(root + suffix if cls else name)
                result[key] = formatted
            else:
                cls = species._registry.get(p)
                suffix = getattr(cls, "suffix", "") if cls else ""
                result[key] = root + suffix if cls else p

        return result


# Lumped Class
@species.register
class lumped(species):
    """
    Lumped species class representing the aggregate behavior of multiple individual species.

    This class is used to group species that are chemically or functionally equivalent
    in a given context (e.g., radicals with similar reactivity, isomers, or chain-length variants).

    The `lumped` object behaves like a single species whose concentration is the sum
    of its components. It inherits from the `species` base class but overrides key attributes.

    Parameters:
        species_list (list[species]): list of `species` instances being lumped.
        color (str or [RGB] values): color to apply for plotting (default=None for automatic assignment).
        linestyle (str, default="-"): Linestyle used in plots.

    Attributes:
        species (list[species]): List of `species` instances being lumped.
        name (str): Default name is "lump_" + joined individual names.
        index (None): Lumped species are not indexed in the mixture.
        concentration (float): Total concentration of included species.
        shortname (str): Default shorthand label ("LM").
        root (None): No root is assigned.
        reactiveFunction (None): No reactive function is assigned.
        freevalence (None): No freevalence is assigned.
        _color: color to apply for plotting.
        _linestyle (str, default="-"): Linestyle used in plots.
        _lumped (True): True for lumped species.

    Usage Example:
        >>> LOOH1 = species.create("L1OOH", concentration=0.1)
        >>> LOOH2 = species.create("L2OOH", concentration=0.2)
        >>> LOOH3 = species.create("L3OOH", concentration=0.3)
        >>> LOOH = lumped([LOOH1, LOOH2, LOOH3])
        >>> print(LOOH.name)
        lump_L1OOH_L2OOH_L3OOH
        >>> print(LOOH.concentration)
        0.6

    Notes:
        - Lumped species are useful for reducing model complexity or matching experimental resolution.
        - They can be passed to plotting, integration, or export functions just like other species.
        - Operations on lumped species do not propagate changes back to the individual species.

    See also:
        - species.__or__: Supports the `|` operator for combining species.
    """

    def __init__(self, species_list,color=None,linestyle="-"):
        self.species = species_list
        self.name = "lump_" + "_".join(sp.name for sp in species_list)
        self.index = None
        self.concentration = sum(sp.concentration for sp in species_list)
        self.shortname = "LM"
        self.root = None
        self.reactiveFunction = None
        self.freevalence = None
        self._color = color
        self._linestyle = linestyle
        self._lumped = True


# %% Specific Chemical Classes derived from species
"""
Specific Reaction Rate Constants
--------------------------------

This section provides structured kinetic data extracted from the two primary sources:

[1] Touffet M., Smith P., Vitrac O., *A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures*,
    Food Research International, 173(1), 2023, 113289. https://doi.org/10.1016/j.foodres.2023.113289

[2] Touffet M., Vitrac O., *Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms*,
    Food Chemistry, 481, 2025, 143952. https://doi.org/10.1016/j.foodchem.2025.143952

Each reaction is represented by a unique fingerprint (e.g., "L1H + HO• → L1• + H2O") and registered via the `reactionRateDB` class.

The database supports:
    - Arrhenius parameters (k₀, Ea, T₀)
    - Diffusion-limited kinetics (Smoluchowski model)
    - Viscosity dependence (η(T))
    - Cross-reaction inference via geometric averaging
    - Consistent unit handling and data export to pandas DataFrame

🔎 Redox Ranking Justification (`_redoxRank`)
--------------------------------------------

The `_redoxRank` attribute is used to consistently order reactants in A + B → C + D schemes.
It reflects both electron donor/acceptor strength and oxidation state:

| Species Class           | Description                           | Oxidation State | Redox Rank | Notes                        |
|------------------------|---------------------------------------|-----------------|------------|------------------------------|
| Carboxylic Acid        | COOH group                            | +3              | 1          | Stable product               |
| Aldehyde/Ketone        | CHO or C=O                            | +2 to +1        | 2          | From beta-scission           |
| Alcohol                | COH                                   | 0               | 3          | From reduction               |
| Alkyl radical (C•)     | Unpaired electron on C                | ~0              | 4          | Initial oxidation product    |
| Alkoxyl radical (RO•)  | Unpaired electron on O                | –1              | 5          | From ROOH decomposition      |
| Hydroperoxide (ROOH)   | OOH group                             | –1              | 6          | Source of radicals           |
| Peroxyl radical (ROO•) | O–O• group                            | ~0              | 7          | Common chain carrier         |
| Dioxygen (O₂)          | Triplet O₂                            | 0               | 8          | Strong oxidant               |

This ranking ensures that symmetric reactions are ordered consistently and
that cross-reactions can be constructed with predictable canonical forms.


🧪 Oxidation State (OS) Estimation and Allylic Context
------------------------------------------------------

The allylic context of RH donors (CH-type) affects both bond dissociation energy and redox potential. This table connects allylicity to effective electron-donating character and thus reactivity.

| Species         | Allylic context   | OS estimate (central C) | Justification                                  |
|-----------------|-------------------|--------------------------|------------------------------------------------|
| `monoAllylicCH` | C=C–CH            | –1                       | Single allylic C–H                             |
| `diAllylicCH`   | C=C–CH–C=C        | ~–0.5 to –1              | More resonance, slightly more electron-rich    |
| `triAllylicCH`  | C=C–CH–C=C–C=C    | ~0                       | High delocalization, lower effective reduction |

This classification is crucial for automated rule generation and mechanistic modeling across large chemical libraries or in generative workflows (LLM-driven).
"""

@species.register
class oxygen(species):
    """Dioxygen (•OO•)"""
    classAliases = ['O2']
    className = "dioxygen (triplet state)"
    defaultName = "oxygen"
    defaultRoot = ""
    suffix = "O2"
    allylic = 0
    reactiveFunction = "•OO•"
    freevalence = 2
    freevalenceHolder = 'O'
    reactWith = ['C•']
    product = ["-"]
    interpretation = """O2 is fixed"""
    oxidationState = 0
    _redoxRank = 8
    _color = "Red"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicCH(species):
    """Aliphatic (CH) on monoallylic site"""
    classAliases = ['C1H', 'R1H', 'L1H', 'P1H']
    className = "aliphatic with monoallylic C-H"
    defaultName = "L1H"
    defaultRoot = "L1"
    suffix = "H"
    allylic = 1
    reactiveFunction = "CH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = ['HO•', 'CO•', 'COO•']
    product = ['monoAllylicC']
    interpretation = """labile H abstraction"""
    oxidationState = -1
    _redoxRank = 3
    _color = "DimGray"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCH(species):
    """Aliphatic (CH) on diallylic site"""
    classAliases = ['C2H', 'R2H', 'L2H', 'P2H']
    className = "aliphatic with diallylic C-H"
    defaultName = "L2H"
    defaultRoot = "L2"
    suffix = "H"
    allylic = 2
    reactiveFunction = "CH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = ['HO•', 'CO•', 'COO•']
    product = ['diAllylicC']
    interpretation = """labile H abstraction"""
    oxidationState = -0.5
    _redoxRank = 2.5
    _color = "DarkGray"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCH(species):
    """Aliphatic (CH) on triallylic site"""
    classAliases = ['C3H', 'R3H', 'L3H', 'P3H']
    className = "aliphatic with triallylic C-H"
    defaultName = "L3H"
    defaultRoot = "L3"
    suffix = "H"
    allylic = 3
    reactiveFunction = "CH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = ['HO•', 'CO•', 'COO•']
    product = ['triAllylicC']
    interpretation = """labile H abstraction"""
    oxidationState = 0
    _redoxRank = 2
    _color = "Silver"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicCOOH(species):
    """Hydroperoxide (COOH) on monoallylic site"""
    classAliases = ['C1OOH', 'R1OOH', 'L1OOH', 'P1OOH']
    className = "hydroperoxide on monoallylic site"
    defaultName = "L1OOH"
    defaultRoot = "L1"
    suffix = "OOH"
    allylic = 1
    reactiveFunction = "COOH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = [None, 'COOH']
    product = [{'monoAllylicCO','peroxyl'},('monoAllylicCO','monoAllylicCOO')]
    interpretation = """monomolecular or bimolecular decomposition"""
    oxidationState = -1
    _redoxRank = 6
    _color = "DarkRed"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCOOH(species):
    """Hydroperoxide (COOH) on diallylic site"""
    classAliases = ['C2OOH', 'R2OOH', 'L2OOH', 'P2OOH']
    className = "hydroperoxide on diallylic site"
    defaultName = "L2OOH"
    defaultRoot = "L2"
    suffix = "OOH"
    allylic = 2
    reactiveFunction = "COOH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = [None, 'COOH']
    product = [{'diAllylicCO','peroxyl'},('diAllylicCO','diAllylicCOO')]
    interpretation = """monomolecular or bimolecular decomposition"""
    oxidationState = -1
    _redoxRank = 6
    _color = "FireBrick"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCOOH(species):
    """Hydroperoxide (COOH) on triallylic site"""
    classAliases = ['C3OOH', 'R3OOH', 'L3OOH', 'P3OOH']
    className = "hydroperoxide on diallylic site"
    defaultName = "L3OOH"
    defaultRoot = "L3"
    suffix = "OOH"
    allylic = 3
    reactiveFunction = "COOH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = [None, 'COOH']
    product = [{'triAllylicCO','peroxyl'},('triAllylicCO','triAllylicCOO')]
    interpretation = """monomolecular or bimolecular decomposition"""
    oxidationState = -1
    _redoxRank = 6
    _color = "Crimson"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicC(species):
    """Alkyl radical (C•) on monoallylic site"""
    classAliases = ['C1', 'R1', 'L1', 'P1']
    className = "alkyl radical on monoallylic site"
    defaultName = "L1•"
    defaultRoot = "L1"
    suffix = "•"
    allylic = 1
    reactiveFunction = "C•"
    freevalence = 1
    freevalenceHolder = 'C'
    reactWith = ['•OO•','C•','COO•']
    product = ['monoAllylicCOO','terminationPolymers','terminationPolymers']
    interpretation = """O2 addition"""
    oxidationState = 0
    _redoxRank = 4
    _color = "DarkBlue"
    _linestyle = ":"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicC(species):
    """Alkyl radical (C•) on diallylic site"""
    classAliases = ['C2', 'R2', 'L2', 'P2']
    className = "alkyl radical on diallylic site"
    defaultName = "L2•"
    defaultRoot = "L2"
    suffix = "•"
    allylic = 2
    reactiveFunction = "C•"
    freevalence = 1
    freevalenceHolder = 'C'
    reactWith = ['•OO•','C•','COO•']
    product = ['diAllylicCOO','terminationPolymers','terminationPolymers']
    interpretation = """O2 addition"""
    oxidationState = 0
    _redoxRank = 4
    _color = "MediumBlue"
    _linestyle = ":"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicC(species):
    """Alkyl radical (C•) on triallylic site"""
    classAliases = ['C3', 'R3', 'L3', 'P3']
    className = "alkyl radical on triallylic site"
    defaultName = "L3•"
    defaultRoot = "L3"
    suffix = "•"
    allylic = 3
    reactiveFunction = "C•"
    freevalence = 1
    freevalenceHolder = 'C'
    reactWith = ['•OO•','C•','COO•']
    product = ['triAllylicCOO','terminationPolymers','terminationPolymers']
    interpretation = """O2 addition"""
    oxidationState = 0
    _redoxRank = 4
    _color = "RoyalBlue"
    _linestyle = ":"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicCO(species):
    """Alkoxyl radical (CO•) on monoallylic site"""
    classAliases = ['C1O', 'R1O', 'L1O', 'P1O']
    className = "alkoxyl radical on monoallylic site"
    defaultName = "L1O•"
    defaultRoot = "L1"
    suffix = "O•"
    allylic = 1
    reactiveFunction = "CO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = [None, 'CH']
    product = ['monoAllylicCeqO', 'monoAllylicCOH']
    interpretation = """beta-scission or reduction to alcohol"""
    oxidationState = -1
    _redoxRank = 5
    _color = "OrangeRed"
    _linestyle = "-."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCO(species):
    """Alkoxyl radical (CO•) on diallylic site"""
    classAliases = ['C2O', 'R2O', 'L2O', 'P2O']
    className = "alkoxyl radical on diallylic site"
    defaultName = "L2O•"
    defaultRoot = "L2"
    suffix = "O•"
    allylic = 2
    reactiveFunction = "CO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = [None, 'CH']
    product = ['diAllylicCeqO', 'diAllylicCOH']
    interpretation = """beta-scission or reduction to alcohol"""
    oxidationState = -1
    _redoxRank = 5
    _color = "Tomato"
    _linestyle = "-."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCO(species):
    """Alkoxyl radical (CO•) on triallylic site"""
    classAliases = ['C3O', 'R3O', 'L3O', 'P3O']
    className = "alkoxyl radical on triallylic site"
    defaultName = "L3O•"
    defaultRoot = "L3"
    suffix = "O•"
    allylic = 3
    reactiveFunction = "CO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = [None, 'CH']
    product = ['triAllylicCeqO', 'triAllylicCOH']
    interpretation = """beta-scission or reduction to alcohol"""
    oxidationState = -1
    _redoxRank = 5
    _color = "Coral"
    _linestyle = "-."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicCOO(species):
    """Hydroperoxyl radical (COO•) on monoallylic site"""
    classAliases = ['C1OO', 'R1OO', 'L1OO', 'P1OO']
    className = "hydroperoxyl radical on monoallylic site"
    defaultName = "L1OO•"
    defaultRoot = "L1"
    suffix = "OO•"
    allylic = 1
    reactiveFunction = "COO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = ['CH', 'COO•','C•']
    product = ['monoAllylicCOOH', 'terminationPolymers','terminationPolymers']
    interpretation = """H abstraction termination"""
    oxidationState = 0
    _redoxRank = 7
    _color = "DeepPink"
    _linestyle = "--"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCOO(species):
    """Hydroperoxyl radical (COO•) on diallylic site"""
    classAliases = ['C2OO', 'R2OO', 'L2OO', 'P2O']
    className = "hydroperoxyl radical on diallylic site"
    defaultName = "L2OO•"
    defaultRoot = "L2"
    suffix = "OO•"
    allylic = 2
    reactiveFunction = "COO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = ['CH', 'COO•','C•']
    product = ['diAllylicCOOH', 'terminationPolymers','terminationPolymers']
    interpretation = """H abstraction termination"""
    oxidationState = 0
    _redoxRank = 7
    _color = "HotPink"
    _linestyle = "--"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCOO(species):
    """Hydroperoxyl radical (COO•) on triallylic site"""
    classAliases = ['C3OO', 'R3OO', 'L3OO', 'P3OO']
    className = "hydroperoxyl radical on triallylic site"
    defaultName = "L3OO•"
    defaultRoot = "L3"
    suffix = "OO•"
    allylic = 3
    reactiveFunction = "COO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = ['CH', 'COO•','C•']
    product = ['triAllylicCOOH', 'terminationPolymers','terminationPolymers']
    interpretation = """H abstraction termination"""
    oxidationState = 0
    _redoxRank = 7
    _color = "LightPink"
    _linestyle = "--"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class terminationPolymers(species):
    """Generic stable polymer termination product"""
    classAliases = ['polymer']
    className = "stable polymer termination product"
    defaultName = "termination product (polymer)"
    defaultRoot = "polymer"
    suffix = "-polymer"
    allylic = None
    reactiveFunction = "-polymer"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product =None
    interpretation = """stable termination products (polymers)"""
    oxidationState = 0
    _redoxRank = None
    _color = "Black"
    _linestyle = "-"

@species.register
class monoAllylicCeqO(species):
    """Ketone (C=O) on monoallylic site (stable)"""
    classAliases = ['C1=O', 'R1=O', 'L1=O', 'P1=O']
    className = "ketone on monoallylic site"
    defaultName = "L1=O"
    defaultRoot = "L1"
    suffix = "=O"
    allylic = 1
    reactiveFunction = "C=O"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable ketone"""
    oxidationState = 1
    _redoxRank = 2
    _color = "Maroon"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCeqO(species):
    """Ketone (C=O) on diallylic site (stable)"""
    classAliases = ['C2=O', 'R2=O', 'L2=O', 'P2=O']
    className = "ketone on diallylic site"
    defaultName = "L2=O"
    defaultRoot = "L2"
    suffix = "=O"
    allylic = 2
    reactiveFunction = "C=O"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable ketone"""
    oxidationState = 1
    _redoxRank = 2
    _color = "SaddleBrown"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCeqO(species):
    """Ketone (C=O) on triallylic site (stable)"""
    classAliases = ['C3=O', 'R3=O', 'L3=O', 'P3=O']
    className = "ketone on triallylic site"
    defaultName = "L3=O"
    defaultRoot = "L3"
    suffix = "=O"
    allylic = 3
    reactiveFunction = "C=O"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable ketone"""
    oxidationState = 1
    _redoxRank = 2
    _color = "Sienna"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


@species.register
class monoAllylicCHO(species):
    """Aldehyde (CHO) on monoallylic site (stable)"""
    classAliases = ['C1=O', 'R1=O', 'L1=O', 'P1=O']
    className = "aldehyde on monoallylic site"
    defaultName = "L1HO"
    defaultRoot = "L1"
    suffix = "=O"
    allylic = 1
    reactiveFunction = "CHO"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable aldehyde"""
    oxidationState = 2
    _redoxRank = 2
    _color = "DarkGoldenrod"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCHO(species):
    """Aldehyde (CHO) on diallylic site (stable)"""
    classAliases = ['C2=O', 'R2=O', 'L2=O', 'P2=O']
    className = "aldehyde on diallylic site"
    defaultName = "L2HO"
    defaultRoot = "L2"
    suffix = "=O"
    allylic = 2
    reactiveFunction = "CHO"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable aldehyde"""
    oxidationState = 2
    _redoxRank = 2
    _color = "Goldenrod"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class triAllylicCHO(species):
    """Aldehyde (CHO) on triallylic site (stable)"""
    classAliases = ['C3=O', 'R3=O', 'L3=O', 'P3=O']
    className = "aldehyde on triallylic site"
    defaultName = "L3HO"
    defaultRoot = "L3"
    suffix = "=O"
    allylic = 3
    reactiveFunction = "CHO"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable aldehyde"""
    oxidationState = 2
    _redoxRank = 2
    _color = "Gold"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class monoAllylicCOH(species):
    """Alcohol (COH) on monoallylic site (stable)"""
    classAliases = ['C1OH', 'R1OH', 'L1OH', 'P1OH']
    className = "alcohol on monoallylic site"
    defaultName = "L1OH"
    defaultRoot = "L1"
    suffix = "OH"
    allylic = 1
    reactiveFunction = "COH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable alcohol"""
    oxidationState = 0
    _redoxRank = 3
    _color = "IndianRed"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class diAllylicCOH(species):
    """Alcohol (COH) on dioallylic site (stable)"""
    classAliases = ['C1OH', 'R1OH', 'L1OH', 'P1OH']
    className = "alcohol on diallylic site"
    defaultName = "L2OH"
    defaultRoot = "L2"
    suffix = "OH"
    allylic = 2
    reactiveFunction = "COH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable alcohol"""
    oxidationState = 0
    _redoxRank = 3
    _color = "Salmon"
    _linestyle = "-"

@species.register
class triAllylicCOH(species):
    """Alcohol (COH) on trioallylic site (stable)"""
    classAliases = ['C1OH', 'R1OH', 'L1OH', 'P1OH']
    className = "alcohol on triallylic site"
    defaultName = "L3OH"
    defaultRoot = "L3"
    suffix = "OH"
    allylic = 3
    reactiveFunction = "COH"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable alcohol"""
    oxidationState = 0
    _redoxRank = 3
    _color = "LightSalmon"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class peroxyl(species):
    """Peroxyl radical (HO•) class"""
    classAliases = ['OH']
    className = "peroxyl radical"
    defaultName = "HO•"
    defaultRoot = ""
    suffix = "HO•"
    allylic = 0
    reactiveFunction = "HO•"
    freevalence = 1
    freevalenceHolder = 'O'
    reactWith = ['CH']
    product = ['H2O']
    interpretation = """labile H abstraction"""
    oxidationState = 0
    _redoxRank = 7
    _color = "DarkOrange"
    _linestyle = "-."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


@species.register
class H2O(species):
    """Water class (stable)"""
    classAliases = ['H2O']
    className = "water"
    defaultName = "H2O"
    defaultRoot = ""
    suffix = "H2O"
    allylic = 0
    reactiveFunction = "None"
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product = None
    interpretation = """stable water"""
    oxidationState = 0
    _redoxRank = 7
    _color = "DeepSkyBlue"
    _linestyle = "-"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

# %% Base Reaction Rate Constants Class
class reactionRateDB:
    """
    Database and container class for reaction rate constants.

    This class enables registration, storage, and retrieval of kinetic parameters
    associated with elementary reactions involved in radical oxidation mechanisms.
    It is designed to support both monomolecular and bimolecular reactions, including
    Arrhenius parameters, solvent-dependent viscosity corrections, and diffusion limits.

    Each entry corresponds to a uniquely defined reaction fingerprint (e.g., "L1OOH -> L1O• + HO•")
    and may include confidence intervals, diffusion-limited behavior, and Smoluchowski model parameters.

    === Attributes ===
    Class Attributes:
        _registry (dict): Maps canonical fingerprints to lists of registered entries.
        _fingerprint_substitutions (list): List of (pattern, replacement) regex pairs
            used to transform fingerprints queried via `.get()` or `.get_latest()`.

    Instance Attributes:
        fingerprint (str): Reaction descriptor in canonical form.
        T0 (float): Reference temperature [°C].
        Tunits (str): Units of T0.
        k0 (float): Pre-exponential factor (Arrhenius rate at T0).
        k0_CI (float or None): Confidence interval on k0.
        kunits (str): Units of k0 (e.g., "m³·mol⁻¹·s⁻¹").
        Ea (float): Activation energy [J/mol].
        Ea_CI (float or None): Confidence interval on Ea.
        Eaunits (str): Units of activation energy.
        source (str): Source or reference for the entry (e.g., literature citation).
        diffusionLimited (bool): Whether rate is capped by diffusion.
        ratio_rgh (float): Radius-of-gyration to hydrodynamic-radius ratio.
        eta0 (float): Reference viscosity [Pa·s].
        eta_Ea (float): Activation energy for viscosity [J/mol].
        eta_T0 (float): Reference temperature for viscosity model [°C].
        etaunits (str): Units for viscosity.
        phi (float): Scaling factor used for cross-reaction inference.

    === Methods ===
    Instance:
        __init__: Creates and registers a new rate constant entry.
        __iter__: Iterator over main fields.
        __str__, __repr__: Human-readable and detailed representations.

    Class:
        get(fingerprint): Return all entries matching the given fingerprint or a substituted one.
        get_latest(fingerprint): Return the latest entry registered for a (possibly substituted) fingerprint.
        to_dataframe(): Export the full database to a pandas DataFrame.
        register_fingerprint_substitution(pattern, replacement): Register a regex substitution to transform
            a queried fingerprint into a canonical one, allowing lookup under different naming schemes.

    === Fingerprint Substitution ===
    If a reaction fingerprint is not found in the registry, a set of regex-based substitutions
    is applied in sequence to find alternative forms (e.g., "L1H" → "R1H"). This enables
    the reuse of reaction rate constants across different naming conventions.

    Example:
        reactionRateDB.register_fingerprint_substitution(r"L([123])", r"R\\1")

        # Register using R-prefixed names
        reactionRateDB(fingerprint="R1H + R1O• -> R1• + R1OH", ...)

        # Lookup via L-prefixed form will succeed:
        reactionRateDB.get_latest("L1H + L1O• -> L1• + L1OH")

    Notes:
        - Only the first matching substitution is applied.
        - Substitution is used only for lookup; registered fingerprints remain unchanged.

    === Usage Example ===
        reactionRateDB(
            fingerprint="L1OOH -> L1O• + HO•",
            T0=140,
            k0=2.11e-4,
            Ea=74.8e3,
            source="Touffet et al. 2023"
        )

    See also:
        - reaction: Uses this class to populate kinetic parameters.
        - mixture.populateReactionRates(): Assigns values to all reactions in a mixture.
    """


    _registry = defaultdict(list)

    def __init__(self, *, fingerprint, T0, k0, k0_CI=None, Ea, Ea_CI=None,
                 Tunits = "°C", kunits = "m³·mol⁻¹·s⁻¹", Eaunits = "J·mol⁻¹", source="",
                 diffusionLimited = False, ratio_rgh = 0.77,
                 eta0 = 1.9e-3, eta_Ea = 17.0e3, eta_T0 = 80, etaunits = "Pa·s",
                 phi = 2.0
                 ):

        self.fingerprint = fingerprint
        self.T0 = T0
        self.Tunits = Tunits
        self.k0 = k0
        self.k0_CI = k0_CI
        self.kunits = kunits
        self.Ea = Ea
        self.Ea_CI = Ea_CI
        self.Eaunits = Eaunits
        self.source = source
        self.diffusionLimited = diffusionLimited
        self.ratio_rgh = ratio_rgh
        self.eta0 = eta0
        self.eta_Ea = eta_Ea
        self.eta_T0 = eta_T0
        self.etaunits = etaunits
        self.phi = phi

        # Register in global registry (class attribute)
        self.__class__._registry[fingerprint].append(self)
        # Substitution rules (class attribute)
        # _fingerprint_substitutions: list[tuple[re.Pattern, str]] = []
        self.__class__._fingerprint_substitutions = []

    @classmethod
    def register_fingerprint_substitution(cls, pattern: str, replacement: str):
        """
        Register a regex-based substitution rule to rewrite reaction fingerprints.
        The first rule that matches will be used.

        Example:
            pattern = r"L([123])"       # matches L1, L2, L3
            replacement = r"R\1"        # replaces with R1, R2, R3
        """
        compiled = re.compile(pattern)
        cls._fingerprint_substitutions.append((compiled, replacement))

    @classmethod
    def _apply_fingerprint_substitution(cls, fingerprint: str) -> str:
        for pattern, repl in cls._fingerprint_substitutions:
            if pattern.search(fingerprint):
                return pattern.sub(repl, fingerprint)
        return fingerprint  # no substitution applied

    def __iter__(self):
        yield from {
            "fingerprint": self.fingerprint,
            "T0": self.T0,
            "Tunits": self.Tunits,
            "k0": self.k0,
            "k0_CI": self.k0_CI,
            "kunits": self.kunits,
            "Ea": self.Ea,
            "Ea_CI": self.Ea_CI,
            "Eaunits": self.Eaunits,
            "source": self.source,
        }.items()

    def __repr__(self):
        props = {
            "reaction": self.fingerprint,
            "T0": f"{self.T0:.1f} [{self.Tunits}]",
            "k0": f"{self.k0:.3g} [{self.kunits}]",
            "k0_CI": f"{self.k0_CI:.3g} [{self.kunits}]" if self.k0_CI else "N/A",
            "Ea": f"{self.Ea:.3g} [{self.Eaunits}]",
            "Ea_CI": f"{self.Ea_CI:.3g} [{self.Eaunits}]" if self.Ea_CI else "N/A",
            "from": self.source if self.source and self.source !="" else "N/A"
        }
        width = max(len(k) for k in props)
        print("\n".join(f"{k.rjust(width)}: {v}" for k, v in props.items()))
        return str(self)

    def __str__(self):
        return f"<{self.fingerprint}> reaction rate constants"

    @classmethod
    def get(cls, fingerprint):
        """Returns all entries matching the fingerprint or its substituted variant."""
        actual = fingerprint
        if fingerprint not in cls._registry:
            actual = cls._apply_fingerprint_substitution(fingerprint)
        return cls._registry.get(actual, [])

    @classmethod
    def get_latest(cls, fingerprint):
        """Returns the last entry added for a fingerprint (with substitution if needed)."""
        entries = cls.get(fingerprint)
        return entries[-1] if entries else None

    @classmethod
    def to_dataframe(cls):
        return pd.DataFrame([
            {
                "reaction": r.fingerprint,
                "T0 (°C)": r.T0,
                "k0 (SI units)": r.k0,
                #"k0_CI (SI units)": r.k0_CI,
                f"Ea ({r.Eaunits})": r.Ea,
                #f"Ea_CI ({r.Eaunits})": r.Ea_CI,
                "eta0 (r.etaunits)": r.eta0,
            }
            for rgroup in cls._registry.values()
            for r in rgroup
        ])


# %% Specific Reaction Rate Constants
"""
Specific Reaction Rate Constants
--------------------------------

This section registers experimental or estimated reaction rate constants using the `reactionRateDB` class. The database enables consistent retrieval, inference, and propagation of kinetic parameters (Arrhenius and Smoluchowski terms) across all supported reactions.

Sources:
    [1] Touffet M., Smith P., Vitrac O., *A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures*,
        Food Research International, 173(1), 2023, 113289. https://doi.org/10.1016/j.foodres.2023.113289

    [2] Touffet M., Vitrac O., *Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms*,
        Food Chemistry, 481, 2025, 143952. https://doi.org/10.1016/j.foodchem.2025.143952

Functionalities:
    - Arrhenius and Smoluchowski kinetic models
    - Inheritance of cross-reactions via geometric mean logic
    - Viscosity dependence of diffusion-limited reactions
    - Documentation of each kinetic entry (units, source)
    - Automatic integration into simulation objects
"""

# R1-3 (monomolecular decomposition of hydroperoxides) in Touffet et al., 2023
kd1 = reactionRateDB(
    fingerprint="L1OOH -> L1O• + HO•",
    T0=140,
    k0=3.0e-4, #2.11e-4,
    k0_CI=1.86e-4,
    kunits = "s⁻¹",
    Ea=89.0e3, #74.8e3,
    Ea_CI=5.58e3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kd2 = reactionRateDB(
    fingerprint="L2OOH -> L2O• + HO•",
    T0=140,
    k0=7.0e-4, #5.93e-4,
    k0_CI=5.23e-4,
    kunits = "s⁻¹",
    Ea=79.0e3, #64.7e3,
    Ea_CI=4.26e3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kd3 = reactionRateDB(
    fingerprint="L3OOH -> L3O• + HO•",
    T0=140,
    k0=9.0e-4, #6.11e-4,
    k0_CI=5.44e-4,
    kunits = "s⁻¹",
    Ea=67.0e3, #51.2e3,
    Ea_CI=8.23e3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)


# R4-12 (H-abstraction from alkoxyl radical, only self-reactions) in Touffet et al., 2023
kb1 = reactionRateDB(
    fingerprint="L1H + L1O• -> L1• + L1OH",
    T0=30,
    k0=3.80e3,
    k0_CI=0.23e3,
    Ea=14.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kb2 = reactionRateDB(
    fingerprint="L2H + L2O• -> L2• + L2OH",
    T0=30,
    k0=8.80e3,
    k0_CI=0.53e3,
    Ea=14.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kb3 = reactionRateDB(
    fingerprint="L3H + L3O• -> L3• + L3OH",
    T0=140,
    k0=13.0e3,
    k0_CI=0.78e3,
    Ea=14.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R13-15 (H-abstraction from peroxyl radical) in Touffet et al., 2023
kbb1 = reactionRateDB(
    fingerprint="L1H + HO• -> L1• + H2O",
    T0=140,
    k0=2.0e6, # guessed from kbb2
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kbb2 = reactionRateDB(
    fingerprint="L2H + HO• -> L2• + H2O",
    T0=140,
    k0=1.0e7,
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kbb3 = reactionRateDB(
    fingerprint="L3H + HO• -> L3• + H2O",
    T0=140,
    k0=1.0e7, # guessed from kbb2
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R16-18 (O2 addition) in Touffet et al., 2023
ka1 = reactionRateDB(
    fingerprint="L1• + O2 -> L1OO•",
    T0=140,
    k0=1.0e6,
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ka2 = reactionRateDB(
    fingerprint="L2• + O2 -> L2OO•",
    T0=140,
    k0=1.0e6,
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ka3 = reactionRateDB(
    fingerprint="L3• + O2 -> L3OO•",
    T0=140,
    k0=1.0e6,
    k0_CI=0,
    Ea=0,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R19-27 (H-abstraction, only self-terms) in Touffet et al., 2023
kp1 = reactionRateDB(
    fingerprint="L1H + L1OO• -> L1• + L1OOH",
    T0=140,
    k0=18.9e-2,
    k0_CI=0,
    Ea=47.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kp2 = reactionRateDB(
    fingerprint="L2H + L2OO• -> L2• + L2OOH",
    T0=140,
    k0=42.2e-2,
    k0_CI=0,
    Ea=28.0e3,
    Ea_CI=None,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)


kp3 = reactionRateDB(
    fingerprint="L3H + L3OO• -> L3• + L3OOH",
    T0=140,
    k0=84.3e-2,
    k0_CI=0,
    Ea=23.8e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R28-30 (Beta-scissions) in Touffet et al., 2023
kB1 = reactionRateDB(
    fingerprint="L1O• -> L1=O",
    T0=140,
    k0=2.28e6, # value kB2
    k0_CI=0,
    kunits = "s⁻¹",
    Ea=47.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kB2 = reactionRateDB(
    fingerprint="L2O• -> L2=O",
    T0=140,
    k0=2.28e6,
    k0_CI=0,
    kunits = "s⁻¹",
    Ea=47.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kB3 = reactionRateDB(
    fingerprint="L3O• -> L3=O",
    T0=140,
    k0=2.28e6, # value kB2
    k0_CI=0,
    kunits = "s⁻¹",
    Ea=47.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R31-36 (termination products via hydropoxyl radicals, self terms) in Touffet et al., 2023
kt1 = reactionRateDB(
    fingerprint="L1OO• + L1OO• -> L1-polymer + L1-polymer",
    T0=140,
    k0=4.00e5,
    k0_CI=0.0,
    Ea=0.0e3,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.9e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kt2 = reactionRateDB(
    fingerprint="L2OO• + L2OO• -> L2-polymer + L2-polymer",
    T0=140,
    k0=1.5e6, # value kB2
    k0_CI=0.0,
    Ea=0.0e3,
    Ea_CI=0,
    diffusionLimited = True,
    eta0 = 1.6e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kt3 = reactionRateDB(
    fingerprint="L3OO• + L3OO• -> L3-polymer + L3-polymer",
    T0=140,
    k0=1.5e7, # value kB2
    k0_CI=0.0,
    Ea=0.0e3,
    Ea_CI=0,
    diffusionLimited = True,
    eta0 = 1.5e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R37-42 (termination products via alkyl radicals, self terms) in Touffet et al., 2023
ktt1 = reactionRateDB(
    fingerprint="L1• + L1• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5, # from ktt2
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.9e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ktt2 = reactionRateDB(
    fingerprint="L2• + L2OO• -> L2-polymer + L2-polymer",
    T0=140,
    k0=1.0e5,
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.6e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ktt3 = reactionRateDB(
    fingerprint="L3• + L3OO• -> L3-polymer + L3-polymer",
    T0=140,
    k0=1.0e5, # from ktt3
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.5e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R43-51 (termination products via alkyl+alcoxyl radicals) in Touffet et al., 2023
kttt1 = reactionRateDB(
    fingerprint="L1• + L1OO• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5, # from kttt2
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.9e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kttt2 = reactionRateDB(
    fingerprint="L2• + L2• -> L2-polymer + L2-polymer",
    T0=140,
    k0=1.0e5,
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.6e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kttt3 = reactionRateDB(
    fingerprint="L3• + L3• -> L3-polymer + L3-polymer",
    T0=140,
    k0=1.0e5, # from kttt3
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    diffusionLimited = True,
    eta0 = 1.5e-3,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)


"""
Entries from Table 1 in
    Maxime Touffet, Olivier Vitrac,
    Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms,
    Food Chemistry, 481, 2025, 143952
    https://doi.org/10.1016/j.foodchem.2025.143952

Note that the values defined here are automatically registered at instantiation
"""
kdbi1 = reactionRateDB(
    fingerprint="L1OOH + L1OOH -> L1O• + L1OO•",
    T0=40,
    k0=4e-10,
    k0_CI=0.0,
    Ea=100e3,
    Ea_CI=0.0,
    source="Table 1 in https://doi.org/10.1016/j.foodchem.2025.143952"
)

kdbi2 = reactionRateDB(
    fingerprint="L2OOH + L2OOH -> L2O• + L2OO•",
    T0=40,
    k0=4.5e-9,
    k0_CI=0.0,
    Ea=110e3,
    Ea_CI=0.0,
    source="Table 1 in https://doi.org/10.1016/j.foodchem.2025.143952"
)

kdbi3 = reactionRateDB(
    fingerprint="L3OOH + L3OOH -> L3O• + L3OO•",
    T0=140,
    k0=14e-9,
    k0_CI=0.0,
    kunits = "s⁻¹",
    Ea=110e3,
    Ea_CI=0.0,
    source="Table 1 in https://doi.org/10.1016/j.foodchem.2025.143952"
)

# %% For testing and debugging
if __name__ == '__main__':

    # database
    ktable = reactionRateDB.to_dataframe()
    print(ktable)

    # Oxidation of a FAME mixture
    oil = mixture()
    oil.add("L1H",concentration=3000)
    oil.add("L2H",concentration=1000)
    oil.add("L3H",concentration=500)
    oil.add("L1OOH",concentration=100)
    oil.add("O2",concentration=10)
    oil.addProducts()
    oil.addReactions()
    oil.populateReactionRates()
    oil.reactions
    polar = oil.lumped_polar_compounds()
    oilmodel = mixtureKinetics(oil) # kinetic model
    oilmodel.solve(10*24*3600,60)
    oilmodel.plot()
    df = oilmodel.results_as_dataframe(["L1H","L2H","L3H","L1OOH","L2OOH","L3OOH"])
    print(df)

    # low-level examples
    O2 = oxygen()
    mono = species.create('C1H', concentration=1.0)
    print(mono)
    mono.isradical
    False
    alk = species.create('C1', concentration=0.01)
    ox = species.create('O2')
    rxn = alk + ox  # placeholder for reaction
    lump = mono | alk
    print(lump.name)
    P1OOH = species.create('P1OOH',concentration=300)

    # mixture
    M = mixture(name="mymixture", description="this is my mixture")
    M.add("P1OOH", concentration=0.2)
    M.add("L2OOH", concentration=0.5)

    print(M["P1OOH"])     # Access by name
    print(M[1])           # Access by index
    print(M.L2OOH)        # Access by attribute
    print(len(M))         # → 2

    # Submixture
    sub = M[[0, 1]]
    print(sub)

    # Deletion
    M["P1OOH"] = []
    print(len(M))         # → 1

    # Replacement
    from copy import deepcopy
    new_sp = deepcopy(M.L2OOH)
    new_sp.name = "L2OOH_new"
    M["L2OOH"] = new_sp

    # Iteration and display
    for sp in M:
        print(sp)

    print(repr(M))
    print(str(M))
