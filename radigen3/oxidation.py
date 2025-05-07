"""
Module: radigen3.oxidation.py

Core module for simulating oxidation reactions in complex mixtures.

This kernel supports the generative construction and simulation of radical-driven oxidation mechanisms in systems such as edible oils, fuels, and polymers. The underlying approach is modular and combinatorial, enabling flexibility across classes of reactive functions and substrates.

Current status: functional prototype under active development.

Typical usage example (oxidation of fatty acid methyl esters, FAMEs):

    oil = mixture()
    oil.add("L1H", concentration=3000)
    oil.add("L1OOH", concentration=100)
    oil.add("O2", concentration=10)
    oil.addProducts()
    oil.addReactions()
    oil.populateReactionRates()

    # ouput
    oil.reactions yield (* means parameterized reaction):

        [*R0: L1H + HO• -> L1• + H2O,
         *R1: L1H + L1O• -> L1• + L1OH,
         *R2: L1H + L1OO• -> L1• + L1OOH,
         *R3: L1OOH -> L1O• + HO•,
         *R4: L1OOH + L1OOH -> L1O• + L1OO•,
         *R5: L1• + O2 -> L1OO•,
         *R6: L1• + L1• -> L1-polymer + L1-polymer,
         *R7: L1• + L1OO• -> L1-polymer + L1-polymer,
         *R8: L1O• -> L1=O,
         *R9: L1OO• + L1OO• -> L1-polymer + L1-polymer]

Implemented features:
    - Representation of reactive species with structural and kinetic attributes (✔)
    - Combinatorial construction of reaction networks (✔)
    - Parameterization of reaction rate constants for self- and cross-reactions (➤ in progress)
    - Diffusion-limited reaction rate modeling (✘ pending)
    - Variable oxygenation and temperature coupling (✘ pending)
    - Stoichiometric matrix construction and kinetic ODE integration (✘ pending)
    - Support for dynamic simulation conditions (✘ future scope)


This work builds upon the following two publications, which respectively address the mechanistic modeling of autoxidation kinetics and the predictive simulation of FAME mixture oxidizability:

[2] Touffet M., Vitrac O., Temperature-dependent kinetics of unsaturated fatty acid methyl esters: Modeling autoxidation mechanisms, Food Chemistry, 481 (2025), 143952. https://doi.org/10.1016/j.foodchem.2025.143952
[1] Touffet M., Smith P., Vitrac O., A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures, Food Research International, 173(1) (2023), 113289. https://doi.org/10.1016/j.foodres.2023.113289


This module is part of the Generative Simulation Initiative.
For questions or contributions, contact: olivier.vitrac@gmail.com

Revision: 2025-05-07
"""


# %% Dependencies
import re, math
import pandas as pd
from collections import defaultdict

# %% Reaction Class
class reaction:
    def __init__(self, *, A, C, B=None, D=None,
                 index=None, k0=1.0, Ea=45.0e3, T0=25.0,
                 Tunits = "°C", kunits = "m³/mol/s", Eaunits = "J/mol",
                 model="default"):
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
        self.index = index
        self.k0 = k0
        self.kunits = kunits
        self.Ea = Ea
        self.Eaunits = Eaunits
        self.T0 = T0
        self.Tunits = Tunits
        self.model = model
        self.assigned = False

    @property
    def isselfroot(self):
        """Returns True if the reaction is bimolecular and involves same roots (e.g., L1OOH+L1OOH)"""
        return self.A.allylic == self.B.allylic if self.B is not None else True

    @property
    def iscrossroot(self):
        """Returns True if the reaction is monomolecular or involves different roots (e.g., L1OOH+L2OOH)"""
        return self.A.allylic != self.B.allylic if self.B is not None else False

    @property
    def ismonomolecular(self):
        """Returns True if the reaction is monomolecular"""
        return self.B is None

    @property
    def isbimolecular(self):
        """Returns True if the reaction is bimolecular"""
        return self.B is not None

    @property
    def k(self):
        """
        Returns a function k(T) giving the Arrhenius rate constant at temperature T (°C).
        """
        def arrhenius(T):
            R = 8.314  # J/mol·K
            T_kelvin = 273.15 + T
            T0_kelvin = 273.15 + self.T0
            exponent = -self.Ea / R * (1/T_kelvin - 1/T0_kelvin)
            return self.k0 * math.exp(exponent)
        return arrhenius

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
            f'k(T0={self.T0:.1f} [{T0_units}])': f"{self.k(self.T0):.3g} [{k0_units}]",
            'Ea': f"{self.Ea:.3g} [{Ea_units}]",
            'model': self.model,
        }
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
            if A_rank > B_rank:
                return B, A, D, C

        return A, B, C, D


# %% Base Mixture Class
class mixture:
    """
    Represents a mixture of chemical species, such as those found in a reacting system.

    Attributes:
        name (str): The name of the mixture.
        description (str): Descriptive text.
        substances (list): List of species instances, each with an assigned index.
    """
    def __init__(self, name="mixture", description=""):
        self.name = name
        self.description = description
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
        core_attrs = {"name", "description", "_substances", "_indexmap", "_namemap"}
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
        rhash = reaction.hash_reagents(A, B)
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

                if len(Aps) == len(Bps) and len(Aps) > 1 and type(A) == type(B):
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
            fp = rxn.fingerprint
            entry = reactionRateDB.get_latest(fp)
            if entry:
                rxn.k0 = entry.k0
                rxn.Ea = entry.Ea
                rxn.T0 = entry.T0
                rxn.kunits = entry.kunits
                rxn.Eaunits = entry.Eaunits
                rxn.Tunits = entry.Tunits
                rxn.model = "values from "+entry.source
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

# %% Base Species Class
class species:
    """
    Base class for chemical species in radical oxidation modeling.

    This class serves as the root of a registry-based class hierarchy.
    Each subclass of `species` corresponds to a specific chemical entity,
    such as a radical, hydroperoxide, or stable product.

    Attributes:
        name (str): Long name of the species.
        shortname (str): Short name (alias or abbreviation).
        root (str): Common root to track species family.
        concentration (float): Current concentration of the species.
        index (int or None): Optional index identifier.

    Class Attributes:
        _registry (dict): Maps class names and aliases to class objects.

    Methods:
        __iter__(): Iterates over the main attributes as key-value pairs.
        __repr__(): Pretty-prints internal state.
        __str__(): Returns a short description.
        __or__(other): Defines `|` as lumping operator (returns a lumped object).
        __add__(other): Placeholder to represent reaction with `+`.
        register(cls): Class decorator to register species and their aliases.
        create(name, **kwargs): Factory method using name or alias.

    Properties:
        isradical: True if the species has free valence.
        isradical_onC: True if radical is centered on C.
        isradical_onO: True if radical is centered on O.
        ishydroperoxide: True if reactiveFunction is "COOH".
        isreactive: True if the species has defined `reactWith`.
        islumped: True if the species is a `lumped` subclass.

    Example:
        >>> mono = species.create('C1H', concentration=1.0)
        >>> print(mono)
        <L1H - index=None> of class monoAllylicCH with conc.=1.0
        >>> mono.isradical
        False
        >>> alk = species.create('C1', concentration=0.01)
        >>> ox = species.create('O2')
        >>> rxn = alk + ox  # placeholder for a reaction
        >>> lump = mono | alk  # lumping operation
        >>> print(lump.name)
        lump_L1H_L1•
    """
    _registry = {}

    def __init__(self, *, root=None, name='', shortname='', concentration=0.0, index=None, _alias_used=None):
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
        return self.freevalence > 0 if hasattr(self, "freevalence") else False

    @property
    def isradical_onC(self):
        """True if radical center is on carbon."""
        return self.freevalenceHolder == "C" if self.freevalence else False

    @property
    def isradical_onO(self):
        """True if radical center is on oxygen."""
        return self.freevalenceHolder == "O" if self.freevalence else False

    @property
    def ishydroperoxide(self):
        """True if reactiveFunction is 'COOH', interpreted as a hydroperoxide."""
        return self.reactiveFunction == "COOH"

    @property
    def isreactive(self):
        """True if the species leads to at least one product."""
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


# %% Lumped Class
@species.register
class lumped(species):
    def __init__(self, species_list):
        self.species = species_list
        self.name = "lump_" + "_".join(sp.name for sp in species_list)
        self.index = None
        self.concentration = sum(sp.concentration for sp in species_list)
        self.shortname = "LM"
        self.root = None



# %% Specific Chemical Classes derived from species

# _redoxRank is used to which reagent is A and B in A + B -> C + D to get a unique representation of each reaction
# enforced rule: _redoxRank(B) > _redoxRank(A)

# | Species Class   | Description                 | Allylicity | Role in Redox      | Suggested `_redoxRank` |
# | --------------- | --------------------------- | ---------- | ------------------ | ---------------------- |
# | `monoAllylicCH` | Primary RH donor (labile H) | 1          | Mild reductant     | 3                      |
# | `diAllylicCH`   | Slightly more reactive      | 2          | Mild reductant     | 2.5                    |
# | `triAllylicCH`  | Very labile H               | 3          | Stronger reductant | 2                      |

# | Species         | Allylic context | OS estimate (central C) | Justification                                  |
# | --------------- | --------------- | ----------------------- | ---------------------------------------------- |
# | `monoAllylicCH` | C=C–CH          | –1                      | Single allylic C–H                             |
# | `diAllylicCH`   | C=C–CH–C=C      | \~–0.5 to –1            | More resonance, slightly more electron-rich    |
# | `triAllylicCH`  | C=C–CH–C=C–C=C  | \~0                     | High delocalization, lower effective reduction |

# | Species Type           | Central Atom | Oxidation State | Redox Rank  |
# | ---------------------- | ------------ | --------------- | ----------- |
# | alkyl radical (C•)     | C            | \~0             | 4           |
# | hydroperoxide (ROOH)   | O            | –1              | 6           |
# | alkoxyl radical (RO•)  | O            | –1              | 5           |
# | peroxyl radical (ROO•) | O            | \~0             | 7           |
# | dioxygen (O₂)          | O            | 0               | 8 (oxidant) |
# | alcohol (ROH)          | C/O          | \~0             | 3           |
# | aldehyde/ketone        | C            | +1/+2           | 2           |
# | carboxylic acid        | C            | +3              | 1           |

"""
    Derived species classes
        reactiveFunction and reactWith determine the possible reactions between substances/species
        None is used for monomolecular reactions (decomposition)
        product lists products as:
            - regular strings → single products.
            - set → unordered products (used for homolytic monomolecular decompositions),
            - tuple → ordered products (used for bimolecular reactions)
        Use name, root, suffix to name substance/species.
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
    reactiveFunction = None
    freevalence = 0
    freevalenceHolder = None
    reactWith = None
    product =None
    interpretation = """stable termination products (polymers)"""
    oxidationState = 0
    _redoxRank = None

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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

@species.register
class peroxyl(species):
    """Peroxyl radical (OH•) class"""
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

# %% Base Reaction Rate Constants Class
class reactionRateDB:
    _registry = defaultdict(list)

    def __init__(self, *, fingerprint, T0, k0, k0_CI=None, Ea, Ea_CI=None,
                 Tunits = "°C", kunits = "m³/mol/s", Eaunits = "J/mol", source=""):
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

        # Register in global registry
        self.__class__._registry[fingerprint].append(self)

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
        """Returns all entries matching the fingerprint."""
        return cls._registry.get(fingerprint, [])

    @classmethod
    def get_latest(cls, fingerprint):
        """Returns the last entry added for a given fingerprint."""
        entries = cls._registry.get(fingerprint, [])
        return entries[-1] if entries else None

    @classmethod
    def to_dataframe(cls):
        return pd.DataFrame([
            {
                "reaction": r.fingerprint,
                "T0 (°C)": r.T0,
                "k0 (SI units)": r.k0,
                "k0_CI (SI units)": r.k0_CI,
                f"Ea ({r.Eaunits})": r.Ea,
                f"Ea_CI ({r.Eaunits})": r.Ea_CI,
            }
            for rgroup in cls._registry.values()
            for r in rgroup
        ])


# %% Specific Reaction Rate Constants
"""
Entries from Table 2 in
    Maxime Touffet, Paul Smith, Olivier Vitrac,
    A comprehensive two-scale model for predicting the oxidizability of fatty acid methyl ester mixtures,
    Food Research International, 173(1),2023,113289
    https://doi.org/10.1016/j.foodres.2023.113289.

Note that the values defined here are automatically registered at instantiation
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
    k0=1.0e7, # guessed from kbb3
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
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kt2 = reactionRateDB(
    fingerprint="L2OO• + L2OO• -> L2-polymer + L2-polymer",
    T0=140,
    k0=1.5e6, # value kB2
    k0_CI=0.0,
    Ea=0.0e3,
    Ea_CI=0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kt3 = reactionRateDB(
    fingerprint="L3OO• + L3OO• -> L3-polymer + L3-polymer",
    T0=140,
    k0=1.5e7, # value kB2
    k0_CI=0.0,
    Ea=0.0e3,
    Ea_CI=0,
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
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ktt2 = reactionRateDB(
    fingerprint="L1• + L1OO• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5,
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

ktt3 = reactionRateDB(
    fingerprint="L2• + L2OO• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5, # from ktt3
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

# R43-51 (termination products via alkyl+alcoxyl radicals) in Touffet et al., 2023
kttt1 = reactionRateDB(
    fingerprint="L3• + L3OO• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5, # from kttt2
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kttt2 = reactionRateDB(
    fingerprint="L1• + L1• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5,
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
    source="Table 2 in https://doi.org/10.1016/j.foodres.2023.113289"
)

kttt3 = reactionRateDB(
    fingerprint="L1• + L1• -> L1-polymer + L1-polymer",
    T0=140,
    k0=1.0e5, # from kttt3
    k0_CI=0.0,
    Ea=0.0,
    Ea_CI=0.0,
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

    # basic FAME oxidation
    oil = mixture()
    oil.add("L1H",concentration=3000)
    oil.add("L1OOH",concentration=100)
    oil.add("O2",concentration=10)
    oil.addProducts()
    oil.addReactions()
    oil.populateReactionRates()
    oil.reactions

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
