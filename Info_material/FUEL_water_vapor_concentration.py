"""
Water vapor concentration calculator (Spyder version)

Estimate mole and mass fraction of water vapor in the products of
lean combustion (φ < 1) for hydrocarbon or oxygenated fuels.
"""

from dataclasses import dataclass
import re

# Constants
MW = {
    "H2O": 18.01528,
    "CO2": 44.0095,
    "O2": 32.0,
    "N2": 28.0134,
}
AIR_RATIO = 3.76  # mol N2 per mol O2

@dataclass(frozen=True)
class Fuel:
    name: str
    formula: str

    def atoms(self):
        """Return (x, y, z) for CxHyOz"""
        m = re.fullmatch(r"C(\d*)H(\d*)(?:O(\d*))?", self.formula, re.IGNORECASE)
        if not m:
            raise ValueError(f"Invalid formula: {self.formula}")
        x = int(m.group(1)) if m.group(1) else 1
        y = int(m.group(2)) if m.group(2) else 0
        z = int(m.group(3)) if m.group(3) else 0
        return x, y, z

    @property
    def a_st(self):
        x, y, z = self.atoms()
        return x + y / 4 - z / 2

    def water_mole_fraction(self, phi):
        lam = 1 / phi
        if lam <= 1:
            raise ValueError("Only valid for lean mixtures (phi < 1)")
        x, y, _ = self.atoms()
        n_H2O = y / 2
        n_CO2 = x
        n_O2 = (lam - 1) * self.a_st
        n_N2 = AIR_RATIO * lam * self.a_st
        total = n_H2O + n_CO2 + n_O2 + n_N2
        return n_H2O / total

    def water_mass_fraction(self, phi):
        lam = 1 / phi
        x, y, _ = self.atoms()
        n_H2O = y / 2
        n_CO2 = x
        n_O2 = (lam - 1) * self.a_st
        n_N2 = AIR_RATIO * lam * self.a_st
        mass = (
            n_H2O * MW["H2O"]
            + n_CO2 * MW["CO2"]
            + n_O2 * MW["O2"]
            + n_N2 * MW["N2"]
        )
        return n_H2O * MW["H2O"] / mass

# --- Fuel presets (editable) ---
fuel_db = {
    "Jet-A": Fuel("Jet-A", "C12H23"),
    "n-heptane": Fuel("n-Heptane", "C7H16"),
    "butanol": Fuel("n-Butanol", "C4H10O"),
    "HEFA-SPK": Fuel("HEFA-SPK", "C12H26"),
    "FT-1": Fuel("FT-1", "C11H24"),
    "FT-2": Fuel("FT-2", "C13H28"),
}

# === USER INPUTS BELOW THIS LINE ===

phi = 0.24
fuel_name = "Jet-A"  # Choose from: Jet-A, n-heptane, butanol, HEFA-SPK, FT-1, FT-2

# === DO NOT MODIFY BELOW THIS LINE ===

fuel = fuel_db[fuel_name]
X_H2O = fuel.water_mole_fraction(phi)
w_H2O = fuel.water_mass_fraction(phi)
lam = 1 / phi

print(f"Fuel      : {fuel.name} ({fuel.formula})")
print(f"Phi (φ)   : {phi}")
print(f"Lambda λ  : {lam:.2f} (excess-air factor)")
print(f"X_H2O     : {X_H2O:.5f} ({X_H2O*100:.2f} %)")
print(f"w_H2O     : {w_H2O:.5f} ({w_H2O*100:.2f} %)")
