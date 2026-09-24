import math
import itertools
try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False

# ============================================================
#  1. Exact Database: Ionic Radii from User Images (Angstroms)
# ============================================================
IONIC_RADII = {
    # Rock Salt / Spinel A-site / Perovskite A-site
    "Mg": 0.72,
    "Ni": 0.69,
    "Co": 0.74,
    "Zn": 0.74,  # VI base (0.88 for IV is noted, but standard calculation uses main val)
    "Cu": 0.73,
    "La": 1.032,
    "Ce": 1.01,  # Ce3+
    "Pr": 0.99,
    "Nd": 0.983,
    "Sm": 0.958,
    "Gd": 0.938,
    "Y":  0.90,
    "Ca": 1.00,
    "Sr": 1.18,
    "Ba": 1.35,
    "Pb": 1.20,
    
    # Spinel B-site / Pyrochlore B-site / Fluorite
    "Al": 0.535,
    "Cr": 0.615,
    "Fe": 0.645, # Fe3+
    "Ga": 0.620,
    "Ti": 0.605,
    "Zr": 0.720,
    "Hf": 0.710,
    "Sn": 0.690,
    "Nb": 0.640,
    "Th": 1.050,
    "Mn": 0.645  # Mn3+
}

# ============================================================
#  2. Structure Configurations & Constants
# ============================================================
STRUCTURES = {
    "1": {
        "name": "Rock Salt",
        "formula": "AO",
        "C_A": 1, "C_B": 0, "C_O": 1,
        "site_A_opts": ["Mg", "Ni", "Co", "Zn", "Cu", "La", "Ce", "Pr", "Nd", "Sm"],
        "site_B_opts": [],
        "n_Omega_A": 1.0,
        "n_Omega_B": 0.0,
    },
    "2": {
        "name": "Spinel",
        "formula": "AB2O4",
        "C_A": 1, "C_B": 2, "C_O": 4,
        "site_A_opts": ["Mg", "Ni", "Co", "Zn", "Cu"],
        "site_B_opts": ["Al", "Cr", "Fe", "Ga", "Ti"],
        "n_Omega_A": 1/24,
        "n_Omega_B": 1/2,
    },
    "3": {
        "name": "Pyrochlore",
        "formula": "A2B2O7",
        "C_A": 2, "C_B": 2, "C_O": 7,
        "site_A_opts": ["Y", "Gd", "La", "Nd", "Sm", "Ce", "Hf"],
        "site_B_opts": ["Fe", "Cr", "Ti", "Zr", "Sn", "Nb"],
        "n_Omega_A": 1/4,
        "n_Omega_B": math.sqrt(2)/12,
    }
}
R = 8.314

# ============================================================
#  3. Physics Calculators
# ============================================================
def calc_site_params(els, fracs):
    if not els: return 0.0, 0.0
    r_bar = sum(fracs[i] * IONIC_RADII[els[i]] for i in range(len(els)))
    delta = math.sqrt(sum(fracs[i] * (1 - IONIC_RADII[els[i]] / r_bar)**2 for i in range(len(els)))) * 100
    return r_bar, delta

def calc_site_entropy(fracs):
    if not fracs: return 0.0
    return sum(x * math.log(x) for x in fracs if x > 0)

# ============================================================
#  4. Generators & Evaluators
# ============================================================
def generate_compositions(step, num_elements):
    if num_elements == 0: return [()]
    if num_elements == 1: return [(1.0,)]
        
    steps = int(round(1.0 / step))
    combos = []
    for c in itertools.combinations_with_replacement(range(steps + 1), num_elements):
        if sum(c) == steps:
            for p in set(itertools.permutations(c)):
                if all(val > 0 for val in p): 
                    combos.append(tuple(val * step for val in p))
    return combos

def evaluate_composition(struct, els_A, fracs_A, els_B, fracs_B):
    r_A, delta_A = calc_site_params(els_A, fracs_A)
    r_B, delta_B = calc_site_params(els_B, fracs_B) if els_B else (0.0, 0.0)

    # 1. Delta_r_star
    delta_r_star = math.sqrt(delta_A**2 + delta_B**2)
    
    # 2. Delta_r_N
    term_A = struct['n_Omega_A'] * (delta_A**2)
    term_B = struct['n_Omega_B'] * (delta_B**2)
    delta_r_N = math.sqrt(term_A + term_B)

    # 3. Delta_S_mix_N
    C_total = struct['C_A'] + struct['C_B'] + struct['C_O']
    S_A_term = (struct['C_A'] / C_total) * calc_site_entropy(fracs_A)
    S_B_term = (struct['C_B'] / C_total) * calc_site_entropy(fracs_B)
    Delta_S_mix_N = -R * (S_A_term + S_B_term)

    # 4. Criteria Check
    passed = False
    if struct['name'] == "Rock Salt":
        if (6.5 <= Delta_S_mix_N <= 8.5) and (2 <= delta_r_N <= 8) and (delta_r_star < 8):
            passed = True
    elif struct['name'] == "Spinel":
        if (4.5 <= Delta_S_mix_N <= 6.75) and (0 <= delta_r_N <= 12.5) and (delta_r_star < 8):
            passed = True
    elif struct['name'] == "Pyrochlore":
        ratio = (r_A / r_B) if r_B > 0 else 0
        if (1.75 <= Delta_S_mix_N <= 4.0) and (0 <= delta_r_N <= 18) and (delta_r_star < 5) and (ratio > 1.46):
            passed = True

    return passed, {"dS_N": Delta_S_mix_N, "dr_N": delta_r_N, "dr_star": delta_r_star, "r_A": r_A, "r_B": r_B}

# ============================================================
#  5. Main Application
# ============================================================
def main():
    print("=" * 60)
    print("   HEO DESIGNER - ACCURATE DATABASE (User Images)")
    print("=" * 60)
    print(" Select Target Structure:")
    print("   1. Rock Salt   (AO)")
    print("   2. Spinel      (AB2O4)")
    print("   3. Pyrochlore  (A2B2O7)")
    
    choice = input("\n Enter number (1/2/3): ").strip()
    if choice not in STRUCTURES:
        print(" ❌ Invalid choice."); return
        
    struct = STRUCTURES[choice]
    print(f"\n 🔸 Selected: {struct['name']} [{struct['formula']}]")
    
    # --- GET SITE A ---
    print(f"\n [SITE A] Allowed: {', '.join(struct['site_A_opts'])}")
    raw_A = input(" Enter elements for Site A (space-separated): ").strip().split()
    els_A = [e.capitalize() for e in raw_A if e.capitalize() in struct['site_A_opts']]
    if not els_A:
        print(" ❌ No valid elements entered."); return

    # --- GET SITE B ---
    els_B = []
    if struct['site_B_opts']:
        print(f"\n [SITE B] Allowed: {', '.join(struct['site_B_opts'])}")
        raw_B = input(" Enter elements for Site B (space-separated): ").strip().split()
        els_B = [e.capitalize() for e in raw_B if e.capitalize() in struct['site_B_opts']]
        if not els_B:
            print(" ❌ No valid elements entered."); return

    # --- Scanning Step Setup ---
    total_elements = len(els_A) + len(els_B)
    step = 0.05
        
    print(f"\n ⚙️ Searching combinations (Step = {step*100:.0f}%)...")
    
    combos_A = generate_compositions(step, len(els_A))
    combos_B = generate_compositions(step, len(els_B))
    
    valid_compositions = []
    total_checked = 0
    
    for fA in combos_A:
        for fB in combos_B:
            total_checked += 1
            passed, params = evaluate_composition(struct, els_A, fA, els_B, fB)
            if passed:
                valid_compositions.append({"fA": fA, "fB": fB, "params": params})

    print(f" 📊 Checked {total_checked} possibilities.")
    
    if not valid_compositions:
        print("\n ❌ No composition found that meets ALL criteria for a HIGH-ENTROPY single phase.")
        print(" 💡 This structure might need more elements to increase entropy, or the size differences are too large.")
        return
        
    print(f" ✅ Found {len(valid_compositions)} valid combinations.")
    
    # --- PRINT RANGES ---
    print("\n" + "="*60)
    print("   CONCENTRATION RANGES FOR SINGLE-PHASE FORMATION")
    print("="*60)
    
    print("\n [ SITE A RANGES ]")
    for i, el in enumerate(els_A):
        fracs = [item["fA"][i] * 100 for item in valid_compositions]
        print(f"   ▶ {el:4}: {min(fracs):5.1f}%  to  {max(fracs):5.1f}%")
        
    if els_B:
        print("\n [ SITE B RANGES ]")
        for i, el in enumerate(els_B):
            fracs = [item["fB"][i] * 100 for item in valid_compositions]
            print(f"   ▶ {el:4}: {min(fracs):5.1f}%  to  {max(fracs):5.1f}%")

    # --- Print an Example ---
    best = valid_compositions[len(valid_compositions)//2]
    print("\n" + "="*60)
    print("   EXAMPLE OF A VALID COMPOSITION")
    print("="*60)
    
    str_A = " ".join([f"{el}{best['fA'][i]*100:.0f}" for i, el in enumerate(els_A)])
    str_B = " ".join([f"{el}{best['fB'][i]*100:.0f}" for i, el in enumerate(els_B)])
    
    if struct['name'] == "Rock Salt": print(f" Composition: ({str_A})O")
    elif struct['name'] == "Spinel": print(f" Composition: ({str_A})({str_B})2O4")
    elif struct['name'] == "Pyrochlore": print(f" Composition: ({str_A})2({str_B})2O7")
        
    p = best["params"]
    print(f"\n   ΔS_mix^N : {p['dS_N']:5.2f} J/mol·K")
    print(f"   δ_r^N    : {p['dr_N']:5.2f}%")
    print(f"   δ_r^*    : {p['dr_star']:5.2f}%")
    if struct['name'] == "Pyrochlore":
        print(f"   r_A / r_B: {p['r_A']/p['r_B']:5.3f}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()