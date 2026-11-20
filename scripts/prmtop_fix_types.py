#!/usr/bin/env python3
"""
prmtop_fix_types.py

Check (and fix) the Amber atom types of a prmtop before it is used by CP2K.

CP2K (FIST and QM/MM with CONN_FILE_FORMAT AMBER) does not use the Amber atom
types the way Amber does. It builds one atomic KIND per (atom type, charge)
pair, and only the names of these kinds are used to tell the atoms apart:

  - the type names are converted to upper case, so types that differ only by
    case (GLYCAM Ho/Oh and ff99SB HO/OH, GAFF ca and ff14SB CA, ...) end up in
    the same group of kinds;
  - when a type has more than one charge, its kinds are named TYPE1, TYPE2,
    ... (one per charge). If the topology also has a type called TYPE1 (e.g.
    the backbone amide H and the aliphatic H1), atoms with different charges
    end up in the same kind, with a single charge. Released versions of CP2K
    (up to 2026.2) do not stop or warn, the electrostatics are simply wrong.
    This was fixed in the development version in October 2026 (cp2k/cp2k#6191);
  - the element of each kind is guessed from its name (trailing digits
    removed) using built-in lists of Amber, CHARMM and GROMOS types and the
    periodic table, unless a &KIND section gives the ELEMENT. Some names are
    guessed wrong (the Amber Ca2+ type C0 becomes carbon) or not at all (Na+,
    Cl-, many GAFF types), and the run stops with "Unknown element for KIND";
  - the AMBER_ATOM_TYPE section is read by splitting at whitespace, so a type
    that fills its 4 characters (e.g. Ca2+) is glued to the next one.

The script reproduces the kind names that CP2K builds, reports these problems
and, unless --dry-run is given, writes a new prmtop in which the offending
types are renamed. New names are taken from the type lists that CP2K already
knows (e.g. H1 -> HZ, Oh -> OA, C0 -> CAL), so that CP2K guesses the right
element without extra &KIND sections. Only the AMBER_ATOM_TYPE labels change:
the force-field parameters are indexed by number in the prmtop, so Amber
programs give exactly the same energy with both files.

Types whose element CP2K cannot guess at all (e.g. Na+ and Cl-) are reported
but not renamed (unless --rename-unknown is given): give them a &KIND section
with ELEMENT in the CP2K input, as in the inputs of the tutorial.

Usage examples:
    python3 prmtop_fix_types.py -p amylase_prep.prmtop --dry-run
    python3 prmtop_fix_types.py -p amylase_prep.prmtop -o amylase.prmtop
    python3 prmtop_fix_types.py -p amylase_prep.prmtop -o amylase.prmtop --map H1:HL Ho:HOL
    python3 prmtop_fix_types.py -p amylase_ee.prmtop --dry-run     # check another topology
"""

import argparse
import sys
from collections import OrderedDict, defaultdict

try:
    import parmed
    from parmed.periodic_table import Element as PT_ELEMENT
except ImportError as exc:  # pragma: no cover
    sys.exit(f"ERROR: {exc}. Install ParmEd with: pip install parmed (it also comes with AmberTools)")

VERSION = "2.0"

# Type -> element maps built into CP2K (src/mm_mapping_library.F). Lookup order: AMBER, CHARMM, GROMOS.
AMBER_MAP = [
    ("BR", "Br"), ("C", "C"), ("CA", "C"), ("CB", "C"), ("CC", "C"), ("CD", "C"), ("CI", "C"), ("CK", "C"),
    ("CP", "C"), ("CM", "C"), ("CS", "C"), ("CN", "C"), ("CQ", "C"), ("CR", "C"), ("CT", "C"), ("CV", "C"),
    ("CW", "C"), ("C*", "C"), ("CX", "C"), ("CY", "C"), ("CZ", "C"), ("CO", "C"), ("TG", "C"), ("2C", "C"),
    ("3C", "C"), ("C0", "Ca"), ("F", "F"), ("H", "H"), ("HC", "H"), ("HA", "H"), ("HO", "H"), ("HS", "H"),
    ("HW", "H"), ("HP", "H"), ("HZ", "H"), ("I", "I"), ("IM", "Cl"), ("IP", "Na"), ("IB", "Na"), ("MG", "Mg"),
    ("N", "N"), ("NA", "N"), ("NB", "N"), ("NC", "N"), ("N*", "N"), ("NT", "N"), ("NY", "N"), ("ND", "N"),
    ("NL", "N"), ("TN", "N"), ("OA", "O"), ("O", "O"), ("OW", "O"), ("OH", "O"), ("OD", "O"), ("OS", "O"),
    ("OP", "O"), ("P", "P"), ("S", "S"), ("SH", "S"), ("CU", "Cu"), ("FE", "Fe"), ("Li", "Li"), ("K", "K"),
    ("Rb", "Rb"), ("Cs", "Cs")]
CHARMM_MAP = [
    ("CA", "C"), ("CAL", "Ca"), ("C", "C"), ("CC", "C"), ("CD", "C"), ("CE", "C"), ("CEL", "C"), ("CES", "Cs"),
    ("CLA", "Cl"), ("CL", "C"), ("CM", "C"), ("CN1A", "C"), ("CN1T", "C"), ("CN3A", "C"), ("CN3B", "C"),
    ("CN3C", "C"), ("CN3D", "C"), ("CN3T", "C"), ("CN5G", "C"), ("CN7B", "C"), ("CN7C", "C"), ("CN7D", "C"),
    ("CN8B", "C"), ("CNA", "C"), ("CN", "C"), ("CNE", "C"), ("CPA", "C"), ("CPB", "C"), ("CP", "C"),
    ("CPH", "C"), ("CPM", "C"), ("CPT", "C"), ("CS", "C"), ("CT", "C"), ("CTL", "C"), ("CY", "C"),
    ("DUM", "H"), ("FE", "Fe"), ("FNA", "F"), ("FN", "F"), ("HA", "H"), ("HAL", "H"), ("HB", "H"),
    ("HC", "H"), ("HCL", "H"), ("HE", "H"), ("HE", "He"), ("HEL", "H"), ("H", "H"), ("HL", "H"),
    ("HN3B", "H"), ("HN3C", "H"), ("HNE", "H"), ("HN", "H"), ("HNP", "H"), ("HOL", "H"), ("HP", "H"),
    ("HR", "H"), ("HS", "H"), ("HT", "H"), ("MG", "Mg"), ("NC", "N"), ("NE", "Ne"), ("NH3L", "N"),
    ("NH", "N"), ("N", "N"), ("NN1C", "N"), ("NN2B", "N"), ("NN2C", "N"), ("NN2G", "N"), ("NN2U", "N"),
    ("NN3A", "N"), ("NN3G", "N"), ("NN3I", "N"), ("NN", "N"), ("NPH", "N"), ("NP", "N"), ("NR", "N"),
    ("NR", "N"), ("NTL", "N"), ("NY", "N"), ("O2L", "O"), ("OBL", "O"), ("OB", "O"), ("OCL", "O"),
    ("OC", "O"), ("OHL", "O"), ("OH", "O"), ("OM", "O"), ("ON1C", "O"), ("ON2b", "O"), ("ON6B", "O"),
    ("ON", "O"), ("O", "O"), ("OSL", "O"), ("OS", "O"), ("OT", "O"), ("P2", "P"), ("PL", "P"), ("POT", "K"),
    ("P", "P"), ("SL", "S"), ("SM", "S"), ("SOD", "Na"), ("S", "S"), ("SS", "S"), ("ZN", "Zn"), ("CN6B", "C")]
GROMOS_MAP = [
    ("H1", "H"), ("CA", "C"), ("HA", "H"), ("SD", "S"), ("OG", "O"), ("CG", "C"), ("HD", "H"), ("SG", "S"),
    ("CZ", "C"), ("OH", "O"), ("C2", "C"), ("H8", "H"), ("O4*", "O"), ("O5*", "O"), ("PA", "P"), ("O1A", "O"),
    ("MG", "Mg"), ("Na", "Na"), ("Cl", "Cl"), ("OW", "O"), ("HW", "H"), ("HW", "H")]
PT_UPPER = {el.upper(): el for el in PT_ELEMENT if el and el != "EP"}
# Names tried first for the most common elements (all are recognized by CP2K); the full lists above follow.
PREFERRED = {
    "H": ["HZ", "HN", "HL", "HR", "HT", "HB", "HOL", "HAL", "HEL", "HNE", "HNP"],
    "C": ["CZ", "CY", "CX", "CV", "CQ", "CI", "CK", "CS", "CO", "CE", "CTL", "CEL", "CNA", "CNE", "CPA"],
    "N": ["NT", "NL", "NY", "ND", "NH", "NP", "NR", "NN", "NC", "NPH", "NTL"],
    "O": ["OA", "OD", "OP", "OB", "OC", "OM", "OT", "ON", "OHL", "OBL", "OCL", "OSL"],
    "S": ["SM", "SL", "SS", "SD", "SG"],
    "Ca": ["CAL"], "Na": ["IP", "SOD"], "Cl": ["IM", "CLA"], "K": ["POT"]}


def parse_args():
    p = argparse.ArgumentParser(
        description="Check the Amber atom types of a prmtop for CP2K and rename the problematic ones.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-p", "--prmtop", required=True, help="Input Amber topology")
    p.add_argument("-o", "--output", default=None, help="Output Amber topology (required unless --dry-run)")
    p.add_argument("--map", nargs="+", default=[], metavar="OLD:NEW",
                   help="Explicit new names, e.g. H1:HL Ho:HOL (overrides the automatic choice)")
    p.add_argument("--keep", nargs="+", default=[], metavar="TYPE",
                   help="Types that must not be renamed, whatever the checks say")
    p.add_argument("--rename-unknown", action="store_true",
                   help="Also rename the types whose element CP2K cannot guess (e.g. Na+ -> IP, Cl- -> IM). "
                        "Remove the matching &KIND sections from the CP2K inputs if you do this")
    p.add_argument("--map-file", default="type_renaming.dat", help="Table of the renamed types (OLD NEW element atoms)")
    p.add_argument("--dry-run", action="store_true", help="Only report the problems and what would be renamed")
    p.add_argument("-q", "--quiet", action="store_true", help="Do not print the list of CP2K kinds per type")
    a = p.parse_args()
    if not a.dry_run and not a.output:
        p.error("-o/--output is required unless --dry-run is given")
    if a.output and a.output == a.prmtop:
        p.error("the output topology must be different from the input one")
    return a


def strip_digits(name):
    return name.rstrip("0123456789")


def cp2k_element(kind_name):
    """Element that CP2K guesses for a kind name of an Amber topology (no &KIND section)."""
    base = strip_digits(kind_name.upper())
    for table in (AMBER_MAP, CHARMM_MAP, GROMOS_MAP):
        for k, el in table:
            if base == k:  # CP2K compares the upper-case name with the list entries as written
                return el
    if len(base) <= 2 and base in PT_UPPER:
        return PT_UPPER[base]
    return None


def element_of(atom):
    return PT_ELEMENT[atom.atomic_number] if atom.atomic_number > 0 else "X"


def analyse(parm):
    """Types, CP2K kind names and the problems found. Returns a dict."""
    charges = parm.parm_data["CHARGE"]
    types = OrderedDict()
    for at in parm.atoms:
        t = types.setdefault(at.type, {"el": element_of(at), "n": 0, "q": set()})
        t["n"] += 1
        t["q"].add(charges[at.idx])
    groups = OrderedDict()
    for name, t in types.items():
        g = groups.setdefault(name.upper(), {"types": [], "q": set()})
        g["types"].append(name)
        g["q"] |= t["q"]
    kinds = OrderedDict()
    owners = defaultdict(set)
    for up, g in groups.items():
        names = [up] if len(g["q"]) == 1 else [f"{up}{k}" for k in range(1, len(g["q"]) + 1)]
        kinds[up] = names
        for kn in names:
            owners[kn].add(up)
    clashes = {kn: sorted(o) for kn, o in owners.items() if len(o) > 1}
    return {"types": types, "groups": groups, "kinds": kinds, "clashes": clashes}


def counter_risk(name, group_uppers):
    """Groups whose counters can produce this type name (TYPE + 1, 2, ...)."""
    up = name.upper()
    out = []
    for i in range(1, len(up)):
        base, digits = up[:i], up[i:]
        if digits.isdigit() and digits[0] != "0" and base in group_uppers and base != up:
            out.append(base)
    return out


def find_problems(info, keep, rename_unknown):
    types, groups = info["types"], info["groups"]
    uppers = set(groups)
    reasons = OrderedDict()
    unknown = OrderedDict()

    def add(t, why):
        if t not in keep:
            reasons.setdefault(t, []).append(why)

    for t in types:
        if len(t) >= 4:
            add(t, "4 characters (glued to the next type when CP2K reads the prmtop)")
        risk = counter_risk(t, uppers)
        if risk:
            add(t, "name can be produced by the kind counter of " + "/".join(risk))
    for up, g in groups.items():
        if len(g["types"]) > 1:
            keeper = up if up in g["types"] else max(g["types"], key=lambda x: types[x]["n"])
            for t in g["types"]:
                if t != keeper:
                    add(t, f"same upper-case name as {keeper}")
    for t, d in types.items():
        guess = cp2k_element(t)
        if guess is None:
            if rename_unknown:
                add(t, "element not recognized by CP2K")
            else:
                unknown[t] = d["el"]
        elif guess != d["el"]:
            add(t, f"CP2K reads it as {guess}, element is {d['el']}")
    return reasons, unknown


def choose_names(info, reasons, explicit):
    types = info["types"]
    taken = {t.upper() for t in types if t not in reasons}
    kept_digit_types = [t for t in types if t not in reasons and t[-1].isdigit()]
    taken |= {n.upper() for n in explicit.values()}
    pool = []
    seen = set()
    for names in PREFERRED.values():
        for ku in names:
            if ku not in seen:
                seen.add(ku)
                pool.append(ku)
    for table in (AMBER_MAP, CHARMM_MAP, GROMOS_MAP):
        for k, _ in table:
            ku = k.upper()
            if ku in seen or len(ku) > 3 or any(c.isdigit() or c in "*+-" for c in ku):
                continue
            seen.add(ku)
            pool.append(ku)
    rename = OrderedDict()
    for t in reasons:
        el = types[t]["el"]
        if t in explicit:
            rename[t] = explicit[t]
            continue
        new = None
        for cand in pool:
            if cand in taken or cp2k_element(cand) != el:
                continue
            # a new group name must not make a kept type (e.g. H5) reachable by its counters
            if any(cand in counter_risk(k, {cand}) for k in kept_digit_types):
                continue
            new = cand
            break
        if new is None:
            sys.exit(f"ERROR: no free CP2K type name for element {el} (type {t}); use --map {t}:NEW "
                     "and a &KIND section with ELEMENT")
        rename[t] = new
        taken.add(new)
    return rename


def report(info, label, quiet):
    types, kinds = info["types"], info["kinds"]
    nkinds = sum(len(v) for v in kinds.values())
    print(f"{label}: {len(types)} atom types, {nkinds} CP2K kinds")
    if not quiet:
        print(f"  {'type':<6s} {'element':<8s} {'atoms':>7s} {'charges':>8s}  CP2K kinds")
        for t, d in types.items():
            names = kinds[t.upper()]
            txt = names[0] if len(names) == 1 else f"{names[0]} ... {names[-1]}"
            print(f"  {t:<6s} {d['el']:<8s} {d['n']:>7d} {len(d['q']):>8d}  {txt}")
    if info["clashes"]:
        txt = ", ".join(f"{k} <- {'/'.join(v)}" for k, v in list(info["clashes"].items())[:10])
        print(f"  kind names shared by different types: {len(info['clashes'])} ({txt})")
    else:
        print("  kind names shared by different types: none")


def main():
    a = parse_args()
    parm = parmed.load_file(a.prmtop)
    if not isinstance(parm, parmed.amber.AmberParm):
        sys.exit("ERROR: the topology is not an Amber prmtop")
    print(f"prmtop_fix_types.py v{VERSION}")

    explicit = OrderedDict()
    for item in a.map:
        if ":" not in item:
            sys.exit(f"ERROR: --map entries must be OLD:NEW, got {item}")
        old, new = item.split(":", 1)
        if len(new) > 3 or new[-1].isdigit():
            sys.exit(f"ERROR: new type name {new} must have at most 3 characters and must not end with a digit")
        explicit[old] = new

    info = analyse(parm)
    report(info, a.prmtop, a.quiet)
    reasons, unknown = find_problems(info, set(a.keep), a.rename_unknown)
    for old in explicit:
        if old not in info["types"]:
            sys.exit(f"ERROR: type {old} (from --map) is not in {a.prmtop}")
        reasons.setdefault(old, ["explicit --map"])

    if unknown:
        print("\nTypes whose element CP2K cannot guess (give them a &KIND section with ELEMENT in &SUBSYS):")
        for t, el in unknown.items():
            names = info["kinds"][t.upper()]
            print(f"  {t:<6s} element {el:<3s} kind name(s): {', '.join(names)}")

    if not reasons:
        print("\nNo atom type needs to be renamed.")
        if not a.dry_run:
            parm.write_parm(a.output)
            print(f"Topology copied to {a.output}")
        return

    rename = choose_names(info, reasons, explicit)
    print(f"\n{'old':<6s} {'new':<5s} {'element':<8s} {'atoms':>7s}  reason")
    for old, new in rename.items():
        d = info["types"][old]
        print(f"{old:<6s} {new:<5s} {d['el']:<8s} {d['n']:>7d}  {'; '.join(reasons[old])}")
    if a.dry_run:
        print(f"\n{len(rename)} type(s) would be renamed (dry run, nothing written).")
        return

    for i, at in enumerate(parm.atoms):
        if at.type in rename:
            at.type = rename[at.type]
            parm.parm_data["AMBER_ATOM_TYPE"][i] = at.type
    parm.write_parm(a.output)

    with open(a.map_file, "w") as fh:
        fh.write(f"# prmtop_fix_types.py v{VERSION}: {a.prmtop} -> {a.output}\n# old new element atoms\n")
        for old, new in rename.items():
            d = info["types"][old]
            fh.write(f"{old} {new} {d['el']} {d['n']}\n")

    new_info = analyse(parmed.load_file(a.output))
    print()
    report(new_info, a.output, True)
    left, _ = find_problems(new_info, set(a.keep), a.rename_unknown)
    if left:
        print("  WARNING: problems left: " + "; ".join(f"{t}: {', '.join(r)}" for t, r in left.items()))
    print(f"\n{len(rename)} type(s) renamed. Topology written to {a.output}; table in {a.map_file}")
    print("Only the type labels changed: an Amber single point gives the same energy with both topologies.")


if __name__ == "__main__":
    main()
