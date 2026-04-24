#!/usr/bin/env python3
"""
02_run_docking.py
=================
Main molecular docking pipeline.

Reads batch_config.csv and for each compound-protein pair:
  1. Converts SMILES → 3D SDF (RDKit ETKDGv3 + MMFF94)
  2. Converts SDF → PDBQT (meeko or obabel fallback)
  3. Converts prepped PDB → receptor PDBQT
  4. Runs AutoDock Vina
  5. Parses and saves results

Usage
-----
  python scripts/02_run_docking.py --config config/batch_config.csv --out results/

Requirements
------------
  pip install rdkit meeko requests scipy gemmi
  AutoDock Vina 1.2.5 on PATH
  UCSF ChimeraX prepped PDB files in results/
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

# ── Colour output ─────────────────────────────────────────────────
G="\033[92m"; Y="\033[93m"; R="\033[91m"; C="\033[96m"; B="\033[1m"; X="\033[0m"
def ok(m):   print(f"{G}  ✓  {m}{X}")
def warn(m): print(f"{Y}  ⚠  {m}{X}")
def die(m):  print(f"{R}  ✗  {m}{X}"); sys.exit(1)
def step(m): print(f"\n{B}{C}▶  {m}{X}")
def info(m): print(f"     {m}")


# ── AutoDock atom type mapping ────────────────────────────────────
TYPE_MAP = {
    "C":"C","N":"N","O":"OA","S":"SA","H":"H","P":"P",
    "F":"F","CL":"Cl","BR":"Br","I":"I",
    "FE":"Fe","ZN":"Zn","MG":"Mg","CA":"Ca","MN":"Mn",
}


def build_receptor_pdbqt(pdb_path: Path, pdbqt_path: Path):
    """Convert prepped PDB to PDBQT with correct AutoDock format."""
    with open(pdb_path) as fh, open(pdbqt_path, "w") as out:
        for line in fh:
            rec = line[:6].strip()
            if rec == "TER":
                out.write("TER\n"); continue
            if rec not in ("ATOM", "HETATM"): continue
            aname = line[12:16].strip()
            if aname.startswith("H"): continue  # merge non-polar H

            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                an = aname.lstrip("0123456789 ")
                element = an[:1].upper()
                if aname in ("CA","CB","CG","CD","CE","CZ","CH",
                             "CD1","CD2","CE1","CE2","CG1","CG2","CG3"):
                    element = "C"
                elif aname in ("ND1","ND2","NE","NE1","NE2","NH1","NH2","NZ"):
                    element = "N"
                elif aname in ("OD1","OD2","OE1","OE2","OG","OG1","OH","OXT"):
                    element = "O"
                elif aname == "SD":
                    element = "S"

            ad_type = TYPE_MAP.get(element, "C")
            base    = line[:66].ljust(66)
            charge  = f"{0.0:10.3f}"
            out.write(f"{base}{charge} {ad_type}\n")


def smiles_to_pdbqt(smiles: str, name: str, out_dir: Path, safe_name: str) -> Path:
    """SMILES → 3D SDF → PDBQT via RDKit + meeko."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        die("RDKit not installed. Run: conda install -c conda-forge rdkit")

    sdf_path  = out_dir / f"{safe_name}.sdf"
    pdbqt_out = out_dir / f"{safe_name}_ligand.pdbqt"

    # Handle selenium — replace Se with S (S-analog approximation)
    if "[Se]" in smiles:
        warn(f"{name}: selenium replaced with sulfur (S-Methylcysteine analog)")
        smiles = smiles.replace("[Se]", "S")

    # SMILES → 3D
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {name}: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)

    # MMFF94 minimization (2000 iterations)
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    if ff:
        ff.Minimize(maxIts=2000)
    else:
        AllChem.UFFGetMoleculeForceField(mol).Minimize(maxIts=2000)

    mol.SetProp("_Name", name)
    Chem.SDWriter(str(sdf_path)).write(mol)

    # SDF → PDBQT (meeko preferred)
    mk = shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")
    if mk:
        r = subprocess.run([mk, "-i", str(sdf_path), "-o", str(pdbqt_out)],
                           capture_output=True, text=True)
        if r.returncode == 0 and pdbqt_out.exists():
            return pdbqt_out
        warn(f"meeko failed: {r.stderr[-200:]}")

    if shutil.which("obabel"):
        subprocess.run(["obabel", str(sdf_path), "-O", str(pdbqt_out),
                        "--gen3d", "-p", "7.4", "--partialcharge", "gasteiger"],
                       capture_output=True, text=True)
        if pdbqt_out.exists():
            return pdbqt_out

    raise RuntimeError(f"Cannot convert {name} to PDBQT. Install meeko or obabel.")


def run_vina(receptor: Path, ligand: Path, out_dir: Path,
             safe_name: str, center: tuple, size: tuple,
             exhaustiveness: int = 8, num_modes: int = 9) -> tuple[Path, str]:
    """Run AutoDock Vina and return (poses_path, stdout)."""
    vina = shutil.which("vina")
    if not vina:
        die("AutoDock Vina not found on PATH.")

    poses_path = out_dir / f"{safe_name}_poses.pdbqt"
    cfg_path   = out_dir / f"{safe_name}_vina.cfg"

    cfg_path.write_text(textwrap.dedent(f"""\
        receptor = {receptor.resolve()}
        ligand   = {ligand.resolve()}
        center_x = {center[0]}
        center_y = {center[1]}
        center_z = {center[2]}
        size_x   = {size[0]}
        size_y   = {size[1]}
        size_z   = {size[2]}
        exhaustiveness = {exhaustiveness}
        num_modes      = {num_modes}
        energy_range   = 3
        out = {poses_path.resolve()}
    """))

    r = subprocess.run([vina, "--config", str(cfg_path)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Vina failed:\n{r.stderr[-500:]}")
    return poses_path, r.stdout + r.stderr


def parse_vina_scores(vina_output: str) -> list[dict]:
    rows = re.findall(
        r"^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)",
        vina_output, re.MULTILINE
    )
    return [{"mode": m, "affinity": float(a),
             "rmsd_lb": float(rl), "rmsd_ub": float(ru)}
            for m, a, rl, ru in rows]


def affinity_label(v: float) -> str:
    if v <= -10: return f"{G}Strong{X}"
    if v <= -7:  return f"{G}Moderate{X}"
    if v <= -5:  return f"{Y}Weak{X}"
    return f"{R}Very weak{X}"


def load_config(config_path: Path, only: list[str] | None) -> list[dict]:
    jobs = []
    with open(config_path, newline="") as f:
        for row in csv.DictReader(f):
            pid = row["protein_id"].strip()
            if only and pid not in only: continue
            if not row["smiles"].strip(): continue
            jobs.append({
                "protein_id":    pid,
                "pdb":           Path(row["pdb_path"].strip()),
                "compound_name": row["compound_name"].strip(),
                "cid":           row.get("cid", "").strip(),
                "smiles":        row["smiles"].strip(),
                "pchembl":       float(row.get("pchembl", 0) or 0),
                "center":        (float(row["center_x"]), float(row["center_y"]), float(row["center_z"])),
                "size":          (float(row["size_x"]),   float(row["size_y"]),   float(row["size_z"])),
            })
    return jobs


def main():
    p = argparse.ArgumentParser(description="Molecular docking pipeline")
    p.add_argument("--config",  required=True, help="Path to batch_config.csv")
    p.add_argument("--out",     required=True, help="Output directory")
    p.add_argument("--only",    nargs="+",     help="Run only these protein IDs")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--num-modes",      type=int, default=9)
    args = p.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs = load_config(Path(args.config), args.only)
    ok(f"Loaded {len(jobs)} docking jobs")

    # Build receptor PDBQTs (once per unique protein)
    step("Building receptor PDBQT files")
    receptor_cache = {}
    seen_proteins = set()
    for job in jobs:
        pid = job["protein_id"]
        if pid in seen_proteins: continue
        seen_proteins.add(pid)
        prepped_pdb = out_dir / f"{pid}_prepped.pdb"
        pdbqt       = out_dir / f"{pid}_receptor.pdbqt"
        if not prepped_pdb.exists():
            warn(f"{pid}_prepped.pdb not found — run 01_protein_prep.cxc first")
            continue
        if not pdbqt.exists():
            build_receptor_pdbqt(prepped_pdb, pdbqt)
        ok(f"{pid} receptor PDBQT ready")
        receptor_cache[pid] = pdbqt

    # Run docking
    step("Running docking")
    all_results  = []
    best_results = []
    failed       = []

    for i, job in enumerate(jobs, 1):
        pid   = job["protein_id"]
        cname = job["compound_name"]
        safe  = re.sub(r"[^a-zA-Z0-9]", "_", f"{pid}_{cname}")
        print(f"\n  [{i}/{len(jobs)}] {pid} × {cname} (pChEMBL {job['pchembl']:.2f})")

        if pid not in receptor_cache:
            print(f"{R}  ✗  Skipping — receptor not prepared{X}")
            failed.append(f"{pid} × {cname}")
            continue

        try:
            ligand_pdbqt = smiles_to_pdbqt(job["smiles"], cname, out_dir, safe)
            t0 = time.time()
            poses_path, vina_out = run_vina(
                receptor=receptor_cache[pid],
                ligand=ligand_pdbqt,
                out_dir=out_dir,
                safe_name=safe,
                center=job["center"],
                size=job["size"],
                exhaustiveness=args.exhaustiveness,
                num_modes=args.num_modes,
            )
            elapsed = time.time() - t0
            scores  = parse_vina_scores(vina_out)
            if not scores:
                raise RuntimeError("No scores parsed")

            best = scores[0]["affinity"]
            ok(f"Best: {best} kcal/mol  {affinity_label(best)}  ({elapsed:.0f}s)")

            for s in scores:
                all_results.append({
                    "protein_id": pid, "compound_name": cname,
                    "cid": job["cid"], "pchembl": job["pchembl"],
                    "mode": s["mode"], "affinity_kcal_mol": s["affinity"],
                    "rmsd_lb": s["rmsd_lb"], "rmsd_ub": s["rmsd_ub"],
                    "poses_file": str(poses_path),
                })
            best_results.append({
                "protein_id": pid, "compound_name": cname,
                "cid": job["cid"], "pchembl": job["pchembl"],
                "affinity_kcal_mol": best, "poses_file": str(poses_path),
            })

        except Exception as e:
            print(f"{R}  ✗  FAILED: {e}{X}")
            failed.append(f"{pid} × {cname}")

    # Save results
    step("Saving results")
    if all_results:
        f1 = out_dir / "combined_docking_results.csv"
        with open(f1, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
            w.writeheader(); w.writerows(all_results)
        ok(f"Full results   → {f1}")

    if best_results:
        f2 = out_dir / "best_poses_summary.csv"
        with open(f2, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(best_results[0].keys()))
            w.writeheader()
            w.writerows(sorted(best_results, key=lambda x: x["affinity_kcal_mol"]))
        ok(f"Best poses     → {f2}")

    # Summary table
    print(f"\n{B}{'─'*65}{X}")
    print(f"  {'Protein':<15} {'Compound':<25} {'pChEMBL':>8} {'Affinity':>10}")
    print(f"  {'─'*63}")
    for r in sorted(best_results, key=lambda x: x["affinity_kcal_mol"]):
        print(f"  {r['protein_id']:<15} {r['compound_name']:<25} "
              f"{r['pchembl']:>8.2f} {r['affinity_kcal_mol']:>10.2f}")

    if failed:
        print(f"\n{R}{B}Failed: {len(failed)} jobs{X}")
        for f in failed: info(f"• {f}")
    else:
        print(f"\n{G}{B}All {len(best_results)} jobs completed!{X}")


if __name__ == "__main__":
    main()
