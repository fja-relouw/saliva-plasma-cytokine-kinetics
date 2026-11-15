"""Tables made in Python (run by run_all.m, or on its own: python python/tables.py).

  table1_demographics.csv     Table 1: participant characteristics, median (IQR)
  below_lloq_by_time.csv      Table S1: % of values below the LLOQ per cytokine and time
                              (pooled placebo, anakinra and tocilizumab groups)
  below_lloq_overall.csv      Methods: overall % below the LLOQ (all four groups, all times)
  physicochemical.csv         Tables S12-S13: size, charge and hydrophobicity of the mature chains

Plasma values count as below the LLOQ when reported at the lower limit of the assay (the
per-cytokine minimum of the plasma sheet); salivary values when the Olink value is below the
LLOQ or missing (failed wells excluded). Only the study day is analysed (Day 0 and its 24-h
sample, coded D0 and D1). Protein sequences: UniProtKB release 2026_04 (sequences.fasta).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "DEDICATE-LPS SALIVA.xlsx"
EXPORT = ROOT / "data" / "DEDICATE-LPS_study_export_20250908.csv"
OUT = ROOT / "results" / "tables"

CYTO = ["TNF", "IL6", "IL8", "IL10", "IL1RA", "GCSF", "IP10", "MCP1", "MIP1a", "MIP1B"]
ASSAY = {"TNF": "TNF", "IL6": "IL6", "CXCL8": "IL8", "IL10": "IL10", "IL1RN": "IL1RA",
         "CSF3": "GCSF", "CXCL10": "IP10", "CCL2": "MCP1", "CCL3": "MIP1a", "CCL4": "MIP1B"}
ALL_GROUPS = ["placebo", "dexamethason", "tocilizumab", "anakinra"]
POOLED = ["placebo", "anakinra", "tocilizumab"]
EXCLUDE = {"IL1RA": "anakinra", "IL6": "tocilizumab"}     # the drug acts on the cytokine itself
DAYS = ["D0", "D1"]
TIMES = [30, 60, 90, 120, 180, 240, 360, 480, 1440]


def fmt(x, d):
    """Midline decimal point and true minus sign."""
    return f"{x:.{d}f}".replace("-", "−").replace(".", "·")


def load_data():
    p = pd.read_excel(XLSX, sheet_name="plasma cytokines ")
    p.columns = ["id", "day", "time", "cyto", "conc", "group"]
    p["conc"] = pd.to_numeric(p["conc"], errors="coerce")
    p["group"] = p["group"].astype(str).str.strip().str.lower()
    p = p[p["cyto"].isin(CYTO) & p["conc"].notna()].copy()
    p["blq"] = p["conc"] == p["cyto"].map(p.groupby("cyto")["conc"].min())
    p = p[p["day"].isin(DAYS)].copy()
    p.loc[p["day"] == "D1", "time"] = 1440
    group_of = p.drop_duplicates("id").set_index("id")["group"]

    s = pd.read_excel(XLSX, sheet_name="Olink flex saliva")
    s = s[s["SampleType"] == "SAMPLE"].copy()
    qv = pd.to_numeric(s["QuantifiedValue"], errors="coerce")
    s["blq"] = qv.isna() | (qv < pd.to_numeric(s["LLOQ"], errors="coerce"))
    s["failed"] = s["SampleID"].map(s.groupby("SampleID")["blq"].sum()) >= 8
    parts = s["SampleID"].str.replace(" ", "_").str.split("_")
    s["id"], s["day"] = parts.str[0], parts.str[1]
    s["time"] = [float(x[2][1:]) * 60 if len(x) == 3 else 1440.0 for x in parts]
    s["cyto"] = s["Assay"].map(ASSAY)
    s["group"] = s["id"].map(group_of)
    s = s[s["cyto"].notna() & s["day"].isin(DAYS)].copy()
    return p, s, group_of


def demographics(group_of, saliva_ids):
    rows = [("Age, years", "Age_screening", 1, 0), ("Weight, kg", "Weight_screening", 1, 1),
            ("Height, cm", "Length_screening", 100, 0), ("Body-mass index, kg/m²", "BMI", 1, 1),
            ("Heart rate, beats per min", "HR_screening", 1, 0),
            ("Systolic blood pressure, mm Hg", "SBP_screening", 1, 0),
            ("Diastolic blood pressure, mm Hg", "DBP_screening", 1, 0)]
    cols = [("Placebo", ["placebo"]), ("Dexamethasone", ["dexamethason"]), ("Tocilizumab", ["tocilizumab"]),
            ("Anakinra", ["anakinra"]), ("All participants", ALL_GROUPS), ("Saliva subgroup", None)]
    d = pd.read_csv(EXPORT, sep=";", encoding="utf-8-sig").copy()
    d["group"] = d["Participant Id"].map(group_of)
    d = d[d["group"].isin(ALL_GROUPS)]
    sel = {name: d[d["Participant Id"].isin(saliva_ids)] if g is None else d[d["group"].isin(g)] for name, g in cols}
    out = []
    for lab, col, scale, dec in rows:
        row = {"Characteristic": lab}
        for name, x in sel.items():
            q1, m, q3 = np.percentile(pd.to_numeric(x[col]) * scale, [25, 50, 75])
            row[f"{name} (n={len(x)})"] = f"{fmt(m, dec)} ({fmt(q1, dec)}–{fmt(q3, dec)})"
        out.append(row)
    row = {"Characteristic": "Sex, male"}
    for name, x in sel.items():
        row[f"{name} (n={len(x)})"] = f"{(x['Gender'] == 0).sum()} ({100 * (x['Gender'] == 0).mean():.0f}%)"
    out.append(row)
    return pd.DataFrame(out)


def below_lloq(p, s):
    rows, overall = [], []
    for matrix, d in (("plasma", p), ("saliva", s)):
        pooled = d[d["group"].isin(POOLED)]
        pooled = pooled[[EXCLUDE.get(c) != g for c, g in zip(pooled["cyto"], pooled["group"])]]
        for k in CYTO:
            dk = pooled[pooled["cyto"] == k]
            row = {"matrix": matrix, "Cytokine": k}
            for t in TIMES:
                dt = dk[dk["time"] == t]
                row[f"t{t}_min"] = 100 * dt["blq"].mean() if len(dt) else float("nan")
            row["n_participants"] = dk["id"].nunique()
            rows.append(row)
        a = d[d["group"].isin(ALL_GROUPS)]
        overall.append({"matrix": matrix, "n_values": len(a), "n_below": int(a["blq"].sum()),
                        "pct_below": 100 * a["blq"].mean()})
    return pd.DataFrame(rows), pd.DataFrame(overall)


def physicochemical():
    # accession: (cytokine, first and last residue of the mature chain, copies in the native form)
    mature = {"P13236": ("MIP1B", 24, 92, 1), "P10147": ("MIP1a", 24, 92, 1), "P13500": ("MCP1", 24, 99, 1),
              "P10145": ("IL8", 28, 99, 1), "P02778": ("IP10", 22, 98, 1), "P22301": ("IL10", 19, 178, 2),
              "P05231": ("IL6", 30, 212, 1), "P01375": ("TNF", 77, 233, 3), "P09919": ("GCSF", 30, 207, 1),
              "P18510": ("IL1RA", 26, 177, 1)}
    pk = {"Nterm": 8.6, "Cterm": 3.6, "C": 8.5, "D": 3.9, "E": 4.1, "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1}  # EMBOSS
    rows = []
    for rec in SeqIO.parse(Path(__file__).resolve().parent / "sequences.fasta", "fasta"):
        acc = rec.id.split("|")[0]
        cyto, a, b, copies = mature[acc]
        seq = str(rec.seq)[a - 1:b]
        pa = ProteinAnalysis(seq)
        mass = pa.molecular_weight()
        pos = 1 / (1 + 10 ** (7.4 - pk["Nterm"])) + sum(seq.count(x) / (1 + 10 ** (7.4 - pk[x])) for x in "KRH")
        neg = 1 / (1 + 10 ** (pk["Cterm"] - 7.4)) + sum(seq.count(x) / (1 + 10 ** (pk[x] - 7.4)) for x in "DEYC")
        aliphatic = 100 * (seq.count("A") + 2.9 * seq.count("V") + 3.9 * (seq.count("I") + seq.count("L"))) / len(seq)
        rows.append({"Cytokine": cyto, "UniProt": acc, "length_aa": len(seq), "monomer_mass_Da": mass,
                     "native_mass_kDa": copies * mass / 1000, "pI": pa.isoelectric_point(),
                     "net_charge_pH7_4": pos - neg, "Rmin_nm": 0.066 * (copies * mass) ** (1 / 3),  # Erickson 2009
                     "GRAVY": pa.gravy(), "aliphatic_index": aliphatic,                             # Ikai 1980
                     "aromaticity": pa.aromaticity(), "n_cys": seq.count("C")})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p, s, group_of = load_data()
    demographics(group_of, set(s["id"])).to_csv(OUT / "table1_demographics.csv", index=False, encoding="utf-8")
    by_time, overall = below_lloq(p, s[~s["failed"]])
    by_time.to_csv(OUT / "below_lloq_by_time.csv", index=False)
    overall.to_csv(OUT / "below_lloq_overall.csv", index=False)
    physicochemical().to_csv(OUT / "physicochemical.csv", index=False)
    print(f"Python tables written to {OUT}")


if __name__ == "__main__":
    main()
