import csv, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Ellipse, FancyArrowPatch, Polygon
import numpy as np

root = sys.argv[1]
out = os.path.join(root, "Fig9_Model_Root"); os.makedirs(out, exist_ok=True)

ann = {r["gene"]: r["protein_name"] for r in csv.DictReader(open(os.path.join(root, "wheat_gene_annotation.tsv"), encoding="utf-8"), delimiter="\t")}
degs = [(r[""], float(r["log2FC"])) for r in csv.DictReader(open(os.path.join(root, "Root_DEGs.csv"), encoding="utf-8"))]
CU = "Cu⁺-ATPase"
FAM = {
    "Metal transp.": r"^Metal transporter$", "ZIP": r"^Zinc transporter", "COPT": r"^Copper transport",
    "KUP/HAK": r"^Potassium transporter", "SWEET": r"^Bidirectional sugar", "ABCG": r"ABC transporter G family",
    "SOD": r"^Superoxide dismutase", "CAT": r"^Catalase", "POD": r"^Peroxidase", "GPX": r"^Glutathione peroxidase",
    "PRX": r"eroxiredoxin", "GST": r"glutathione transferase|Glutathione S-transferase|^GST ",
    "GR": r"^Glutathione reductase", "MT": r"etallothionein",
    "NAS": r"icotianamine synthase", "DMAS": r"eoxymugineic", "YSL": r"YSL|Metal-nicotianamine",
    CU: r"^P-type Cu", "HMA-domain": r"^HMA domain", "ABCB": r"ABC transporter B family|MDR-like",
    "OPT": r"ligopeptide transporter",
    "ABCC": r"ABC transporter C family", "V-PPase": r"diphosphatase", "V-ATPase": r"V-type proton ATPase",
    "MATE": r"DETOXIFICATION", "TIP": r"onoplast intrinsic|Delta tonoplast", "VIT": r"^Vacuolar iron",
    "PAL": r"henylalanine ammonia", "C4H": r"cinnamate 4-monooxygenase", "4CL": r"4-coumarate",
    "CCR": r"innamoyl-CoA reductase", "CAD": r"innamyl-alcohol", "LAC": r"accase",
    "WRKY": r"WRKY", "NAC": r"NAC domain", "bHLH": r"BHLH|bHLH", "AP2/ERF": r"AP2/ERF", "bZIP": r"^BZIP|bZIP",
    "MYB": r"MYB|myb",
}
V = {}
with open(os.path.join(out, "fig9_family_values.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["family", "regex_on_UniProt_name", "n_DEGs", "n_up", "n_down", "mean_log2FC", "genes"])
    for k, p in FAM.items():
        v = [(g, l) for g, l in degs if re.search(p, ann.get(g, ""))]
        L = [l for _, l in v]
        if L:
            up, dn = sum(x > 0 for x in L), sum(x < 0 for x in L)
            V[k] = (sum(L) / len(L), len(L), up, dn)
            w.writerow([k, p, len(L), up, dn, round(sum(L) / len(L), 2), ";".join(g for g, _ in v)])
HUB = {"PRX", "GPX", "GST", "GR", "ZIP", "COPT", CU}
GS = {x["Group"]: (int(x["Up"]), int(x["DEGs"]), float(x["Percent_up"]))
      for x in csv.DictReader(open(os.path.join(root, "Fig7_GeneGroups_Root_main", "Root_gene_group_summary.csv")))}

RED, BLUE, MIX, GOLD = "#C0392B", "#2E6DA4", "#7B5EA7", "#E0A100"
fig = plt.figure(figsize=(16, 8.2), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 160); ax.set_ylim(0, 82); ax.axis("off")


def rbox(x, y, w, h, fc, ec, lw=1.2, r=0.8, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z))


def arrow_txt(mean, n, up, dn):
    if min(up, dn) / n >= 0.3:
        return "±", MIX, f"{abs(mean):.1f}" if False else f"{mean:+.1f}"
    return ("↑", RED, f"{mean:+.1f}") if mean > 0 else ("↓", BLUE, f"{mean:+.1f}")


def gene(x, y, key, w=9.6, h=4.6, fs=8.8, label=None, z=5):

    rbox(x, y, w, h, "white", "#555555", lw=0.8, r=0.35, z=z)
    ax.text(x + w / 2, y + h * 0.70, label or key, ha="center", va="center", fontsize=fs, fontweight="bold", zorder=z + 1)
    if key in V:
        m, n, up, dn = V[key]; a, c, s = arrow_txt(m, n, up, dn)
        ax.text(x + w / 2, y + h * 0.28, f"{a}{s} ({n})", ha="center", va="center", fontsize=fs - 0.4,
                color=c, fontweight="bold", zorder=z + 1)
    else:
        ax.text(x + w / 2, y + h * 0.28, "n.s.", ha="center", va="center", fontsize=fs - 0.4, color="grey", zorder=z + 1)
    if key in HUB:
        ax.text(x + w - 0.2, y + h - 0.1, "★", ha="center", va="center", fontsize=fs + 3.5, color=GOLD, zorder=z + 2)


def mem_transporter(x, key, label=None):
    w, h = 7.4, 8.2; y = 61.2
    rbox(x, y, w, h, "white", "#444444", lw=1.0, r=0.4, z=4)
    lab = label or key
    ax.text(x + w / 2, y + h - 2.3, lab, ha="center", va="center", fontsize=8.0 if "\n" in lab else 9.5,
            fontweight="bold", linespacing=1.0, zorder=5)
    if key in V:
        m, n, up, dn = V[key]; a, c, s = arrow_txt(m, n, up, dn)
        ax.text(x + w / 2, y + 3.0, f"{a}{s}", ha="center", va="center", fontsize=10.5, color=c, fontweight="bold", zorder=5)
        ax.text(x + w / 2, y + 1.0, f"({n} DEGs)", ha="center", va="center", fontsize=7.3, color="#444444", zorder=5)
    if key in HUB:
        ax.text(x + w - 0.3, y + h + 0.1, "★", ha="center", va="center", fontsize=14, color=GOLD, zorder=6)


def cd(x, y, r=1.15, z=6):
    ax.add_patch(Circle((x, y), r, fc="#9B1B1B", ec="white", lw=0.8, zorder=z))
    ax.text(x, y, "Cd", ha="center", va="center", fontsize=7.5, color="white", fontweight="bold", zorder=z + 1)


def arr(p1, p2, color="#2E7D4F", ls="-", lw=1.4, rad=0.0, z=3):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=12, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, zorder=z))


DASH = (0, (3.2, 2.2))
CD_C, TF_C, FLUX_C = "#1F3C99", "#5E3A86", "#2E7D4F"


def route(pts, color, dashed=False, lw=1.5, z=3):
    xs, ys = zip(*pts)
    ax.plot(xs[:-1] + (xs[-1],), ys, color=color, lw=lw, ls=DASH if dashed else "-",
            solid_capstyle="butt", dash_capstyle="butt", zorder=z)
    (x0, y0), (x1, y1) = pts[-2], pts[-1]
    L = max(abs(x1 - x0), abs(y1 - y0))
    f = min(1.0, 1.2 / L) if L else 1.0
    ax.add_patch(FancyArrowPatch((x1 - (x1 - x0) * f, y1 - (y1 - y0) * f), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=13, color=color, lw=lw, shrinkA=0, shrinkB=0, zorder=z + 1))


def title2(x, y, t1, t2, c1, c2=None, fs=9.8):
    ax.text(x, y, t1, ha="center", va="center", fontsize=fs, fontweight="bold", color=c1)
    if t2:
        ax.text(x, y - 1.6, t2, ha="center", va="center", fontsize=7.9, style="italic", color=c2 or c1)



rbox(0.3, 70, 159.4, 11.7, "#F6DCDC", "#C98A8A", lw=1.2, r=1.2)
ax.text(1.6, 79.2, "APOPLAST (root exterior / cell-wall space)", fontsize=15, fontweight="bold", color="#8E1B1B", va="center")
ax.text(1.6, 75.9, "textile wastewater / Cd²⁺ exposure", fontsize=11.5, color="#8E1B1B", va="center")
for (x, y) in [(64, 78.8), (60, 73.2), (72, 75.4), (81, 79.2), (89, 73.4), (98, 77.6), (107, 73.0)]:
    cd(x, y)
rbox(0.3, 0.3, 159.4, 60.6, "#EEF5F0", "#7FAE8F", lw=1.2, r=1.2)
ax.text(2.0, 57.6, "CYTOSOL", fontsize=15, fontweight="bold", color="#1E5B3A", va="center")


ax.add_patch(Rectangle((0.3, 61.8), 159.4, 7.2, fc="#F2C27B", ec="none", zorder=2))
xs = np.arange(0.9, 159.7, 1.05)
for yy in (68.4, 62.4):
    ax.scatter(xs, np.full_like(xs, yy), s=16, c="#D9822B", edgecolors="#9A5A17", linewidths=0.4, zorder=3)
for x0 in xs:
    ax.plot([x0, x0], [63.0, 67.8], color="#C9933F", lw=0.5, zorder=2)
rbox(1.2, 66.0, 21.5, 3.0, "white", "#8E1B1B", lw=1.0, r=0.3, z=4)
ax.text(11.95, 67.5, "PLASMA MEMBRANE", ha="center", va="center", fontsize=10, fontweight="bold", color="#8E1B1B", zorder=5)
mem = [("Metal transp.", "Metal transp.\n(NRAMP-type)"), ("ZIP", None), ("COPT", None), ("KUP/HAK", None),
       ("ABCG", "ABCG\n(efflux)"), ("SWEET", None)]
for i, (k, lab) in enumerate(mem):
    mem_transporter(46 + i * 9.6, k, lab)

star = [(18 + (4.6 if i % 2 == 0 else 2.6) * np.cos(t), 53 + (4.6 if i % 2 == 0 else 2.6) * np.sin(t))
        for i, t in enumerate(np.linspace(0, 2 * np.pi, 25)[:-1])]
ax.add_patch(Polygon(star, closed=True, fc="#7E3FB0", ec="#4B1F70", lw=1, zorder=4))
ax.text(18, 53.6, "ROS", ha="center", va="center", fontsize=10, color="white", fontweight="bold", zorder=5)
ax.text(18, 51.8, "burst", ha="center", va="center", fontsize=8.5, color="white", zorder=5)
ax.text(7.6, 47.4, "Cd²⁺ → Fenton-type chemistry → O₂⁻•, H₂O₂, •OH",
        ha="left", va="center", fontsize=8, style="italic", color="#333333")
up, n, pct = GS["Oxidative stress"]
title2(19.5, 45.4, "Antioxidant enzymes", f"oxidative-stress DEGs: {up}/{n} up ({pct:.0f}%)", "#1E5B3A", RED)
for i, k in enumerate(["SOD", "CAT", "POD", "GPX", "PRX"]):
    gene(2.0 + i * 7.15, 37.6, k, w=6.8, fs=8.4)
ax.text(7.2, 36.4, "→ lower ROS, H₂O₂, MDA (as in the biochemical assays)",
        ha="left", va="center", fontsize=7.8, style="italic", color="#333333")
route([(13.4, 53.0), (5.4, 53.0), (5.4, 42.2)], RED)

rbox(1.5, 21.6, 37.5, 13.8, "#F4FAF5", "#6E9C7C", lw=1.0)
title2(20.25, 34.0, "Glutathione-based Cd chelation & detoxification", "top GO/KEGG terms; GST, GPX, GR = PPI hub genes",
       "#1E5B3A", "#2F6B45", fs=9.5)
rbox(2.4, 26.6, 5.0, 4.6, "white", "#555555", lw=0.8, r=0.35, z=5)
ax.text(4.9, 28.9, "GSH", ha="center", va="center", fontsize=8.8, fontweight="bold", zorder=6)
gene(8.1, 26.6, "GST", w=9.8); gene(18.5, 26.6, "GR", w=9.2); gene(28.3, 26.6, "MT", w=9.6)
ax.text(20.25, 24.8, "GST conjugates GSH → Cd–thiol complexes; GR regenerates GSH",
        ha="center", va="center", fontsize=7.7, style="italic", color="#333333")
ax.text(20.25, 23.0, "metallothioneins (MT) bind cytosolic Cd²⁺",
        ha="center", va="center", fontsize=7.7, style="italic", color="#333333")
route([(20.25, 21.6), (20.25, 20.7)], FLUX_C)

rbox(1.5, 7.3, 37.5, 13.4, "#E5EEF8", "#7E9FC5", lw=1.0)
title2(20.25, 19.5, "Vacuolar Cd sequestration & storage", "tonoplast pumps load Cd–GSH / Cd–PC complexes",
       "#2B4F7A", fs=9.5)
for i, k in enumerate(["ABCC", "V-PPase", "V-ATPase"]):
    gene(2.4 + i * 12.1, 12.5, k, w=11.6)
for i, k in enumerate(["MATE", "TIP", "VIT"]):
    gene(2.4 + i * 12.1, 7.7, k, w=11.6)

rbox(43, 34.5, 47, 22.5, "#E6EFF9", "#7E9FC5", lw=1.0)
title2(66.5, 55.2, "Nicotianamine (NA) chelation", "NAS = strongest induction among all Cd-response families", "#2B4F7A", fs=10)
for (x, t) in [(46, "Met"), (54, "SAM"), (62, "NA")]:
    rbox(x, 47.0, 5.6, 3.0, "white", "#555555", lw=0.8, r=0.35, z=5)
    ax.text(x + 2.8, 48.5, t, ha="center", va="center", fontsize=9, fontweight="bold", zorder=6)
arr((51.7, 48.5), (53.9, 48.5), color="#333333", lw=1.1); arr((59.7, 48.5), (61.9, 48.5), color="#333333", lw=1.1)
ax.text(60.8, 50.8, "↑NAS", ha="center", va="center", fontsize=9, color=RED, fontweight="bold")
ax.add_patch(Circle((64.8, 41.6), 2.3, fc="#8CC9A0", ec="#2E7D4F", lw=1.0, zorder=5))
ax.text(64.8, 41.6, "Cd-NA", ha="center", va="center", fontsize=7.5, fontweight="bold", zorder=6)
arr((64.8, 46.9), (64.8, 44.0), color="#2E7D4F", lw=1.1)
gene(71.0, 46.6, "NAS", w=16.5, label="NAS family"); gene(71.0, 41.0, "DMAS", w=16.5)
ax.text(66.5, 36.3, "NA binds free Cd²⁺ → phloem loading → moved away from photosynthetic tissues",
        ha="center", va="center", fontsize=7.6, style="italic", color="#2B4F7A")

rbox(43, 21.0, 47, 12.4, "#E6EFF9", "#7E9FC5", lw=1.0)
up, n, pct = GS["Transporters"]
title2(66.5, 31.8, "Metal efflux & membrane transporters", f"transporter DEGs: {up}/{n} up ({pct:.0f}%)", "#2B4F7A", RED, fs=10)
for i, k in enumerate([CU, "HMA-domain", "ABCB", "OPT"]):
    gene(44.0 + i * 11.4, 22.6, k, w=10.8)

rbox(43, 7.3, 47, 12.6, "#F4FAF5", "#6E9C7C", lw=1.0)
title2(66.5, 18.4, "Phloem: long-distance transport", None, "#1E5B3A", fs=9.6)
gene(50, 12.2, "YSL", w=13)
ax.add_patch(Circle((73, 14.5), 2.1, fc="#8CC9A0", ec="#2E7D4F", lw=1.0, zorder=5))
ax.text(73, 14.5, "Cd-NA", ha="center", va="center", fontsize=7.2, fontweight="bold", zorder=6)
ax.text(66.5, 9.2, "↓ Cd reaching grain / leaves (proposed)", ha="center", va="center", fontsize=8.2, color=BLUE, fontweight="bold")
route([(87.5, 43.3), (91.6, 43.3), (91.6, 14.5), (90.0, 14.5)], FLUX_C)
route([(63.0, 14.5), (70.9, 14.5)], FLUX_C)

rbox(93, 1.0, 66, 15.8, "#F1E9F7", "#A88BC4", lw=1.0)
title2(126, 15.2, "Root endodermal cell wall: phenylpropanoid → lignin barrier",
       "lignin deposition proposed to restrict apoplastic Cd²⁺ entry into the stele", "#5B2E83", fs=9.8)
for i, k in enumerate(["PAL", "C4H", "4CL", "CCR", "CAD", "POD", "LAC"]):
    gene(94.2 + i * 9.25, 7.4, k, w=8.9, fs=8.6)
for row in range(2):
    for j in range(10):
        x0 = 94.6 + j * 6.3 + (3.1 if row else 0)
        if x0 + 5.6 > 158.4:
            continue
        ax.add_patch(FancyBboxPatch((x0, 1.9 + row * 2.6), 5.6, 2.0, boxstyle="round,pad=0,rounding_size=0.3",
                                    fc="#C9A8E0", ec="#9C79BC", lw=0.6, zorder=3))

ax.add_patch(Ellipse((133, 39.5), 44, 37, fc="#8E63B8", ec="#5E3A86", lw=2, zorder=2))
ax.add_patch(Ellipse((133, 40), 16, 13, fc="#6A4294", ec="none", zorder=3))
rbox(119, 57.0, 28, 3.6, "#EDE3F6", "#5E3A86", lw=1.2, r=0.6, z=4)
ax.text(133, 58.8, "Nucleus", ha="center", va="center", fontsize=13, fontweight="bold", color="#3F2266", zorder=5)
up, n, pct = GS["Transcription factors"]
ax.text(133, 54.2, "Transcriptional control", ha="center", va="center", fontsize=9.5, fontweight="bold",
        style="italic", color="white", zorder=5)
rbox(121.5, 50.0, 23, 2.8, "white", RED, lw=1.0, r=0.4, z=5)
ax.text(133, 51.4, f"TF DEGs: {up}/{n} up ({pct:.0f}%)", ha="center", va="center", fontsize=8.6, color=RED,
        fontweight="bold", zorder=6)
tf_pos = {"WRKY": (114.5, 43.6), "NAC": (139.5, 43.6), "bHLH": (114.5, 36.0), "AP2/ERF": (139.5, 36.0),
          "bZIP": (119.5, 28.4), "MYB": (134.5, 28.4)}
for k, (x, y) in tf_pos.items():
    gene(x, y, k, w=12.0, fs=9.0, z=6)
route([(115.8, 51.0), (90.0, 51.0)], TF_C, dashed=True)
route([(116.8, 27.0), (90.0, 27.0)], TF_C, dashed=True)
route([(126.0, 21.9), (126.0, 16.8)], TF_C, dashed=True)


for xc in (51.5, 59.3):
    route([(xc, 61.2), (xc, 57.0)], CD_C, dashed=True)
route([(47.6, 61.2), (47.6, 59.3), (41.2, 59.3), (41.2, 28.9), (39.0, 28.9)], CD_C, dashed=True)
rbox(1.5, 0.4, 89, 6.5, "white", "#888888", lw=0.8, r=0.5)
leg = [("↑ +x.x (n)", RED, "up-regulated: mean log2FC (no. of DEGs)"),
       ("↓ −x.x (n)", BLUE, "down-regulated"),
       ("± x.x (n)", MIX, "mixed (≥30% of DEGs in each direction)"),
       ("★", GOLD, "hub gene, stress-response PPI network (Fig. 8D)")]
for i, (sym, c, txt) in enumerate(leg):
    xx = 2.5 + (i % 2) * 40.0; yy = 5.6 - (i // 2) * 2.1
    ax.text(xx, yy, sym, color=c, fontsize=9.5, fontweight="bold", va="center")
    ax.text(xx + 8.5, yy, txt, fontsize=8.4, va="center", color="#222222")
cd(83.0, 5.6, r=0.9); ax.text(84.3, 5.6, "Cd²⁺", fontsize=8.4, va="center")
ax.add_patch(Circle((83.0, 3.5), 0.9, fc="#8CC9A0", ec="#2E7D4F", lw=0.8, zorder=5))
ax.text(84.3, 3.5, "Cd-NA", fontsize=8.4, va="center")
ax.text(46.0, 1.35, "Values: reanalysed public RNA-seq (wild-type roots, Cd vs control). FeO(B)-NP effects on these pathways were not assessed.",
        ha="center", va="center", fontsize=7.3, style="italic", color="#555555")

for ext in ("png", "tif"):
    kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}
    fig.savefig(os.path.join(out, f"Fig9_root_cell_model.{ext}"), dpi=600, **kw)
print("saved; families with values:", len(V))
