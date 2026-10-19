import csv, os, re, sys

src, out = sys.argv[1], sys.argv[2]
v11 = re.compile(r"(TraesCS[0-9A-Z]+02G[0-9]+)\.\d+")

gene2go, ann, gene2up = {}, {}, {}
with open(src, encoding="utf-8") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        genes = set(v11.findall(row["EnsemblPlants"] or ""))
        gos = [g.strip() for g in (row["Gene Ontology IDs"] or "").split(";") if g.strip()]
        for g in genes:
            gene2go.setdefault(g, set()).update(gos)
            gene2up.setdefault(g, set()).add(row["Entry"])
            name = re.sub(r"\s*\(EC [^)]*\)", "", row["Protein names"]).split(" (")[0].strip()
            rec = (row["Reviewed"] == "reviewed", row["Entry"], name)
            if g not in ann or rec[0] and not ann[g][0]:
                ann[g] = rec

with open(os.path.join(out, "wheat_gene2go.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["gene", "go"])
    for g in sorted(gene2go):
        for t in sorted(gene2go[g]):
            w.writerow([g, t])
with open(os.path.join(out, "wheat_gene_annotation.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["gene", "uniprot", "reviewed", "protein_name"])
    for g in sorted(ann):
        w.writerow([g, ann[g][1], "yes" if ann[g][0] else "no", ann[g][2]])

with open(os.path.join(out, "wheat_gene2uniprot.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["gene", "uniprot"])
    for g in sorted(gene2up):
        for u in sorted(gene2up[g]):
            w.writerow([g, u])

print("genes with UniProt entry:", len(ann), "| genes with GO:", sum(1 for g in gene2go if gene2go[g]),
      "| gene-GO pairs:", sum(len(v) for v in gene2go.values()))
