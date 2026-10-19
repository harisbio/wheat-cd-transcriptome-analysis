
import csv, os, sys, time, urllib.parse, urllib.request

root, tis, gfolder = sys.argv[1], sys.argv[2], sys.argv[3]
score = sys.argv[4] if len(sys.argv) > 4 else "700"
API = "https://string-db.org/api/tsv/"
CALLER = "wheat_Cd_reanalysis"
out = os.path.join(root, "ppi", tis); os.makedirs(out, exist_ok=True)


def post(method, params, tries=5):
    data = urllib.parse.urlencode(params).encode()
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(API + method, data=data), timeout=600) as r:
                return r.read().decode("utf-8")
        except Exception as e:
            err = e; time.sleep(5 * (k + 1))
    raise SystemExit(f"STRING {method} failed: {err}")


degs = [r for r in csv.DictReader(open(os.path.join(root, f"{tis}_DEGs.csv")))]
deg_ids = [r[""] for r in degs]
g2u = {}
for r in csv.DictReader(open(os.path.join(root, "wheat_gene2uniprot.tsv")), delimiter="\t"):
    g2u.setdefault(r["gene"], []).append(r["uniprot"])

# 1) UniProt -> STRING
accs = sorted({a for g in deg_ids for a in g2u.get(g, [])})
acc2str = {}
for i in range(0, len(accs), 1000):
    txt = post("get_string_ids", {"identifiers": "\r".join(accs[i:i + 1000]), "species": "4565",
                                  "limit": "1", "echo_query": "1", "caller_identity": CALLER})
    for r in csv.DictReader(txt.splitlines(), delimiter="\t"):
        if "stringId" in r:
            acc2str[r["queryItem"]] = r["stringId"]
    time.sleep(1)

gene2str = {}
for g in deg_ids:
    for a in g2u.get(g, []):
        if a in acc2str:
            gene2str[g] = acc2str[a]; break
with open(os.path.join(out, "gene2string.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["gene", "string_id"])
    for g, s in gene2str.items():
        w.writerow([g, s])
print(f"{tis}: {len(deg_ids)} DEGs, {len(gene2str)} mapped to STRING")


# 2) networks
def network(genes, name):
    ids = sorted({gene2str[g] for g in genes if g in gene2str})
    txt = post("network", {"identifiers": "\r".join(ids), "species": "4565",
                           "required_score": score, "caller_identity": CALLER})
    path = os.path.join(out, f"network_{name}.tsv")
    open(path, "w", encoding="utf-8").write(txt)
    n_edges = max(0, len(txt.strip().splitlines()) - 1)
    print(f"  {name}: {len(ids)} proteins submitted, {n_edges} interactions (score >= {score})")
    time.sleep(1)

def network_blocks(genes, name, block=1000):
    """STRING's API accepts at most 2000 proteins per network call. Split the proteins into
    blocks and query every pair of blocks, so that every protein pair is covered; then de-duplicate."""
    ids = sorted({gene2str[g] for g in genes if g in gene2str})
    blocks = [ids[i:i + block] for i in range(0, len(ids), block)]
    header, edges = None, {}
    pairs = [(i, j) for i in range(len(blocks)) for j in range(i + 1, len(blocks))] or [(0, 0)]
    for n, (i, j) in enumerate(pairs, 1):
        sub = blocks[i] + (blocks[j] if j != i else [])
        txt = post("network", {"identifiers": "\r".join(sub), "species": "4565",
                               "required_score": score, "caller_identity": CALLER})
        lines = txt.strip().splitlines()
        header = header or lines[0]
        for ln in lines[1:]:
            f = ln.split("\t"); key = tuple(sorted((f[0], f[1])))
            edges.setdefault(key, ln)
        print(f"    block pair {n}/{len(pairs)}: {len(edges)} unique interactions so far", flush=True)
        time.sleep(1)
    with open(os.path.join(out, f"network_{name}.tsv"), "w", encoding="utf-8") as fh:
        fh.write(header + "\n" + "\n".join(edges.values()) + "\n")
    print(f"  {name}: {len(ids)} proteins, {len(edges)} interactions (score >= {score})")

network_blocks(deg_ids, "all_DEGs")
for grp in ["Oxidative_stress", "Transcription_factors", "Transporters"]:
    rows = list(csv.DictReader(open(os.path.join(root, gfolder, f"{tis}_{grp}_DEGs.csv"))))
    network([r["gene"] for r in rows], grp)
