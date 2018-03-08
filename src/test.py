import json

with open("conv_507F3FF70D871519818695426.txt", "r") as f:
    text = f.read()

rows = text.split("\n")[1:]
needed = {}
for row in rows:
    try:
        gene, entrez_id, species, _ = row.split("\t")
        if species == "Homo sapiens":
            needed[gene] = entrez_id
    except:
        pass
with open("gene-di_conversion.json", "w") as out:
    json.dump(needed, out)
