ACCESSIONS = {
    "A_fumigatus": "GCF_000002655.1",
    "A_flavus":    "GCF_009017415.1",
    "A_niger":     "GCF_000002855.4",
    "A_oryzae":    "GCF_000184455.2",
    "A_nidulans":  "GCF_000011425.1",
    "A_terreus":   "GCF_000149615.1",
}

for species, accession in ACCESSIONS.items():
    print(f"{species}: {accession}")