#!/usr/bin/env python3

import ujson
import sys

def main():
    id_file = sys.argv[1]
    db_file = sys.argv[2]
    query = {}
    q_order = []
    with open(id_file) as infile:
        for line in infile:
            someid = line.strip()
            if someid:
                q_order.append(someid)
                query[someid] = []
    uniprot = ujson.load(open(db_file))
    for entry in uniprot:
        uniprot_id, accessions, seq_len, features, pdb = entry
        is_query = False
        hit_ids = []
        for ids in accessions:
            is_hit = [anyid in query for anyid in ids]
            if any(is_hit):
                hit_ids.append(ids[is_hit.index(True)])
                is_query = True
        if not is_query:
            continue
        subcellular = []
        for subcell in features["Subcellular"]:
            subcellular.append(subcell[2])
        for hit_id in hit_ids:
            query[hit_id].append(subcellular)
    for query_id in q_order:
        sys.stdout.write("{}\t{}\n".format(
            query_id, "|".join([":".join(subcell) for subcell in query[query_id]])))

if __name__ == "__main__":
    sys.exit(main())
