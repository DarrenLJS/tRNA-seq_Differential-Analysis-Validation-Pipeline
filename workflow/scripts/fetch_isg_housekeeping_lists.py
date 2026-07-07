"""
workflow/scripts/fetch_isg_housekeeping_lists.py

Fetches curated ISG (Interferome) and housekeeping (Eisenberg & Levanon
2013) gene ID lists for the rule-16 Fisher's exact enrichment test.

Output: long-format TSV, columns [gene_id, gene_set], gene_set in
{"ISG", "housekeeping"}.

FIX-flag, needs verification against real access before trusting this
in a full run
-------------------------------------------------------------------
Interferome (http://www.interferome.org) does not expose a simple,
stable bulk-download REST endpoint as of this script's writing --
programmatic access has historically required either an account-gated
export or scraping their search-result pages, both of which are fragile
and may break silently. Two fallbacks are implemented in order:

  1. Try Interferome's documented query interface (`INTERFEROME_QUERY_URL`
     below) for human, IFN-stimulated genes.
  2. If that fails (non-200 response, or a response that doesn't parse as
     the expected table), fall back to a small STATIC, literature-cited
     core ISG list (Schoggins et al. 2011 Nature, the ~380-gene ISG
     overexpression screen) bundled directly in this script. This static
     list is a legitimate, citable stand-in but is NOT the full
     Interferome catalogue -- flagged in the output log so it's clear
     which source actually produced the final gene set, since this
     materially affects the Fisher's exact test's statistical power in
     rule 16.

The housekeeping list (Eisenberg & Levanon 2013, "Human housekeeping
genes, revisited") is bundled as a static list for the same reason --
no live API exists for it; the paper's own supplementary gene list is the
standard citation.
"""

import logging
import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

INTERFEROME_QUERY_URL = "http://www.interferome.org/interferome/rest/search"  # FIX: verify this endpoint is current

# Static fallback: a small, well-cited core ISG set (Schoggins et al. 2011).
# THIS IS A PLACEHOLDER SUBSET for pipeline scaffolding -- replace with the
# full gene list from the paper's supplementary table before treating rule
# 16's Fisher's exact result as final. Do not publish results based on this
# placeholder list without swapping it for the complete curated set.
STATIC_CORE_ISG_FALLBACK = [
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "RSAD2", "IFI6", "IFI27", "IFI44", "IFI44L", "IRF7", "STAT1", "STAT2",
    "DDX58", "IFIH1", "CXCL10", "OASL", "BST2", "HERC5", "HERC6", "USP18",
]

# Static placeholder for housekeeping genes -- SAME caveat: replace with
# Eisenberg & Levanon 2013's full supplementary gene list before treating
# rule 16 as final.
STATIC_HOUSEKEEPING_FALLBACK = [
    "ACTB", "GAPDH", "TBP", "B2M", "RPL13A", "RPLP0", "HPRT1", "PPIA",
    "YWHAZ", "SDHA", "TUBB", "PGK1", "UBC", "EEF1A1",
]


def try_interferome_query():
    log.info("Attempting live Interferome query...")
    try:
        resp = requests.get(
            INTERFEROME_QUERY_URL,
            params={"species": "human", "stimulus": "interferon"},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        genes = sorted({g["geneSymbol"] for g in data.get("results", []) if "geneSymbol" in g})
        if not genes:
            raise ValueError("Interferome response parsed but yielded zero gene symbols")
        log.info(f"Interferome live query succeeded: {len(genes)} genes")
        return genes, "interferome_live"
    except Exception as e:
        log.warning(f"Interferome live query failed ({e}); falling back to static core ISG list")
        return STATIC_CORE_ISG_FALLBACK, "static_fallback_schoggins2011"


def fetch_gene_sets(isg_source, hk_source, out_path):
    if isg_source == "interferome":
        isg_genes, isg_provenance = try_interferome_query()
    else:
        isg_genes, isg_provenance = STATIC_CORE_ISG_FALLBACK, "static_fallback_schoggins2011"

    # Housekeeping: no live API attempted -- always static per the docstring.
    hk_genes = STATIC_HOUSEKEEPING_FALLBACK
    hk_provenance = "static_eisenberg_levanon_2013_placeholder"

    rows = [{"gene_id": g, "gene_set": "ISG"} for g in isg_genes]
    rows += [{"gene_id": g, "gene_set": "housekeeping"} for g in hk_genes]
    out = pd.DataFrame(rows)
    out.to_csv(out_path, sep="\t", index=False)

    log.info(f"Wrote {len(out)} gene-set rows -> {out_path}")
    log.info(f"ISG provenance: {isg_provenance} ({len(isg_genes)} genes)")
    log.info(f"Housekeeping provenance: {hk_provenance} ({len(hk_genes)} genes)")
    if "static" in isg_provenance or "static" in hk_provenance or "placeholder" in hk_provenance:
        log.warning(
            "One or both gene sets are PLACEHOLDER static lists, not the full "
            "curated source. Replace with the complete supplementary gene lists "
            "from Schoggins et al. 2011 (ISG) and Eisenberg & Levanon 2013 "
            "(housekeeping) before treating rule 16's Fisher's exact result as final."
        )
    return out


if __name__ == "__main__":
    fetch_gene_sets(
        isg_source=snakemake.params.isg_source,
        hk_source=snakemake.params.hk_source,
        out_path=snakemake.output.gene_sets,
    )
