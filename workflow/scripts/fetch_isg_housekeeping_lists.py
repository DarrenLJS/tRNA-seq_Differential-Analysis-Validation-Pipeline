"""
workflow/scripts/fetch_isg_housekeeping_lists.py

Builds the ISG (Interferome v2.0, Rusinova et al. 2013, NAR) and
housekeeping (Eisenberg & Levanon 2013, Trends Genet) gene ID lists for
the rule-16 Fisher's exact enrichment test, from locally-staged export
files -- NOT a live API call.

WHY THIS IS A LOCAL-FILE STEP, NOT A LIVE FETCH
-------------------------------------------------
  - ISG list, RECOMMENDED SOURCE (CHANGED 2026-07-17): Interferome v2.0
    was found to be unreliable to reach in practice (see below), so the
    default source is now MSigDB's HALLMARK_INTERFERON_ALPHA_RESPONSE
    gene set (Liberzon et al. 2015, Cell Systems, doi:10.1016/j.cels.
    2015.12.004) -- 97 human genes up-regulated in response to alpha
    interferon, curated by the Broad Institute from real IFN-alpha-vs-
    control expression datasets (GSE31019, GSE31193, GSE43723), CC-BY-4.0
    licensed, and downloadable with NO login and NO browser session at:
    https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_INTERFERON_ALPHA_RESPONSE&fileType=grp
    Confirmed live and fetchable 2026-07-17. This is Type I IFN-alpha
    response rather than Interferome's IFN-beta-specific query -- alpha
    and beta both signal through the same ISGF3/ISRE pathway with heavily
    overlapping downstream programs, so this is a defensible substitute,
    but it IS a different curation methodology/gene count than what the
    proposal originally specified (Interferome, IFN-beta-filtered) --
    name this substitution explicitly in the methods section rather than
    presenting it as if it were the original Interferome query.
  - Interferome v2.0 (the originally proposed source) is a browser
    query/download tool, not an API -- CHECKED 2026-07-17: both
    interferome.its.monash.edu.au/interferome/ and
    interferome.org/interferome/ resolve to real, populated pages
    (Search, Database Statistics, How To Cite, a User Login page) in
    search-result snippets, so the DB itself does not appear down, but it
    was still unreachable in practice (confirmed by the user after
    retrying) -- if you want to keep using Interferome specifically
    instead of the MSigDB substitute above, retry later, try both
    domains, or check whether an account is required for search/export.
  - Eisenberg & Levanon (2013) "Human housekeeping genes, revisited" is
    NOT distributed as a journal supplementary file -- CONFIRMED
    2026-07-17: the paper's own Table 1 is only a small ~12-gene
    calibration subset ("genes proposed for calibration"), not the full
    housekeeping list. The full list of 3804 genes is hosted directly by
    the authors at http://www.tau.ac.il/~elieis/HKG/HK_genes.txt
    (confirmed live, no login), tab-separated gene_symbol + RefSeq
    accession, with NO header row in the raw file -- prepend a header
    line "gene_symbol\trefseq_id" after downloading it so this script's
    default first-column reader picks up the right column name.
    ALTERNATIVE, ALSO CONFIRMED WORKING 2026-07-17: MSigDB's
    HSIAO_HOUSEKEEPING_GENES gene set (Hsiao et al. 2001, Physiol
    Genomics, doi:10.1152/physiolgenomics.00025.2001) -- 395 genes
    identified as expressed across 19 normal human tissues, no login
    required: https://www.gsea-msigdb.org/gsea/msigdb/human/
    download_geneset.jsp?geneSetName=HSIAO_HOUSEKEEPING_GENES&fileType=grp
    This is an OLDER, SMALLER, DIFFERENT-METHODOLOGY list than Eisenberg
    & Levanon (2001 microarray-based presence/absence across 19 tissues,
    vs. 2013 RNA-seq-based expression-stability across a larger panel) --
    if the tau.ac.il source is reachable at all, prefer it as the closer
    match to what the proposal specified; use Hsiao only if tau.ac.il is
    also blocked, and name the substitution explicitly in the methods
    section either way.
Silently substituting a small hardcoded gene list when a live query fails
would mislabel the output as "Interferome"/"Eisenberg & Levanon" derived
when it isn't -- misleading for a methods section. This script instead
fails loudly with instructions until the real files are staged locally.

EXPECTED INPUT FILES (see config stage2_references.isg_list_path /
housekeeping_list_path)
------------------------------------------------------------------
  - ISG list: plain text/CSV/TSV with one gene symbol per row (whatever
    column Interferome's export uses -- pass the column name via
    stage2_references.isg_gene_col if it isn't the first column).
  - Housekeeping list: same, for Eisenberg & Levanon's supplementary
    table (pass housekeeping_gene_col if needed).

OUTPUT
------
Long-format TSV, columns [gene_id, gene_symbol, gene_set], gene_set in
{"ISG", "housekeeping"}.

GENE ID FIX (2026-07-18)
-------------------------
Both source lists are gene-SYMBOL keyed (MSigDB HALLMARK_INTERFERON_
ALPHA_RESPONSE and the Eisenberg & Levanon housekeeping list both use
gene_symbol as their identifier), but this previously flowed straight
into an output column named "gene_id" -- silently degenerating rule 16's
Fisher's exact test (the join against Ensembl-ID-keyed
gene_translation_scores left gene_set as NA for every row). Symbols are
now mapped to Ensembl gene IDs via references.ensembl_gtf (same GTF
build_codon_usage_table.py / fetch_watson_polysome_data.py already
parse). Unmapped/ambiguous symbols are dropped and counted/logged.
"""

import logging
import os
import gzip
import re

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def build_symbol_to_ensembl_map(gtf_path):
    """Parse a gzipped Ensembl GTF's 'gene' feature lines into a
    gene_name (symbol) -> gene_id (Ensembl, version stripped) map.
    Ambiguous symbols (same gene_name -> >1 distinct gene_id) are
    dropped rather than resolved arbitrarily.
    """
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')

    symbol_to_ids = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            id_match = gene_id_re.search(attrs)
            name_match = gene_name_re.search(attrs)
            if not id_match or not name_match:
                continue
            gene_id = id_match.group(1).split(".")[0]
            gene_name = name_match.group(1)
            symbol_to_ids.setdefault(gene_name, set()).add(gene_id)

    ambiguous = {sym: ids for sym, ids in symbol_to_ids.items() if len(ids) > 1}
    if ambiguous:
        log.warning(
            f"{len(ambiguous)} gene symbol(s) in the GTF map to more than one "
            "Ensembl gene_id -- dropped from the symbol->Ensembl map."
        )

    symbol_to_ensembl = {
        sym: next(iter(ids)) for sym, ids in symbol_to_ids.items() if len(ids) == 1
    }
    log.info(
        f"Built symbol->Ensembl map from {gtf_path}: "
        f"{len(symbol_to_ensembl)} unambiguous symbols "
        f"({len(ambiguous)} ambiguous symbols dropped)"
    )
    return symbol_to_ensembl


def _read_gene_list(path, gene_col, label, instructions):
    if not path or not os.path.exists(path):
        raise RuntimeError(
            f"{label} gene list not found at '{path}'.\n{instructions}"
        )
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    if gene_col:
        if gene_col not in df.columns:
            raise RuntimeError(
                f"{label} file '{path}' has columns {list(df.columns)}, "
                f"but configured gene column '{gene_col}' isn't among them. "
                "Fix stage2_references.isg_gene_col / housekeeping_gene_col."
            )
        genes = df[gene_col]
    else:
        genes = df.iloc[:, 0]
    genes = sorted(set(genes.dropna().astype(str).str.strip()))
    log.info(f"{label}: loaded {len(genes)} genes from {path}")
    return genes


ISG_INSTRUCTIONS = (
    "RECOMMENDED: download MSigDB's HALLMARK_INTERFERON_ALPHA_RESPONSE gene "
    "set (Liberzon et al. 2015, Cell Systems, doi:10.1016/j.cels.2015.12.004) "
    "-- 97 genes, no login required, confirmed live 2026-07-17: "
    "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
    "?geneSetName=HALLMARK_INTERFERON_ALPHA_RESPONSE&fileType=grp -- save "
    "it, strip the first two header/comment lines, add back a single "
    "'gene_symbol' header line, and stage it at the path given in config "
    "stage2_references.isg_list_path. Note in the methods section that this "
    "replaces the originally planned Interferome (IFN-beta-specific) query "
    "with a Type I IFN-alpha response set, since Interferome was found to "
    "be unreachable. ALTERNATIVE: if you'd rather keep using Interferome "
    "v2.0 (Rusinova et al. 2013, doi:10.1093/nar/gks1215) specifically, try "
    "http://interferome.its.monash.edu.au/interferome/ or "
    "http://interferome.org/interferome/ (both confirmed live 2026-07-17 "
    "in search results; a login page exists, so export may need a free "
    "account) -> Search -> Species: Human, Interferon Type: Type I, "
    "Sub-type: IFN-beta -> export -> stage it at the same config path."
)

HOUSEKEEPING_INSTRUCTIONS = (
    "PRIMARY: this is NOT a journal supplementary file -- the paper's own "
    "Table 1 is only a small ~12-gene calibration subset. Download the "
    "full 3804-gene list directly from the authors (Eisenberg & Levanon "
    "2013, 'Human housekeeping genes, revisited', Trends Genet "
    "29(10):569-574, doi:10.1016/j.tig.2013.05.010) at "
    "http://www.tau.ac.il/~elieis/HKG/HK_genes.txt (confirmed live, no "
    "login, 2026-07-17). The raw file has NO header row (tab-separated "
    "gene_symbol, RefSeq accession) -- prepend a header line "
    "'gene_symbol\\trefseq_id' before staging it at the path given in "
    "config stage2_references.housekeeping_list_path, or this script's "
    "default first-column reader will silently drop the first gene. "
    "ALTERNATIVE if tau.ac.il is unreachable: MSigDB's "
    "HSIAO_HOUSEKEEPING_GENES set (Hsiao et al. 2001, 395 genes, no "
    "login, confirmed live 2026-07-17): "
    "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
    "?geneSetName=HSIAO_HOUSEKEEPING_GENES&fileType=grp -- older/smaller/"
    "different methodology than Eisenberg & Levanon, name the "
    "substitution explicitly in the methods section if used."
)


def fetch_gene_sets(
    isg_list_path, housekeeping_list_path,
    isg_gene_col, housekeeping_gene_col,
    isg_source, housekeeping_source,
    gtf_path, out_path,
):
    isg_genes = _read_gene_list(
        isg_list_path, isg_gene_col, f"ISG ({isg_source})", ISG_INSTRUCTIONS
    )
    hk_genes = _read_gene_list(
        housekeeping_list_path, housekeeping_gene_col,
        f"Housekeeping ({housekeeping_source})", HOUSEKEEPING_INSTRUCTIONS,
    )

    symbol_to_ensembl = build_symbol_to_ensembl_map(gtf_path)

    def _map_genes(symbols, label):
        mapped = []
        n_unmapped = 0
        for sym in symbols:
            ens = symbol_to_ensembl.get(sym)
            if ens is None:
                n_unmapped += 1
                continue
            mapped.append((ens, sym))
        if n_unmapped:
            log.warning(
                f"{label}: {n_unmapped} of {len(symbols)} gene symbols did not map "
                "to an Ensembl gene_id via the GTF (not found, or ambiguous) -- "
                "dropped."
            )
        return mapped

    isg_mapped = _map_genes(isg_genes, f"ISG ({isg_source})")
    hk_mapped = _map_genes(hk_genes, f"Housekeeping ({housekeeping_source})")

    rows = [
        {"gene_id": ens, "gene_symbol": sym, "gene_set": "ISG"}
        for ens, sym in isg_mapped
    ]
    rows += [
        {"gene_id": ens, "gene_symbol": sym, "gene_set": "housekeeping"}
        for ens, sym in hk_mapped
    ]
    out = pd.DataFrame(rows)
    out.to_csv(out_path, sep="\t", index=False)

    log.info(f"Wrote {len(out)} gene-set rows -> {out_path}")
    log.info(f"ISG: {len(isg_mapped)} genes mapped to Ensembl (source: {isg_source})")
    log.info(
        f"Housekeeping: {len(hk_mapped)} genes mapped to Ensembl "
        f"(source: {housekeeping_source})"
    )
    if isg_source != "interferome":
        log.warning(
            f"ISG source is '{isg_source}', not the originally proposed Interferome "
            "v2.0 IFN-beta query -- note this substitution explicitly in the methods "
            "section (see this script's docstring for why)."
        )
    if housekeeping_source != "eisenberg_levanon_2013":
        log.warning(
            f"Housekeeping source is '{housekeeping_source}', not the originally "
            "proposed Eisenberg & Levanon (2013) list -- note this substitution "
            "explicitly in the methods section (see this script's docstring for why)."
        )
    return out


if __name__ == "__main__":
    fetch_gene_sets(
        isg_list_path=snakemake.params.isg_list_path,
        housekeeping_list_path=snakemake.params.housekeeping_list_path,
        isg_gene_col=snakemake.params.isg_gene_col,
        housekeeping_gene_col=snakemake.params.housekeeping_gene_col,
        isg_source=snakemake.params.isg_source,
        housekeeping_source=snakemake.params.housekeeping_source,
        gtf_path=snakemake.input.gtf,
        out_path=snakemake.output.gene_sets,
    )
