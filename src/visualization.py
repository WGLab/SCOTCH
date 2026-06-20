from preprocessing import load_pickle, find_parent_isoform
import os
import shutil
import numpy as np
import pandas as pd
import pysam
import re
import subprocess
import warnings
from collections import defaultdict
from scipy.io import mmread

_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


# ---------------------------------------------------------------------------
# 1. Gene info & isoform selection
# ---------------------------------------------------------------------------

def _extract_tid(attr_str):
    """Return transcript_id from a GTF attribute string, if present."""
    match = _TRANSCRIPT_ID_RE.search(attr_str)
    return match.group(1) if match else None


def _extract_isoform_from_var(col_name, gene_name):
    """Return isoform name when a count-matrix column belongs exactly to gene_name."""
    prefix = f'{gene_name}_'
    if not str(col_name).startswith(prefix):
        return None
    return str(col_name)[len(prefix):]

def _find_filtered_gtf(target):
    """Locate the SCOTCH filtered GTF, tolerating the _N_P suffix variants."""
    ref_dir = os.path.join(target, 'reference')
    # Try exact name first (older runs)
    exact = os.path.join(ref_dir, 'SCOTCH_updated_annotation_filtered.gtf')
    if os.path.exists(exact):
        return exact
    # Look for pattern SCOTCH_updated_annotation_filtered_*.gtf
    candidates = [f for f in os.listdir(ref_dir)
                  if f.startswith('SCOTCH_updated_annotation_filtered') and f.endswith('.gtf')]
    if candidates:
        return os.path.join(ref_dir, sorted(candidates)[0])
    # Fall back to the unfiltered annotation
    fallback = os.path.join(ref_dir, 'SCOTCH_updated_annotation.gtf')
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(f'No SCOTCH annotation GTF found in {ref_dir}')


def _find_compatible_csv(target, gene_name):
    """Locate the compatible matrix CSV for a gene."""
    cm_dir = os.path.join(target, 'compatible_matrix')
    candidates = [f for f in os.listdir(cm_dir) if f.startswith(gene_name + '_')]
    if not candidates:
        return None
    return os.path.join(cm_dir, candidates[0])


def sub_gtf(gene_name, gtf_path):
    """Filter GTF to rows belonging to *gene_name*."""
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,
                       dtype={0: str})
    gtf.columns = ['seqname', 'source', 'feature', 'start', 'end',
                    'score', 'strand', 'frame', 'attribute']
    filtered = gtf[gtf['attribute'].str.contains(f'gene_name "{gene_name}"', regex=False)].copy()
    if filtered.empty:
        raise ValueError(f'Gene "{gene_name}" not found in {gtf_path}')
    transcript_ids = filtered['attribute'].apply(_extract_tid).dropna().unique().tolist()
    gene_chr = filtered['seqname'].iloc[0]
    gene_start = int(filtered['start'].min())
    gene_end = int(filtered['end'].max())
    gene_strand = filtered['strand'].iloc[0]
    return filtered, transcript_ids, gene_chr, gene_start, gene_end, gene_strand


def get_gene_info(gene_name, target, novel_pct=0.1):
    """Return gene coordinates, filtered GTF DataFrame, and known/novel isoform lists.

    Reuses the isoform-selection logic from the old visualisation code: keep
    all non-novel isoforms and novel isoforms whose expression fraction
    exceeds *novel_pct*.
    """
    gtf_path = _find_filtered_gtf(target)
    filtered_gtf, transcript_ids, gene_chr, gene_start, gene_end, gene_strand = sub_gtf(gene_name, gtf_path)

    # --- Identify major novel isoforms via count matrix ---
    count_dir = os.path.join(target, 'count_matrix')
    mtx_files = [f for f in os.listdir(count_dir)
                 if 'adata_transcript' in f and f.endswith('.mtx')]
    pkl_files = [f for f in os.listdir(count_dir)
                 if 'adata_transcript' in f and f.endswith('.pickle')]
    known_isoform_names = [t for t in transcript_ids if not str(t).startswith('novel')]
    novel_isoform_names = []

    if mtx_files and pkl_files:
        mtx = mmread(os.path.join(count_dir, mtx_files[0])).tocsc()
        pkl = load_pickle(os.path.join(count_dir, pkl_files[0]))
        gene_col_idx = []
        gene_isoforms = []
        for i, col in enumerate(pkl['var']):
            isoform_name = _extract_isoform_from_var(col, gene_name)
            if isoform_name is not None:
                gene_col_idx.append(i)
                gene_isoforms.append(isoform_name)
        if gene_col_idx:
            sub = mtx[:, gene_col_idx]
            col_sums = np.asarray(sub.sum(axis=0)).ravel()
            total = float(col_sums.sum())
            if total > 0:
                col_pct = col_sums / total
                over_thresh = np.flatnonzero(col_pct > novel_pct).tolist()
                novel_isoform_names = [
                    gene_isoforms[i]
                    for i in over_thresh
                    if gene_isoforms[i].startswith('novel')
                ]

    # Filter the GTF to only the selected isoforms
    selected = set(known_isoform_names + novel_isoform_names)
    tid_series = filtered_gtf['attribute'].apply(_extract_tid)
    # Keep gene-level rows (no transcript_id) plus rows with selected transcripts
    keep_mask = tid_series.isna() | tid_series.isin(selected)
    filtered_gtf = filtered_gtf[keep_mask].copy()

    return (gene_chr, gene_start, gene_end, gene_strand,
            filtered_gtf, known_isoform_names, novel_isoform_names)


# ---------------------------------------------------------------------------
# 2. Cell type map
# ---------------------------------------------------------------------------

def load_cell_type_map(cell_type_file):
    """Read a two-column CSV/TSV (barcode, cell_type). Returns dict or None."""
    if cell_type_file is None:
        return None
    sep = '\t' if cell_type_file.endswith('.tsv') else ','
    df = pd.read_csv(cell_type_file, sep=sep, header=0)
    df.columns = df.columns.str.strip()
    col0, col1 = df.columns[0], df.columns[1]
    return dict(zip(df[col0].astype(str), df[col1].astype(str)))


# ---------------------------------------------------------------------------
# 3. Classify reads into groups
# ---------------------------------------------------------------------------

def _build_cbumi_sets(target, gene_name, known_isoform_names, novel_isoform_names, gene_strand):
    """From the compatible matrix CSV, build CBUMI sets for known / each novel isoform."""
    csv_path = _find_compatible_csv(target, gene_name)
    if csv_path is None:
        return set(), {}
    df = pd.read_csv(csv_path)
    cols = df.columns.tolist()
    cols[0] = 'CBUMI'
    df.columns = cols

    # Known
    known_cols = [c for c in known_isoform_names if c in df.columns]
    cbumi_known = set()
    if known_cols:
        mask = df[known_cols].gt(0).any(axis=1)
        cbumi_known = set(df.loc[mask, 'CBUMI'])

    # Novel — group child isoforms under each major parent
    novel_cbumi = {}
    novel_cols_df = df.filter(like='novelIsoform')
    if not novel_cols_df.empty and novel_isoform_names:
        child_ids = [int(c.split('_')[1]) for c in novel_cols_df.columns]
        for nname in novel_isoform_names:
            nid = int(nname.split('_')[1])
            col_indices = [i for i, cid in enumerate(child_ids)
                           if find_parent_isoform([nid], cid, gene_strand) is not None]
            if col_indices:
                mask = novel_cols_df.iloc[:, col_indices].gt(0).any(axis=1)
                novel_cbumi[nname] = set(df.loc[mask, 'CBUMI'])

    return cbumi_known, novel_cbumi


def _resolve_bam_path(bam_path, gene_chr):
    """Handle BAM file vs BAM directory (per-chromosome)."""
    if os.path.isfile(bam_path):
        return bam_path
    candidates = [f for f in os.listdir(bam_path)
                  if f.endswith('.bam') and '.' + gene_chr + '.' in f]
    if candidates:
        return os.path.join(bam_path, candidates[0])
    raise FileNotFoundError(
        f'No BAM file found for {gene_chr} in directory {bam_path}')


def classify_reads(bam_path, gene_chr, gene_start, gene_end,
                   cbumi_known, novel_cbumi_dict,
                   cell_type_map, cell_types,
                   separate_known_novel, barcode_cell, barcode_umi,
                   sample_name):
    """Fetch reads from BAM and assign each to one or more group labels.

    Returns dict {group_label: [pysam.AlignedSegment, ...]}.
    """
    bam_file = _resolve_bam_path(bam_path, gene_chr)
    reads_by_group = defaultdict(list)
    all_novel_cbumis = set()
    for s in novel_cbumi_dict.values():
        all_novel_cbumis |= s
    use_tag_classification = barcode_cell is not None and barcode_umi is not None

    ct_filter = None
    if cell_type_map and cell_types and cell_types != ['bulk']:
        if cell_types != ['all']:
            ct_filter = set(cell_types)

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for idx, read in enumerate(bam.fetch(gene_chr, gene_start, gene_end)):
            cb = ub = None
            if use_tag_classification:
                if read.has_tag(barcode_cell):
                    cb = read.get_tag(barcode_cell)
                if read.has_tag(barcode_umi):
                    ub = read.get_tag(barcode_umi)

            cbumi = (
                f'{cb}_{ub}'
                if cb is not None and ub is not None
                else (read.query_name or f'untagged_read_{idx}')
            )

            # Determine cell type group only when a cell barcode is available.
            ct_label = None
            if cell_type_map and cb is not None:
                ct_label = cell_type_map.get(str(cb))
                if ct_label is None:
                    continue  # barcode not in annotation
                if ct_filter and ct_label not in ct_filter:
                    continue

            # Skip CB/UMI classification entirely for bulk/untagged reads.
            if not use_tag_classification or cb is None or ub is None:
                reads_by_group[sample_name].append(read)
                continue

            is_known = cbumi in cbumi_known
            is_novel = cbumi in all_novel_cbumis
            base = ct_label if ct_label else sample_name

            if separate_known_novel:
                if is_known:
                    reads_by_group[f'{base}::known'].append(read)
                elif is_novel:
                    reads_by_group[f'{base}::novel'].append(read)
                else:
                    reads_by_group[f'{base}::other'].append(read)
            else:
                reads_by_group[base].append(read)

    return dict(reads_by_group)


# ---------------------------------------------------------------------------
# 4. Coverage extraction
# ---------------------------------------------------------------------------

def extract_coverage(reads_by_group, gene_start, gene_end, bin_size,
                     overlay_known_novel=False):
    """Compute per-position depth per group from classified reads.

    Returns a DataFrame: pos, depth, group [, status].
    """
    length = gene_end - gene_start + 1
    rows = []

    for group, reads in reads_by_group.items():
        cov = np.zeros(length, dtype=np.int32)
        for read in reads:
            for start, end in read.get_blocks():
                s = max(start, gene_start) - gene_start
                e = min(end, gene_end + 1) - gene_start
                if s < e:
                    cov[s:e] += 1

        # Bin
        if bin_size > 1:
            n_bins = (length + bin_size - 1) // bin_size
            for bi in range(n_bins):
                s = bi * bin_size
                e = min(s + bin_size, length)
                pos = gene_start + s
                depth = float(np.mean(cov[s:e]))
                rows.append((pos, depth, group))
        else:
            nz = np.nonzero(cov)[0]
            rows.extend((gene_start + int(i), int(cov[i]), group) for i in nz)

    df = pd.DataFrame(rows, columns=['pos', 'depth', 'group'])

    if overlay_known_novel:
        # Split group "SampleA::known" → group="SampleA", status="known"
        parts = df['group'].str.split('::', n=1, expand=True)
        if parts.shape[1] == 2:
            df['group'] = parts[0]
            df['status'] = parts[1]
        else:
            df['status'] = 'all'

    return df


# ---------------------------------------------------------------------------
# 5. Junction extraction
# ---------------------------------------------------------------------------

def _get_exon_boundaries(filtered_gtf):
    """Extract a set of (exon_end, next_exon_start) pairs from GTF exon rows
    to identify annotated splice junctions."""
    exon_rows = filtered_gtf[filtered_gtf['feature'] == 'exon'].copy()
    boundaries = set()
    if exon_rows.empty:
        return boundaries
    # Get unique (exon_start, exon_end) per transcript
    tid_col = exon_rows['attribute'].apply(_extract_tid)
    exon_rows = exon_rows.assign(tid=tid_col)
    for _, grp in exon_rows.groupby('tid'):
        exons = grp[['start', 'end']].sort_values('start').values
        for i in range(len(exons) - 1):
            # Junction = (end_of_exon_i, start_of_exon_i+1)
            boundaries.add((int(exons[i, 1]), int(exons[i + 1, 0])))
    return boundaries


def extract_junctions(reads_by_group, exon_boundaries,
                      junction_min_reads, junction_annotated_only):
    """Parse CIGAR 'N' ops to find splice junctions.

    Returns DataFrame: start, end, count, group, is_annotated.
    """
    junction_counts = defaultdict(lambda: defaultdict(int))
    annotated_junctions = {
        (eb_end + delta_end, eb_start + delta_start)
        for eb_end, eb_start in exon_boundaries
        for delta_end in range(-2, 3)
        for delta_start in range(-2, 3)
    }

    for group, reads in reads_by_group.items():
        for read in reads:
            if read.cigartuples is None:
                continue
            ref_pos = read.reference_start
            for op, length in read.cigartuples:
                if op == 3:  # N = skipped region (intron)
                    junc_start = ref_pos
                    junc_end = ref_pos + length
                    junction_counts[group][(junc_start, junc_end)] += 1
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                    ref_pos += length

    rows = []
    for group, juncs in junction_counts.items():
        for (js, je), count in juncs.items():
            if count < junction_min_reads:
                continue
            annotated = (js, je) in annotated_junctions
            if junction_annotated_only and not annotated:
                continue
            rows.append((js, je, count, group, annotated))

    return pd.DataFrame(rows, columns=['start', 'end', 'count', 'group', 'is_annotated'])


# ---------------------------------------------------------------------------
# 6. Gene model
# ---------------------------------------------------------------------------

def build_gene_model(filtered_gtf, novel_isoform_names):
    """Convert GTF exon rows into a tidy DataFrame of exon/intron features.

    Returns DataFrame: transcript_id, feature, start, end, strand, is_novel.
    """
    exon_rows = filtered_gtf[filtered_gtf['feature'] == 'exon'].copy()
    if exon_rows.empty:
        return pd.DataFrame(columns=['transcript_id', 'feature', 'start', 'end',
                                      'strand', 'is_novel'])
    tid_col = exon_rows['attribute'].apply(_extract_tid)
    exon_rows = exon_rows.assign(tid=tid_col)
    exon_rows = exon_rows.dropna(subset=['tid'])

    novel_set = set(novel_isoform_names)
    rows = []
    for tid, grp in exon_rows.groupby('tid'):
        strand = grp['strand'].iloc[0]
        is_novel = tid in novel_set
        exons = grp[['start', 'end']].sort_values('start').values
        for i, (es, ee) in enumerate(exons):
            rows.append((tid, 'exon', int(es), int(ee), strand, is_novel))
            if i < len(exons) - 1:
                rows.append((tid, 'intron', int(ee), int(exons[i + 1, 0]),
                             strand, is_novel))

    return pd.DataFrame(rows, columns=['transcript_id', 'feature', 'start', 'end',
                                        'strand', 'is_novel'])


# ---------------------------------------------------------------------------
# 7. Write intermediate files
# ---------------------------------------------------------------------------

def write_intermediate_files(output_dir, gene_name, gene_chr, gene_start,
                             gene_end, gene_strand,
                             coverage_df, junctions_df, gene_model_df):
    """Write TSV files consumed by R/scotch_plot.R.

    Returns dict mapping file type to absolute path.
    """
    vis_dir = os.path.join(output_dir, 'visualization', gene_name)
    os.makedirs(vis_dir, exist_ok=True)

    paths = {}
    cov_path = os.path.join(vis_dir, 'coverage.tsv')
    coverage_df.to_csv(cov_path, sep='\t', index=False)
    paths['coverage'] = cov_path

    junc_path = os.path.join(vis_dir, 'junctions.tsv')
    junctions_df.to_csv(junc_path, sep='\t', index=False)
    paths['junctions'] = junc_path

    gm_path = os.path.join(vis_dir, 'gene_model.tsv')
    gene_model_df.to_csv(gm_path, sep='\t', index=False)
    paths['gene_model'] = gm_path

    meta_path = os.path.join(vis_dir, 'metadata.tsv')
    pd.DataFrame([{
        'gene_name': gene_name,
        'gene_chr': gene_chr,
        'gene_start': gene_start,
        'gene_end': gene_end,
        'gene_strand': gene_strand,
    }]).to_csv(meta_path, sep='\t', index=False)
    paths['metadata'] = meta_path

    return paths


# ---------------------------------------------------------------------------
# 8. Optional sub-BAM saving
# ---------------------------------------------------------------------------

def save_subbams(reads_by_group, bam_template_path, output_dir, gene_name,
                 gene_chr):
    """Write sorted + indexed sub-BAM per group."""
    vis_dir = os.path.join(output_dir, 'visualization', gene_name)
    os.makedirs(vis_dir, exist_ok=True)
    template_bam = pysam.AlignmentFile(_resolve_bam_path(bam_template_path, gene_chr), 'rb')

    for group, reads in reads_by_group.items():
        safe_name = group.replace('::', '_').replace(' ', '_')
        unsorted = os.path.join(vis_dir, f'{gene_name}_{safe_name}.bam')
        sorted_path = os.path.join(vis_dir, f'{gene_name}_{safe_name}_sorted.bam')
        with pysam.AlignmentFile(unsorted, 'wb', template=template_bam) as out:
            for read in reads:
                out.write(read)
        subprocess.run(['samtools', 'sort', '-o', sorted_path, unsorted],
                       check=True)
        subprocess.run(['samtools', 'index', sorted_path], check=True)
        os.remove(unsorted)

    template_bam.close()


# ---------------------------------------------------------------------------
# 9. Optional sub-GTF saving
# ---------------------------------------------------------------------------

def save_subgtf(filtered_gtf, output_dir, gene_name):
    """Write filtered GTF for the gene."""
    vis_dir = os.path.join(output_dir, 'visualization', gene_name)
    os.makedirs(vis_dir, exist_ok=True)
    out_path = os.path.join(vis_dir, f'{gene_name}_SCOTCH_filtered.gtf')
    filtered_gtf.to_csv(out_path, sep='\t', index=False, header=False, quoting=3)
    return out_path


# ---------------------------------------------------------------------------
# 10. Call R plotting script
# ---------------------------------------------------------------------------

def call_rscript(intermediate_files, output_format, width, height,
                 annotation_scale, overlay_known_novel, separate_known_novel):
    """Call R/scotch_plot.R via subprocess."""
    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'R', 'scotch_plot.R'
    )
    if not os.path.isfile(r_script):
        print(f'Warning: R script not found at {r_script}. '
              'Intermediate files have been saved; you can run the R script manually.')
        return

    if shutil.which('Rscript') is None:
        print('Warning: Rscript not found on PATH. '
              'Intermediate files have been saved; you can run the R script manually.')
        return

    vis_dir = os.path.dirname(intermediate_files['coverage'])
    gene_name = os.path.basename(vis_dir)
    output_path = os.path.join(vis_dir, f'{gene_name}_plot.{output_format}')

    cmd = [
        'Rscript', r_script,
        '--coverage', intermediate_files['coverage'],
        '--junctions', intermediate_files['junctions'],
        '--gene_model', intermediate_files['gene_model'],
        '--metadata', intermediate_files['metadata'],
        '--output', output_path,
        '--format', output_format,
        '--width', str(width),
        '--height', str(height),
        '--annotation_scale', str(annotation_scale),
    ]
    if overlay_known_novel:
        cmd.append('--overlay')
    if separate_known_novel:
        cmd.append('--separate')

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'R script failed (exit {result.returncode}):\n{result.stderr}')
    else:
        print(f'Plot saved to: {output_path}')


# ---------------------------------------------------------------------------
# 11. Main entry point
# ---------------------------------------------------------------------------

def visualization(gene, bam_files, targets, sample_names,
                  novel_pct, output_dir,
                  cell_type_file, cell_types,
                  separate_known_novel, overlay_known_novel,
                  save_bam, save_gtf,
                  junction_min_reads, junction_annotated_only,
                  bin_size, output_format,
                  width, height, annotation_scale,
                  barcode_cell, barcode_umi, **_ignored):
    """Orchestrate data extraction and R plotting for a single gene."""

    if sample_names is None:
        sample_names = [f'sample{i}' for i in range(len(targets))]
    if not (len(bam_files) == len(targets) == len(sample_names)):
        raise ValueError(
            'bam_files, targets, and sample_names must have the same length: '
            f'got {len(bam_files)}, {len(targets)}, and {len(sample_names)}'
        )

    cell_type_map = load_cell_type_map(cell_type_file)

    # Collect coverage / junction data across samples
    all_coverage = []
    all_junctions = []
    gene_chr = gene_start = gene_end = gene_strand = None
    filtered_gtf = None
    known_isoforms = []
    novel_isoforms = []

    for bam_path, target, sname in zip(bam_files, targets, sample_names):
        # 1. Gene info
        info = get_gene_info(gene, target, novel_pct)
        curr_chr, curr_start, curr_end, curr_strand = info[:4]
        curr_filtered_gtf = info[4]
        curr_known_isoforms = info[5]
        curr_novel_isoforms = info[6]

        if gene_chr is None:
            gene_chr, gene_start, gene_end, gene_strand = (
                curr_chr, curr_start, curr_end, curr_strand
            )
            filtered_gtf = curr_filtered_gtf
            known_isoforms = curr_known_isoforms
            novel_isoforms = curr_novel_isoforms
        else:
            current_meta = (curr_chr, curr_start, curr_end, curr_strand)
            reference_meta = (gene_chr, gene_start, gene_end, gene_strand)
            if current_meta != reference_meta:
                warnings.warn(
                    f'Gene metadata mismatch for target {target}: '
                    f'{current_meta} != {reference_meta}. '
                    'Using metadata from the first target.'
                )

        # 2. Build CBUMI sets
        cbumi_known, novel_cbumi = _build_cbumi_sets(
            target, gene, curr_known_isoforms, curr_novel_isoforms, curr_strand)

        # 3. Classify reads
        reads_by_group = classify_reads(
            bam_path, curr_chr, curr_start, curr_end,
            cbumi_known, novel_cbumi,
            cell_type_map, cell_types,
            separate_known_novel or overlay_known_novel,
            barcode_cell, barcode_umi, sname)

        # 4. Coverage
        cov_df = extract_coverage(reads_by_group, curr_start, curr_end,
                                  bin_size, overlay_known_novel)
        all_coverage.append(cov_df)

        # 5. Junctions
        exon_bounds = _get_exon_boundaries(curr_filtered_gtf)
        junc_df = extract_junctions(reads_by_group, exon_bounds,
                                    junction_min_reads, junction_annotated_only)
        all_junctions.append(junc_df)

        # Optional BAM saving
        if save_bam:
            save_subbams(reads_by_group, bam_path, output_dir, gene, curr_chr)

    # Merge across samples
    coverage_merged = pd.concat(all_coverage, ignore_index=True)
    junctions_merged = pd.concat(all_junctions, ignore_index=True)

    # Gene model from the first target's filtered GTF.
    gene_model_df = build_gene_model(filtered_gtf, novel_isoforms)

    # Write intermediate files
    paths = write_intermediate_files(
        output_dir, gene, gene_chr, gene_start, gene_end, gene_strand,
        coverage_merged, junctions_merged, gene_model_df)

    # Optional GTF saving
    if save_gtf and filtered_gtf is not None:
        save_subgtf(filtered_gtf, output_dir, gene)

    # Call R
    call_rscript(paths, output_format, width, height, annotation_scale,
                 overlay_known_novel, separate_known_novel)
