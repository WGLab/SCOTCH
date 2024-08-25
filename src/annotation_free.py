from collections import defaultdict
import pysam




def get_continuous_blocks(data):
    blocks = []
    if len(data) == 0:
        return blocks
    if isinstance(data, dict):
        sorted_positions = sorted(data.keys())
    elif isinstance(data, list):
        sorted_positions = sorted(data)
    else:
        raise TypeError("Input must be a dictionary or a list.")
    current_start = sorted_positions[0]
    # for a single base, pysam uses 0-based
    # we will use the same 0-based half open system as in bed file: e.g. (0,3): 1st,2nd,3rd bases
    current_end = sorted_positions[0] + 1
    for pos in sorted_positions[1:]:
        if pos == current_end:
            current_end = pos + 1
        else:
            blocks.append((current_start, current_end))
            current_start = pos
            current_end = pos + 1
    blocks.append((current_start, current_end))
    return blocks

def count_and_sort_tuple_elements(tuple_list):
    freq_dict = defaultdict(int)
    for a, b in tuple_list:
        freq_dict[a] += 1
        freq_dict[b] += 1
    sorted_freq_dict = dict(sorted(freq_dict.items()))
    return sorted_freq_dict


def merge_boundaries_by_evidence(boundaries, merge_distance=10):
    # Sort boundaries by value (evidence) in descending order
    sorted_boundaries = sorted(boundaries.items(), key=lambda item: item[1], reverse=True)
    merged_boundaries = {}
    while sorted_boundaries:
        # Take the boundary with the highest evidence
        base_position, base_value = sorted_boundaries.pop(0)
        merged_value = base_value
        to_remove = []
        # Check nearby positions within the merge_distance
        for i, (position, value) in enumerate(sorted_boundaries):
            if abs(position - base_position) <= merge_distance:
                merged_value += value
                to_remove.append(i)
        # Remove merged positions from sorted_boundaries
        for index in sorted(to_remove, reverse=True):
            sorted_boundaries.pop(index)
        # Add the merged boundary to the result
        merged_boundaries[base_position] = merged_value
    return merged_boundaries


def get_splicejuction_from_read(read):
    junctions = []
    ref_pos = read.reference_start
    for cigartype, cigarlen in read.cigartuples:
        if cigartype == 3:  # 'N' operation indicates a splice junction
            junction_start = ref_pos
            junction_end = ref_pos + cigarlen
            junctions.append((junction_start, junction_end))
        if cigartype in (0, 2, 3):
            ref_pos += cigarlen
    return junctions #this returns introns e.g. (103, 105) is 0 based containing 2 bases



###main function
##TODO: update existing annotation, generate annotation, no repeat read bam file
def get_non_overlapping_exons(bam_file, chrom, gene_start, gene_end, coverage_threshold=20):
    bam = pysam.AlignmentFile(bam_file, "rb")
    exons, coverage, junctions = [], {}, []
    #get read coverage
    for read in bam.fetch(chrom, gene_start, gene_end):
        if read.reference_start >= gene_start and read.reference_end < gene_end:
            junctions += get_splicejuction_from_read(read)
            for read_start, read_end in read.get_blocks():
                for pos in range(read_start, read_end):
                    if pos in coverage:
                        coverage[pos] += 1
                    else:
                        coverage[pos] = 1
    #only keep positions that meet the coverage threshold
    coverage = {pos: cov for pos, cov in coverage.items() if cov > coverage_threshold}
    coverage_blocks = get_continuous_blocks(coverage)
    coverage_blocks = [(a,b) for a, b in coverage_blocks if b-a >= 20]#delete less than 20bp exons
    if len(coverage_blocks)==0:
        return exons
    #use splicing junctions to find sub-exons
    boundaries = count_and_sort_tuple_elements(junctions)
    filtered_boundaries = [{} for _ in coverage_blocks]
    for boundary, freq in boundaries.items():
        for i, (start, end) in enumerate(coverage_blocks):
            if start < boundary < end:
                filtered_boundaries[i][boundary] = freq
                break
    Filtered_Boundaries = []
    for filtered_boundaries_ in filtered_boundaries:
        merged_boundaries = merge_boundaries_by_evidence(filtered_boundaries_, merge_distance=10)
        sorted_boundaries = dict(sorted(merged_boundaries.items()))
        filtered_sorted_boundaries = {k: v for k, v in sorted_boundaries.items() if v > 20}
        Filtered_Boundaries.append(filtered_sorted_boundaries)
    for i in range(len(coverage_blocks)):
        start, end = coverage_blocks[i]
        positions = [start]+list(Filtered_Boundaries[i].keys())+[end]
        for ii in range(len(positions)-1):
            exons.append((positions[ii], positions[ii]))
    bam.close()
    return exons


def get_genes_from_bam(bam_file, coverage_threshold = 5, min_region_size=50):
    bam = pysam.AlignmentFile(bam_file, "rb")
    coverage = defaultdict(lambda: defaultdict(int))
    chromosomes = bam.references
    for chrom in chromosomes:
        for read in bam.fetch(chrom):
            if not read.is_unmapped:
                for pos in range(read.reference_start, read.reference_end):
                    coverage[chrom][pos] += 1
    bam.close()
    genes = {}
    for chrom, cov_dict in coverage.items():
        cov_dict = {pos: cov for pos, cov in cov_dict.items() if cov > coverage_threshold}
        coverage_blocks = get_continuous_blocks(cov_dict)
        genes[chrom] = [(a, b) for a, b in coverage_blocks if b - a > min_region_size]
    return genes
