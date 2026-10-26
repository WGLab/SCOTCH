"""GTF attribute access and SCOTCH's internal/external transcript IDs."""
import re


def get_gtf_attribute(attributes, name):
    if not isinstance(attributes, str):
        return None
    match = re.search(r'(?:^|;)\s*' + re.escape(name) + r'\s+"([^"]*)"', attributes)
    return match.group(1) if match else None


def get_scotch_transcript_id(attributes):
    """Use the internal ID for matrix/PKL matching, with legacy GTF fallback."""
    return (get_gtf_attribute(attributes, 'scotch_transcript_id')
            or get_gtf_attribute(attributes, 'transcript_id'))


def transcript_attributes(gene_id, isoform, novel=False):
    if novel:
        return (f'transcript_id "{gene_id}_{isoform}"; '
                f'transcript_name "{isoform}"; scotch_transcript_id "{isoform}"')
    return f'transcript_id "{isoform}"; transcript_name "{isoform}"'


def qualify_novel_transcript(attributes):
    """Upgrade existing SCOTCH novel records without changing other attributes."""
    isoform = get_scotch_transcript_id(attributes)
    gene_id = get_gtf_attribute(attributes, 'gene_id')
    if not gene_id or not isoform or not isoform.startswith('novelIsoform_'):
        return attributes
    external_id = f'{gene_id}_{isoform}'
    attributes = re.sub(r'((?:^|;)\s*transcript_id\s+")[^"]*(")',
                        lambda match: match.group(1) + external_id + match.group(2), attributes)
    if get_gtf_attribute(attributes, 'scotch_transcript_id') is None:
        attributes = attributes.rstrip().rstrip(';') + f'; scotch_transcript_id "{isoform}";'
    return attributes
