"""Test GTF export/consumers without importing unrelated BAM/plotting packages."""
import ast
import csv
import logging
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest

import pandas as pd
from joblib import Parallel, delayed

SRC = Path(__file__).resolve().parents[1] / 'src'
sys.path.insert(0, str(SRC))
from gtf_utils import (get_gtf_attribute, get_scotch_transcript_id,
                       transcript_attributes, qualify_novel_transcript)
import reference


def load_functions(filename, names):
    tree = ast.parse((SRC / filename).read_text())
    nodes = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name in names]
    namespace = dict(pd=pd, os=os, csv=csv, Parallel=Parallel, delayed=delayed,
                     get_gtf_attribute=get_gtf_attribute,
                     get_scotch_transcript_id=get_scotch_transcript_id,
                     transcript_attributes=transcript_attributes,
                     qualify_novel_transcript=qualify_novel_transcript)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SRC / filename), 'exec'), namespace)
    return namespace


class GtfTranscriptTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        helpers = load_functions('preprocessing.py', ['merge_exons', 'get_transcript_id', 'add_novel_isoform_to_gtf'])
        self.helpers = helpers
        self.export = load_functions('compatible.py', ['convert_to_gtf'])
        self.export['merge_exons'] = helpers['merge_exons']

    def export_fixture(self, existing=None):
        data = {}
        for gene_id in ('G1', 'G2'):
            data[gene_id] = [[{
                'geneID': gene_id, 'geneName': 'SAME_NAME', 'geneChr': 'chr1',
                'geneStart': 100, 'geneEnd': 250, 'geneStrand': '+',
            }, [(100, 150), (200, 250)], {'novelIsoform_3': [0, 1], f'KNOWN_{gene_id}': [0]}]]
        path = self.root / 'export.gtf'
        self.export['convert_to_gtf'](data, str(path), gtf_df=existing)
        df = pd.read_csv(path, sep='\t', header=None, names=[
            'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
        return path, df

    def test_export_unique_across_same_named_genes_and_roundtrip(self):
        path, df = self.export_fixture()
        transcript_rows = df[df.feature == 'transcript']
        ids = transcript_rows.attribute.map(lambda a: get_gtf_attribute(a, 'transcript_id')).tolist()
        self.assertCountEqual(ids, ['G1_novelIsoform_3', 'G2_novelIsoform_3', 'KNOWN_G1', 'KNOWN_G2'])
        for _, row in df[df.feature != 'gene'].iterrows():
            gene = get_gtf_attribute(row.attribute, 'gene_id')
            internal = get_scotch_transcript_id(row.attribute)
            if internal == 'novelIsoform_3':
                self.assertEqual(get_gtf_attribute(row.attribute, 'transcript_id'), f'{gene}_{internal}')
            self.assertEqual(reference.parse_info(row.attribute)['transcript_id'], internal)
        # Exporting an already exported GTF must not duplicate transcripts or prefixes.
        _, again = self.export_fixture(existing=df)
        pd.testing.assert_frame_equal(df, again)

    def test_legacy_novel_records_are_upgraded_on_export(self):
        _, df = self.export_fixture()
        for i, attrs in df.attribute.items():
            internal = get_gtf_attribute(attrs, 'scotch_transcript_id')
            if internal:
                gene = get_gtf_attribute(attrs, 'gene_id')
                df.at[i, 'attribute'] = (attrs.replace(f'transcript_id "{gene}_{internal}"',
                                                     f'transcript_id "{internal}"')
                                        .replace(f'; scotch_transcript_id "{internal}"', ''))
        _, upgraded = self.export_fixture(existing=df)
        self.assertEqual(len(upgraded), len(df))
        ids = upgraded[upgraded.feature == 'transcript'].attribute.map(
            lambda a: get_gtf_attribute(a, 'transcript_id'))
        self.assertEqual(ids.nunique(), 4)

    def test_attribute_order_and_exact_keys(self):
        for attrs in (
            'transcript_id "G1_novelIsoform_3"; scotch_transcript_id "novelIsoform_3";',
            'scotch_transcript_id "novelIsoform_3"; transcript_id "G1_novelIsoform_3";',
        ):
            self.assertEqual(get_gtf_attribute(attrs, 'transcript_id'), 'G1_novelIsoform_3')
            self.assertEqual(self.helpers['get_transcript_id'](attrs), 'novelIsoform_3')
        self.assertIsNone(get_gtf_attribute('other_transcript_id "WRONG";', 'transcript_id'))
        self.assertEqual(get_scotch_transcript_id('transcript_id "ENST1";'), 'ENST1')
        self.assertEqual(get_scotch_transcript_id('transcript_id "novelIsoform_3";'), 'novelIsoform_3')

    def test_count_filter_uses_internal_id_and_preserves_external_id(self):
        path, _ = self.export_fixture()
        directory = self.root / 'reference'
        directory.mkdir()
        (directory / 'SCOTCH_updated_annotation.gtf').write_text(path.read_text())
        functions = load_functions('count_matrix.py', ['_extract_attribute', 'filter_gtf'])
        instance = types.SimpleNamespace(target=[str(self.root)], novel_read_n_str='1',
                                        novel_read_pct_str='0', logger=logging.getLogger(__name__),
                                        novel_isoform_del_dict={'SAME_NAME': ['novelIsoform_3']})
        instance._extract_attribute = types.MethodType(functions['_extract_attribute'], instance)
        functions['filter_gtf'](instance)
        output = (directory / 'SCOTCH_updated_annotation_filtered_1_0.gtf').read_text()
        self.assertNotIn('novelIsoform_3', output)
        self.assertIn('transcript_id "KNOWN_G1"', output)
        self.assertIn('transcript_id "KNOWN_G2"', output)

    def test_visualization_matches_internal_novel_names(self):
        _, df = self.export_fixture()
        functions = load_functions('visualization.py', ['_extract_tid', 'build_gene_model'])
        model = functions['build_gene_model'](df, ['novelIsoform_3'])
        novel = model[model.is_novel]
        self.assertFalse(novel.empty)
        self.assertEqual(set(novel.transcript_id), {'novelIsoform_3'})

    def test_add_novel_isoform_writes_both_ids(self):
        path, _ = self.export_fixture()
        output = self.root / 'added.gtf'
        self.helpers['add_novel_isoform_to_gtf'](
            str(path), 'G1', 'SAME_NAME', [([(101, 150)], 'SAME_NAME-novelIsoformGroup03')], str(output))
        text = output.read_text()
        self.assertIn('transcript_id "G1_SAME_NAME-novelIsoformGroup03"', text)
        self.assertIn('scotch_transcript_id "SAME_NAME-novelIsoformGroup03"', text)


if __name__ == '__main__':
    unittest.main()
