"""Exercise real TSV/PKL summary I/O without the unrelated BAM dependencies."""
import ast
import logging
import os
from pathlib import Path
import pickle
import re
import tempfile
import unittest

from joblib import Parallel, delayed
import numpy as np
import pandas as pd


def load_summary_functions():
    source = Path(__file__).resolve().parents[1] / 'src' / 'compatible.py'
    names = {
        'prepare_read_selection_df', 'build_read_selection_df',
        'build_read_selection_mapping', 'read_auxiliary_mapping_file',
        'is_per_gene_auxiliary_file', 'get_gene_name_from_auxiliary_filename',
        'summarise_auxiliary', '_summarise_auxiliary_sample',
    }
    constants = {'READ_MAPPING_COLUMNS', 'AUXILLARY_MAPPING_SUFFIX', 'MERGED_MAPPING_FILENAME'}
    nodes = [node for node in ast.parse(source.read_text()).body
             if (isinstance(node, ast.FunctionDef) and node.name in names)
             or (isinstance(node, ast.Assign) and any(
                 isinstance(target, ast.Name) and target.id in constants for target in node.targets))]
    namespace = dict(os=os, pickle=pickle, re=re, pd=pd, np=np,
                     Parallel=Parallel, delayed=delayed)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), namespace)
    return namespace['summarise_auxiliary']


class ParseSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.summarise = load_summary_functions()
        self.logger = logging.getLogger('parse-summary-tests')

    def write_mapping(self, target, gene='GENE', isoform='TX1', folder='auxiliary', read='shared_read'):
        directory = target / folder
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f'{gene}_ID1_read_isoform_exon_mapping.tsv'
        pd.DataFrame([{
            'Read': read, 'Isoform': isoform, 'Exon Index': '0',
            'Exon Coordinates': '100,200', 'Cell': 'CELL', 'Umi': 'UMI',
            'CBUMI': 'CELL_UMI', 'geneName': gene, 'geneID': 'ID1',
            'geneChr': 'chr1', 'MappingScore': 1,
        }]).to_csv(path, sep='\t', index=False)
        return path

    def summary(self, target=None, **kwargs):
        self.summarise(str(target or self.root), logger=self.logger, **kwargs)

    def assert_mapping(self, target, gene='GENE', isoform='TX1'):
        folder = target / 'auxiliary'
        df = pd.read_csv(folder / 'all_read_isoform_exon_mapping.tsv', sep='\t')
        self.assertEqual(df['geneName'].tolist(), [gene])
        self.assertEqual(df['Isoform'].tolist(), [isoform])
        self.assertEqual(df['Keep'].tolist(), [1])
        with (folder / 'read_selection.pkl').open('rb') as handle:
            self.assertEqual(pickle.load(handle), {'CELL_UMI': (f'{gene}_ID1', isoform)})

    def test_parse_samples_are_separate_and_rerun_is_safe(self):
        a, b = self.root / 'samples' / 'A', self.root / 'samples' / 'B'
        first = self.write_mapping(a, isoform='TX_A')
        second = self.write_mapping(b, isoform='TX_B', folder='auxillary')
        (self.root / 'samples' / 'notes.txt').write_text('not a sample')
        # A top-level auxiliary directory must not hide the nested samples.
        top = self.write_mapping(self.root, isoform='TOP')
        self.summary(platform='parse-ont')
        self.assert_mapping(a, isoform='TX_A')
        self.assert_mapping(b, isoform='TX_B')
        self.assertFalse(first.exists())
        self.assertFalse(second.exists())
        self.assertTrue(top.exists())
        before = {p: p.read_bytes() for p in (a / 'auxiliary').iterdir()}
        self.summary(platform='parse-ont')
        self.assertEqual(before, {p: p.read_bytes() for p in (a / 'auxiliary').iterdir()})
        self.assert_mapping(b, isoform='TX_B')

    def test_parse_without_top_level_auxiliary(self):
        sample = self.root / 'samples' / '1'
        self.write_mapping(sample)
        self.summary(platform='parse-ont')
        self.assert_mapping(sample)
        self.assertFalse((self.root / 'auxiliary').exists())

    def test_subset_update_reaches_each_sample(self):
        samples = [self.root / 'samples' / name for name in ('A', 'B')]
        for sample in samples:
            self.write_mapping(sample)
        self.summary(platform='parse-ont')
        for sample in samples:
            self.write_mapping(sample, isoform='TX_UPDATED')
        self.summary(platform='parse-ont', gene_subset=['GENE'])
        for sample in samples:
            self.assert_mapping(sample, isoform='TX_UPDATED')

    def test_10x_uses_top_level_even_with_samples_directory(self):
        self.write_mapping(self.root)
        nested = self.write_mapping(self.root / 'samples' / 'A')
        self.summary(platform='10x-ont')
        self.assert_mapping(self.root)
        self.assertTrue(nested.exists())

    def test_direct_parse_sample_target(self):
        sample = self.root / 'samples' / 'A'
        self.write_mapping(sample)
        self.summary(target=sample, platform='parse-ont')
        self.assert_mapping(sample)


if __name__ == '__main__':
    unittest.main()
