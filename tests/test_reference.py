"""Run with: python -m unittest discover -s tests -v"""
import itertools
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import reference


class StopBeforeExternalTools(Exception):
    pass


class ReferenceTests(unittest.TestCase):
    def test_attribute_order_does_not_change_bed_output(self):
        for type_key in ('gene_type', 'gene_biotype'):
            values = {'gene_id': 'G1', type_key: 'protein_coding',
                      'gene_name': 'TRBV19', 'transcript_id': 'T1'}
            for order in itertools.permutations(values):
                for kinds in (('gene', 'exon'), ('exon', 'gene')):
                    with self.subTest(type_key=type_key, order=order, kinds=kinds):
                        with tempfile.TemporaryDirectory() as temp:
                            path = Path(temp) / 'reference.gtf'
                            rows = []
                            for kind in kinds:
                                attrs = ' '.join(
                                    f'{key} "{values[key]}";' for key in order
                                    if kind == 'exon' or key != 'transcript_id')
                                rows.append(f'7\thavana\t{kind}\t100\t200\t.\t+\t.\t{attrs}\n')
                            path.write_text(''.join(rows))
                            # Check the actual files handed to bedtools/bedops.
                            with patch.object(reference.os, 'system', side_effect=StopBeforeExternalTools):
                                with self.assertRaises(StopBeforeExternalTools):
                                    reference.generate_reference_df(str(path))
                            self.assertEqual((Path(temp) / 'gene_gtf.bed').read_text(),
                                             '7\t99\t200\tG1\tTRBV19\tprotein_coding\t+\n')
                            self.assertEqual((Path(temp) / 'exon_gtf.bed').read_text(),
                                             '7\t99\t200\tG1\tTRBV19\tT1\n')

    def test_explicit_gene_type_takes_precedence_over_alias(self):
        for attributes in (
            'gene_type "lncRNA"; gene_biotype "protein_coding";',
            'gene_biotype "protein_coding"; gene_type "lncRNA";',
        ):
            self.assertEqual(reference.parse_info(attributes), {'gene_type': 'lncRNA'})


if __name__ == '__main__':
    unittest.main()
