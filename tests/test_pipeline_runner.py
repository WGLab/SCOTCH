import ast
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

SRC = Path(__file__).resolve().parents[1] / 'src'
sys.path.insert(0, str(SRC))
import run_pipeline as runner


class PipelineRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='scotch runner ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ('reads.bam', 'reads.bam.bai', 'genes.gtf', 'genome.fa', 'genome.fa.fai'):
            (self.root / name).touch()
        self.config_path = self.root / 'my_run.sh'
        self.config_path.write_text('''BAM_PATHS=("reads.bam")
TARGET_DIRS=("output")
REFERENCE_GTF="genes.gtf"
REFERENCE_FASTA="genome.fa"
PLATFORM="10x-ont"
BULK=false
RUNNER="local"
BARCODE_CELL="CB"
BARCODE_UMI="UB"
CPUS=2
COMPATIBLE_JOBS=3
MAX_CONCURRENT_JOBS=2
''')

    def config(self, extra=''):
        with self.config_path.open('a') as handle:
            handle.write(extra)
        return runner.load_config(self.config_path)

    def test_defaults_relative_paths_and_cpu_budget(self):
        config = self.config()
        runner.validate_inputs(config)
        self.assertEqual(config['BAM_PATHS'], [str((self.root / 'reads.bam').resolve())])
        self.assertEqual(config['ANNOTATION_CPUS'], 2)
        self.assertEqual(config['SUMMARY_CPUS'], 2)
        self.assertEqual(config['COMPATIBLE_CPUS'], 1)
        self.assertTrue(config['UPDATE_GTF'])
        self.assertFalse(config['SAVE_CSV'])
        self.assertEqual(config['REFERENCE_PKL'], '')

    def test_invalid_boolean_and_cpu_overcommit_rejected(self):
        with self.assertRaisesRegex(ValueError, 'BULK'):
            self.config('BULK=maybe\n')
        with self.assertRaisesRegex(ValueError, 'budget'):
            self.config('BULK=false\nCOMPATIBLE_CPUS=3\n')

    def test_empty_output_path_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'empty path'):
            self.config('TARGET_DIRS=("")\n')

    def test_data_and_runner_choices_cannot_be_left_blank(self):
        for key in ('BULK', 'PLATFORM', 'RUNNER'):
            with self.subTest(key=key):
                original = self.config_path.read_text()
                with self.assertRaisesRegex(ValueError, key + ' must be explicitly set'):
                    self.config(f'{key}=""\n')
                self.config_path.write_text(original)

    def test_single_cell_tags_required_bulk_and_parse_can_leave_empty(self):
        with self.assertRaisesRegex(ValueError, 'BARCODE_CELL'):
            self.config('BARCODE_CELL=""\nBARCODE_UMI=""\n')
        for extra in ('BULK=true\n', 'BULK=false\nPLATFORM=parse-ont\n'):
            config = self.config(extra)
            command = runner.stage_command(config, 'annotation')
            self.assertNotIn('--barcode_cell', command)
            self.assertNotIn('--barcode_umi', command)

    def test_multiple_bams_keep_input_output_order(self):
        config = self.config('BAM_PATHS=("sample_a.bam" "sample_b.bam")\n'
                             'TARGET_DIRS=("output_a" "output_b")\n')
        command = runner.stage_command(config, 'compatible')
        i = command.index('--bam')
        self.assertEqual(command[i+1:i+3], config['BAM_PATHS'])
        i = command.index('--target')
        self.assertEqual(command[i+1:i+3], config['TARGET_DIRS'])
        for target in map(Path, config['TARGET_DIRS']):
            (target / 'reference').mkdir(parents=True)
            (target / 'reference/metageneStructureInformation.pkl').touch()
        (Path(config['TARGET_DIRS'][0]) / 'reference/geneStructureInformation.pkl').touch()
        runner.check_outputs(config, 'annotation')

    def test_missing_index_and_duplicate_targets_rejected(self):
        config = self.config()
        (self.root / 'genome.fa.fai').unlink()
        with self.assertRaisesRegex(ValueError, 'FASTA index'):
            runner.validate_inputs(config)
        with self.assertRaisesRegex(ValueError, 'distinct'):
            self.config('BAM_PATHS=("reads.bam" "reads.bam")\nTARGET_DIRS=("output" "output")\n')

    def test_commands_use_real_cli_options(self):
        import argparse
        parser = argparse.ArgumentParser()
        for node in ast.parse((SRC / 'main_preprocessing.py').read_text()).body:
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                f = node.value.func
                if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) and f.value.id == 'parser' and f.attr == 'add_argument':
                    exec(compile(ast.Module(body=[node], type_ignores=[]), '<cli>', 'exec'), {'parser': parser})
        config = self.config('BULK=true\nSAVE_CSV=true\nGROUP_NOVEL=false\nUPDATE_GTF=false\n')
        for stage in runner.STAGES:
            parsed = parser.parse_args(runner.stage_command(config, stage, 2)[2:])
            self.assertTrue(parsed.bulk)
            if stage == 'compatible':
                self.assertEqual((parsed.total_jobs, parsed.job_index), (3, 2))
            if stage == 'count':
                self.assertTrue(parsed.save_csv)
                self.assertFalse(parsed.group_novel)

    def test_dry_run_does_not_create_outputs(self):
        result = subprocess.run([sys.executable, str(SRC / 'run_pipeline.py'),
                                 str(self.config_path), '--dry-run'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--job_index 2', result.stdout)
        self.assertFalse((self.root / 'output').exists())

    def prepare_fake_pipeline(self, fail_index=-1, missing_annotation=False):
        # Execute the real orchestrator; substitute only the expensive SCOTCH CLI.
        shim = self.root / 'fake_python'
        shim.write_text(f'#!{sys.executable}\n' + '''import os, sys, time
from pathlib import Path
if Path(sys.argv[1]).name == 'run_pipeline.py':
    os.execv(sys.executable, [sys.executable] + sys.argv[1:])
args = sys.argv[2:]
def value(key):
    return args[args.index(key) + 1]
target = Path(value('--target'))
target.mkdir(parents=True, exist_ok=True)
stage = value('--task')
with (target / 'events').open('a') as f:
    f.write(stage + '\\n')
ref = target / 'reference'
ref.mkdir(exist_ok=True)
def touch(name):
    (ref / name).touch()
if stage == 'annotation':
    touch('geneStructureInformation.pkl')
    if not MISSING_ANNOTATION:
        touch('metageneStructureInformation.pkl')
elif stage == 'compatible matrix':
    index = int(value('--job_index'))
    if index == FAIL_INDEX:
        sys.exit(7)
    time.sleep(0.05)
    touch(f'metageneStructureInformationwNovel_{index}.pkl')
    touch(f'metageneStructureInformationwNovel_{index}.gtf')
elif stage == 'summary':
    assert len(list(ref.glob('metageneStructureInformationwNovel_*.pkl'))) == 3
    touch('metageneStructureInformationwNovel.pkl')
    touch('SCOTCH_updated_annotation.gtf')
    aux = target / 'auxiliary'
    aux.mkdir(exist_ok=True)
    (aux / 'read_selection.pkl').touch()
    (aux / 'all_read_isoform_exon_mapping.tsv').touch()
elif stage == 'count matrix':
    assert (target / 'auxiliary/read_selection.pkl').exists()
    touch('SCOTCH_updated_annotation_filtered_0_0.gtf')
    out = target / 'count_matrix'
    out.mkdir(exist_ok=True)
    for level in ('gene', 'transcript'):
        for ext in ('mtx', 'pickle'):
            (out / f'adata_{level}_0_0.{ext}').write_text('output')
'''.replace('FAIL_INDEX', str(fail_index)).replace('MISSING_ANNOTATION', repr(missing_annotation)))
        shim.chmod(0o755)
        config = self.config('PYTHON=' + shlex.quote(str(shim)) + '\n')
        directory = self.root / 'run logs'
        directory.mkdir()
        (directory / 'config.json').write_text(json.dumps(config))
        return config, directory

    def test_local_pipeline_runs_all_partitions_before_summary(self):
        config, directory = self.prepare_fake_pipeline()
        runner.run_local(config, directory)
        events = (self.root / 'output/events').read_text().splitlines()
        self.assertEqual(events, ['annotation'] + ['compatible matrix'] * 3 + ['summary', 'count matrix'])
        self.assertEqual(len(list(directory.glob('compatible_*.log'))), 3)

    def test_failed_partition_stops_summary_and_count(self):
        config, directory = self.prepare_fake_pipeline(fail_index=1)
        with self.assertRaisesRegex(RuntimeError, 'compatible_1 failed'):
            runner.run_local(config, directory)
        events = (self.root / 'output/events').read_text().splitlines()
        self.assertNotIn('summary', events)
        self.assertNotIn('count matrix', events)

    def test_missing_outputs_fail_even_if_cli_exits_zero(self):
        config, directory = self.prepare_fake_pipeline(missing_annotation=True)
        with self.assertRaisesRegex(RuntimeError, 'annotation failed'):
            runner.run_local(config, directory)
        self.assertEqual((self.root / 'output/events').read_text(), 'annotation\n')

    def cleanup_fixture(self, parse=False):
        config = self.config('KEEP_INTERMEDIATE=false\n')
        if parse:
            config['PLATFORM'] = 'parse-ont'
        target = Path(config['TARGET_DIRS'][0])
        ref = target / 'reference'
        ref.mkdir(parents=True)
        for name in ('metageneStructureInformationwNovel.pkl', 'SCOTCH_updated_annotation.gtf',
                     'SCOTCH_updated_annotation_filtered_0_0.gtf'):
            (ref / name).write_text('reference')
        samples = [target / 'samples/a', target / 'samples/b'] if parse else [target]
        disposable = []
        preserved = []
        for sample in samples:
            for directory in ('auxiliary', 'count_matrix', 'compatible_matrix'):
                (sample / directory).mkdir(parents=True, exist_ok=True)
            for name in ('read_selection.pkl', 'all_read_isoform_exon_mapping.tsv'):
                preserved.append(sample / 'auxiliary' / name)
            for level in ('gene', 'transcript'):
                for ext in ('mtx', 'pickle'):
                    preserved.append(sample / 'count_matrix' / f'adata_{level}_0_0.{ext}')
            disposable.append(sample / 'compatible_matrix/GENE_ENSG1.csv')
            preserved.append(sample / 'compatible_matrix/log.txt')
        (target / 'bam').mkdir()
        (ref / 'files_jobs').mkdir()
        disposable += [target / 'bam/bam.Info.pkl', target / 'bam/bam.Info2.sqlite',
                       ref / 'files_jobs/metageneStructureInformationwNovel_0.gtf']
        preserved += [target / 'bam/reads.bam', ref / 'files_jobs/unknown.pkl']
        for path in disposable + preserved:
            path.write_text('test data')
        run_dir = self.root / 'cleanup run'
        run_dir.mkdir()
        return config, run_dir, disposable, preserved

    def test_cleanup_parse_preserves_final_outputs_and_unknown_files(self):
        config, directory, disposable, preserved = self.cleanup_fixture(parse=True)
        runner.cleanup_intermediates(config, directory)
        self.assertTrue(all(not p.exists() for p in disposable))
        self.assertTrue(all(p.exists() for p in preserved))
        report = json.loads((directory / 'cleanup.json').read_text())
        self.assertEqual(report['status'], 'complete')
        self.assertEqual(report['bytes_removed'], 9 * len(disposable))

    def test_cleanup_requires_every_sample_output_before_deleting(self):
        config, directory, disposable, preserved = self.cleanup_fixture(parse=True)
        next(p for p in preserved if '/samples/b/count_matrix/' in str(p)).write_text('')
        with self.assertRaisesRegex(RuntimeError, 'nonempty output'):
            runner.cleanup_intermediates(config, directory)
        self.assertTrue(all(p.exists() for p in disposable))

    def test_cleanup_default_keep_and_symlink_and_input_protection(self):
        config, directory, disposable, _ = self.cleanup_fixture()
        config['KEEP_INTERMEDIATE'] = True
        runner.cleanup_intermediates(config, directory)
        self.assertTrue(all(p.exists() for p in disposable))
        config['KEEP_INTERMEDIATE'] = False
        config['REFERENCE_PKL'] = str(disposable[1])
        link = disposable[0].parent / 'LINK_ENSG2.csv'
        source = self.root / 'outside.csv'
        source.write_text('keep')
        link.symlink_to(source)
        runner.cleanup_intermediates(config, directory)
        self.assertTrue(disposable[1].exists())
        self.assertTrue(link.is_symlink())
        self.assertEqual(source.read_text(), 'keep')

    def test_local_cleanup_runs_after_count(self):
        config, directory = self.prepare_fake_pipeline()
        config['KEEP_INTERMEDIATE'] = False
        (directory / 'config.json').write_text(json.dumps(config))
        runner.run_local(config, directory)
        self.assertEqual(json.loads((directory / 'cleanup.json').read_text())['status'], 'complete')

    def test_failed_pipeline_never_runs_cleanup(self):
        config, directory = self.prepare_fake_pipeline(fail_index=1)
        config['KEEP_INTERMEDIATE'] = False
        (directory / 'config.json').write_text(json.dumps(config))
        with self.assertRaises(RuntimeError):
            runner.run_local(config, directory)
        self.assertFalse((directory / 'cleanup.json').exists())
        self.assertFalse((directory / 'cleanup.log').exists())

    def test_slurm_cleanup_depends_on_count(self):
        config = self.config('RUNNER=slurm\nKEEP_INTERMEDIATE=false\n')
        directory = self.root / 'slurm'
        directory.mkdir()
        calls = []
        def submit(command, **kwargs):
            calls.append(command)
            return types.SimpleNamespace(stdout=str(100 + len(calls)))
        with patch.object(runner.shutil, 'which', return_value='/usr/bin/sbatch'), \
                patch.object(runner.subprocess, 'run', side_effect=submit):
            runner.submit_slurm(config, directory)
        self.assertEqual(len(calls), 5)
        self.assertIn('afterok:104', calls[-1])
        self.assertIn('scotch_cleanup', calls[-1])

    def test_slurm_array_dependencies_and_frozen_configuration(self):
        config = self.config('RUNNER=slurm\nPARTITION=compute\n')
        directory = self.root / 'slurm logs'
        directory.mkdir()
        calls = []
        def submit(command, **kwargs):
            calls.append(command)
            return types.SimpleNamespace(stdout=f'{100 + len(calls)};cluster\n')
        with patch.object(runner.shutil, 'which', return_value='/usr/bin/sbatch'), \
                patch.object(runner.subprocess, 'run', side_effect=submit):
            runner.submit_slurm(config, directory)
        self.assertNotIn('--dependency', calls[0])
        self.assertIn('0-2%2', calls[1])
        for i in range(1, 4):
            self.assertIn(f'afterok:{100+i}', calls[i])
        self.assertEqual(json.loads((directory / 'jobs.json').read_text())['count'], '104')
        script = (directory / 'compatible.sh').read_text()
        self.assertIn('SLURM_ARRAY_TASK_ID', script)
        self.assertIn('config.json', script)
        for path in directory.glob('*.sh'):
            subprocess.run(['bash', '-n', str(path)], check=True)

    def test_summary_exceptions_are_not_reported_as_success(self):
        tree = ast.parse((SRC / 'main_preprocessing.py').read_text())
        node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'run_summary')
        from unittest.mock import Mock
        cp = Mock()
        cp.summarise_auxiliary.side_effect = RuntimeError('test error')
        namespace = dict(args=types.SimpleNamespace(target=['out'], bulk=False, platform='parse-ont', workers=2),
                         cp=cp, setup_logger=Mock(return_value=(Mock(), 'summary.log')), copy_log_to_targets=Mock())
        exec(compile(ast.Module(body=[node], type_ignores=[]), '<summary>', 'exec'), namespace)
        with self.assertRaisesRegex(RuntimeError, 'Summary failed'):
            namespace['run_summary']()
        self.assertEqual(cp.summarise_auxiliary.call_args.kwargs['workers'], 2)

    def test_copy_logs_preserves_stage_and_partition_filename(self):
        source = self.root / 'compatible_2.log'
        source.write_text('partition 2')
        target = self.root / 'second_target'
        target.mkdir()
        (target / 'annotation.log').write_text('annotation')
        tree = ast.parse((SRC / 'main_preprocessing.py').read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'copy_log_to_targets')
        import shutil
        namespace = dict(os=os, shutil=shutil)
        exec(compile(ast.Module(body=[node], type_ignores=[]), '<logging>', 'exec'), namespace)
        namespace['copy_log_to_targets'](str(source), [str(self.root), str(target)])
        self.assertEqual((target / 'compatible_2.log').read_text(), 'partition 2')
        self.assertEqual((target / 'annotation.log').read_text(), 'annotation')


if __name__ == '__main__':
    unittest.main()
