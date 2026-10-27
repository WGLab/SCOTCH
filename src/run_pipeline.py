"""Config-based local/SLURM scheduling of the existing SCOTCH CLI (stdlib only)."""
import argparse
from contextlib import suppress
import datetime
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
STAGES = ('annotation', 'compatible', 'summary', 'count')
ARRAYS = ('BAM_PATHS', 'TARGET_DIRS')
BOOLEANS = ('BULK', 'UPDATE_GTF', 'SAVE_MEM_ANN', 'SAVE_MEM', 'GROUP_NOVEL', 'SAVE_CSV')
INTEGERS = ('CPUS', 'COMPATIBLE_JOBS', 'MAX_CONCURRENT_JOBS', 'ANNOTATION_CPUS',
            'COMPATIBLE_CPUS', 'SUMMARY_CPUS', 'COUNT_CPUS',
            'NOVEL_DISCOVERY_MAX_READS', 'NOVEL_DISCOVERY_MAX_ROUNDS', 'NOVEL_READ_N')
SCALARS = ('REFERENCE_GTF', 'REFERENCE_FASTA', 'REFERENCE_PKL', 'PLATFORM',
           'BARCODE_CELL', 'BARCODE_UMI', 'BUILD', 'RUNNER', 'PYTHON',
           'PARTITION', 'ACCOUNT', 'QOS', 'NOVEL_READ_PCT') + BOOLEANS + INTEGERS + tuple(
               f'{stage.upper()}_{resource}' for stage in STAGES for resource in ('MEM', 'TIME'))


def load_config(path):
    """Source trusted Bash config and read explicit fields as NUL-delimited data."""
    path = Path(path).expanduser().resolve(strict=True)
    script = 'set -eo pipefail\nsource "$1" >&2\nsource "$2" >&2\n'
    for key in SCALARS:
        script += f'printf "%s\\0" "${{{key}}}"\n'
    for key in ARRAYS:
        script += f'printf "%s\\0" "${{#{key}[@]}}" "${{{key}[@]}}"\n'
    result = subprocess.run(['bash', '-c', script, 'scotch-config',
                             str(ROOT / 'config_template.sh'), str(path)],
                            check=True, stdout=subprocess.PIPE)
    values = iter(result.stdout.decode().split('\0')[:-1])
    config = {key: next(values) for key in SCALARS}
    for key in ARRAYS:
        config[key] = [next(values) for _ in range(int(next(values)))]
    for key in ('PLATFORM', 'BULK', 'RUNNER'):
        if not config[key].strip():
            raise ValueError(f'{key} must be explicitly set in your config')
    for key in BOOLEANS:
        value = config[key].lower()
        if value not in ('true', 'false'):
            raise ValueError(f'{key} must be true or false')
        config[key] = value == 'true'
    for key in INTEGERS:
        config[key] = int(config[key])
        if config[key] < 0:
            raise ValueError(f'{key} must be nonnegative')
    for key in ('CPUS', 'COMPATIBLE_JOBS', 'MAX_CONCURRENT_JOBS', 'COMPATIBLE_CPUS'):
        if config[key] < 1:
            raise ValueError(f'{key} must be at least 1')
    config['NOVEL_READ_PCT'] = float(config['NOVEL_READ_PCT'])
    if not 0 <= config['NOVEL_READ_PCT'] <= 1:
        raise ValueError('NOVEL_READ_PCT must be between 0 and 1')
    if config['RUNNER'] not in ('local', 'slurm'):
        raise ValueError('RUNNER must be local or slurm')
    if config['PLATFORM'] not in ('10x-ont', '10x-pacbio', 'parse-ont'):
        raise ValueError('Unsupported PLATFORM')
    if config['BULK'] and config['PLATFORM'] == 'parse-ont':
        raise ValueError('Parse mode requires single-cell data (BULK=false)')
    if not config['BULK'] and config['PLATFORM'] != 'parse-ont':
        for key in ('BARCODE_CELL', 'BARCODE_UMI'):
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9]', config[key]):
                raise ValueError(f'{key} must explicitly name the two-character BAM tag for single-cell data')
    def resolve(value):
        p = Path(value).expanduser()
        return str((path.parent / p).resolve())
    for key in ARRAYS:
        if any(not value.strip() for value in config[key]):
            raise ValueError(f'{key} contains an empty path')
        config[key] = [resolve(value) for value in config[key]]
    for key in ('REFERENCE_GTF', 'REFERENCE_FASTA', 'REFERENCE_PKL'):
        if config[key]:
            config[key] = resolve(config[key])
    python = config['PYTHON'] or sys.executable
    if '/' in python:
        python = resolve(python)
    config['PYTHON'] = shutil.which(python)
    if not config['PYTHON']:
        raise ValueError(f'Python executable not found: {python}')
    if not config['BAM_PATHS'] or len(config['BAM_PATHS']) != len(config['TARGET_DIRS']):
        raise ValueError('BAM_PATHS and TARGET_DIRS must be nonempty and have matching lengths')
    if config['PLATFORM'] == 'parse-ont' and len(config['BAM_PATHS']) != 1:
        raise ValueError('Parse mode requires exactly one BAM input and one target')
    targets = [Path(p) for p in config['TARGET_DIRS']]
    for i, target in enumerate(targets):
        if target.exists() and not target.is_dir():
            raise ValueError(f'Output target is not a directory: {target}')
        if any(target == other or target in other.parents or other in target.parents
               for other in targets[:i]):
            raise ValueError('Output target directories must be distinct and non-nested')
    for stage in STAGES:
        cpu_key = f'{stage.upper()}_CPUS'
        config[cpu_key] = config[cpu_key] or config['CPUS']
        if config['RUNNER'] == 'local' and config[cpu_key] > config['CPUS']:
            raise ValueError(f'{cpu_key} exceeds the local CPUS budget')
        if not re.fullmatch(r'\d+(?:[KMGTP])?', config[f'{stage.upper()}_MEM'], re.I):
            raise ValueError(f'{stage.upper()}_MEM must be a SLURM memory value such as 32G')
        if not re.fullmatch(r'(?:\d+-)?\d+:\d{2}:\d{2}', config[f'{stage.upper()}_TIME']):
            raise ValueError(f'{stage.upper()}_TIME must be HH:MM:SS or D-HH:MM:SS')
    config['CONFIG_PATH'] = str(path)
    return config


def validate_inputs(config):
    for key in ('REFERENCE_GTF', 'REFERENCE_FASTA'):
        if not config[key] or not Path(config[key]).is_file():
            raise ValueError(f'{key} is not an existing file: {config[key]}')
    if not Path(config['REFERENCE_FASTA'] + '.fai').is_file():
        raise ValueError('Missing genome FASTA index (.fai); run samtools faidx first')
    if config['REFERENCE_PKL'] and not Path(config['REFERENCE_PKL']).is_file():
        raise ValueError('REFERENCE_PKL does not exist')
    for entry in config['BAM_PATHS']:
        path = Path(entry)
        bams = sorted(path.glob('*.bam')) if path.is_dir() else [path]
        if not bams:
            raise ValueError(f'No BAM files in {path}')
        for bam in bams:
            if not bam.is_file():
                raise ValueError(f'BAM does not exist: {bam}')
            if not any(p.is_file() for p in (Path(str(bam) + '.bai'), bam.with_suffix('.bai'),
                                             Path(str(bam) + '.csi'), bam.with_suffix('.csi'))):
                raise ValueError(f'Missing BAM index: {bam}')


def stage_command(config, stage, index=0):
    task = {'compatible': 'compatible matrix', 'count': 'count matrix'}.get(stage, stage)
    command = [config['PYTHON'], str(ROOT / 'src/main_preprocessing.py'), '--task', task,
               '--platform', config['PLATFORM'], '--target', *config['TARGET_DIRS'],
               '--workers', str(config[f'{stage.upper()}_CPUS'])]
    if config['BULK']:
        command.append('--bulk')
    if stage in ('annotation', 'compatible'):
        command += ['--bam', *config['BAM_PATHS'], '--reference', config['REFERENCE_GTF'],
                    '--reference_pkl', config['REFERENCE_PKL'] or 'None',
                    '--reference_genome_fasta', config['REFERENCE_FASTA']]
        for key, option in (('BARCODE_CELL', '--barcode_cell'), ('BARCODE_UMI', '--barcode_umi'),
                            ('BUILD', '--build')):
            if key.startswith('BARCODE_') and (config['BULK'] or config['PLATFORM'] == 'parse-ont'):
                continue
            if config[key]:
                command += [option, config[key]]
    if stage == 'annotation':
        command.append('--update_gtf' if config['UPDATE_GTF'] else '--update_gtf_off')
        if config['SAVE_MEM_ANN']:
            command.append('--save_mem_ann')
    if stage == 'compatible':
        command += ['--total_jobs', str(config['COMPATIBLE_JOBS']), '--job_index', str(index),
                    '--novel_discovery_max_reads', str(config['NOVEL_DISCOVERY_MAX_READS']),
                    '--novel_discovery_max_rounds', str(config['NOVEL_DISCOVERY_MAX_ROUNDS']),
                    '--save_mem' if config['SAVE_MEM'] else '--save_mem_off']
    if stage == 'count':
        command += ['--novel_read_n', str(config['NOVEL_READ_N']),
                    '--novel_read_pct', str(config['NOVEL_READ_PCT']), '--save_mtx',
                    '--group_novel' if config['GROUP_NOVEL'] else '--group_novel_off']
        if config['SAVE_CSV']:
            command.append('--save_csv')
    return command


def check_outputs(config, stage, index=0):
    for sample_index, target in enumerate(map(Path, config['TARGET_DIRS'])):
        reference = target / 'reference'
        if stage == 'annotation':
            # The annotation step copies only the shared metagene pickle to later targets.
            required = [reference / 'metageneStructureInformation.pkl']
            if sample_index == 0:
                required.append(reference / 'geneStructureInformation.pkl')
        elif stage == 'compatible':
            required = [reference / f'metageneStructureInformationwNovel_{index}.{ext}' for ext in ('pkl', 'gtf')]
        elif stage == 'summary':
            required = [reference / 'metageneStructureInformationwNovel.pkl',
                        reference / 'SCOTCH_updated_annotation.gtf']
            samples = [target]
            if config['PLATFORM'] == 'parse-ont':
                samples_dir = target / 'samples'
                samples = sorted(p for p in samples_dir.iterdir() if p.is_dir()) if samples_dir.is_dir() else []
                if not samples:
                    raise RuntimeError(f'No Parse sample directories were generated under {target}')
            for sample in samples:
                required += [sample / 'auxiliary' / name for name in
                             ('all_read_isoform_exon_mapping.tsv', 'read_selection.pkl')]
        else:
            return
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise RuntimeError(f'{stage} did not produce expected outputs: ' + ', '.join(missing))


def stage_environment():
    env = os.environ.copy()
    # Parallelism is controlled by stage workers and partition concurrency.
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[name] = '1'
    env['PYTHONUNBUFFERED'] = '1'
    return env


def run_stage(config, stage, index=0):
    command = stage_command(config, stage, index)
    print(shlex.join(command), flush=True)
    subprocess.run(command, check=True, env=stage_environment())
    check_outputs(config, stage, index)


def run_local(config, run_dir):
    active = []
    concurrency = min(config['MAX_CONCURRENT_JOBS'], config['CPUS'] // config['COMPATIBLE_CPUS'])
    def launch(stage, index=0):
        name = f'{stage}_{index}' if stage == 'compatible' else stage
        log = (run_dir / f'{name}.log').open('w')
        command = [config['PYTHON'], str(Path(__file__).resolve()), '--manifest',
                   str(run_dir / 'config.json'), '--stage', stage, '--index', str(index)]
        print(f'Running {name}; log: {log.name}', flush=True)
        try:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
        except BaseException:
            log.close()
            raise
        active.append((process, log, name))
    def reap():
        for process, log, name in active[:]:
            status = process.poll()
            if status is not None:
                log.close()
                active.remove((process, log, name))
                if status:
                    raise RuntimeError(f'{name} failed (exit {status}); see {log.name}')
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    previous_sigterm = signal.signal(signal.SIGTERM, interrupted)
    try:
        for stage in STAGES:
            jobs = range(config['COMPATIBLE_JOBS']) if stage == 'compatible' else range(1)
            limit = concurrency if stage == 'compatible' else 1
            for index in jobs:
                while len(active) >= limit:
                    reap()
                    time.sleep(0.1)
                reap()
                launch(stage, index)
            while active:
                reap()
                time.sleep(0.1)
    finally:
        for process, _, _ in active:
            if process.poll() is None:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGTERM)
        for process, log, _ in active:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            log.close()
        signal.signal(signal.SIGTERM, previous_sigterm)


def slurm_script(config, run_dir, stage):
    command = [config['PYTHON'], str(Path(__file__).resolve()), '--manifest',
               str(run_dir / 'config.json'), '--stage', stage]
    text = '#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.join(command)
    if stage == 'compatible':
        text += ' --index "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"'
    return text + '\n'


def sbatch_command(config, run_dir, stage, dependency=None):
    prefix = stage.upper()
    command = ['sbatch', '--parsable', '--job-name', f'scotch_{stage}',
               '--ntasks', '1',
               '--cpus-per-task', str(config[f'{prefix}_CPUS']),
               '--mem', config[f'{prefix}_MEM'], '--time', config[f'{prefix}_TIME'],
               '--output', str(run_dir / f'{stage}_%A_%a.log')]
    for key in ('PARTITION', 'ACCOUNT', 'QOS'):
        if config[key]:
            command += ['--' + key.lower(), config[key]]
    if stage == 'compatible':
        command += ['--array', f'0-{config["COMPATIBLE_JOBS"] - 1}%{config["MAX_CONCURRENT_JOBS"]}']
    if dependency:
        command += ['--dependency', 'afterok:' + dependency]
    command.append(str(run_dir / f'{stage}.sh'))
    return command


def submit_slurm(config, run_dir):
    if not shutil.which('sbatch'):
        raise ValueError('sbatch not found; run on your SLURM cluster or use RUNNER=local')
    submitted = {}
    dependency = None
    try:
        for stage in STAGES:
            (run_dir / f'{stage}.sh').write_text(slurm_script(config, run_dir, stage))
            result = subprocess.run(sbatch_command(config, run_dir, stage, dependency),
                                    check=True, text=True, stdout=subprocess.PIPE)
            job_id = result.stdout.strip().split(';')[0]
            if not job_id.isdigit():
                raise RuntimeError(f'Unexpected sbatch response: {result.stdout!r}')
            submitted[stage] = job_id
            (run_dir / 'jobs.json').write_text(json.dumps(submitted, indent=2) + '\n')
            print(f'Submitted {stage}: {job_id}', flush=True)
            dependency = job_id
    except BaseException:
        if submitted:
            print('Submission stopped. Previously submitted jobs: ' + ', '.join(submitted.values())
                  + '. Inspect with squeue; cancel with scancel if needed.', file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', help='Bash config copied from config_template.sh')
    parser.add_argument('--dry-run', action='store_true', help='Validate config and print commands without running/submitting')
    parser.add_argument('--manifest', help=argparse.SUPPRESS)
    parser.add_argument('--stage', choices=STAGES, help=argparse.SUPPRESS)
    parser.add_argument('--index', type=int, default=0, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.manifest:
        if not args.stage:
            parser.error('--manifest requires --stage')
        config = json.loads(Path(args.manifest).read_text())
        if args.stage == 'compatible' and not 0 <= args.index < config['COMPATIBLE_JOBS']:
            parser.error('Compatible partition index is out of range')
        run_stage(config, args.stage, args.index)
        return
    if not args.config:
        parser.error('Provide a config file, e.g. bash scotch.sh my_run.sh')
    config = load_config(args.config)
    validate_inputs(config)
    if args.dry_run:
        for stage in STAGES:
            indices = range(config['COMPATIBLE_JOBS']) if stage == 'compatible' else range(1)
            for index in indices:
                print(shlex.join(stage_command(config, stage, index)))
        if config['RUNNER'] == 'slurm':
            dependency = None
            for stage in STAGES:
                print(shlex.join(sbatch_command(config, Path('<run-directory>'), stage, dependency)))
                dependency = '<' + stage + '-job-id>'
        return
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run_dir = Path(config['TARGET_DIRS'][0]) / '.scotch_runs' / timestamp
    run_dir.mkdir(parents=True)
    (run_dir / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    print(f'Run settings and logs: {run_dir}', flush=True)
    if config['RUNNER'] == 'slurm':
        submit_slurm(config, run_dir)
    else:
        run_local(config, run_dir)
        print('SCOTCH preprocessing completed.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f'SCOTCH: {exc}', file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
