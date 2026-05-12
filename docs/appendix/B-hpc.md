# Appendix B — HPC and SLURM Basics

Almost every calculation in this book is, in practice, run on a
high-performance computing (HPC) cluster. The wall-clock cost of a
serious DFT calculation, a million-atom MD trajectory, or the
training of a foundation MLIP makes a laptop or workstation
impractical for anything beyond the smallest examples. This
appendix is a working primer on the cluster environment most
readers will encounter: a SLURM-managed Linux cluster with login
nodes, compute nodes, modules, and tiered storage.

It is *not* a tutorial on Linux shell usage; we assume the reader
is comfortable with `ls`, `cd`, `cp`, `mv`, `chmod`, `grep`, `find`,
and writing simple shell scripts. The focus here is on the
HPC-specific conventions that catch new users.

---

## B.1 The cluster anatomy

A typical academic cluster contains three classes of machine.

**Login nodes.** A small number ($\sim 1$–$8$) of nodes that users
SSH into. They are *not* for running computations. Their job is to
let you edit files, submit jobs, transfer data, and check status.
Running heavy work on a login node — a slow Python plotting script
might be tolerated; a DFT calculation will not be — will quickly
attract administrative attention. The standard convention is that
*anything that takes more than a few seconds of CPU or appreciable
memory* belongs on a compute node.

**Compute nodes.** The bulk of the cluster. Each is a multi-socket
machine (often $2 \times 32$- or $2 \times 64$-core CPUs), perhaps
with several GPUs (NVIDIA H100, A100, or older V100/P100), with
$256$ GB–$2$ TB of RAM, connected to its neighbours by a low-latency
fabric (InfiniBand, Slingshot). Compute nodes are *reserved through
the scheduler*, not accessed directly. They typically have *no
internet access*, which has implications for code installation.

**Storage nodes.** Behind the scenes; typically a parallel filesystem
(Lustre, BeeGFS, GPFS). The user encounters this filesystem through
several *mount points* with different policies — see Section B.4.

The scheduler itself runs on dedicated management nodes. The user
interacts with it via the SLURM command-line tools.

---

## B.2 SLURM essentials

Most academic clusters run SLURM (Simple Linux Utility for Resource
Management). The handful of commands below cover $95$% of routine
use.

### Submitting jobs

```bash
sbatch myjob.sh
```

Submits a batch job described by the script `myjob.sh`. The script
contains both SLURM directives (`#SBATCH ...`) and the commands to
run. Returns a job ID, e.g. `Submitted batch job 12345678`.

### Monitoring jobs

```bash
squeue -u $USER            # your jobs only
squeue -u $USER -t RUNNING # only those currently running
squeue -j 12345678          # a specific job
```

The output columns name the job, the partition, the state (`PD` =
pending, `R` = running, `CF` = configuring, `CG` = completing), and
the resources allocated. `sinfo` shows partition-level availability.

### Cancelling jobs

```bash
scancel 12345678              # cancel by ID
scancel -u $USER              # cancel all your jobs
scancel -u $USER -t PENDING   # only cancel pending ones
```

### Interactive sessions

For debugging, an *interactive* shell on a compute node is
invaluable.

```bash
salloc -p debug -t 01:00:00 --gres=gpu:1 --mem=64G
```

This requests one GPU and $64$ GB of RAM on the `debug` partition
for one hour, and drops you into a shell on the compute node. The
session ends when the time expires or you `exit`.

### Job history

```bash
sacct -u $USER --starttime=2026-01-01 --format=JobID,JobName,State,Elapsed,MaxRSS
```

Useful when investigating why a job died (`OUT_OF_MEMORY`,
`TIMEOUT`, `FAILED`, `CANCELLED`, `COMPLETED`).

---

## B.3 Sample job scripts

The directives in the header determine resource allocation, runtime
limit and output paths. The rest of the script is plain bash.

### A basic CPU job

```bash
#!/bin/bash
#SBATCH --job-name=qe_relax
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module purge
module load quantum-espresso/7.3 openmpi/4.1

cd "${SLURM_SUBMIT_DIR}"
mpirun -n 32 pw.x < relax.in > relax.out
```

Key points: `%x` expands to the job name, `%j` to the job ID, and
the `logs/` directory must exist before submission. `module purge`
followed by explicit loads is a discipline against contaminating
environments.

### A GPU job

```bash
#!/bin/bash
#SBATCH --job-name=mace_train
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

module purge
module load cuda/12.4 cudnn/9.0
source ~/venv/mace/bin/activate

cd "${SLURM_SUBMIT_DIR}"
python train_mace.py --config config.yaml
```

Note `--gres=gpu:a100:2`: this requests two GPUs of type A100. Many
sites expose multiple GPU types and the type matters both for
performance and for queue length.

### An MPI job spanning multiple nodes

```bash
#!/bin/bash
#SBATCH --job-name=vasp_relax_large
#SBATCH --partition=cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

module purge
module load vasp/6.4 intel-mpi/2021

cd "${SLURM_SUBMIT_DIR}"
mpirun -n ${SLURM_NTASKS} vasp_std
```

`--mem=0` means "all memory on each node". For multi-node MPI work,
`SLURM_NTASKS` is the total task count (here $4\times 64 = 256$),
and is the natural number to pass to `mpirun`.

### An array job

When the same workflow must run on many independent inputs, an
*array job* is the right structure.

```bash
#!/bin/bash
#SBATCH --job-name=screen
#SBATCH --partition=cpu
#SBATCH --array=1-200%10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/screen-%A_%a.out

set -euo pipefail

module purge
module load python/3.11

CASE_ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
WORK="${SLURM_SUBMIT_DIR}/cases/${CASE_ID}"

cd "${WORK}"
python run_one_case.py
```

`--array=1-200%10` requests $200$ tasks but limits concurrent
execution to $10$. The `%A_%a` placeholders in the log filename give
the array job ID and task ID.

---

## B.4 Modules

Most clusters expose installed software through environment modules
(`Lmod` or `environment-modules`). These are not Python modules but
shell-level helpers that prepend directories to `PATH`,
`LD_LIBRARY_PATH`, etc.

```bash
module avail                  # list available modules
module avail python           # list versions matching "python"
module load python/3.11       # load a specific version
module load python            # load the default version (avoid; pin explicitly)
module list                   # show currently loaded modules
module unload python/3.10     # unload
module purge                  # unload everything
```

Two strong recommendations.

**Always `module purge` at the start of a job script.** The login
shell may inherit a different set of modules to those required by
the job. Explicit purging followed by explicit loads makes the
environment reproducible.

**Pin the version.** `module load gcc` loads whatever the
administrators have set as the default, and that default will change
over time. `module load gcc/12.3` is reproducible.

For Python in particular, the right pattern is usually to load only
the Python interpreter from the module system and to manage all
further dependencies through a project-specific virtual environment.

```bash
module load python/3.11
python -m venv ~/venv/proj
source ~/venv/proj/bin/activate
pip install -r requirements.txt
```

The `venv` lives on your home filesystem and is not affected when
the cluster software stack is updated.

---

## B.5 Storage tiers

Cluster filesystems are organised into tiers with very different
policies. The exact names vary by site but the structure is nearly
universal.

| Mount | Typical name | Size | Speed | Backed up? | Purged? |
|---|---|---|---|---|---|
| Home | `/home/$USER` | $\sim 50$ GB | Slow | Yes | No |
| Scratch | `/scratch/$USER` | $\sim 10$ TB | Fast | No | Yes (days–weeks) |
| Project | `/project/$GROUP` | $\sim 1$ TB | Moderate | Yes | No |
| Archive | `/archive/$GROUP` | Effectively unlimited | Slow | Yes | No |

The discipline:

- **Home** is for source code, scripts, dotfiles, virtual
  environments. Quota is tight.
- **Scratch** is for the *working* files of running jobs: input
  geometries, trajectory files, intermediate output. Fast but
  unbacked. Files older than a few weeks are typically purged
  automatically.
- **Project** is for shared data that the whole group needs to keep:
  curated training sets, published results.
- **Archive** is for long-term storage of large datasets you might
  want again in five years but not next week.

Two common mistakes:

1. *Running jobs from `/home`*. The home filesystem cannot sustain
   the I/O of a parallel job. The result is a wallclock-killing
   bottleneck and an unhappy administrator. Always set the working
   directory of a job to scratch.

2. *Leaving large data on scratch indefinitely*. The purge will
   delete it eventually, and you will discover this when the
   directory is gone, not before. Move outputs to project or archive
   as soon as the job completes.

---

## B.6 File transfer

Most data exchange with a cluster happens through one of three
protocols.

### `rsync`

The workhorse for incremental, resumable transfers. To copy a local
directory to the cluster:

```bash
rsync -avzP ./local_data/ user@cluster:/scratch/user/data/
```

Flags: `-a` archive (preserves permissions, times, symlinks), `-v`
verbose, `-z` compress in transit, `-P` show progress and allow
resumption. Trailing slash on the source directory matters: with the
slash, the *contents* of `local_data/` are copied; without, the
directory itself is created on the remote.

### `scp` and `sftp`

`scp` is a simple one-shot copy:

```bash
scp myfile.tar.gz user@cluster:/scratch/user/
```

`sftp` is an interactive FTP-like session. Both are fine for
small or one-off transfers but lack rsync's incrementality.

### Globus

For very large datasets ($> 100$ GB), or transfers between two
clusters, Globus is the standard solution. It is a managed service
that handles checksums, retries, and parallelism. The user
interaction is through a web interface; the bulk transfer runs in
the background. Most institutions have a Globus endpoint configured;
ask your local support team.

---

## B.7 Job-script common mistakes

A non-exhaustive but representative list, all of which the editors
of this book have personally committed.

**Requesting too many CPUs.** A naive `--ntasks-per-node=128` on a
$64$-core node will queue forever (the request is unsatisfiable) or
be silently re-interpreted in a way you do not want. Match the
request to the actual hardware (`scontrol show node <nodename>` will
tell you the cores).

**Requesting too much memory.** `--mem=1T` on a node with
$256$ GB of RAM will queue forever. Always check the partition's
memory cap with `sinfo --Format=NodeList,Memory`.

**Forgetting `module purge`.** Inherited modules from the login
environment may cause obscure linking errors in the job. Purge and
re-load.

**Hard-coded paths to your home directory.** Job scripts checked
into a shared repository should use environment variables
(`${SLURM_SUBMIT_DIR}`, `${SCRATCH}`, `${HOME}`) rather than
absolute paths.

**Missing output directories.** SLURM does not create the
directories named by `--output`/`--error`; the job fails silently if
they do not exist. Either `mkdir -p logs` before submission or
ensure the directory is present in your repository.

**Hidden internet access requirements.** A compute node typically
has no internet. A job that runs `pip install`, `wget`, `huggingface
download` will fail. Install dependencies and download data sets on
the login node (or via a dedicated data-transfer node), not in the
job itself.

**Not cleaning scratch.** When a job series finishes, archive the
useful files and delete the rest. The cluster filesystem is a shared
resource; leaving terabytes of stale trajectories is anti-social.

**Underestimating walltime.** `--time=02:00:00` looked like enough,
the job needs $2:30:00$, the run terminates at the limit with a
partially written checkpoint. Estimate generously and use SLURM's
job-step preemption hooks (`SIGTERM` $5$ minutes before kill) to
trigger graceful saves.

---

## B.8 A reproducibility checklist

For each computational study in this book the editors recommend the
following minimal record of the cluster environment used.

1. The exact `sbatch` command and the script body, version-controlled
   alongside the input data.
2. The output of `module list` at the start of the job, captured to
   the log.
3. The output of `pip freeze` (for Python work) or the equivalent
   compiler / library versions (for compiled codes), captured to the
   log.
4. The hostname and partition (`hostname`, `${SLURM_JOB_PARTITION}`,
   `${SLURM_JOB_NODELIST}`) — useful when investigating
   hardware-dependent anomalies.
5. The exact input files, including any random seeds, in a state that
   can be re-run without modification.

A typical preamble that implements all five:

```bash
echo "host: $(hostname)"
echo "partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "nodes: ${SLURM_JOB_NODELIST:-unknown}"
echo "---- modules ----"
module list 2>&1
echo "---- python env ----"
pip freeze 2>&1 | head -100
echo "---- end preamble ----"
```

This is two pages of log, recoverable forever, that turns
"my old calculation does not reproduce" from a research crisis into
a routine diff.

---

The next appendix consolidates the terminology used throughout the
book.
