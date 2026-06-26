# Aspera CLI Setup and Configuration Guide

This document describes how to configure the IBM Aspera CLI (`ascli`) within the `mdd` Mamba environment to download FASTQ cohort datasets from the European Nucleotide Archive (ENA) at high speed.

---

## 📌 Context & Dependency Workaround

The `aspera-cli` package from Bioconda requires Ruby >= 3.x. However, Ruby 3.x removed the standard `webrick` library from its core, resulting in a crash (`LoadError: cannot load such file -- webrick`) when calling `ascli`.

To preserve the environment integrity and portability without modifying the global system packages:
1. We install `webrick` locally into the Mamba environment using its local `gem` manager.
2. We configure the `GEM_PATH` environment variable so the Ruby interpreter knows to search the environment's `share/rubygems` folder.

---

## 🛠️ Step-by-Step Installation

### Step 1: Install `aspera-cli` (if not already installed)
Ensure your environment is active and install the package:
```bash
mamba install -y -c bioconda aspera-cli
```

### Step 2: Install the Missing Ruby Dependencies
Since Ruby version 3.4.0+ removed `webrick` and `bigdecimal` from its default packages, install them inside the environment's Gem library:
```bash
gem install webrick bigdecimal
```


### Step 3: Export the `GEM_PATH`
Since Mamba's Ruby installation uses a non-standard gem home directory, export `GEM_PATH` in your terminal session to make the installed gems loadable:
```bash
export GEM_PATH="$CONDA_PREFIX/share/rubygems"
```

> [!TIP]
> You can add `export GEM_PATH="$CONDA_PREFIX/share/rubygems"` directly into your SLURM download batch scripts or `.bashrc` so you do not have to run it manually every time.

---

## ⚙️ Aspera Configuration & ENA Integration

### Step 4: Install the Standalone `ascp` Binary
`ascli` is a high-level manager that requires the compiled `ascp` transfer engine to do the actual data movement. Install the SDK binary into the local folder:
```bash
ascli config ascp install
```
This downloads and extracts the compiled `ascp` binary to `~/.aspera/sdk/ascp`.

### Step 5: Configure the ENA Server Preset
We configure a preset named `era` so we don't have to specify connection details or keys for each file. 

*Note: In `aspera-cli v4.20.0` (installed via Bioconda), the module namespace is `Aspera::Ascp::Installation` rather than `Fasp::Installation`, and we retrieve the bypass key path via `.path(:ssh_private_rsa)`.*

We set the transfer rate limit to **400 Mbps (50 MB/s)** using the `--ts` (transfer specification) flag. This is essential in HPC environments to prevent Aspera from aggressively saturating the network and triggering firewall/IDS UDP packet drop policies.

Run the following configuration command:
```bash
ascli conf preset update era \
  --url=ssh://fasp.sra.ebi.ac.uk:33001 \
  --username=era-fasp \
  --ssh-keys=@ruby:Aspera::Ascp::Installation.instance.path\(:ssh_private_rsa\) \
  --ts=@json:'{"target_rate_kbps":400000}'
```
---

## 🚀 Downloading Files

Once configured, files can be downloaded from ENA at high speeds using the `-Pera` preset option.

### Interactive Terminal (standard usage)

In an interactive terminal, `ascli` detects a TTY and renders the full progress bar (%, speed, ETA) automatically. No wrapper needed:

```bash
ascli -Pera server download <remote_path> --to-folder=<local_destination>
```

> [!IMPORTANT]
> The destination folder argument MUST use the equals sign (`--to-folder=/path/to/dir`). Using a space (e.g. `--to-folder /path/to/dir`) will cause `ascli` to raise `OptionParser::MissingArgument` and fail.


### Slurm / Log File Redirection & Failed Progress Logging Attempts

When stdout/stderr are redirected to a log file (as in Slurm batch jobs), the underlying Ruby progress bar library (`ruby-progressbar`) detects that there is no active TTY (interactive terminal) and silently disables all progress/speed/ETA output. 

Several methods were attempted to force pseudo-TTY allocation to get informative progress logs, but each was ultimately rejected due to operational issues:

1. **`script` Command Wrapper**
   * **Attempt**: Wrapping the download command in the standard Unix `script` tool (e.g. `script -q -c "ascli..."`).
   * **Failure Reason**: `script` blocks indefinitely between files in a loop. It waits for an EOF on the PTY master, which only occurs once all processes holding the PTY slave file descriptor close it. If `ascp` (or any subprocess) briefly keeps a file descriptor reference open after `ascli` exits, `script` hangs, preventing the loop from advancing to the next file.

2. **`python3 pty.spawn()` standard library tool**
   * **Attempt**: Using inline Python to invoke `pty.spawn(['ascli', ...])`.
   * **Failure Reason**: Re-spawning Python processes on every loop iteration adds overhead. Furthermore, the command formatting is highly verbose, sensitive to shell parameter expansion, and prone to argument parsing mismatch with `ascli` option structures (such as the space vs equals sign for `--to-folder`).

3. **`expect` Package (`unbuffer` utility)**
   * **Attempt**: Wrapping the download command in `unbuffer ascli ...`.
   * **Failure Reason**: While this successfully bypassed the TTY detection and rendered live progress bars in the log files, it introduced severe instability. When a network/UDP data transfer timeout occurred (`Session data transfer timeout`), the pseudo-TTY wrapper prevented `ascp` and `ascli` from exiting cleanly. The processes hung indefinitely, causing the Slurm loop to freeze and keep the job active without making any progress.

**Final Decision: Direct Standard Execution**

For maximum reliability, predictability, and clean error handling in production Slurm environments, the standard `ascli` command must be run directly without any TTY wrappers. While this omits the progress-bar meter, it guarantees that timeouts and failures are raised immediately, allowing the script to proceed or exit correctly:

```bash
ascli -Pera server download <remote_path> --to-folder=<local_destination>
```

### `scancel` / Process Termination

`ascp` (the underlying Aspera binary at `~/.aspera/sdk/ascp`) **ignores SIGTERM** by default. This means `scancel` alone will not kill a running transfer — the `ascp` process lingers as an orphan.

The fix is a `SIGTERM` trap in the batch script that explicitly kills all child processes:

```bash
cleanup() {
    pkill -P $$ 2>/dev/null || true
    exit 1
}
trap cleanup SIGTERM SIGINT
```

This is already implemented in [../_hpc/download.sbatch](../_hpc/download.sbatch).
