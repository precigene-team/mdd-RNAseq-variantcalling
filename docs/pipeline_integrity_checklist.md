<!-- [[docs/pipeline_integrity_checklist.md]] -->
<!-- Goal: Checklist to guarantee pipeline data integrity and prevent batch/version inconsistencies during hybrid local-remote processing. -->

# Pipeline Data Integrity Checklist (Hybrid Local-HPC Processing)

This checklist provides a set of validation steps to ensure data consistency, software compatibility, and pipeline integrity when processing samples concurrently on a local workstation and the Slurm HPC cluster.

---

## 💻 1. Environment & Software Alignment (Pre-Run Verification)

To prevent batch effects and formatting incompatibilities:
*   [ ] **GATK Version Check:** Run `gatk --version` on both the server and local machine. Ensure the version is exactly **v4.6.2.0** on both.
*   [ ] **Java Version Check:** Verify that Java JRE/JDK matches (recommended OpenJDK 17) to prevent serialization errors.
*   [ ] **Reference Genome Hash Check:** Verify that the local genome FASTA matches the server genome FASTA exactly:
    *   Server: `md5sum _data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
    *   Local: `md5sum _data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
    *   *(Hashes must match exactly).*

---

## 🔄 2. Bidirectional Execution Coordination

To avoid race conditions and double-processing:
*   [ ] **Reverse Queue Traversal:** Ensure the local script is configured to process the cohort in reverse order (`sort -r`), while the HPC server processes in normal alphabetical order.
*   [ ] **Active Run Audit:** Run `squeue -u phongdinh` on the server before starting the local worker. Confirm which indices are running to verify that they are at opposite ends of the queue.
*   [ ] **Lock Verification:** Confirm that the local worker successfully creates a `${sample_name}.local_lock` file in `_data/logs/split_n_cigar/` on the server before transferring the BAM file.

---

## 🔍 3. Data Integrity & Integrity Checks (Post-Processing)

To ensure files are complete, uncorrupted, and properly integrated:
*   [ ] **Non-Zero Size Check:** Never rely on simple file existence checks. Verify that output files are non-empty:
    ```bash
    [ -s "_data/variant/split_n_cigar/sample.split.bam" ] && [ -s "_data/variant/split_n_cigar/sample.split.bai" ]
    ```
*   [ ] **BAM Header Validation:** Use `samtools quickcheck` to verify that the split BAM file has a valid EOF marker and header:
    ```bash
    samtools quickcheck _data/variant/split_n_cigar/sample.split.bam
    ```
*   [ ] **Index Verification:** Check that the uploaded `.split.bai` index file is valid and recognized:
    ```bash
    samtools idxstats _data/variant/split_n_cigar/sample.split.bam > /dev/null
    ```
*   [ ] **Stale Lock Cleanup:** In case of a local worker crash or script abort, proactively log onto the server and delete any orphaned `.local_lock` files to prevent pending Slurm array slots from getting permanently skipped.

---

## 🛠️ 4. HPC Pipeline Recovery & Integration

*   [ ] **Check skipped count:** Verify the server logs to confirm that finished samples uploaded by the local workstation are properly recognized and skipped by the Slurm jobs:
    ```
    [Skip] sample_name (Split CIGAR BAM already exists)
    ```
*   [ ] **HaplotypeCaller alignment:** Ensure that once Step 3 completes, HaplotypeCaller (Step 4) begins with 4 CPUs per task and accesses the correct split BAM paths on the server.
