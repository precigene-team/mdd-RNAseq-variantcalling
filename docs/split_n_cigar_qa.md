<!-- [[docs/split_n_cigar_qa.md]] -->
<!-- Goal: Document the Q&A conversation about SplitNCigar and HaplotypeCaller progress and runtimes. -->

# Q&A: SplitNCigar and HaplotypeCaller Execution Progress

**Date:** July 5, 2026  
**Time:** 00:03 ICT  

---

### **Q1: Why is SplitNCigar taking 7+ hours, and why has not a single job completed?**
* **A:** 
  1. **Massive Data Scale:** The active jobs are processing control/MDD `_input` BAM files, which are very large (**15.89 GB to 24.34 GB** in size, containing **1.1+ billion reads** each). 
  2. **Baseline Benchmark:** GATK `SplitNCigarReads` is single-threaded and I/O bound. A benchmark of a completed large sample (`female_ctrl_1_input`, 20.89 GB) showed it takes **~13.7 hours** to run.
  3. **Progress Verification:** The logs indicate the active jobs are currently on **Chromosome 8**, which means they are roughly **50% complete** and on track to finish overnight (13–15 hours total).
  4. **Queue Concurrency Limit:** The job array has a concurrency limit of 5 (`VC_ARRAY_CONCURRENCY=5`). Since these 5 large tasks started at the exact same time, they are occupying all 5 available slots, leaving other tasks queued in pending state.

---

### **Q2: Is there a faster way to run SplitNCigarReads?**
* **A:** Yes, there are three primary ways to accelerate it:
  1. **Increase Concurrency (`VC_ARRAY_CONCURRENCY`):** Increase the job array concurrency limit (e.g., to 10 or 12) in the orchestrator script to run more samples in parallel. This is the simplest change, but is limited by the server's CPU/RAM and disk write capacity.
  2. **Use GATK4 Spark Mode (`SplitNCigarReadsSpark`):** Run the Spark-based version of the tool with the `--spark-master local[N]` argument to allow multi-threading (e.g., 4 or 8 CPUs) within a single sample, speeding it up by 3–4x.
  3. **Scatter-Gather Parallelization by Chromosome:** Split the genome into intervals, run the jobs concurrently by chromosome using `-L`, and merge the output BAMs afterward. This is the standard GATK best practice, but adds significant script complexity.

---

### **Q3: Is HaplotypeCaller (Step 4) going to take this long, or will it be the same/worse?**
* **A:** HaplotypeCaller will be **significantly faster** (estimated at **3 to 6 hours** per large sample, which is 2x to 3x faster than SplitNCigar) due to:
  1. **Exon-only target:** In RNA-seq, HaplotypeCaller only processes transcribed regions (exons), ignoring introns and intergenic regions (~97% of the genome).
  2. **Low I/O Overhead:** Unlike SplitNCigar which splits reads and writes out a new, massive **20+ GB BAM file**, HaplotypeCaller only writes out a small, compressed **VCF file** (a few megabytes), avoiding the disk I/O bottleneck.
