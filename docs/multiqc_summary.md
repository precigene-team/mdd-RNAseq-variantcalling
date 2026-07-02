# MultiQC Upstream Diagnostics Summary Report
**MDD PFC Cohort (78 Libraries: 39 Inputs vs. 39 IP/m6A Samples)**

This report summarizes the raw metrics extracted from the aggregated MultiQC workflow report (`llms-full.txt`).

---

## 📊 Summary Statistics Table

| Metric | Sample Type | Mean | Min | Max | Sample Count |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Total Reads (M)** | IP (m6A) | 58.05M | 33.41M | 124.66M | 39 |
| | Input | 66.02M | 19.97M | 97.05M | 39 |
| **Uniquely Aligned Reads (%)** | IP (m6A) | 93.25% | 91.21% | 94.64% | 39 |
| | Input | 52.60% | 20.73% | 78.50% | 39 |
| **Multimapped Reads (%)** | IP (m6A) | 2.49% | 1.05% | 5.18% | 39 |
| | Input | 29.03% | 8.05% | 59.63% | 39 |
| **Raw Duplication Rate (%)** | IP (m6A) | 64.52% | 52.24% | 77.00% | 78 |
| | Input | 57.55% | 24.92% | 79.94% | 78 |
| **featureCounts Assignment Rate (%)** | IP (m6A) | 68.24% | 58.48% | 76.01% | 39 |
| | Input | 8.30% | 1.77% | 20.47% | 39 |

---

## 🔍 Key Findings & Diagnostics

### 1. High Multi-Mapping Rates in Input Samples
* **Observation**: Uniquely aligned reads for Input samples average **52.60%**, with some samples dropping as low as **20.73%**. Conversely, the IP samples show highly efficient alignment, averaging **93.25%** uniquely aligned.
* **Explanation**: The Input fraction captures baseline total RNA, which is dominated by ribosomal RNA (rRNA) that escaped pre-mapping depletion. These rRNA reads map to the hundreds of repetitive rDNA locus units across the genome, resulting in an average multi-mapping rate of **29.03%** for Inputs, compared to just **2.49%** for IP.

### 2. Duplication Profile Differences
* **Observation**: Raw duplication levels are significantly higher in the IP (m6A) fraction, averaging **64.52%** compared to the Input fraction's **57.55%**.
* **Explanation**: The IP fraction selects specifically for methylated transcripts ($m^6A$), which drastically reduces library complexity. This bottlenecking effect naturally leads to higher sequence duplication levels during PCR amplification, which is normal and expected for MeRIP-seq libraries.

### 3. featureCounts Assignment Rates
* **Observation**: featureCounts assignment rates are low for Inputs (average **8.30%**), while IP samples show high feature assignment (average **68.24%**).
* **Explanation**: The low assignment in Inputs is a consequence of the high multi-mapping rate (ribosomal sequences) described above. Since featureCounts assigns only uniquely mapping reads by default, the multi-mapping alignments are discarded, artificially lowering the raw assignment rate calculation.
