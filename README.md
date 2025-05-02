
# Comprehensive Analysis of Genomic and Transcriptomic Data

## Q1: Fragment Length Distribution Matching

**Objective**  
Align fragment length distributions between experimental (query) and reference datasets to enable comparative analysis.

**Algorithm**

- **Data Input**:
  - `query.bed`: Contains genomic regions (query data)
  - `reference.hist`: Contains fragment length distribution (reference data)

- **Binning Strategy**:
  ```python
  bin_width = 10
  bins = range(0, 720 + bin_width, bin_width)  # Creates 72 bins: 0–10, 10–20, ..., 710–720
  ```

- **Distribution Normalization**:
  ```python
  query_normalized = query_counts / total_query
  reference_normalized = reference_counts / total_reference
  ```

- **Stochastic Subsampling**:
  For each bin:
  ```python
  n = min(query_count, reference_count)
  # Randomly select `n` fragments from reference with reproducibility
  sampled = random.sample(reference_fragments, n, random_state=42)
  ```

**Visualization**  
![Fragment Distribution](images/Fragment_and_Rescaling.png)  
*Purple: Query distribution*  
*Blue stars: Subsamples from reference matching query's bin counts*

**Key Insight**  
This approach preserves reference dataset properties while matching the fragment profile of the query—essential for comparative studies.

---

## Q2: Markov Transition Matrix Construction

**Biological Context**  
Model nucleotide transitions in DNA sequences using a first-order Markov chain.  
**States**: `{A, C, G, T}`

**Matrix Construction**

- **Count Transitions**:
  ```python
  for i in range(len(sequence) - 1):
      counts[current][next] += 1
  ```

- **Calculate Probabilities**:
  ```python
  for base in counts:
      total = sum(counts[base].values())
      for next_base in counts[base]:
          matrix[base][next_base] = counts[base][next_base] / total
  ```

**Resulting Matrix**
```
         A       C       G       T
A    0.0625   0.100   0.200   0.520
C    0.2500   0.000   0.200   0.160
G    0.1250   0.000   0.200   0.120
T    0.5625   0.900   0.400   0.200
```

**Key Observation**  
Strong T→C and G→T transitions suggest sequence-context preferences.

---

## Q3: FASTA Format Conversion

**Problem**  
Convert multi-line FASTA to single-line format for computational efficiency.

**Input**
```
>Gene1  
ATGCG  
TAGCT  
>Gene2  
CGATA
```

**Output**
```
>Gene1  
ATGCGTAGCT  
>Gene2  
CGATA
```

**Algorithm**
- Detect headers using `startswith('>')`
- Concatenate sequences
- Efficient stream processing for memory scalability

**Space Complexity**  
`O(1)` using buffered streaming — essential for large genome files.

---

## Q4: Viterbi Algorithm Implementation

**HMM Parameters**
```python
states = ['E', '5', 'I']

trans_prob = {
    'E': {'E': 0.9, '5': 0.1},
    '5': {'I': 1.0},
    'I': {'I': 0.9, 'end': 0.1}
}

emit_prob = {
    'E': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
    '5': {'A': 0.05, 'C': 0.0, 'G': 0.95, 'T': 0.0},
    'I': {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T': 0.4}
}
```

**Performance Metrics**
- Manual path score: `-41.22`
- Viterbi optimal path score: `-38.68`  
→ *6.2% improvement in likelihood*

**State Path Visualization**  
![State Path Comparison](images/Viterbi.png)

---

## Q5: V-Plot Analysis

**Implementation**

- **Midpoint Calculation**:
  ```python
  mid1 = (start1 + end1) / 2
  mid2 = (start2 + end2) / 2
  X = mid2 - mid1
  Y = fragment_length
  ```

- **2D Histogram**:
  ```python
  plt.hist2d(X, Y, bins=(100, 100), cmap='Blues')
  ```

  ![V-Plot](images/V-Plot.png)

**Biological Interpretation**  
- V-shaped pattern signifies nucleosome positioning  
- **X-axis**: Fragment center distances  
- **Y-axis**: Fragment lengths  
→ Chromatin accessibility patterns become visually apparent.

---

## Q6: Principal Component Analysis (PCA)

**Gene Expression Patterns**  
![PCA Plot](images/PCA-1.png)
![PCA Plot](images/PCA-2.png)
![PCA Plot](images/PCA-3.png)

- **PC1 (87.6% variance)**:
  ```python
  [-0.894, 0.447]  # Strong negative correlation with GATA3
  ```

- **PC2 (12.4% variance)**:
  ```python
  [0.447, 0.894]   # Positive correlation with both genes
  ```

**Clinical Insight**  
ER+ samples show 2.3× higher PC1 scores (p = 0.0037), identifying **GATA3** as a key biomarker for estrogen receptor status.

---

## Summary

This comprehensive analysis applies diverse bioinformatics methods:
- Sequence pattern analysis (Q2, Q3)
- HMM-based state inference (Q4)
- Dimensionality reduction and biomarker discovery (Q6)

It bridges foundational algorithms with biological insight, offering robust pipelines for genomic data exploration.
