# BEE-102-Spring-2025-Assignment
This is a Computational Biology Assignment in our course BEE-102 Introduction to Computational Biology

# Comprehensive Analysis of Genomic and Transcriptomic Data

## Q1: Fragment Length Distribution Matching  
### Objective  
Align fragment length distributions between experimental (query) and reference datasets to enable comparative analysis.  

### Algorithm  
**Data Input**:  
- Query data: `query.bed` (genomic regions)  
- Reference data: `reference.hist` (fragment length distribution)  

**Binning Strategy**:  
bin_width = 10
bins = range(0, 720+bin_width, bin_width) # 72 bins (0-10, 10-20,...710-720)


**Distribution Normalization**:  
query_normalized = query_counts / total_query
reference_normalized = reference_counts / total_reference


**Stochastic Subsampling**:  
For each bin:  
- `n = min(query_count, reference_count)`  
- Randomly select `n` fragments from reference using `random_state=42` for reproducibility  

### Visualization  
![Fragment Distribution](attachment:imagerence_distribution.png)  
- **Purple**: Query distribution  
- **Blue stars**: Subsamples from reference matching query's bin counts  

### Key Insight  
Maintains reference dataset characteristics while matching query's fragment profile, crucial for comparative analyses.

---

## Q2: Markov Transition Matrix Construction  
### Biological Context  
Model nucleotide transitions in DNA sequences using a first-order Markov chain with states: `{A, C, G, T}`.  

### Matrix Construction Process  
1. **Count Transitions**:  
for i in range(len(sequence)-1):
counts[current][next] += 1


2. **Calculate Probabilities**:  
for base in counts:
total = sum(counts[base].values())
for next_base in counts[base]:
matrix[base][next_base] = counts[base][next_base]/total


### Resulting Matrix  
|       | A     | C     | G     | T     |
|-------|-------|-------|-------|-------|
| **A** | 0.0625| 0.100 | 0.200 | 0.520 |
| **C** | 0.2500| 0.000 | 0.200 | 0.160 |
| **G** | 0.1250| 0.000 | 0.200 | 0.120 |
| **T** | 0.5625| 0.900 | 0.400 | 0.200 |

### Key Observation  
Strong **T→C** and **G→T** transitions suggest specific sequence context preferences.

---

## Q3: FASTA Format Conversion  
### Problem Solved  
Convert multi-line FASTA to single-line format for computational efficiency:  

**Input**:  
"Gene1
ATGCG
TAGCT
Gene2
CGATA"


**Output**:  
"Gene1
ATGCGTAGCT
Gene2
CGATA"


### Algorithm  
- Header detection with `startswith('>')`  
- Sequence concatenation  
- **Space Complexity**: `O(1)` through stream processing (efficient for large genomes).  

---

## Q4: Viterbi Algorithm Implementation  
### Hidden Markov Model Parameters  
states = ['E', '5', 'I']
trans_prob = {
'E': {'E': 0.9, '5': 0.1},
'5': {'I': 1.0},
'I': {'I': 0.9, 'end': 0.1}
}
emit_prob = {
'E': {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
'5': {'A':0.05, 'C':0.0, 'G':0.95, 'T':0.0},
'I': {'A':0.4, 'C':0.1, 'G':0.1, 'T':0.4}
}


### Performance Metrics  
- Given path score: `-41.22`  
- Viterbi optimal path score: `-38.68`  
- **6.2% likelihood improvement** over manual annotation  

### State Path Visualization  
![State Path Comparison](attachment:viterbi_path.png)  

---

## Q5: V-Plot Analysis  
### Technical Implementation  
1. **Midpoint Calculation**:  
mid1 = (start1 + end1)/2
mid2 = (start2 + end2)/2
X = mid2 - mid1 # Distance between fragment centers
Y = fragment_length # Fragment length


2. **2D Histogram**:  
plt.hist2d(X, Y, bins=(100,100), cmap='Blues')


### Biological Interpretation  
**V-shaped pattern** indicates nucleosome positioning:  
- **X-axis**: Distance between fragment centers  
- **Y-axis**: Fragment length  
- Density reveals chromatin accessibility patterns  

---

## Q6: Principal Component Analysis  
### Gene Expression Patterns  
![PCA Plot](attachment:pca_plot.png)  

- **PC1 (87.6% variance)**: `[-0.894, 0.447]` (strong negative correlation with *GATA3*)  
- **PC2 (12.4% variance)**: `[0.447, 0.894]` (positive correlation with both genes)  

### Clinical Correlation  
- **ER+ samples** show **2.3× higher PC1 scores** (*p=0.0037*), suggesting *GATA3* as a key biomarker for estrogen receptor status.  

---

This analysis integrates bioinformatics techniques-from sequence analysis to machine learning-to reveal genomic regulation insights and clinical biomarkers.

