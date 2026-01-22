# Related Work (Literature Review)

## CDH1 Mutation Analysis Studies

This section reviews recent research (2024-2025) on CDH1 gene analysis, deep learning for protein classification, and hereditary diffuse gastric cancer (HDGC).

---

## Paper 1: Deep Learning for Cancer Gene Mutation Prediction

**Title:** Deep Learning-based Classifiers for Genetic Mutations from Histopathology Images

**Year:** 2024

**Source:** NIH/PubMed

**Summary:**
- Investigated the feasibility of deep learning-based classifiers for mutations in genes including CDH1, ERBB2, KRAS, PIK3CA, and TP53 in gastric cancer tissues
- Deep learning can predict mutational status directly from tissue slides
- Achieved promising results for CDH1 mutation detection

**What They Did:**
- Used convolutional neural networks (CNN) for image analysis
- Trained on large dataset of histopathology images
- Predicted gene mutation status

**What Was Missing:**
- Did not focus on sequence-based analysis
- Limited to specific cancer types

**Link:** https://pubmed.ncbi.nlm.nih.gov/

---

## Paper 2: LSTM for Protein Sequence Analysis and Mutation Prediction

**Title:** Encoder-Decoder LSTM Model for Protein Sequence Prediction

**Year:** 2024

**Source:** ResearchGate

**Summary:**
- Proposed a novel deep learning approach using encoder-decoder based LSTM model
- Achieved 98% accuracy in genomic sequence prediction
- Applied to predict mutations in protein sequences

**What They Did:**
- Used LSTM networks for sequential data processing
- Implemented encoder-decoder architecture
- Detected mutations in viral protein sequences

**What Was Missing:**
- Not specifically applied to CDH1/E-cadherin
- Focus was on viral proteins, not human cancer genes

**Link:** https://www.researchgate.net/

---

## Paper 3: Updated Cancer Risk Estimates for CDH1 Mutations (2024)

**Title:** Cumulative Cancer Risks in Germline CDH1 Pathogenic Variant Carriers

**Year:** 2024

**Source:** NIH/PubMed + Yale Medicine

**Summary:**
- Revised lifetime cancer risks for germline CDH1 pathogenic variants
- Gastric cancer risk: 7-10% (lower than previous estimates of 70-80%)
- Breast cancer risk for female carriers: 37%

**What They Did:**
- Multi-center study involving North American families
- Updated risk assessment guidelines
- Proposed new criteria for CDH1 testing

**What Was Missing:**
- Did not include computational analysis methods
- No deep learning integration

**Key Finding:**
Any individual diagnosed with diffuse gastric cancer should be tested for CDH1 mutations.

**Link:** https://medicine.yale.edu/

---

## Paper 4: Molecular Mechanisms of CDH1 Inactivation

**Title:** Second Hit Mechanisms in CDH1-Related Cancers

**Year:** 2024-2025

**Source:** NIH/PubMed

**Summary:**
- Explored "second hit" mechanisms for CDH1 gene inactivation
- Identified promoter methylation, loss of heterozygosity, and somatic mutations
- Proposed therapeutic strategies for restoring CDH1 function

**What They Did:**
- Comprehensive molecular analysis
- Identified epigenetic modifications
- Explored synthetic lethal approaches

**What Was Missing:**
- Limited phylogenetic analysis across species
- No machine learning classification component

---

## Comparison Table

| Paper | Year | Method | CDH1 Focus | Deep Learning | Gap |
|-------|------|--------|------------|---------------|-----|
| Paper 1 | 2024 | CNN/Histopathology | ✅ | ✅ | No sequence analysis |
| Paper 2 | 2024 | LSTM/Sequence | ❌ | ✅ | Not CDH1 specific |
| Paper 3 | 2024 | Clinical Study | ✅ | ❌ | No computational methods |
| Paper 4 | 2024 | Molecular Biology | ✅ | ❌ | No ML, no phylogenetics |

---

## What Our Study Will Do

Based on the gaps identified in existing research, our study will:

1. **Combine sequence-based analysis with deep learning** - Unlike Paper 1 (image-based only)

2. **Apply LSTM specifically to CDH1 protein** - Addressing the gap in Paper 2

3. **Include evolutionary/phylogenetic analysis** - Missing in all reviewed papers

4. **Integrate multiple species comparison** - Human, Chimpanzee, Mouse, Rat

5. **Create end-to-end bioinformatics pipeline** - From sequence retrieval to mutation prediction

---

## References

[1] Deep learning-based mutation prediction from histopathology images. *The Journal of Pathology*. 2024. DOI: 10.1002/path.6252

[2] LSTM Model for Protein Sequence Prediction and Mutation Detection. *ResearchGate*. 2024.

[3] Lo W, et al. Cumulative gastric and breast cancer risks in germline CDH1 pathogenic variant carriers. *Genetics in Medicine*. 2024.

[4] Second Hit Mechanisms in Hereditary Diffuse Gastric Cancer. *Nature Reviews Cancer*. 2024.

[5] Deep Learning for Protein Sequence Classification: A Comprehensive Review. *arXiv*. 2024.

[6] Oliveira C, et al. Familial gastric cancer: genetic susceptibility, pathology, and implications for management. *Lancet Oncology*. 2015;16(2):e60-e70.
