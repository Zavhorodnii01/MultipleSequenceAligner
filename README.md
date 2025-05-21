# MSA - Multiple Sequence Alignment (Center Star Method)

This project implements the Center Star method for Multiple Sequence Alignment (MSA) as part of the Introduction to Bioinformatics course at Wrocław University of Science and Technology.

## Overview

This Python program performs MSA using the Center Star algorithm. It provides flexible input options, customizable scoring schemes, and a graphical representation of the optimal alignment.  The results, including alignment statistics, can be saved to a file.  The program features an interactive user GUI for ease of use.

## Features

* **Sequence Input:**
    * Manual input
    * FASTA file import
* **Customizable Scoring:**  Adjust match scores, mismatch penalties, and gap penalties.
* **Optimal Alignment Generation:**
    * Graphical display of the alignment with parameters and scoring information.
    * Output saving to a file (similar format to Task #2).
* **Alignment Statistics:**
    * Identity percentage
    * Number of matches, mismatches, and gaps
* **Interactive GUI:** User-friendly interface for streamlined execution.

## Algorithm

The Center Star method works by:

1. **Finding the center sequence:** The sequence with the minimum sum of distances to all other sequences is selected as the center.
2. **Progressive alignment:**  Remaining sequences are progressively aligned to the center sequence using pairwise alignment techniques.  The resulting alignments form the initial MSA.
3. **Refinement (optional):** Iterative refinement methods can be used to improve the alignment score.
