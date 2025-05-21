import numpy as np

class MultipleSequenceAligner:
    """
        A class for performing multiple sequence alignment using a Center-Start-Method.
    """
    def __init__(self, sequences, match_score, mismatch_score, gap_penalty, match, substitution, gap):
        """
        Initializes the aligner with input sequences and scoring parameters.

        Args:
            sequences (list): List of (name, sequence) tuples.
            match_score (int): Score for a character match.
            mismatch_score (int): Score for a character mismatch.
            gap_penalty (int): Penalty for a gap.
            match (int): Score used in final alignment match scoring.
            substitution (int): Score used for mismatch in final scoring.
            gap (int): Penalty used for gap in final scoring.
        """
        self.sequences = sequences
        self.__match_score = match_score
        self.__mismatch_score = mismatch_score
        self.__gap_penalty = gap_penalty
        self.__match = match
        self.__substitution = substitution
        self.__gap = gap
        self.__matrices = self._fill_all_matrices()
        self.__central_sequence = self._find_central_sequence()
        self.__alignments_with_cs = self._align_sequences_along_with_cs()
        self.__merged_cs = self._merge_central_sequence()
        self.__final_alignments = self._compute_final_alignments()

    def _fill_all_matrices(self):
        """
        Generates pairwise alignment matrices for all sequence pairs.

        Returns:
            list: List of tuples (seq1, seq2, matrix).
        """
        matrices = []
        num_sequences = len(self.sequences)
        for i in range(num_sequences):
            for j in range(num_sequences):
                if j > i:
                    matrices.append((self.sequences[i], self.sequences[j],
                                     self._fill_matrix(self.sequences[i][1], self.sequences[j][1])))

        return matrices

    def _find_central_sequence(self):
        """
        Identifies the sequence with the highest cumulative alignment score.

        Returns:
            tuple: The central sequence (name, sequence).
        """
        all_scores = []
        for sequence in self.sequences:
            score = 0
            for se1, seq2, matrix in self.__matrices:
                if sequence == se1 or sequence == seq2:
                    score += self._get_matrix_score(matrix)
            all_scores.append((sequence, score))

        central_sequence = max(all_scores, key=lambda x: x[1])[0]

        return central_sequence

    def _align_sequences_along_with_cs(self):
        """
        Aligns each sequence to the central sequence.

        Returns:
            list: List of aligned sequence pairs.
        """
        alignments = []
        for seq1,seq2,matrix in self.__matrices:
            if self.__central_sequence == seq1:
                alignments.append(self._align_two_sequences(seq1,seq2,matrix))

            elif self.__central_sequence == seq2:
                (name1, align1), (name2, align2) = self._align_two_sequences(seq1,seq2,matrix)
                alignments.append(((name2, align2), (name1, align1)))


        return alignments

    def _merge_central_sequence(self):
        """
        Merges the aligned central sequence with consistent gaps.

        Returns:
            tuple: Central sequence with gaps merged across all alignments.
        """
        alignments = [alignment[0][1] for alignment in self.__alignments_with_cs]
        alignments = [list("".join(seq)) for seq in alignments]
        i = 0
        merged_cs = []

        while True:
            # Stop if all sequences are passed
            if all(i >= len(seq) for seq in alignments):
                break

            # Get current column in case a sequence is shorter, uses gap
            current_col = [seq[i] if i < len(seq) else '-' for seq in alignments]

            # If any sequence has a gap it inserts a gap into all that don't
            if '-' in current_col:
                for seq in alignments:
                    if i >= len(seq) or seq[i] != '-':
                        seq.insert(i, '-')
                merged_cs.append('-')
                i += 1
            else:
                merged_cs.append(alignments[0][i])
                i += 1

        seq_name = self.__alignments_with_cs[0][0][0]

        return (seq_name,merged_cs)

    def _compute_final_alignments(self):
        """
        Inserts missing gaps and return final aligned sequences.

        Returns:
            list: Fully aligned sequences.
        """
        final_alignments = []
        alignments = [alignment[1] for alignment in self.__alignments_with_cs]
        gaps_indexes = [i for i, char in enumerate(self.__merged_cs[1]) if char == '-']

        # Appends the central sequence
        final_alignments.append(self.__merged_cs)

        # Inserts gaps where needed
        for alignment in alignments:
            new_alignment = alignment[1]
            if len(new_alignment) < len(self.__merged_cs[1]):
                for gap_index in gaps_indexes:
                    if len(new_alignment) < len(self.__merged_cs[1]):
                        new_alignment.insert(gap_index, '-')
                final_alignments.append((alignment[0], new_alignment))
            else:
                final_alignments.append(alignment)

        return final_alignments

    def get_final_alignments(self):
        """
        Gets the list of final aligned sequences.

        Returns:
            ist: Aligned sequences with names.
        """
        return self.__final_alignments

    def _fill_matrix(self, first_seq, second_seq):
        """
        Constructs the scoring matrix for a single optimal path using dynamic programming.

        Returns:
            numpy.ndarray: The filled scoring matrix with backtracking pointers (int, (int, int)).
        """
        matrix = np.empty((len(first_seq) + 1, len(second_seq) + 1), dtype=object)
        first_num = 0

        # Initialize first row
        for i in range(len(second_seq) + 1):
            matrix[0][i] = (first_num, (0, i - 1))
            first_num += self.__gap_penalty

        # Initialize first column
        first_num = 0
        for i in range(len(first_seq) + 1):
            matrix[i][0] = (first_num, (i - 1, 0))
            first_num += self.__gap_penalty

        # Fill rest of the matrix
        for i in range(1, len(first_seq) + 1):
            for j in range(1, len(second_seq) + 1):
                vertical = matrix[i - 1][j][0] + self.__gap_penalty
                horizontal = matrix[i][j - 1][0] + self.__gap_penalty
                if first_seq[i - 1] == second_seq[j - 1]:
                    diagonal = matrix[i - 1][j - 1][0] + self.__match_score
                else:
                    diagonal = matrix[i - 1][j - 1][0] + self.__mismatch_score

                if vertical == max(vertical, horizontal, diagonal):
                    matrix[i][j] = (vertical, (i - 1, j))
                elif horizontal == max(horizontal, vertical, diagonal):
                    matrix[i][j] = (horizontal, (i, j - 1))
                else:
                    matrix[i][j] = (diagonal, (i - 1, j - 1))

        return matrix

    def _align_two_sequences(self, first_seq_inp, second_seq_inp, matrix):
        """
        Retrieve the optimal alignment path and aligned sequences.

        Args:
            first_seq_inp (tuple): (name, sequence_string)
            second_seq_inp (tuple): (name, sequence_string)
            matrix (list): Scoring matrix with traceback info.

        Returns:
            tuple: A tuple containing two aligned sequences in tuple form:
                   ((name1, aligned_sequence1),(name2, aligned_sequence2)
        """


        name1, seq1 = first_seq_inp
        name2, seq2 = second_seq_inp

        optimal_path = []
        align1 = []
        align2 = []
        i = len(matrix) - 1
        j = len(matrix[0]) - 1

        # Traceback path
        while matrix[i][j][1] != (0, 0):
            optimal_path.append((i, j))
            i, j = matrix[i][j][1]
        optimal_path.append((i, j))
        optimal_path.append((0, 0))
        optimal_path = optimal_path[::-1]

        # Create aligned sequences
        for i in range(len(optimal_path) - 1):
            if optimal_path[i][0] == optimal_path[i + 1][0]:
                align1.append("-")
                align2.append(seq2[optimal_path[i + 1][1] - 1])
            elif optimal_path[i][1] == optimal_path[i + 1][1]:
                align1.append(seq1[optimal_path[i + 1][0] - 1])
                align2.append("-")
            else:
                align1.append(seq1[optimal_path[i + 1][0] - 1])
                align2.append(seq2[optimal_path[i + 1][1] - 1])

        return ((name1, align1), (name2, align2))

    def _get_matrix_score(self, matrix):
        """
        Retrieve the highest alignment score from the matrix.

        Returns:
            str: The alignment score as a string.
        """
        return matrix[len(matrix) - 1][len(matrix[0]) - 1][0]

    def _get_matrices(self):
        """
        Get the scoring matrix for a single optimal alignment path with backtracking pointers for each cell in the matrix

        Returns:
            numpy.ndarray: The filled scoring matrix with backtracking pointers.
        """
        return self.__matrices

    def get_statistics(self):
        """
        Calculates statistics of the final alignment (match, mismatch, gap, identity%).

        Returns:
            dict: Statistics of alignment.
        """
        num_match = 0
        num_mismatch = 0
        num_gap = 0
        total_columns = len(self.__final_alignments[0][1])  # Assuming all sequences have the same length

        for i in range(total_columns):
            column_nucleotides = [seq[i] for _, seq in self.__final_alignments]

            if all(nucl == column_nucleotides[0] for nucl in column_nucleotides):
                num_match += 1

            elif any(nucl == '-' for nucl in column_nucleotides):
                num_gap += 1
            else:
                num_mismatch += 1

        identity_percent = (num_match / total_columns) * 100

        return {
            "match": num_match,
            "mismatch": num_mismatch,
            "gap": num_gap,
            "identity_percent": round(identity_percent, 2)
        }

    def get_score(self):
        """
        Computes the final alignment score for all sequence pairs.

        Returns:
            int: Total alignment score.
        """
        sum = 0
        for k in range(len(self.__final_alignments)):
            for l in range(len(self.__final_alignments)):
                if l > k:
                    length = len(self.__final_alignments[0][1])
                    seq1 = self.__final_alignments[k][1]
                    seq2 = self.__final_alignments[l][1]
                    for i in range(length):
                        if seq1[i] == "-" and seq2[i] == "-":
                            sum += 0
                        elif seq1[i] == "-" or seq2[i] == "-":
                            sum += self.__gap
                        elif seq1[i] == seq2[i]:
                            sum += self.__match
                        else:
                            sum += self.__substitution
        return sum