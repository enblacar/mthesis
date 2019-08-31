#  __           _       _      ___   ___            _ _                                       _ _      _                 
# / _\ ___ _ __(_)_ __ | |_   / _ \ / _ \ _   _ __ (_) | ___ _   _ _ __    _ __  _ __ ___  __| (_) ___| |_ ___  _ __ ___ 
# \ \ / __| '__| | '_ \| __| | | | | (_) (_) | '_ \| | |/ _ \ | | | '_ \  | '_ \| '__/ _ \/ _` | |/ __| __/ _ \| '__/ __|
# _\ \ (__| |  | | |_) | |_  | |_| |\__, |_  | |_) | | |  __/ |_| | |_) | | |_) | | |  __/ (_| | | (__| || (_) | |  \__ \
# \__/\___|_|  |_| .__/ \__|  \___/   /_/(_) | .__/|_|_|\___|\__,_| .__/  | .__/|_|  \___|\__,_|_|\___|\__\___/|_|  |___/
#                |_|                         |_|                  |_|     |_|                                            

# Script: 09_pileup_predictors.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: The objective of this script is to generate a set of pileup predictors for ML models. They are based on the genomic context of each read, 
# and take into account the MAPQ score and the QUAL of each nucleotide, as well as the coverage of each nucleotide and its consensus towards the variant.

import sys
import os
from math import log10
import pandas as pd
import numpy as np
import copy

class GeneratePredictors:
    def __init__(self, sam, pileup_range):
        self.sam = sam
        self.sam_name = self.sam.split("/")[-1].split(".")[0]
        self.pileup_range = int(pileup_range)
        self.process_SAM()


    def load_matrix(self, first, last, tsv):
        """This function load the TSV file that contains all the reads necessary for the generation 
        of the pileup into memory as a Dataframe. This is done by generating a temporary TSV file that
        is later on removed."""

        # Filter from the TSV file the reads that can be piled-up and store them in a temporary file. 
        os.system(f"""awk '$1 >= {first} && $1 <= {last} {{print $0}}' {tsv} > ./data/tmp/range_{self.sam_name}.tsv""")
        # Load the reads into memory using a DataFrame.
        matrix = pd.read_csv(f"""./data/tmp/range_{self.sam_name}.tsv""", sep = "\t", names = ["POS", "SEQ", "QUAL", "MAPQ", "CIGAR"])
        # Remove the temporary file. 
        os.system(f"""rm ./data/tmp/range_{self.sam_name}.tsv""")
        return matrix

   def process_cigar(self, cigar):
        """This function transform a given CIGAR string into a list of CIGAR operations."""

        # Placeholder for the operations in the CIGAR string.
        list_operations = []
        # Placeholder for each individual operation. 
        string = ""
        if isinstance(cigar, float):
            print(cigar)
        for letter in cigar:
            # First, the numbers come in for the CIGAR.
            try:
                test = int(letter)
                string += letter
            # Once we get to the operation, we add it to the list.
            except:
                string += letter
                list_operations.append(string)
                # Clear the string for the new operation.
                string = ""
        count = 0
        for operation in list_operations:
            num = operation[:-1]
            count += int(num)
        final_list = []
        if count > 100:
            remaining = count - 100
            for operation in list_operations[::-1]:
                if remaining == 0:
                    final_list.append(operation)
                    continue
                num = int(operation[:-1])
                if num > remaining:
                    num -= remaining
                    remaining = 0
                    string = f"{num}{operation[-1]}"
                    final_list.append(string)
                elif num == remaining:
                    remaining = 0
                elif num < remaining:
                    remaining -= num
        else:
            final_list_ordered = list_operations
        final_list_ordered = final_list[::-1]

        return final_list_ordered


    def modify_SEQ(self, seq_in, cigar_list_in):
        """This function transform a given sequence based on a cigar list."""
        seq = seq_in[:] # Make a copy.
        cigar_list = cigar_list_in[:]
        # Placeholder for the new sequence.
        new_seq = ""
        for item in cigar_list:
            # Number of operations.
            num = int(item[:-1])
            # Operation.
            letter = item[-1]
            if letter == "M" and num == len(seq_in):
                return seq_in
            if True:
                # Matches or mismatches.
                if letter in ["M", "X"]:
                    new_seq += seq[:num]
                    seq = seq[num:]

                # Hard-clips or skipped regions.
                elif letter in ["H", "N"]:
                    seq = seq[num:]
                    new_seq += num * " "
                # Deletions.
                elif letter == "D":
                    seq = seq[num:]
                    new_seq += num * "~"
                # Paddings, insertions, soft-clips.
                elif letter in ["P", "I", "S"]:
                    seq = seq[num:]
                # Sequence match.
                elif letter == "=":
                    new_seq = seq

        return new_seq

    def process_SEQ(self, seq_in, original_pos, read_pos, original_end, read_end):
        """This function, given the beginning and end coordinates of the template and target read, 
        trims the targe read so it can form a pileup."""
        if original_pos > read_end or original_end < read_pos:
            return "Remove"
        # By having the positions, we can address which case are we dealing with, and transform the sequence accordingly.
        if original_pos == read_pos and original_end == read_end:
            where = "=="
            cutoff = 0
            seq = seq_in
        elif original_pos > read_pos and original_end == read_end:
            where = "+="
            cutoff = original_pos - read_pos
            seq = seq_in[cutoff:]
        elif original_pos > read_pos and original_end > read_end:
            where = "++"
            cutoff_init = original_pos - read_pos
            cutoff_end = original_end - read_end
            seq = seq_in[cutoff_init:] + " " * cutoff_end
        elif original_pos == read_pos and original_end > read_end:
            where = "=+"
            cutoff_end = original_end - read_end
            seq = seq_in + " " * cutoff_end
        elif original_pos == read_pos and original_end < read_end:
            where = "=-"
            cutoff_end = original_end - read_end # it will be negative
            seq = seq_in[:cutoff_end]
        elif original_pos < read_pos and original_end == read_end:
            where = "-="
            cutoff_init = read_pos - original_pos
            seq = cutoff_init * " " + seq_in
        elif original_pos < read_pos and original_end < read_end:
            where = "--"
            cutoff_init = read_pos - original_pos
            cutoff_end = original_end - read_end
            seq = cutoff_init * " " + seq_in[:cutoff_end]
        elif original_pos > read_pos and original_end < read_end:
            where = "+-"
            cutoff_init = original_pos - read_pos
            cutoff_end = original_end - read_end
            seq = seq_in[cutoff_init:cutoff_end]
        elif original_pos < read_pos and original_end > read_end:
            where = "-+"
            cutoff_init = read_pos - original_pos
            cutoff_end = original_end - read_end
            seq = cutoff_init * " " + seq_in + cutoff_end * " "
        return seq

    def transform_SEQ(self, seq_in, cigar_list, original_pos, read_pos, original_end):
        """This function acts as a wrapper for other functions. It takes a sequence and 
        returns the sequence for the pileup."""
        # Transform the sequence according to the cigar so its end can be computed.
        seq_middle = self.modify_SEQ(seq_in, cigar_list)
        # Compute the read end once the sequence has been processed. 
        read_end = read_pos + len(seq_middle) - 1
        # Obtain the sequence for the pileup. 
        seq_out = self.process_SEQ(seq_middle, original_pos, read_pos, original_end, read_end)
        return seq_out

    def calculate_predictors(self, matrix, seq, mapq):

        # Transform ASCII to Q:
        def ASCII_to_Q(x): return ord(x) - 33
        # Return the probability of error
        def Q_to_prob(x): return 10 ** (-0.1 * x)
        # Wrapper
        def transform_Phred(x): return 1 - Q_to_prob(ASCII_to_Q(x)) if x != "¬" else "¬"
        # Obtain the pileup for the SEQ column.
        seq_matrix = pd.DataFrame(list(item) for item in matrix["PILEUP_SEQ"]) # Transform into a list of lists
        # Obtain the pileup for the QUAL column.
        qual_matrix = pd.DataFrame(list(item) for item in  matrix["PILEUP_QUAL"])
        for index in range(len(qual_matrix.columns)):
            qual_matrix[index] = qual_matrix[index].apply(lambda x: transform_Phred(x))
        # Generate a list for the consensus and coverate predictors.
        consensus_list = []
        coverage_list = []
        for i in range(len(seq_matrix.columns)):
            c = seq_matrix[i]
            coverage = c[c != " "].shape[0]
            coverage_list.append(coverage)
            mask_A = c == "A"
            mask_C = c == "C"
            mask_T = c == "T"
            mask_G = c == "G"
            mask_N = c == "N"
            c_A, c_C, c_T, c_G, c_N = c[mask_A], c[mask_C], c[mask_T], c[mask_G], c[mask_N]
            c_A, c_C, c_T, c_G, c_N = c_A.apply(lambda x: 1), c_C.apply(lambda x: 1), c_T.apply(lambda x: 1), c_G.apply(lambda x: 1), c_N.apply(lambda x: 1)
            p_A, p_C, p_T, p_G, p_N = qual_matrix[i][mask_A], qual_matrix[i][mask_C], qual_matrix[i][mask_T], qual_matrix[i][mask_G], qual_matrix[i][mask_N]
            m_A, m_C, m_T, m_G, m_N = c_A.multiply(p_A).sum(), c_C.multiply(p_C).sum(), c_T.multiply(p_T).sum(), c_G.multiply(p_G).sum(), c_N.multiply(p_N).sum()
            total = m_A + m_C + m_T + m_G + m_N
            if total == 0 and seq[i] == "~":
                consensus_list.append(1)
            elif total == 0 and seq[i] != "~": # Check this in the future. If the actual label is a nt, it shouldn't have a total of 0. This is a patch, not a solution.
                consensus_list.append(1)
            else:
                consensus = {
                    "A" : m_A / total,
                    "T" : m_T / total,
                    "C" : m_C / total,
                    "G" : m_G / total,
                    "N" : m_N / total
                }
                label = seq[i]
                if label != "~":
                    consensus_list.append(consensus[seq[i]])
                else:
                    consensus_list.append(1)

        predictors = [str(log10(coverage_list[i]) * consensus_list[i]) for i in range(len(consensus_list))]
        # Insertions in the middle will generate missing bases: add non_informative values for the remaining bases.
        non_informative = str(log10(30) * 1)
        while len(predictors) != 100:
            if len(predictors) > 100:
                print("Error in predictors length.")
                sys.exit()
            predictors.append(non_informative)
        # Add MAPQ as 101th predictor.
        predictors.append(str(mapq))
        # Transform the list into a line so it can be written into a file.
        predictors_line = "\t".join(predictors)
        return predictors_line


    def process_SAM(self):
        # Generate a results folder. 
        counter = 0
        with open(self.sam, "rt") as inp, open(f"""./data/machine_learning/predictors/{self.sam_name}.predictors""", "wt") as out:
            for line_ in inp:
                counter += 1
                # Get rid of the newline character. 
                line_noend = line_[:-1]
                # Obtain each SAM field.
                line = line_noend.split("\t")
                # POS.
                pos = int(line[3])
                # MAPQ.
                mapq = int(line[4])
                # CIGAR.
                cigar = line[5]
                # SEQ.
                seq = line[9]
                modified = False # Keep track if we have trimmed 1 base for 101bp lectures.
                if len(seq) < 100:
                    continue
                # Trim 101bp reads.

                elif len(seq) > 100:
                    seq = seq[:100]
                    modified = True
                # Obtain a CIGAR list for the template read.
                cigar_list = self.process_cigar(cigar)
                # Get rid of reads with soft or hard clips or initial insertions.
                if cigar_list[0][-1] in ["S", "H", "I"]:
                    continue
                # Obtain the modified template sequence according to the cigar.
                t_seq = self.modify_SEQ(seq, cigar_list)
                # Calculate the end position.

                pos_end = pos + len(t_seq) - 1
                # Calculate the range of the pileup: 
                range_ = self.pileup_range - 1
                # Get the corresponding TSV file. 
                file = f"./chunks/{(pos // 1000000) * 1000000}.tsv"
                # Load the file into memory.
                matrix_in = self.load_matrix(pos - range_, pos + range_, file)
                matrix = matrix_in[matrix_in["CIGAR"] != np.nan] # THIS SHOULDN'T HAPPEN EITHER BUT IT HAPPENS IN SOME CASES, CHECK IN THE FUTURE
                print(counter)
                # Add a column with the list of the cigar operations:
                matrix["CIGAR_LIST"] = matrix["CIGAR"].apply(self.process_cigar)
                # Deduct 1 if the sequence has been modified.
                #matrix["CIGAR_MODIFIED"] = matrix["CIGAR_LIST"].apply(lambda x: "".join(x))
                # Trim all 101bp reads into 100bp.
                matrix["SEQ"] = matrix["SEQ"].apply(lambda x: x[:100] if len(x) > 100 else x)
                # Trim all 101bp qual into 100bp.
                matrix["QUAL"] = matrix["QUAL"].apply(lambda x: x[:100] if len(x) > 100 else x)
                # Transform the position into integers.
                matrix["POS"] = matrix["POS"].apply(lambda x: int(x))
                # Generate the pileup for SEQ.
                matrix["PILEUP_SEQ"] = np.vectorize(self.transform_SEQ)(matrix["SEQ"], matrix["CIGAR_LIST"], pos, matrix["POS"], pos_end)
                # Generate the pileup for QUAL.
                matrix["PILEUP_QUAL"] = np.vectorize(self.transform_SEQ)(matrix["QUAL"], matrix["CIGAR_LIST"], pos, matrix["POS"], pos_end)
                # Filter -1 reads.
                filtered_matrix = matrix[(matrix["PILEUP_SEQ"] != "Remove") & (matrix["PILEUP_QUAL"] != "Remove")]
                # Compute the predictors.
                predictors_line = self.calculate_predictors(filtered_matrix, t_seq, mapq)
                # Write the predictors into the output file. 
                out.write(predictors_line + "\n")

def Main():
    # Get the SAM file.
    sam = sys.argv[1]
    GeneratePredictors(sam, pileup_range = 100)

Main()
