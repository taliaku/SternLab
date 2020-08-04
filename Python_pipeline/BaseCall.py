#! python/python-anaconda3.2019.7

import argparse
import os


def count_num_times_read_matches(ReadLines, Lines, mode, counter=0, plus_counter=1, minus_counter=2):
    READ_IDS_VALUES = {}
    RowNum = 0
    while RowNum < Lines:
        read_record = ReadLines[RowNum].strip()
        read_record_split = read_record.split("\t")
        if mode in ["sr", "SR", "SeqtoRef"]:
            read_id = read_record_split[0].strip()
        elif mode in ["ReftoSeq", "RS", "rs"]:
            read_id = read_record_split[1].strip()
        else:
            raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")
        if read_id not in READ_IDS_VALUES:
            READ_IDS_VALUES[read_id] = [0, 0, 0,
                                        " "]  # READ_IDS_VALUES[read_id] = [read_id_counter, plus_counter, minus_counter, read_id_quality_line]
        READ_IDS_VALUES[read_id][counter] += 1
        strand = read_record_split[6].strip()
        if strand == "plus":
            READ_IDS_VALUES[read_id][plus_counter] += 1
        elif strand == "minus":
            READ_IDS_VALUES[read_id][minus_counter] += 1
        RowNum += 1

    return READ_IDS_VALUES


def reads_id_quality(READ_IDS_VALUES, blast_FilePath, quality_line=3):
    qual_FilePath = blast_FilePath.split(".fasta")[0] + ".qual"
    try:
        with open(qual_FilePath, 'rt') as qual_file:
            ReadLines = qual_file.readlines()
    except:
        raise Exception("Cannot open qual file\n")

    if len(ReadLines) % 2 != 0:
        raise Exception("Unexpected error, number of lines in qual file " + qual_FilePath + " does not divide by two\n")

    RowNum = 0
    Lines = len(ReadLines)
    while RowNum < Lines:
        read_id = ReadLines[RowNum].strip()
        if not read_id.startswith("@"):
            raise Exception("Unexpected error, missing @. Cannot identify read id in first line of read id\n")
        read_id = read_id[1:len(read_id)]
        RowNum += 1
        if read_id in READ_IDS_VALUES:
            read_id_quality_line = ReadLines[RowNum].strip()
            READ_IDS_VALUES[read_id][
                quality_line] = read_id_quality_line  # READ_IDS_VALUES[read_id] = [read_id_counter, plus_counter, minus_counter, read_id_quality_line]
        RowNum += 1

    return READ_IDS_VALUES


def create_ref_seq(ref_FilePath):
    try:
        with open(ref_FilePath, 'rt') as ref_file:
            ReadLines = ref_file.readlines()
    except:
        raise Exception("Cannot open ref file " + ref_FilePath + "\n")

    Lines = len(ReadLines)
    if Lines < 2:
        raise Exception("Unexpected error, empty or missing lines in ref file " + ref_FilePath + "\n")

    if ReadLines[0].startswith(">"):
        ref_genome = " "
        for i in range(1, Lines):
            if i == 1:
                ref_genome = ReadLines[i].strip().upper()
            else:
                ref_genome += ReadLines[i].strip().upper()
    else:
        raise Exception("first line in ref fasta file " + ref_FilePath + " does not start with >\n")

    REF_GENOME = {}
    ref_genome_length = len(ref_genome)
    for i in range(ref_genome_length):
        if ref_genome[i] in ['A', 'C', 'G', 'T']:
            REF_GENOME[i + 1] = ref_genome[i]
        else:
            raise Exception(
                "Found a non valid DNA letter [{}] in position [{}] of the reference genome".format(ref_genome[i],
                                                                                                    i + 1))

    return REF_GENOME


def get_read_id_qscore(blast_FilePath, read_id, Quality_line, ASCII_Q_SCORE):
    Quality_line_length = len(Quality_line)
    READ_ID_QSCORE = {}
    for read_counter in range(Quality_line_length):
        ASCII_sign = Quality_line[read_counter]
        q_score = ASCII_Q_SCORE.get(ASCII_sign)
        try:
            q_score = int(q_score)
            READ_ID_QSCORE[read_counter + 1] = int(q_score)
        except:
            raise Exception(
                "Unexpected error. Q_score value of ASCII sign " + ASCII_sign + " in read is " + read_id + " in " + blast_FilePath + " is not an integer\n")

    return READ_ID_QSCORE


def get_qscore_dict(pipeline_dir):
    ASCII_FilePath = pipeline_dir + "/ascii_table_processed.txt"
    try:
        with open(ASCII_FilePath, 'rt') as ASCII_file:
            ReadLines = ASCII_file.readlines()
    except:
        raise Exception("Cannot open ascii txt file " + ASCII_FilePath + "\n")

    ASCII_Q_SCORE = {}
    RowNum = 0
    Lines = len(ReadLines)
    while RowNum < Lines:
        Line = ReadLines[RowNum].strip()
        Split_line = Line.split(" ")
        ASCII_sign = Split_line[2].strip()
        try:
            q_score = int(Split_line[0].strip())
        except:
            raise Exception(
                "Unexpected error. Cannot convert ASCII sign to integer q-score with file " + ASCII_FilePath + "\n")
        ASCII_Q_SCORE[ASCII_sign] = q_score
        RowNum += 1

    return ASCII_Q_SCORE


def double_position_counter(read_counter, READ_ID_DOUBLE_POSITION_COUNTER, ref_counter, counter=0, mapping_positions=1):
    if read_counter not in READ_ID_DOUBLE_POSITION_COUNTER:
        READ_ID_DOUBLE_POSITION_COUNTER[read_counter] = [0,
                                                         []]  # READ_ID_DOUBLE_POSITION_COUNTER[read_counter] = [read_position_counter,[mapping_positions]]
    READ_ID_DOUBLE_POSITION_COUNTER[read_counter][counter] += 1
    READ_ID_DOUBLE_POSITION_COUNTER[read_counter][mapping_positions].append(ref_counter)

    return READ_ID_DOUBLE_POSITION_COUNTER


def create_read_id_seq_qscore(read_id, read_start, read_end, ref_start, ref_end, strand, aln, READ_ID_QSCORE,
                              REF_GENOME, \
                              READ_ID_BASE_CALL_COUNTER, READ_ID_DOUBLE_POSITION_COUNTER, q_score, mode, counter=0,
                              q_score_sum=1, mutation_positions=2):
    if mode in ["sr", "SR", "SeqtoRef"]:
        mode = "SeqtoRef"
        Seq_pair_position = 0
        Ref_pair_position = 1
        read_counter_direction = 1
        if strand == "minus":
            ref_counter_direction = -1
        else:  # strand == "plus"
            ref_counter_direction = 1

    elif mode in ["ReftoSeq", "RS", "rs"]:
        mode = "ReftoSeq"
        Seq_pair_position = 1
        Ref_pair_position = 0
        ref_counter_direction = 1
        if strand == "minus":
            read_counter_direction = -1
        else:  # strand == "plus":
            read_counter_direction = 1

    else:
        raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

    REVERSE_COMPLEMENT_BASES = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    try:
        ref_counter = int(ref_start)
        read_counter = int(read_start)
    except:
        raise Exception("Unexpected error. ref or read start position is not an integer value\n")

    aln_length = len(aln)

    i = 0
    while i < aln_length:
        match = " "
        while aln[i].isdigit():
            if match == " ":
                match = aln[i]
            else:
                match += aln[i]
            i += 1
            if i == aln_length:
                break
        match = int(match)
        for j in range(match):
            Base = REF_GENOME.get(ref_counter)
            if ref_counter not in READ_ID_BASE_CALL_COUNTER:
                READ_ID_BASE_CALL_COUNTER[ref_counter] = {}
            if Base not in READ_ID_BASE_CALL_COUNTER[ref_counter]:
                READ_ID_BASE_CALL_COUNTER[ref_counter][Base] = [0, 0,
                                                                []]  # READ_ID_BASE_CALL_COUNTER[ref_counter][Base] = [counter, q_score_sum, [mutation_positions]]
            READ_ID_BASE_CALL_COUNTER[ref_counter][Base][counter] += 1
            READ_ID_BASE_CALL_COUNTER[ref_counter][Base][q_score_sum] += READ_ID_QSCORE.get(read_counter)
            READ_ID_DOUBLE_POSITION_COUNTER = double_position_counter(read_counter, READ_ID_DOUBLE_POSITION_COUNTER,
                                                                      ref_counter)
            ref_counter += ref_counter_direction
            read_counter += read_counter_direction

        if i < aln_length:
            mutation = " "
            while not aln[i].isdigit():
                if mutation == " ":
                    mutation = aln[i]
                else:
                    mutation += aln[i]
                i += 1  # if i = len(aln) sys.exit("Unexpected error in aln match in read id " + read_id + "\n")
            mutation_length = len(mutation)
            if mutation_length % 2 != 0:
                raise Exception("Unexpected error, mutation seq " + mutation + " does not divide by two\n")
            k = 0
            gap = 0
            deletion = 0
            while k < mutation_length:
                Base = mutation[k + Seq_pair_position]
                if (mode == "SeqtoRef") and (strand == "minus") and (Base != "-") and (
                        Base != "N"):  # NextSeq sometimes gives N base calling. Pipeline ignores N base later in calculate_read_id_contribution.
                    Base = REVERSE_COMPLEMENT_BASES.get(Base)

                if mutation[k + Ref_pair_position] in ['A', 'C', 'G', 'T']:  # mutation or deletion, NOT insertion
                    if ref_counter not in READ_ID_BASE_CALL_COUNTER:
                        READ_ID_BASE_CALL_COUNTER[ref_counter] = {}
                    if Base not in READ_ID_BASE_CALL_COUNTER[ref_counter]:
                        READ_ID_BASE_CALL_COUNTER[ref_counter][Base] = [0, 0,
                                                                        []]  # READ_ID_BASE_CALL_COUNTER[ref_counter][Base] = [counter, q_score_sum, [mutation_positions]]
                    READ_ID_BASE_CALL_COUNTER[ref_counter][Base][counter] += 1
                    if mutation[k + Seq_pair_position] == "-":  # deletion
                        if (mode == "ReftoSeq") and (strand == "minus"):
                            y = k
                            deletion = 0
                            while mutation[y + Seq_pair_position] == "-":
                                deletion += 1
                                y += 2
                                if y == len(mutation):
                                    break
                        else:
                            deletion += 1
                        if (mode == "ReftoSeq") and (strand == "minus"):
                            read_deletion_counter = read_counter
                        else:  # ((mode == "ReftoSeq") and (strand == "plus")) or (mode == "SeqtoRef")
                            read_deletion_counter = read_counter - 1
                        deletion_counter = str(read_deletion_counter) + 'X' + str(deletion)
                        READ_ID_BASE_CALL_COUNTER[ref_counter][Base][
                            q_score_sum] += q_score  # assigning min q-score for deletions
                        READ_ID_BASE_CALL_COUNTER[ref_counter][Base][mutation_positions].append(deletion_counter)
                        READ_ID_DOUBLE_POSITION_COUNTER = double_position_counter(deletion_counter,
                                                                                  READ_ID_DOUBLE_POSITION_COUNTER,
                                                                                  ref_counter)
                    else:  # mutation
                        READ_ID_BASE_CALL_COUNTER[ref_counter][Base][q_score_sum] += READ_ID_QSCORE.get(read_counter)
                        READ_ID_BASE_CALL_COUNTER[ref_counter][Base][mutation_positions].append(read_counter)
                        READ_ID_DOUBLE_POSITION_COUNTER = double_position_counter(read_counter,
                                                                                  READ_ID_DOUBLE_POSITION_COUNTER,
                                                                                  ref_counter)
                        read_counter += read_counter_direction
                    ref_counter += ref_counter_direction
                    gap = 0

                elif mutation[k + Ref_pair_position] == "-":  # insertion
                    if (mode == "SeqtoRef") and (strand == "minus"):
                        x = k
                        gap = 0
                        while mutation[x + Ref_pair_position] == "-":
                            gap += 1
                            x += 2
                            if x == len(mutation):
                                break
                    else:
                        gap += 1
                    if ((mode == "SeqtoRef") and (strand == "plus")) or (mode == "ReftoSeq"):
                        ref_gap_counter = ref_counter - 1
                    else:  # (mode == "SeqtoRef") and (strand == "minus")
                        ref_gap_counter = ref_counter
                    gap_counter = str(ref_gap_counter) + "." + str(gap)
                    if gap_counter not in READ_ID_BASE_CALL_COUNTER:
                        READ_ID_BASE_CALL_COUNTER[gap_counter] = {}
                    if Base not in READ_ID_BASE_CALL_COUNTER[gap_counter]:
                        READ_ID_BASE_CALL_COUNTER[gap_counter][Base] = [0, 0,
                                                                        []]  # READ_ID_BASE_CALL_COUNTER[gap_counter][Base] = [counter, q_score_sum, [mutation_positions]]
                    READ_ID_BASE_CALL_COUNTER[gap_counter][Base][counter] += 1
                    READ_ID_BASE_CALL_COUNTER[gap_counter][Base][q_score_sum] += READ_ID_QSCORE.get(read_counter)
                    READ_ID_BASE_CALL_COUNTER[gap_counter][Base][mutation_positions].append(read_counter)
                    READ_ID_DOUBLE_POSITION_COUNTER = double_position_counter(read_counter,
                                                                              READ_ID_DOUBLE_POSITION_COUNTER,
                                                                              ref_counter)
                    read_counter += read_counter_direction
                    deletion = 0

                else:
                    raise Exception("Unexpected error, base in reference genome is not A,C,G or T\n")

                k += 2

    try:
        ref_end = int(ref_end)
        read_end = int(read_end)
    except:
        raise Exception("Unexpected error. ref or read end position is not an integer value\n")

    if ref_counter - ref_counter_direction != ref_end or read_counter - read_counter_direction != read_end:
        raise Exception("Unexpected error. end of read_id base calling does not match blast results " + read_id + "\n")

    return READ_ID_BASE_CALL_COUNTER, READ_ID_DOUBLE_POSITION_COUNTER


def create_supporting_data_files(blast_FilePath):
    root_FilePath = os.path.splitext(blast_FilePath)[0]
    root_FilePath = root_FilePath.replace(".fasta", "")
    good_reads_file = root_FilePath + ".good_reads"
    good_mutations_file = root_FilePath + ".good_mutations"
    stats_file = root_FilePath + ".freqs.stats"
    Non_contributing_file = root_FilePath + ".NonContributing"
    frequencies_file = root_FilePath + ".freqs"

    return good_reads_file, good_mutations_file, stats_file, frequencies_file, Non_contributing_file


def remove_multiple_mapping(READ_ID_DOUBLE_POSITION_COUNTER, READ_ID_BASE_CALL_COUNTER, double_mapping_counter,
                            counter=0, mapping_positions=1):
    print("remove_multiple_mapping!!!")
    for read_position in READ_ID_DOUBLE_POSITION_COUNTER:
        count = READ_ID_DOUBLE_POSITION_COUNTER[read_position][counter]
        if count > 1:
            double_mapping_counter += 1
            double_mapping_positions_list = READ_ID_DOUBLE_POSITION_COUNTER[read_position][mapping_positions]
            for ref_position in double_mapping_positions_list:
                if ref_position in READ_ID_BASE_CALL_COUNTER:
                    del READ_ID_BASE_CALL_COUNTER[ref_position]

    return READ_ID_DOUBLE_POSITION_COUNTER, READ_ID_BASE_CALL_COUNTER, double_mapping_counter


def calculate_read_id_contribution(READ_ID_BASE_CALL_COUNTER, TOTAL_BASE_CALL_COUNTER, q_score, num_of_repeats,
                                   REF_GENOME, good_mutations_file, \
                                   Non_contributing_file, read_id, Total_base_counter_per_read_id, counter=0,
                                   q_score_sum=1, mutation_positions=2):
    BASE_POSITION = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 5}
    with open(good_mutations_file, 'at') as good_reads_mutations:
        with open(Non_contributing_file, 'at') as bad_reads_mutations:
            for ref_position in READ_ID_BASE_CALL_COUNTER:
                for base in READ_ID_BASE_CALL_COUNTER[ref_position]:
                    if base != 'N':  # NextSeq sometimes gives N base calling. Pipeline ignores N base.
                        read_id_base_call_counter = READ_ID_BASE_CALL_COUNTER[ref_position][base][counter]
                        read_id_base_call_qscore = READ_ID_BASE_CALL_COUNTER[ref_position][base][q_score_sum]
                        average_q_score = read_id_base_call_qscore / read_id_base_call_counter

                        if (average_q_score >= q_score) and (read_id_base_call_counter >= num_of_repeats) and \
                                ((num_of_repeats > 1) or (
                                        num_of_repeats == 1 and len(READ_ID_BASE_CALL_COUNTER[ref_position]) == 1)):

                            if ref_position not in TOTAL_BASE_CALL_COUNTER:
                                TOTAL_BASE_CALL_COUNTER[ref_position] = [0, 0, 0, 0, 0,
                                                                         0]  # TOTAL_BASE_CALL_COUNTER[ref_position] = [Total_counter, A_counter, C_counter, G_counter, T_counter, -_counter]
                            Base_position = BASE_POSITION.get(base)
                            TOTAL_BASE_CALL_COUNTER[ref_position][Base_position] += 1
                            TOTAL_BASE_CALL_COUNTER[ref_position][counter] += 1
                            Total_base_counter_per_read_id += 1
                            if len(READ_ID_BASE_CALL_COUNTER[ref_position][base][mutation_positions]) > 0:
                                mutation_positions_list = sorted(
                                    READ_ID_BASE_CALL_COUNTER[ref_position][base][mutation_positions])
                                try:
                                    good_reads_mutations.write(
                                        str(ref_position) + "\t" + read_id + "\t" + base + "\t" + str(
                                            mutation_positions_list) + "\n")
                                except:
                                    raise Exception(
                                        "Unexpected error, cannot write into file " + good_mutations_file + "\n")

                        else:
                            if ref_position in REF_GENOME:
                                ref_base = REF_GENOME.get(ref_position)
                            else:
                                ref_base = "-"  # insertion
                            try:
                                bad_reads_mutations.write("Base Calling or mutation at ref position " + str(
                                    ref_position) + " " + ref_base + ">" + base + " with " + \
                                                          str(read_id_base_call_counter) + " repeat/s and Q " + str(
                                    average_q_score) + " was removed from " + read_id + "\n")
                            except:
                                raise Exception(
                                    "Unexpected error, cannot write into file " + Non_contributing_file + "\n")

    return READ_ID_BASE_CALL_COUNTER, TOTAL_BASE_CALL_COUNTER, Total_base_counter_per_read_id


def summarize_total_contribution_statistics_1(CONTRIBUTING_STATS, Contributing_reads_counter,
                                              Contributing_bases_counter, Non_contributing_reads_counter, \
                                              Total_base_counter_per_read_id, match_counter, good_reads_file, read_id,
                                              Quality_line, repeat_read_counter=0, repeat_base_counter=1):
    if Total_base_counter_per_read_id > 0:
        with open(good_reads_file, 'at') as good_reads:
            try:
                good_reads.write("@" + read_id + "\n")
            except:
                raise Exception("Unexpected error, cannot write into file " + good_reads_file + "\n")

        Contributing_reads_counter += 1
        Contributing_bases_counter += Total_base_counter_per_read_id

        if match_counter not in CONTRIBUTING_STATS:
            CONTRIBUTING_STATS[match_counter] = [0,
                                                 0]  # CONTRIBUTING_STATS[match_counter] = [repeat_read_counter, repeat_base_counter]

        CONTRIBUTING_STATS[match_counter][repeat_read_counter] += 1
        CONTRIBUTING_STATS[match_counter][repeat_base_counter] += Total_base_counter_per_read_id

    else:
        Non_contributing_reads_counter += 1  # Counts any read id with zero contribution

    return CONTRIBUTING_STATS, Contributing_reads_counter, Contributing_bases_counter, Non_contributing_reads_counter


def calculate_frequencies(TOTAL_BASE_CALL_COUNTER, REF_GENOME, frequencies_file, counter=0):
    # TOTAL_BASE_CALL_COUNTER[ref_position] = [Total_counter, A_counter, C_counter, G_counter, T_counter, -_counter]
    BASE_POSITION = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: '-'}
    with open(frequencies_file, 'wt') as freqs:
        for ref_position in TOTAL_BASE_CALL_COUNTER:
            for Base_position in [1, 2, 3, 4, 5]:
                base = BASE_POSITION.get(Base_position)
                base_counter = TOTAL_BASE_CALL_COUNTER[ref_position][Base_position]
                total_counter = TOTAL_BASE_CALL_COUNTER[ref_position][counter]
                frequency = round(base_counter / total_counter, 6)
                base_counter = TOTAL_BASE_CALL_COUNTER[ref_position][Base_position]
                if ref_position in REF_GENOME:
                    ref_base = REF_GENOME.get(ref_position)  # mutation or deletion
                else:
                    ref_base = "-"  # insertion
                coverage = TOTAL_BASE_CALL_COUNTER[ref_position][counter]
                try:
                    freqs.write(
                        str(ref_position) + "\t" + ref_base + "\t" + base + "\t" + str(base_counter) + "\t" + str(
                            coverage) + "\t" + str(frequency) + "\n")
                except:
                    raise Exception("Unexpected error, cannot write into file " + frequencies_file + "\n")

    os.system("sort -V " + frequencies_file + " -o " + frequencies_file + "\n")
    os.system("sed -i $'1i ref_position\tref_base\tbase\tbase_counter\tcoverage\tfrequency\n' " + frequencies_file)


def summarize_total_contribution_statistics_2(CONTRIBUTING_STATS, stats_file, double_mapping_counter,
                                              Contributing_reads_counter, Contributing_bases_counter, \
                                              Non_contributing_reads_counter, MATCH_STATISTICS, repeat_read_counter=0,
                                              repeat_base_counter=1):
    with open(stats_file, 'at') as stats:
        try:
            stats.write(
                "\nNumber of reads contributing to frequency counts:\t" + str(Contributing_reads_counter) + "\n")
            stats.write(
                "\nNumber of bases contributing to frequency counts:\t" + str(Contributing_bases_counter) + "\n")
        except:
            raise Exception("Unexpected error, cannot write into file " + stats_file + "\n")

        try:
            stats.write("\nContribution of reads and bases by repeat:\n")
            for repeat in sorted(CONTRIBUTING_STATS.keys()):
                reads_counter = CONTRIBUTING_STATS[repeat][repeat_read_counter]
                bases_counter = CONTRIBUTING_STATS[repeat][repeat_base_counter]
                stats.write(str(repeat) + " repeats, " + str(reads_counter) + " reads called and " + str(
                    bases_counter) + " bases called\n")
            del CONTRIBUTING_STATS
        except:
            raise Exception("Unexpected error, cannot write into file " + stats_file + "\n")

        try:
            stats.write("\nNumber of bases removed due to double mapping:\t" + str(double_mapping_counter) + "\n")
            stats.write("\nNumber of non contributing reads:\t" + str(Non_contributing_reads_counter) + "\n")
        except:
            raise Exception("Unexpected error, cannot write into file " + stats_file + "\n")

        try:
            stats.write("\nDistribution of matches by blast:\n")
            for match_counter in sorted(MATCH_STATISTICS.keys()):
                stats.write(str(match_counter) + "\t" + str(MATCH_STATISTICS[match_counter]) + "\n")
        except:
            raise Exception("Unexpected error, cannot write into file " + stats_file + "\n")


def get_record(ReadLines, RowNum, mode):
    read_record = ReadLines[RowNum].strip()
    read_record_split = read_record.split("\t")
    if len(read_record_split) != 9:
        raise Exception("Unexpected error. Read record in blast file does not have 9 arguments\n")
    if mode in ["sr", "SR", "SeqtoRef"]:
        read_id = read_record_split[0].strip()
    elif mode in ["ReftoSeq", "RS", "rs"]:
        read_id = read_record_split[1].strip()
    else:
        raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

    return read_record_split, read_id


def BaseCall(pipeline_dir, blast_FilePath, ref_FilePath, num_of_repeats, q_score, mode, Protocol, counter=0,
             plus_counter=1, minus_counter=2, quality_line=3, please_remove_multiple_mapping=True):
    try:
        with open(blast_FilePath, 'rt') as read_records:
            blast_lines = read_records.readlines()
    except:
        raise Exception("Cannot open/read file " + blast_FilePath + "\n")

    # initialize parameters
    total_base_call_counter = {}
    contributing_stats = {}
    match_statistics = {}
    contributing_reads_counter = 0
    contributing_bases_counter = 0
    double_mapping_counter = 0
    non_contributing_reads_counter = 0
    good_reads_file, good_mutations_file, stats_file, frequencies_file, Non_contributing_file = create_supporting_data_files(
        blast_FilePath)

    # load some stuff to memory
    ref_genome = create_ref_seq(ref_FilePath)
    ascii_q_score = get_qscore_dict(pipeline_dir)
    read_ids_values = count_num_times_read_matches(blast_lines, len(blast_lines), mode)
    read_ids_values = reads_id_quality(read_ids_values, blast_FilePath)

    blast_line_counter = 0
    while blast_line_counter < len(blast_lines):
        read_record_split, read_id = get_record(blast_lines, blast_line_counter, mode)
        number_of_matches = read_ids_values[read_id][counter]
        number_of_plus_matches = read_ids_values[read_id][plus_counter]
        number_of_minus_matches = read_ids_values[read_id][minus_counter]
        total_base_counter_per_read_id = 0
        match_counter = 1
        if number_of_matches not in match_statistics:
            match_statistics[number_of_matches] = 0
        match_statistics[number_of_matches] += 1
        always_enter_this_if = True if please_remove_multiple_mapping == "N" else False  # TODO: remove this hack when done.
        if always_enter_this_if or (number_of_matches >= num_of_repeats) and (
                (Protocol in ["C", "c", "circular"]) or (
                Protocol in ["L", "l", "linear"] and num_of_repeats == 2 and number_of_plus_matches == 1 and
                number_of_minus_matches == 1
                )
                or (Protocol in ["L", "l", "linear"] and num_of_repeats == 1 and number_of_plus_matches < 2 and
                    number_of_minus_matches < 2)
        ):
            if read_id in read_ids_values:
                Quality_line = read_ids_values[read_id][quality_line]
            else:
                raise Exception("read id " + read_id + " was not found in qual file\n")

            read_id_qscore = get_read_id_qscore(blast_FilePath, read_id, Quality_line, ascii_q_score)
            # For each read_id in blast file, create a dictionary READ_ID_BASE_CALL_COUNTER that will get ref position
            # base count and q-score including mutations and / or gaps.
            read_id_base_call_counter = {}
            # For each read_id, create a dictionary READ_ID_DOUBLE_POSITION_COUNTER that will get how many times each
            # position in the read was counted. Multiple mapped positions will be removed from base calling.
            read_id_double_position_counter = {}

            while match_counter <= number_of_matches:
                if mode in ["sr", "SR", "SeqtoRef"]:
                    read_start = read_record_split[2].strip()
                    read_end = read_record_split[3].strip()
                    ref_start = read_record_split[4].strip()
                    ref_end = read_record_split[5].strip()
                elif mode in ["ReftoSeq", "RS", "rs"]:
                    ref_start = read_record_split[2].strip()
                    ref_end = read_record_split[3].strip()
                    read_start = read_record_split[4].strip()
                    read_end = read_record_split[5].strip()
                else:
                    raise Exception(
                        "Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

                strand = read_record_split[6].strip()
                # length = read_record_split[7].strip()
                aln = read_record_split[8].strip()

                read_id_base_call_counter, read_id_double_position_counter = \
                    create_read_id_seq_qscore(read_id, read_start, read_end, ref_start, ref_end, strand, aln,
                                              read_id_qscore, \
                                              ref_genome, read_id_base_call_counter, read_id_double_position_counter,
                                              q_score, mode)

                match_counter += 1
                blast_line_counter += 1
                if match_counter <= number_of_matches:  # and blast_line_counter < Lines
                    read_record_split, next_read_id = get_record(blast_lines, blast_line_counter, mode)
                    if next_read_id != read_id:
                        raise Exception(
                            "Unexpected error, next_read_id " + next_read_id + " does not match read_id " + read_id + "\n")

            # For each read_id remove positions that were mapped more than once from contributing base calls
            # if please_remove_multiple_mapping == "Y":  #TODO: this doesn't seem to do anything....
            read_id_double_position_counter, read_id_base_call_counter, double_mapping_counter = \
                remove_multiple_mapping(read_id_double_position_counter, read_id_base_call_counter,
                                        double_mapping_counter)
            del read_id_double_position_counter

            # For each read_id in the blast file calculate contribution based on q-score in READ_ID_BASE_CALL_COUNTER.
            # Summarize results in total_base_call_counter.
            read_id_base_call_counter, total_base_call_counter, total_base_counter_per_read_id = \
                calculate_read_id_contribution(read_id_base_call_counter, total_base_call_counter, q_score,
                                               num_of_repeats, ref_genome, \
                                               good_mutations_file, Non_contributing_file, read_id,
                                               total_base_counter_per_read_id)
            del read_id_base_call_counter

            # Create contributing statistics files, part1
            contributing_stats, contributing_reads_counter, contributing_bases_counter, non_contributing_reads_counter = \
                summarize_total_contribution_statistics_1(contributing_stats, contributing_reads_counter,
                                                          contributing_bases_counter, non_contributing_reads_counter,
                                                          total_base_counter_per_read_id, match_counter - 1, \
                                                          good_reads_file, read_id, Quality_line)

        else:
            with open(Non_contributing_file, 'at') as bad_reads_mutations:  # Record any non contributing reads
                while match_counter <= number_of_matches:
                    try:
                        bad_reads_mutations.write("\t".join(read_record_split) + "\n")
                    except:
                        raise Exception("Unexpected error, cannot write into file " + Non_contributing_file + "\n")
                    match_counter += 1
                    blast_line_counter += 1
                    if match_counter <= number_of_matches:
                        read_record_split, next_read_id = get_record(blast_lines, blast_line_counter, mode)
                        if next_read_id != read_id:
                            raise Exception(
                                "Unexpected error, next_read_id " + next_read_id + " does not match read_id " + read_id + "\n")

            non_contributing_reads_counter += 1  # Counts any non contributing read

    del ascii_q_score

    os.system("sort " + Non_contributing_file + " -o " + Non_contributing_file + "\n")

    # Calculate frequencies and write to file
    calculate_frequencies(total_base_call_counter, ref_genome, frequencies_file)

    del total_base_call_counter
    del ref_genome

    # Create contributing statistics files, part2
    summarize_total_contribution_statistics_2(contributing_stats, stats_file, double_mapping_counter,
                                              contributing_reads_counter, contributing_bases_counter,
                                              non_contributing_reads_counter, match_statistics)


def main(args):
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))

    blast_FilePath = args.file_path
    if not (os.path.isfile(blast_FilePath) and os.path.splitext(blast_FilePath)[1] == '.blast'):
        raise Exception(
            "Unexpected error, " + blast_FilePath + " does not exist, is not a file or or is not a blast file\n")

    ref_FilePath = args.ref
    if not (os.path.isfile(ref_FilePath) and os.path.splitext(ref_FilePath)[1] == '.fasta'):
        raise Exception("Unexpected error, " + ref_FilePath + " does not exist, is not a file or is not a fasta file\n")

    min_num_repeats = args.repeats
    if min_num_repeats < 1:
        raise Exception("Unexpected error, min number of repeats is less than 1\n")
    if min_num_repeats < 2:
        print("\nWarning. Running base calling with min number of repeats less than 2\n")
    if min_num_repeats > 2:
        print("\nWarning. Running pipeline with min number of repeats bigger than 2\n")

    q_score = args.q_score
    if q_score < 0 or q_score > 40:
        raise Exception("Unexpected error, q-score value " + str(q_score) + " is not valid\n")
    if q_score < 16:
        print("Warning, running pipeline with q-score value of " + str(q_score) + "\n")

    mode = args.blast_mode
    if mode in ["ReftoSeq", "RS", "rs"]:
        mode = "ReftoSeq"
    elif mode in ["sr", "SR", "SeqtoRef"]:
        mode = "SeqtoRef"
    else:
        raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

    Protocol = args.protocol
    if Protocol == None:
        Protocol = "linear"
    else:
        if Protocol not in ["L", "l", "linear", "C", "c", "circular"]:
            raise Exception(
                "Unexpected error, for linear library prep protocol type 'linear' or 'L', for circular library prep protocol type 'circular' or 'C'\n")
    please_remove_multiple_mapping = args.please_remove_multiple_mapping
    BaseCall(pipeline_dir, blast_FilePath, ref_FilePath, min_num_repeats, q_score, mode, Protocol,
             please_remove_multiple_mapping=please_remove_multiple_mapping)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path", type=str, help="input blast file", required=True)
    parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference seq file (fasta)", required=True)
    parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
    parser.add_argument("-x", "--repeats", type=int, help="number of repeats, default=2", required=False, default=2)
    parser.add_argument("-m", "--blast_mode", type=str,
                        help="mode for blast, for Seq to Ref blast type SR, sr or SeqtoRef, for Ref to Seq blast type RS, rs or ReftoSeq, default = 'SeqtoRef'",
                        required=False, default="SeqtoRef")
    parser.add_argument("-p", "--protocol", type=str,
                        help="Library prep protocol is linear = 'L', 'l' or 'linear', or circular = 'C', 'c' or 'circular'. Default='linear'",
                        required=False, default="linear")
    parser.add_argument("-pr", "--please_remove_multiple_mapping", default=True)
    args = parser.parse_args()
    main(args)
