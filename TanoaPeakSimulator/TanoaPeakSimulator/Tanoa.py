import secrets
import random


def main():
    location_length = 1000
    bp_length = 100
    bp_amount = 200

    genome_file = open("dmel-all-chromosome-r6.38.fasta", "r")
    sequences = genome_file.read()
    sequences_arr = sequences.split('>')
    sequence_object_array = []
    for sequence in sequences_arr:
        sequence_split = sequence.splitlines()
        if len(sequence_split) != 0:
            header = sequence_split[0]
            body = ""
            for x in range(1, len(sequence_split)):
                body += sequence_split[x]
            sequence_object_array.append(Chromosome(header, body))

    # for chromosome in sequence_object_array:
    #     for x in range(2000):
    #         get_random_base_pair_string(location_length, bp_length, chromosome.sequence,
    #                                     chromosome.bp_segment_arr)
    # for chromosome in sequence_object_array:
    #     chromosome.small_bp_segment_list = get_random_base_pair_string_list(chromosome.bp_segment_arr, bp_length,
    #                                                                         bp_amount, chromosome.sequence)
    chromosome = sequence_object_array[1]
    for x in range(100):
        get_random_base_pair_string(location_length, bp_length, chromosome.sequence, chromosome.bp_segment_arr)

    chromosome.small_bp_segment_list = get_random_base_pair_string_list(chromosome.bp_segment_arr, bp_length,
                                                                        bp_amount, chromosome.sequence)

    generate_output_files(chromosome)


def get_random_base_pair_string(sequence_length, bp_length, sequence, array_of_segments):
    if len(sequence) > sequence_length + 101:  # Ignore sequences that don't have enough length to work...
        random_start_location = random.randrange(99, len(sequence) - (sequence_length + 1))

        mini_sequence = sequence[random_start_location: (random_start_location + sequence_length)]
        temp_segment = BasePairSegment(random_start_location, mini_sequence, sequence_length, bp_length)

        array_of_segments.append(temp_segment)


def get_random_base_pair_string_list(list_of_bp_segments, bp_length, bp_amount, original_sequence):
    temp_bp_list = []
    for x in list_of_bp_segments:
        for y in range(bp_amount):
            random_start_location = random.randrange(x.overlap_start, x.end)
            temp_bp_list.append(original_sequence[random_start_location: random_start_location + bp_length])
    return temp_bp_list


def generate_output_files(chromosome):
    output_to_file = ""
    for x in chromosome.bp_segment_arr:
        output_to_file += chromosome.chromosome_name + ":" + str(x.start_location) + "-" + str(x.end) + "\n"
    output_file = open("sample_regions.txt", 'w')
    output_file.write(output_to_file)
    output_file.close()

    random_regions = ""
    for x in chromosome.small_bp_segment_list:
        random_regions += x + "\n"
    samples_file = open("small_base_pair_samples.txt", 'w')
    samples_file.write(random_regions)
    samples_file.close()


class BasePairSegment:
    def __init__(self, start_int, sequence_string, sequence_length, bp_length):
        self.start_location = start_int
        self.sequence = sequence_string
        self.end = start_int + (sequence_length - 1)
        self.overlap_start = self.start_location - (bp_length - 1)


class Chromosome:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.bp_segment_arr = []
        self.small_bp_segment_list = []
        self.chromosome_name = header.split()[0]


main()
