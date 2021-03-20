import secrets
import random


def main():
    base_pair_strings_with_locations = []
    location_length = 1000
    kmer_length = 100

    genome_file = open("dmel-all-chromosome-r6.38.fasta", "r")
    header = genome_file.readline()
    sequence_string = genome_file.read().replace("\n", "")
    print(header)
    print(sequence_string)

    for x in range(2000):
        get_random_base_pair_string(location_length, kmer_length, sequence_string, base_pair_strings_with_locations)

    kmer_list = get_random_base_pair_string_kmer_list(base_pair_strings_with_locations, kmer_length, sequence_string)


def get_random_base_pair_string(sequence_length, kmer_length, sequence, array_of_segments):
    random_start_location = random.randrange(99, len(sequence) - (sequence_length + 1))

    mini_sequence = sequence[random_start_location: (random_start_location + sequence_length)]
    temp_segment = BasePairSegment(random_start_location, mini_sequence, sequence_length, kmer_length)

    array_of_segments.append(temp_segment)


def get_random_base_pair_string_kmer_list(list_of_bp_segments, kmer_length, original_sequence):
    temp_kmer_list = []
    for x in list_of_bp_segments:
        for y in range(kmer_length):
            random_start_location = random.randrange(x.overlap_start, x.end)
            temp_kmer_list.append(original_sequence[random_start_location: random_start_location + kmer_length])
    return temp_kmer_list


class BasePairSegment:
    def __init__(self, start_int, sequence_string, sequence_length, kmer_length):
        self.start_location = start_int
        self.sequence = sequence_string
        self.end = start_int + (sequence_length - 1)
        self.overlap_start = self.start_location - (kmer_length - 1)

class Chromosome:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

main()
