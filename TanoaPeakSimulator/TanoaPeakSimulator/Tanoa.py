import secrets
import random


def main():
    location_length = 1000
    bp_length = 100
    bp_amount = 200
    # This should be in a function with the file as a parameter for when this gets turned into
    genome_file = open("dmel-all-chromosome-r6.38.fasta", "r")  # code which takes arguments
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

    # Code for running all of the chromosomes. Commented out for now since I don't have enough memory...
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


# Gets a psudo-random string of a given amount of characters from the main sequence. This will move to taking
# the buffers on the sides and saving that as well so when sampled later we don't have to return to the
# Very large string. That's inefficient.
# Parmeters: sequence_length: Length of the big sample that you want to take
# bp_length: the length of a smaller base pair that you will sample from the broader region
#   this is setup so that the buffer can be dynamic based on the size of the reads you want to match
#   to the original chromosome are
# Sequence: The actual character sequence, generally attached to a chromosome like 2L or 2R
# array_of_segments: array of BasePairSegments (samples) from the original chromosome.
#   This array is stored within a chromosome as bp_sequence_array and is modifying that object.
def get_random_base_pair_string(sequence_length, bp_length, sequence, array_of_segments):
    if len(sequence) > sequence_length + 101:  # Ignore sequences that don't have enough length to work...
        random_start_location = random.randrange(99, len(sequence) - (sequence_length + 1))

        # gets a random sample from the random start location to the sequence length. This eventually should be
        # random_start_location + sequence_length and + buffer... This is a simple substring
        mini_sequence = sequence[random_start_location: (random_start_location + sequence_length)]
        # makes a BasePairSegment object from the mini_sequence string and associated data
        temp_segment = BasePairSegment(random_start_location, mini_sequence, sequence_length, bp_length)
        # stores that in the passed array from chromosome
        array_of_segments.append(temp_segment)


# Gets smaller reads from larger samples and returns them as a list
# Parameters: list_of_bp_segments: the larger samples, stored as an array in a chromosome object
# bp_length: the length of one single read. Generally between 50-100bps in length
# bp_amount: how many of these smaller reads you want to generate per larger sample
# original_sequence: for now, this is the whole original sequence. This parameter will be replaced
#   by the chromosome and the chromosome will be looped through here, with it storing the sample with
#   it's bp_length buffers so that we can randomly sample that string. I think that will improve
#   performance significantly
def get_random_base_pair_string_list(list_of_bp_segments, bp_length, bp_amount, original_sequence):
    temp_bp_list = []
    for x in list_of_bp_segments: # loops through the list of BasePairSegment objects
        for y in range(bp_amount):
            random_start_location = random.randrange(x.overlap_start, x.end)  # gets a random location to get a read
            # from. These are the smaller bp_length reads.
            temp_bp_list.append(original_sequence[random_start_location: random_start_location + bp_length])
            # adds them to this temp list, and returns this temp list. Which will be stored in a chromosome
            # as its small_bp_segment_list. I should rename all of those to reads_list and the other thing to
            # samples list...
    return temp_bp_list


# generates output files from imported data. These are formatted as a sample_regions.txt file which
# hold the larger samples (where they start and end in the sequence for a chromosome) and also
# a small_base_pair_samples.txt which is just a list of all of the reads from larger samples. This output
# will change to be fasta files and be formatted in that same format
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


# BasePairSegment class
# this class stores a single sample
class BasePairSegment:
    def __init__(self, start_int, sequence_string, sequence_length, bp_length):
        self.start_location = start_int
        self.sequence = sequence_string
        self.end = start_int + (sequence_length - 1)
        self.overlap_start = self.start_location - (bp_length - 1)


# chromosome class
# this class stores the entire sequence and header for an entire chromosome (2L, 2R...)
# this class also stores the longer samples for the chromosome as an array of BasePairSegment objects, and also
# a list which holds the reads from the samples that it has stored as those objects.
class Chromosome:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.bp_segment_arr = []
        self.small_bp_segment_list = []
        self.chromosome_name = header.split()[0]


main()
