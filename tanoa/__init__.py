import random

def get_random_base_pair_string(enhancer_region_length, read_length, sequence, array_of_segments):
    """
    Gets a psudo-random string of a given amount of characters from the main sequence. This will move to taking
    the buffers on the sides and saving that as well so when sampled later we don't have to return to the
    Very large string. That's inefficient.

    :param enhancer_region_length: Length of the big sample that you want to take
    :param read_length: the length of a smaller base pair that you will sample from the broader region this is setup so
    that the buffer can be dynamic based on the size of the reads you want to match to the original chromosome are
    :param sequence: The actual character sequence, generally attached to a chromosome like 2L or 2R
    :param array_of_segments: array of EnhancerRegions (samples) from the original chromosome.
    This array is stored within a chromosome as bp_sequence_array and is modifying that object.
    :return: no returns. modifies variables stored in objects...
    """
    if len(sequence) > (enhancer_region_length + read_length + 1):  # Ignore sequences that don't have enough length
        # to read from
        min_pull_location = read_length - 1  # minimum location at which a region
        # can be taken from to prevent out of bounds
        max_pull_location = len(sequence) - (enhancer_region_length + read_length)
        random_start_location = random.randrange(min_pull_location, max_pull_location)

        # gets a random sample from the random start location to the sequence length. This eventually should be
        # random_start_location + enhancer_region_length and + buffer... This is a simple substring
        mini_sequence = sequence[random_start_location: (random_start_location + enhancer_region_length)]
        # makes a EnhancerRegion object from the mini_sequence string and associated data
        temp_segment = EnhancerRegion(random_start_location, mini_sequence, enhancer_region_length, read_length)
        # stores that in the passed array from chromosome
        array_of_segments.append(temp_segment)


def get_random_base_pair_string_list(list_of_bp_segments, read_length, num_reads,
                                     original_sequence, chromosome, read_number, read_counts_dict):
    """
    Gets smaller reads from larger samples and returns them as a list
    :param read_counts_dict: dictionary of position read and how many times that position has been read
    :param list_of_bp_segments: the larger samples, stored as an array in a chromosome object
    :param read_length: the length of one single read. Generally between 50-100bps in length
    :param num_reads: how many of these smaller reads you want to generate per larger sample
    :param original_sequence: this is the whole original sequence
    :param chromosome: the chromosome name, like 2L, 2R or similar
    :param read_number: an incrementing number from zero which keeps track of which read is which
    :return:
    """
    temp_bp_list = []
    for x in list_of_bp_segments:  # loops through the list of EnhancerRegion objects
        for y in range(num_reads):
            random_start_location = random.randrange(x.overlap_start, x.end)  # gets a random location to get a read
            # from. These are the smaller read_length reads.
            sequence = original_sequence[random_start_location: random_start_location + read_length]
            read_number += 1
            temp_bp_list.append(SingleRead(chromosome, random_start_location, random_start_location + read_length,
                                           sequence, read_number))
            for z in range(random_start_location, random_start_location + read_length):
                if z in read_counts_dict:
                    read_counts_dict[z] += 1
                else:
                    read_counts_dict[z] = 1
            # adds them to this temp list, and returns this temp list. Which will be stored in a chromosome
            # as its reads_list. I should rename all of those to reads_list and the other thing to
            # samples list...
    return temp_bp_list, read_number


def generate_output_files(chromosome, output_file, samples_file, read_counts_file, read_counts_dict):
    """
    generates output files from imported data. These are formatted as a sample_regions.txt file which
    hold the larger samples (where they start and end in the sequence for a chromosome) and also
    a small_base_pair_samples.txt which is just a list of all of the reads from larger samples.
    :param read_counts_file: the file in which to write the read counts to
    :param read_counts_dict: Dictionary with how many times a single location has been read
    :param chromosome: the chromosome object which stores a chromosome and its sequence
    :param output_file: the file in which to write the enhancer sequences to
    :param samples_file: the file to write the reads to (the shorter 50-100bp reads)
    :return: no return
    """
    output_to_file = ""
    for x in chromosome.enhancer_list:
        output_to_file += chromosome.chromosome_name[1:] + ":" + str(x.start_location) + "-" + str(x.end) + "\n"
    output_file.write(output_to_file)

    random_regions = ""
    for x in chromosome.reads_list:
        random_regions += x.output + "\n" + x.sequence + "\n"
    samples_file.write(random_regions)

    read_counts_output = ""
    for x in read_counts_dict.keys():
        read_counts_output += chromosome.chromosome_name[1:] + "," + str(x) + "," + str(read_counts_dict[x]) + "\n"
    read_counts_file.write(read_counts_output)


class EnhancerRegion:
    """
    EnhancerRegion class
    this class stores a single sample
    """

    def __init__(self, start_int, sequence_string, enhancer_region_length, read_length):
        self.start_location = start_int
        self.sequence = sequence_string
        self.end = start_int + (enhancer_region_length - 1)
        self.overlap_start = self.start_location - (read_length - 1)


class SingleRead:
    def __init__(self, chromosome, start, end, sequence, read_number):
        self.chromosome = chromosome.replace(">", "")
        self.read_number = ">" + str(read_number) + " "
        self.output = self.read_number + self.chromosome + ":" + str(start) + "-" + str(end) + ""
        self.sequence = sequence


class Chromosome:
    """
    chromosome class
    this class stores the entire sequence and header for an entire chromosome (2L, 2R...)
    this class also stores the longer samples for the chromosome as an array of EnhancerRegion objects, and also
    a list which holds the reads from the samples that it has stored as those objects.
    """

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.enhancer_list = []
        self.reads_list = []
        self.chromosome_name = header.split()[0]


