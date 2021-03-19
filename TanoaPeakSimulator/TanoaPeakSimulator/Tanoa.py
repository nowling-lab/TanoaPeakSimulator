import secrets


def main():
    genome_file = open("dmel-all-chromosome-r6.38.fasta", "r")
    sequence_string = genome_file.read()
    random_string = get_random_base_pair_string(1000, sequence_string)
    print(random_string)


def get_random_base_pair_string(length, sequence):
    random_start_location = secrets.randbelow(len(sequence) - (length + 1))
    return sequence[random_start_location: (random_start_location + length)]


main()
