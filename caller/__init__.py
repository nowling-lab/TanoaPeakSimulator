
from os import pathsep, write


def call_peaks_calculus(file_name, outputdir, calculation_style, write_depths):
    """
    Takes a depth file and calculates where peaks are within that file
    :param file_name: the path to that depth file
    """
    depths = {}
    chromosome = "null"
    output_file = open(outputdir + "/" + "tanoa_called_peaks.peaks", 'w')
    if write_depths:
        compressed_depths_output = open(outputdir + "/" + "compressed_depths.depths", 'w')
    with open(file_name, "r") as file:
        for line in file:
            line.strip()
            split_line = line.split()
            ch = split_line[0].strip()
            if chromosome == "null":
                chromosome = ch
            if chromosome == ch:
                depths[split_line[1].strip()] = split_line[2].strip()
            else:
                #do processing here
                print("Compressing " + chromosome)
                compressed_depths = compress(depths)
                depths = {}
                if calculation_style == "Calculus":
                    print("Finding critical points")
                    criticals = find_critical(compressed_depths)
                    print("Finding inflection points")
                    inflections = find_inflection(compressed_depths)
                    print("Scrubbing extranious criticals")
                    delete_minimum_criticals(criticals, inflections, compressed_depths)
                    print("removing redundant inflections")
                    remove_redundant_inflections(criticals, inflections)
                    write_peaks_to_file(criticals, inflections, output_file, chromosome)
                elif calculation_style == "Conceptual":
                    print("Calling peaks")
                    peak_list = call_peaks_conceptual(compressed_depths)
                    inflections = []
                    write_peaks_to_file(peak_list, inflections, output_file, chromosome)
                if write_depths:
                    write_depths_to_file(compressed_depths_output, chromosome, compressed_depths)
                print("Peaks in " + chromosome + " called and written to file")
                #####
                chromosome = ch
                depths[split_line[1]] = split_line[2]

def compress(depths):
    depths_compressed = {}
    keys = list(depths.keys())
    start_index = 0
    start = depths[keys[start_index]]
    for x in range(len(keys) - 1):
        end = depths[keys[x + 1]]
        if start != end:
            key =  keys[start_index] + "-" + keys[x] 
            depths_compressed[key] = start
            start = end
            start_index = x + 1

    return depths_compressed

def find_critical(depths):
    """
    Finds critical points using the first derivative using numerical differentation
    :param depths: the depths from the depths file as a dictionary
    """
    index = 0
    critical_list = []
    keys = list(depths.keys())
    for x in keys:
        if index != 0 and index <= len(keys) - 2:
            before = depths[keys[index-1]]
            after = depths[keys[index+1]]
            slope = (int(after) - int(before))/2
            # uses the formula f(x+h) - (f(x-2))/2h so that the slope (first derivate) is cnetered on the 
            # intended point
            if slope <= 0.2 and int(depths[x]) > 1:
                critical_list.append(x)
        index += 1
    return critical_list

def find_inflection(depths):
    """
    Finds inflection points using the second derivatve using numerical differentation
    :param depths: the depths from the depths file as a dictionary
    """
    index = 0
    inflection_list = []
    keys = list(depths.keys())
    for x in keys:
        if index >= 2 and index <= len(keys) - 3:
            left_before =  depths[keys[index - 2]]
            left_after = depths[keys[index - 1]]
            left_slope = (int(left_after) - int(left_before))/2
           
            right_before = depths[keys[index + 1]]
            right_after = depths[keys[index + 2]]
            right_slope = (int(right_after) - int(right_before))/2

            second_derivative = (int(right_slope) - int(left_slope))/2
            # uses the same formula as in find_critical, but uses two first order derivatives to calculate
            # the slope/second order derivative, which would be the inflection points
            if second_derivative <= 0:
                inflection_list.append(x)
        index += 1
    return inflection_list

def delete_minimum_criticals(criticals, inflections, compresssed_depths):
    criticals_to_remove = []
    inflections_to_remove = set()
    for x in range(len(criticals) - 1):
        for y in range(len(inflections) - 1):
            if inflections[y] > criticals[x]:
                if compresssed_depths[inflections[y]] > compresssed_depths[criticals[x]]:
                    criticals_to_remove.append(criticals[x])
                    inflections_to_remove.add(inflections[y])
                    inflections_to_remove.add(inflections[y-1])
                    break
    for x in criticals_to_remove:
        criticals.remove(x)
    for x in inflections_to_remove:
        inflections.remove(x)

def remove_redundant_inflections(criticals, inflections):
    inflections_to_remove = set()
    safe_inflections = set()
    for x in range(len(inflections) - 2):
        rising = inflections[x]
        falling = inflections[x+1]
        found_critical = False
        for y in criticals:
            if rising < y < falling:
                found_critical = True
                safe_inflections.add(rising)
                safe_inflections.add(falling)
        if not found_critical:
            inflections_to_remove.add(rising)
            inflections_to_remove.add(falling)

    for x in safe_inflections:
        if x in inflections_to_remove:
            inflections_to_remove.remove(x)
    for x in inflections_to_remove:
        inflections.remove(x)

def write_peaks_to_file(criticals, inflections, output_file, chromosome):
    output_str = ">" + chromosome + " Peaks: \n"
    for x in range( len(criticals) - 1):
        output_str = output_str + str(x + 1) + ": " + criticals[x] + "\n"
    output_file.write(output_str)

def write_depths_to_file(output, chromosome, compressed_depths):
    output.write(">" + chromosome + " Depths: \n")
    keys = list(compressed_depths.keys())
    for key in keys:
        output.write(key + ": " + compressed_depths[key] + "\n")
        

###########################################
# Below here is for the Conceptual Approach
###########################################

def call_peaks_conceptual(compressed_depths):
    keys = list(compressed_depths.keys())
    peak_list = []
    prev_depth = 0
    is_falling = False
    for index, key in enumerate(keys):
        current_depth = int(compressed_depths[key])
        if prev_depth > current_depth:
            is_max = verify_max_steps(compressed_depths, keys, 5, keys[index - 1], index - 1)
            if is_max:
                peak_list.append(keys[index - 1])
            is_falling = True
        elif current_depth > prev_depth:
            is_falling = False
        prev_depth = current_depth
    return peak_list

def verify_max_bps(compressed_depths, keys, window_length, peak_key, key_index):
    half_window = int(window_length/2)
    split_peak = peak_key.split("-")
    left_of_peak = split_peak[0]
    right_of_peak = split_peak[1]
    middle_peak = int(right_of_peak) - int(left_of_peak)
    left_target = middle_peak - half_window
    right_target = middle_peak + half_window
    window_keys = generate_window_bps(keys, key_index, left_target, right_target)
    return check_max(window_keys, compressed_depths, key_index, peak_key)
    
def verify_max_steps(compressed_depths, keys, steps, peak_key, key_index):
    window_keys = generate_window_steps(keys, key_index, steps)
    return check_max(window_keys, compressed_depths, peak_key)

def check_max(window_keys, compressed_depths, peak_key):
    potential_max = compressed_depths[peak_key]
    is_max = True
    for key in window_keys:
        if compressed_depths[key] > potential_max:
            is_max = False
    return is_max

def generate_window_steps(keys, key_index, steps):
    window_keys = []
    for x in range(1, steps):
        if key_index - x >= 0:
            window_keys.append(keys[key_index - x])
    for x in range(1, steps):
        if key_index + x < len(keys):
            window_keys.append(keys[key_index + x])
    return window_keys

def generate_window_bps(keys, key_index, left_target, right_target):
    window_keys = []
    # Left side first
    left_window_index = find_target(keys, key_index, left_target, -1)
    # Right side next
    right_window_index = find_target(keys, key_index, right_target, 1)
    if right_window_index < len(keys) - 2:
        right_window_index += 1
    for index in range(left_window_index, right_window_index):
        window_keys.append(keys[index])
    return window_keys

def find_target(keys, key_index, target, right_left_modifier):
    found_left_index = False
    start_index = key_index
    while not found_left_index:
        left_answer = keys[start_index]
        split_answer = left_answer.split("-")
        if target <= int(split_answer[0]) or int(split_answer[1]):
            found_left_index = True
        start_index += right_left_modifier
        if start_index == 0:
            return start_index
        elif start_index == len(keys) - 1:
            return len(keys) - 1

    return start_index
