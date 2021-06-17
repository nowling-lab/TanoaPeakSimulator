
def call_peaks(file_name, outputdir):
    """
    Takes a depth file and calculates where peaks are within that file
    :param file_name: the path to that depth file
    """
    depths = {}
    chromosome = "null"
    output_file = open(outputdir + "/" + "tanoa_called_peaks.peaks", 'w')
    with open(file_name, "r") as file:
        for line in file:
            line.strip()
            split_line = line.split(",")
            ch = split_line[0].strip()
            if chromosome == "null":
                chromosome = ch
            if chromosome == ch:
                depths[split_line[1].strip()] = split_line[2].strip()
            else:
                #do processing here
                compressed_depths = compress(depths)
                criticals = find_critical(compressed_depths)
                inflections = find_inflection(compressed_depths)
                delete_minimum_criticals(criticals, inflections, compressed_depths)
                remove_redundant_inflections(criticals, inflections)
                write_peaks_to_file(criticals, inflections, output_file, chromosome)
                #####
                chromosome = ch
                depths = {}
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
            if slope == 0 and int(depths[x]) > 1:
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
            if second_derivative == 0:
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