
def call_peaks(file_name):
    """
    Takes a depth file and calculates where peaks are within that file
    :param file_name: the path to that depth file
    """
    depths = {}
    chromosome = "null"
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
                print(criticals)
                print()
                print(inflections)
                #####
                return
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
            if slope == 0:
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

