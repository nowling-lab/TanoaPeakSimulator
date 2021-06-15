
def call_peaks(file_name):
    """
    Takes a depth file and calculates where peaks are within that file
    :param file_name: the path to that depth file
    """
    depths = {}
    chromosome = "null"
    with open(file_name, "r") as file:
        for line in file:
            split_line = line.split(",")
            ch = split_line[0]
            if chromosome == ch:
                depths[split_line[1]] = split_line[2]
            else:
                #do processing here
                criticals = find_critical(depths)
                inflections = find_inflection(depths)
                print(criticals, inflections)
                #####
                chromosome = ch
                depths = {}
                depths[split_line[1]] = split_line[2]

def find_critical(depths):
    """
    Finds critical points using the first derivative using numerical differentation
    :param depths: the depths from the depths file as a dictionary
    """
    index = 0
    critical_list = []
    for x in depths.keys():
        if index != 0 or index != len(depths.keys()) - 1:
            before = depths.keys()[index-1]
            after = depths.keys()[index+1]
            slope = (after - before)/2
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
    keys = depths.keys()
    for x in keys:
        if index <= 1 or index >= len(keys) - 2:
            left_before =  keys[index - 2]
            left_after = keys[index - 1]
            left_slope = (left_after - left_before)/2
           
            right_before = keys[index + 1]
            right_after = keys[index + 2]
            right_slope = (right_after - right_before)/2

            second_derivative = (right_slope - left_slope)/2
            # uses the same formula as in find_critical, but uses two first order derivatives to calculate
            # the slope/second order derivative, which would be the inflection points
            if second_derivative == 0:
                inflection_list.append(x)
        index += 1
    return inflection_list
