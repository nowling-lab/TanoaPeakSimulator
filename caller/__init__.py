
from os import pathsep, write
from posixpath import split


def read_and_call(file_name, outputdir, write_depths):
    """Takes a depth file and calculates where peaks are within that file

    Args:
        file_name (String): The name of the depth file
        outputdir (String): The directory to output to
        write_depths (Bool): Bool which decides if to write the compressed depth file to file or not
    """
    depths = {}
    chromosome = None
    output_file = open(outputdir + "/" + "tanoa_called_peaks.peaks", 'w')
    if write_depths:
        compressed_depths_output = open(outputdir + "/" + "compressed_depths.depths", 'w')
    else:
        compressed_depths_output = None

    with open(file_name, "r") as file:
        for line in file:
            line.strip()
            split_line = line.split()
            ch = split_line[0].strip()
            if chromosome == None:
                chromosome = ch
                print("Reading in ", ch)
            if chromosome == ch: #Standard chromosome check to run software 1 ch at a time
                depths[split_line[1].strip()] = split_line[2].strip()
            else:
                #do processing here
                print("Compressing " + chromosome)
                compressed_depths = compress(depths) #Compress the depths to make later algorithms work better
                depths = {}                           #and faster
                prossess_and_write_peaks(chromosome, compressed_depths, output_file, write_depths, compressed_depths_output)
                #####
                print("Reading in ", ch)
                chromosome = ch
                depths[split_line[1]] = split_line[2]
        print("Compressing " + chromosome)
        compressed_depths = compress(depths) #Compress the depths to make later algorithms work better
        depths = {} 
        prossess_and_write_peaks(chromosome, depths, output_file, write_depths, compressed_depths_output)
        #This gets the last chromosome in the file 

    output_file.close
    if write_depths:
        compressed_depths_output.close

def prossess_and_write_peaks(chromosome, compressed_depths, output_file, write_depths, compressed_depths_output):
    """High level method for grouping together the method calls for calling and writing peaks

    Args:
        chromosome (Str): The name of the chromosome currently being called
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format
        output_file (File): The file in which to write out for the peaks
        write_depths (Bool): Bool which decides if to write the compressed depth file to file or not
        compressed_depths_output (File): File to output if the option of writing compressed depths was given
    """
    print("Calling peaks")
    peak_list = call_peaks(compressed_depths) #calls peaks
    print("Cleaning peaks") 
    peak_list = clean_peaks(peak_list, compressed_depths)
    write_peaks_to_file(peak_list, output_file, chromosome) #writes them
    if write_depths:
        write_depths_to_file(compressed_depths_output, chromosome, compressed_depths)
    print(str(len(peak_list)) + " Peaks in " + chromosome + " called and written to file")

def compress(depths):
    """Compress takes a depth file and compresses the depths into "steps": The zones of consistent depth 
    registered as one singular depth for the sake of finding peaks

    Args:
        depths (dictionary): The original depth file in position: depth format

    Returns:
        dictionary: The dictionary of compressed depth files in start_pos-end_pos: depth of step
    """
    depths_compressed = {}
    keys = list(depths.keys())
    start_key = keys[0]
    start = depths[keys[0]]
    jump = False
    #TODO Use touples instead of strings for keys... Sorta large undertaking but eh, it'll be fine
    for index, key in enumerate(keys):
        end = depths[key]
        if index > 0:
            int_prev = int(keys[index - 1])
            int_end = int(key)
            if int_end - int_prev != 1:
                jump = True
        if start != end or jump: #Finds all 
            key_out =  str(start_key) + "-" + str(keys[index - 1]) 
            depths_compressed[key_out] = start
            start = end
            start_key = key
            jump = False

    return depths_compressed

def clean_peaks(peak_list, compressed_depths):
    """Removes, extends and registers new peaks for all peaks

    Args:
        peak_list (list): A list of all called peaks from call_peaks(compressed_depths)
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format

    Returns:
        list: A new peak_list, which removes, combines and lengthens peaks to clean up the list
    """
    keys = list(compressed_depths.keys())
    new_peaks = []
    peak_black_dict = {}
    start_index = 0
    for index, peak in enumerate(peak_list):
        if peak not in peak_black_dict:
            start_index, background_index_left, background_index_right = get_background(peak, compressed_depths, keys, start_index)
            if start_index == keys[background_index_right] and keys[background_index_left]:
                new_peaks.append(peak)
            else:
                temp_peak = register_new_peak(background_index_left, background_index_right, peak_list, index, keys, peak_black_dict)
                new_peaks.append(temp_peak)

    return new_peaks

def register_new_peak(background_index_left, background_index_right, peak_list, peak_index, keys, peak_black_dict):
    """Takes a peak, and registers a new peak by extending, combining or removing extranious current nearby peaks

    Args:
        background_index_left (int): index in compressed depths of background left to a peak
        background_index_right (int): index in compressed depths of background right to a peak
        peak_list (list): A list of all called peaks from call_peaks(compressed_depths)
        peak_index (int): Index of given peak in compressed depths dict
        keys (List): The keys for the compressed depths dictionary
        peak_black_dict (dictionary): A list of blacklisted peaks (peaks which were deleted, or in a combination)

    Returns:
        str: A peak in start - end format that is the result of registration
    """
    # Get original peak, left step and right step
    left_end = keys[background_index_left]
    right_end = keys[background_index_right]

    # combine all peaks within the boundary and blacklist necessary peaks
    # get furthest left peak
    found_left = False
    index_copy = peak_index
    left_limit = ""
    while not found_left:
        compare = compare(peak_list[index_copy], left_end)
        if compare == -1:
            peak_black_dict[peak_list[index_copy]] = True
            index_copy -= 1
            if index_copy <= 0:
                break
        elif compare == 0:
            break
        else:
            index_copy += 1
            left_limit = peak_list[index_copy]
            found_left = True

    # get furthest right peak
    found_right = False
    right_limit = ""
    index_copy = peak_index
    while not found_right:
        compare = compare(peak_list[index_copy], right_end)
        if compare == -1:
            peak_black_dict[peak_list[index_copy]] = True
            index_copy += 1
            if index_copy == len(peak_list):
                break
        elif compare == 0:
            break
        else:
            index_copy -= 1
            right_limit = peak_list[index_copy]
            found_right = True

    # use the found variables to combine and make a new peak
    if peak_list[peak_index] in left_limit and right_limit:
        found_left == False
        found_right == False

    peak = ""
    if found_left and found_right:
        left_start = left_limit.split("-")[0]
        right_end = right_limit.split("-")[1]
        peak = left_start + "-" + right_end
    elif found_left:
        left_start = left_limit.split("-")[0]
        right_end =  peak_list[peak_index].split("-")[1]
        peak = left_start + "-" + right_end
    elif found_right:
        left_start = peak_list[peak_index].split("-")[0]
        right_end =  right_limit.split("-")[1]
        peak = left_start + "-" + right_end
    else:
        # steps_between = background_index_right - background_index_left
        # quarter_peak = steps_between//3
        # left_peak = keys[background_index_left + quarter_peak]
        # right_peak = keys[background_index_right - quarter_peak]
        # left_start = left_peak.split("-")[0]
        # right_end = right_peak.split("-")[1]
        left_start = left_end.split("-")[0]
        right_side = right_end.split("-")[1]
        peak = left_start + "-" + right_side
        #peak = peak_list[peak_index]
    # return new combined peak if applicible
    peak_black_dict[peak_list[peak_index]] = True
    return peak

def compare(peak1, peak2):
    """Compares two given peaks

    Args:
        peak1 (string): The first peak to be compared
        peak2 (string): The second peak to be compared

    Returns:
        int: -1 for smaller, 0 for same, 1 for larger (relative to peak 1)
    """
    peak1_start = int(peak1.split("-")[0])
    peak2_start = int(peak2.split("-")[0])
    if peak1_start == peak2_start:
        return 0
    elif peak1_start > peak2_start:
        return 1
    else:
        return -1
    
def get_background(peak, compressed_depths, keys, start_index):
    """Finds the index of a given peak within the compressed depths file, finds the background left and right
        of a given peak

    Args:
        peak (str): The peak which the finding background algroithm is centered on
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format
        keys (List): The keys for the compressed depths dictionary
        start_index (int): Index which to start searching for a peak within the keys list from, given from the previous found peak's index

    Returns:
        [type]: [description]
    """
    peak_index = keys.index(peak, start_index)
    temp_index_left = peak_index
    temp_index_right = peak_index
    background_index_left = 0
    background_index_right = 0
    # Find peak in compressed_depths file
    # Find Left background
    found_background = False
    while not found_background:
        depth = int(compressed_depths[keys[temp_index_left]])
        if depth <= 3:
            found_background = True
            background_index_left = temp_index_left
        elif temp_index_left > 0:
            temp_index_left -= 1
        else:
            background_index_left = 0
    # Find Right Background
    found_background = False
    while not found_background:
        depth = int(compressed_depths[keys[temp_index_right]])
        if depth <= 3: #TODO: CHECK BACKGROUND SIZE. 3 GIVES LESS PEAKS THAN 4. MUST CHECK ACCURACY AGAIN
            found_background = True
            background_index_right = temp_index_right
        elif temp_index_right < len(compressed_depths):
            temp_index_right += 1
        else:
            background_index_right = len(compressed_depths) - 1

    return peak_index, background_index_left, background_index_right


def write_peaks_to_file(peaks, output_file, chromosome):
    """Wries a given list of peaks to file

    Args:
        peaks (List): List of peaks to be written to file
        output_file (File): The file in which to write out for the peaks
        chromosome (str): The name of the chromosome currently being called
    """
    for peak in peaks:
        split_peak = peak.split("-")
        begin_peak = split_peak[0]
        end_peak = split_peak[1]
        output_file.write(chromosome + " " + begin_peak + " " + end_peak + "\n")

def write_depths_to_file(output, chromosome, compressed_depths):
    """Wries depths to file if the Write_Depths field was set to True

    Args:
        output (File): The file object which is used to write to file
        chromosome (String): The name of the chromosome who's depths are being added to the file
        chromosome (str): The name of the chromosome currently being called
    """
    keys = list(compressed_depths.keys())
    for key in keys:
        split_keys = key.split("-")
        output.write(chromosome + " " + split_keys[0] + " " + split_keys[1] + " " + compressed_depths[key] + "\n")
        

def call_peaks(compressed_depths):
    """Searches for and adds found peaks to a list

    Args:
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format

    Returns:
        [List]: A list of found peaks
    """
    keys = list(compressed_depths.keys())
    peak_list = []
    prev_depth = 0
    for index, key in enumerate(keys):
        current_depth = int(compressed_depths[key])
        if current_depth > prev_depth: #Finds locations where the slope is positive...
            is_max = verify_max(compressed_depths, keys, 3, keys[index], index)
            if is_max: #Checks those for being a maximum. What if I did no checking?
                peak_list.append(keys[index])
        prev_depth = current_depth
    return peak_list
    
def verify_max(compressed_depths, keys, steps, peak_key, key_index):
    """Given a peak, checks if it is a local maximum

    Args:
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format
        keys (List): The keys for the compressed depths dictionary
        steps (int): How many steps left and write will be searched to verify a local maximum
        peak_key (string): The key to the given peak being checked
        key_index (int): The index in the keys list of the peak_key key

    Returns:
        [boolean]: A true or false statement to if the given point is a maximum or not
    """
    window_keys = generate_window(keys, key_index, steps)
    return check_max(window_keys, compressed_depths, peak_key)

def check_max(window_keys, compressed_depths, peak_key):
    """Given a list of keys, verifies that the given point(peak_key) corrisponds with the maximum
    value present in depths which those keys corrispond to

    Args:
        window_keys (List): A list of keys x steps away from a given peak
        compressed_depths (Dictionary): A dictionary of compressed depths in location: depth format
        peak_key (string): The key to the given peak being checked

    Returns:
        [bool]: A value which verifies if the given point of interest is a maximum or not
    """
    potential_max = int(compressed_depths[peak_key])
    is_max = True

    if potential_max <= 4: #TODO: CHANGED FROM 4 TO 3 AS I MOVED BACKGROUND TO 3. VERIFY AGAIN...
        return False

    max_equals = 0

    for key in window_keys:
        key_depth = int(compressed_depths[key])
        if key_depth > potential_max:
            is_max = False
        elif key_depth == potential_max:
            max_equals += 1
    if max_equals >= 2:
        is_max = False

    return is_max
    
def generate_window(keys, key_index, steps):
    """Generates a list (window) steps distance to the left and right of a given peak 

    Args:
        keys (List): The keys for the compressed depths dictionary
        key_index (int): The index in the keys list of the peak_key key
        steps (int): How many steps left and write will be searched to verify a local maximum

    Returns:
        [list]: A list of keys that includes all keys steps away from a given location
    """
    window_keys = []
    for x in range(1, steps + 1):
        if key_index - x >= 0:
            window_keys.append(keys[key_index - x])
    for x in range(1, steps + 1):
        if key_index + x < len(keys):
            window_keys.append(keys[key_index + x])
    return window_keys
