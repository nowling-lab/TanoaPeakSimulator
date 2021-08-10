import os

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
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
        output_file (File): The file in which to write out for the peaks
        write_depths (Bool): Bool which decides if to write the compressed depth file to file or not
        compressed_depths_output (File): File to output if the option of writing compressed depths was given
    """

    print("Calling peaks")
    smoothed_compressed = smooth_compressed(compressed_depths, 1, False)
    peak_list = find_maximums(smoothed_compressed) #calls peaks
    smoothed_compressed = smooth_compressed(smoothed_compressed, 1, True) 
    print("Cleaning peaks")
    peak_list = clean_peaks(peak_list, smoothed_compressed, compressed_depths)
    print("Connecting nearby peaks")
    peak_list = connect_peaks(peak_list, compressed_depths)
    print("Writing to file")
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
    start = int(depths[keys[0]])
    jump = False
    for index, key in enumerate(keys):
        end = int(depths[key])
        if index > 0:
            int_prev = int(keys[index - 1])
            int_end = int(key)
            if int_end - int_prev != 1:
                jump = True
        if start != end or jump: #Finds all 
            key_out = int(start_key), int(keys[index -1])
            depths_compressed[key_out] = int(start)
            start = end
            start_key = key
            jump = False

    return depths_compressed

def smooth_compressed(compressed_depths, window, to_round):
    """Takes a compressed depth file and averages each position to window(int) left and right of the point to created a smoothed compressed depth file

    Args:
        compressed_depths (dict): a position_start, position_end : read_depth formated dictionary of depths 
        window (int): The width, left and right, of a given point to formulate a window in which to take an average from to form a smooth compressed point
        to_round (bool): A boolean value which chooses if the smooth compressed should be rounded to the nearest int or not

    Returns:
        dict: A smooth compressed dictionary of position_start, position_end : read_depth format
    """
    half_window = window
    keys = list(compressed_depths.keys())
    new_compressed = {}
    for index, key in enumerate(keys):
        array_to_avg = []
        if index >= half_window:
            for x in range(index-half_window, index):
                array_to_avg.append(compressed_depths[keys[x]])
        else:
            for x in range(0, index):
                array_to_avg.append(compressed_depths[keys[x]])

        if index + half_window < len(keys):
            for x in range(index, index+half_window + 1):
                array_to_avg.append(compressed_depths[keys[x]])
        else:
            for x in range(index, len(keys)):
                array_to_avg.append(compressed_depths[keys[x]])
        
        if to_round:
            average = sum(array_to_avg)/len(array_to_avg)
            average = round(average)
        else:
            average = sum(array_to_avg)/len(array_to_avg)
        new_compressed[key] = average
    return new_compressed
          
def connect_peaks(peak_list, compressed_depths):
    """Connect_peaks takes a list of peaks, and compressed depths and attempts to combine nearby peaks which should be combined, and not combine nearby peaks that shouldn't be combined. This gets rid of a significant amount of the cases of multiple peaks on a single peak

    Args:
        peak_list (list): A list of all called peaks from call_peaks(compressed_depths)
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 

    Returns:
        new_peak_list: A new peak list, which is the result of combining nearby peaks
    """
    search_start = 0
    new_peak_list = []
    new_peak = None
    connect_peaks = False
    for index, peak in enumerate(peak_list):
        if index < len(peak_list) - 1:
            next_peak_start = peak_list[index + 1][0]
            compressed_index = find_in_compressed(peak[1], compressed_depths, search_start)
            search_start = compressed_index
            compressed_index_next = find_in_compressed(next_peak_start, compressed_depths, search_start)

            distance_between_peak_ends = compressed_index_next - compressed_index
            current_peak_width = find_peak_width(peak, compressed_depths)
            next_peak_width = find_peak_width(peak_list[index + 1], compressed_depths)

            too_wide = False
            if current_peak_width and next_peak_width > 900:
                too_wide = True

            found_peak = False
            if peak[0] == 246676:
                print(peak, too_wide)
                print(peak_list[index + 1])
                found_peak = True
            if peak[0] == 247549:
                found_peak = True

            valid_distance = 0 < distance_between_peak_ends <= 2

            if valid_distance and not too_wide:
                if connect_peaks:
                    new_peak = new_peak_list[-1][0], peak_list[index + 1][1]
                    new_peak_list.pop()
                    new_peak_list.append(new_peak)
                else:
                    new_peak = peak[0], peak_list[index + 1][1]
                    new_peak_list.append(new_peak)
                connect_peaks = True
                if found_peak:
                    print(peak, new_peak)
                    print()
                    print(compressed_index_next, compressed_index)
                    found_peak = False
            else:
                connect_peaks = False
                if len(new_peak_list) >= 1 and check_overlap(peak, new_peak_list[-1]) != 0:
                    new_peak_list.append(peak)
                elif len(new_peak_list) == 0:
                    new_peak_list.append(peak)
                elif check_overlap(peak, new_peak_list[-1]) == 0:
                    if peak[0] == new_peak_list[-1][1]:
                        new_peak_list.append(peak)
                else:
                    print(peak, "?")

    if not connect_peaks:
        new_peak_list.append(peak_list[-1])
    return new_peak_list
            
def find_max_in_peak(peak, compressed_depths):
    """Given a given peak, this method finds the maximum depth within that peak

    Args:
        peak (tupol): A peak as defined as a start, end tuple
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 

    Returns:
        int: the maximum depth within the range of values which constitutes a peak
    """
    keys = list(compressed_depths.keys())

    compressed_index = find_in_compressed(peak[0], compressed_depths, 0)
    compressed_index_end = find_in_compressed(peak[1], compressed_depths, compressed_index)

    peak_arr = []
    for x in range(compressed_index, compressed_index_end + 1):
        peak_arr.append(compressed_depths[keys[x]])
    
    return max(peak_arr)
    

def find_peak_width(peak, compressed_depths):
    """Given a peak, this method finds the width in positions of that peak

    Args:
        peak (tupol): A peak as defined as a start, end tuple
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths

    Returns:
        int: The width of a given peak in positions 
    """
    keys = list(compressed_depths.keys())
    
    compressed_index = find_in_compressed(peak[0], compressed_depths, 0)
    compressed_index_end = find_in_compressed(peak[1], compressed_depths, compressed_index)

    start = keys[compressed_index]
    end = keys[compressed_index_end]

    return end[1] - start[0]


def find_in_compressed(peak_start, compressed_depths, search_start_index):
    """Finds the location of the start of a given peak in the compressed depth file. This method also can be called multiple times in a row, with an increasing search start index to reduce runtime of this method

    Args:
        peak_start (int): the position of the start of the given peak that would be found
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
        search_start_index (int): the index in which to begin the start. This is generally the last found index from the previous search to decrease runtime

    Returns:
        [type]: [description]
    """
    keys = list(compressed_depths.keys())
    for index in range(search_start_index, len(keys)):
        key_start, key_end = keys[index]
        if key_start == peak_start or key_end == peak_start:
            return index
    return -1

def clean_peaks(peak_list, smoothed_compressed, compressed_depths):
    """Removes, extends and registers new peaks for all peaks

    Args:
        peak_list (list): A list of all called peaks from call_peaks(compressed_depths)
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
    Returns:
        list: A new peak_list, which removes, combines and lengthens peaks to clean up the list
    """
    keys = list(compressed_depths.keys())
    new_peaks = []
    start_index = 0

    for peak in peak_list:
        start_index, background_index_left, background_index_right = get_background(peak, compressed_depths, smoothed_compressed, keys, start_index)
        temp_peak = register_new_peak(background_index_left, background_index_right, keys)
        new_peaks.append(temp_peak)

    return new_peaks

def register_new_peak(background_index_left, background_index_right, keys):
    """Registers a new peak by creating a peak from the background values given

    Args:
        background_index_left (int): index in compressed depths of background left to a peak
        background_index_right (int): index in compressed depths of background right to a peak
        keys (List): The keys for the compressed depths dictionary

    Returns:
        str: A peak in start - end format that is the result of registration
    """
    
    if background_index_right - background_index_left >= 5:
        peak = keys[background_index_left + 1][0], keys[background_index_right - 1][1]
    else:
        peak = keys[background_index_left][0], keys[background_index_right][1]

    return peak

def compare_peaks(peak1, peak2):
    """Compares two given peaks

    Args:
        peak1 (string): The first peak to be compared
        peak2 (string): The second peak to be compared

    Returns:
        int: -1 for smaller, 0 for same, 1 for larger (relative to peak 1)
    """
    peak1_start = peak1[0]
    peak2_start = peak2[0]
    if peak1_start == peak2_start:
        return 0
    elif peak1_start > peak2_start:
        return 1
    else:
        return -1

def check_overlap(peak1, peak2):
    """check_overlap checks if a peak overlaps, if they do, it returns 0, if not, it returns 1 if the first given peak is greater than the second, or -1 if otherwise

    Args:
        peak1 (tuple): start, end tuple which constitutes a peak
        peak2 (tuple): start, end tuple which constitutes a peak

    Returns:
        [type]: [description]
    """
    peak1_start = peak1[0]
    peak2_start = peak2[0]
    peak1_end = peak1[1]
    peak2_end = peak2[1]
    if peak2_start <= peak1_start <= peak2_end or peak1_start <= peak2_start <= peak1_end:
        return 0
    elif peak1_start > peak2_start:
        return 1
    else:
        return -1

def get_background(peak, compressed_depths, smoothed_compressed, keys, start_index):
    """Finds the index of a given peak within the compressed depths file, finds the background left and right
        of a given peak

    Args:
        peak (str): The peak which the finding background algroithm is centered on
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
        keys (List): The keys for the compressed depths dictionary
        start_index (int): Index which to start searching for a peak within the keys list from, given from the previous found peak's index

    Returns:
        [type]: [description]
    """
    depths_compressed = smoothed_compressed

    peak_index = keys.index(peak, start_index)
    temp_index_left = peak_index
    temp_index_right = peak_index
    background_index_left = None
    background_index_right = None
    # Find peak in compressed_depths file
    # Find Left background
    found_background = False
    background = 5

    depth = compressed_depths[keys[peak_index]]
    if depth < 15:
        depths_compressed = compressed_depths
    elif depth >= 15:
        background = 8
        depths_compressed = compressed_depths
    elif depth >= 50:
        background = 8
        
    while not found_background:
        depth = depths_compressed[keys[temp_index_left]]
        if depth <= background:
            found_background = True
            background_index_left = temp_index_left
        elif temp_index_left > 0:
            temp_index_left -= 1
        else:
            background_index_left = 0
            found_background = True

    # Find Right Background
    found_background = False
    while not found_background:
        depth = depths_compressed[keys[temp_index_right]]
        if depth <= background: 
            found_background = True
            background_index_right = temp_index_right
        elif temp_index_right < len(depths_compressed):
            temp_index_right += 1
        else:
            background_index_right = len(depths_compressed) - 1
            found_background = True

    return peak_index, background_index_left, background_index_right


def write_peaks_to_file(peaks, output_file, chromosome):
    """Wries a given list of peaks to file

    Args:
        peaks (List): List of peaks to be written to file
        output_file (File): The file in which to write out for the peaks
        chromosome (str): The name of the chromosome currently being called
    """
    for peak in peaks:
        output_file.write(chromosome + " " + str(peak[0]) + " " + str(peak[1]) + "\n")

def write_depths_to_file(output, chromosome, compressed_depths):
    """Wries depths to file if the Write_Depths field was set to True

    Args:
        output (File): The file object which is used to write to file
        chromosome (String): The name of the chromosome who's depths are being added to the file
        chromosome (str): The name of the chromosome currently being called
    """
    keys = list(compressed_depths.keys())
    for key in keys:
        output.write(chromosome + " " + str(key) + " " + str(round(compressed_depths[key], 3)) + "\n")
        

def find_maximums(compressed_depths):
    """Searches for and adds found peaks to a list

    Args:
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 

    Returns:
        [List]: A list of found peaks
    """
    keys = list(compressed_depths.keys())
    local_maximums = []
    prev_depth = 0
    for index, key in enumerate(keys):
        current_depth = compressed_depths[key]
        if current_depth > prev_depth: #Finds locations where the slope is positive...
            is_max = verify_max(compressed_depths, keys, 3, keys[index], index)
            if is_max: #Checks those for being a maximum. What if I did no checking?
                local_maximums.append(keys[index]) #TODO: Rename stuff
        prev_depth = current_depth

    return local_maximums
    
def verify_max(compressed_depths, keys, steps, peak_key, key_index):
    """Given a peak, checks if it is a local maximum

    Args:
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
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
        compressed_depths (Dictionary): a position_start, position_end : read_depth formated dictionary of depths 
        peak_key (string): The key to the given peak being checked

    Returns:
        [bool]: A value which verifies if the given point of interest is a maximum or not
    """
    potential_max = compressed_depths[peak_key]
    is_max = True

    if potential_max <= 4: 
        return False

    max_equals = 0

    for key in window_keys:
        key_depth = compressed_depths[key]
        if key_depth > potential_max:
            is_max = False
        # elif key_depth == potential_max:
        #     max_equals += 1

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
