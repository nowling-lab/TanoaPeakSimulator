
def read_depth(file_name):
    depths = []
    with open(file_name, "r") as file:
        depths.append(int(file.read(1)))
    for x in range(1, 50):
        print(depths[x])

