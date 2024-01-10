import re
import os

path = "datasets/wt"

for i in os.listdir(path):
    dir = i
    subdir_files = []
    for j in os.listdir(path + "/" + i):
        if j[0] != "H":
            # print(j)
            continue
        for k in os.listdir(path + "/" + i + "/" + j):
            if k[:14] == "hydrogensAdded":
                # Set up the input and output file paths (can be the same)
                input_file_path = path + "/" + i + "/" + j + "/" + k

                with open(input_file_path, "r") as file:
                    lines = file.readlines()  # Read all lines in the file into a list

                with open(
                    input_file_path, "w"
                ) as file:  # Re-open the file in write mode
                    for line in lines:
                        # If the line has "60" followed by a letter, skip writing it back to the file
                        if re.search("186[a-zA-Z]", line):
                            continue
                        # Write the line back to the file
                        file.write(line)
