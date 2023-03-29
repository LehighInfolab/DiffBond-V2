import getopt
import os
import math
import sys
import argparse
import shutil
import os


def parse_files(path):
    counter = 0
    print("### READING FOLDERS ###")
    for i in os.listdir(path):
        dir = i
        subdir_files = []
        for j in os.listdir(path + "/" + i):
            for k in os.listdir(path + "/" + i + "/" + j):
                if k[:5] == "final":
                    subdir_files.append(path + "/" + i + "/" + j + "/" + k)
        os.system(
            "python DiffBond_v2.py -i "
            + subdir_files[0]
            + " "
            + subdir_files[1]
            + " -m a g -o "
            + dir
        )

        # if counter == 2:
        #     return
        # counter = counter + 1


def main():
    parse_files("1brs_muts")


if __name__ == "__main__":
    main()
