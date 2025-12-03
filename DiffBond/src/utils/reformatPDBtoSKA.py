import getopt
import math
import sys
import optparse

######################################################################################
#### Reformats PDB file so that
######################################################################################


def parsePDB(PDB):
    f = open(PDB, "r")
    l = f.readlines()

    f.close()

    lines = []

    # 	X,Y,Z	coordinates	found	in lines[7],lines[8],lines[9]
    count = 0
    for x in range(len(l)):
        if l[x][0] == "A":
            lines.append(l[x].split())
            count = count + 1

    count = 0
    max = 0
    for x in range(len(lines)):
        lines[x][4] = "A"
        if int(lines[x][5]) > int(max):
            max = int(lines[x][5])
        else:
            while int(lines[x][5]) < int(max):
                lines[x][5] = int(lines[x][5]) + 100
            max = lines[x][5]

    return lines


def writePDB(PDB, simlist):
    PDB = PDB.split(".")
    f = open(PDB[0] + "_edited.pdb", "w+")
    for i in simlist:
        i[0] = i[0].ljust(7)
        i[1] = i[1].rjust(6) + "  "
        i[2] = i[2].ljust(5)
        i[3] = i[3].ljust(5)
        i[4] = i[4].rjust(2)
        i[5] = str(i[5]).rjust(5) + "    "
        i[6] = str("%8.3f" % (float(i[6]))).rjust(9)
        i[7] = str("%8.3f" % (float(i[7]))).rjust(9)
        i[8] = str("%8.3f" % (float(i[8]))).rjust(9)
        i[9] = i[9].rjust(13)
        i[10] = i[10].rjust(13)
        f.write("".join(map(str, i)) + "\n")
    f.close
    print("pdb reformatted")


def main():
    if len(sys.argv) != 2:
        print("Usage: ")
        sys.exit()

    pdbFile = sys.argv[1]

    atoms = parsePDB(pdbFile)
    writePDB(pdbFile, atoms)


if __name__ == "__main__":
    main()
