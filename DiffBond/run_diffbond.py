import getopt
import os
import math
import sys
import argparse
import shutil
import os
import csv

def make_output_dir(pdb, path):
	try:
		os.mkdir("Results/" + str(path) + "/" + str(pdb))
	except Exception as error:
		pass

def parse_skempi():
	## Make a skempi directory for current base
	# base = path.split("/")[-1]
	for i in range(8):
		try:
			os.mkdir("Results/base-" + str(i))
			print("Making new base folder here: Results/base-" + str(i))
		except OSError as error:
			print(error)

	## Open skempi csv file and start reading from base1 and split into folders from PDB and subfolders for index
	temp_dir = ""
	idx_pdb_dict = {}
	base_idx = 0
	with open('../../SKEMPI_dataset_download/skempi_v2.csv', 'r') as csv_file:
		reader = csv.reader(csv_file, delimiter=";")
		for index, row in enumerate(reader):
			## Skip header
			if index == 0:
				continue
			
			## Get the pdb from each row - keep track of index as well
			pdb = row[0].split("_")[0]
			if pdb != temp_dir:
				temp_dir = pdb
				make_output_dir(pdb, "base-" + str(base_idx))
	
			if index % 1000 == 0:
				base_idx += 1
				make_output_dir(pdb, "base-" + str(base_idx))

			## Make a folder for each index in the associated PDB
			try:
				formatted_index = f"{(index):05d}"
				# print(formatted_index)
				os.mkdir("Results/" + "base-" + str(base_idx) + "/" + pdb + "/" + formatted_index)
			except OSError as error:
				# print(error)
				pass

			## Keep track of idx to pdb mapping
			idx_pdb_dict[formatted_index] = "base-" + str(base_idx) + "/" + pdb + "/" + formatted_index

			# if count>1000:
			# 	break
			# else:
			# 	count += 1

	return idx_pdb_dict


def progressbar(it, prefix="", size=60, out=sys.stdout):  # Python3.3+
	"""Prints out a progress bar in stdout in for loops

	Args:
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	it (_type_): _description_
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	prefix (str, optional): _description_. Defaults to "".
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	size (int, optional): _description_. Defaults to 60.
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	out (_type_, optional): _description_. Defaults to sys.stdout.

	Yields:
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	_type_: _description_
	"""
	count = len(it)

	def show(j):
		x = int(size * j / count)
		print(
			"{}[{}{}] {}/{}".format(prefix, "#" * x, "." * (size - x), j, count),
			end="\r",
			file=out,
			flush=True,
		)

	show(0)
	for i, item in enumerate(it):
		yield item
		show(i + 1)
	print("\n", flush=True, file=out)

def calculate_diffbond(file1, file2, output):
	## Calculates diffbond using os.system() command
	os.system(
			"python DiffBond_v2.py -i " + file1 + " " + file2 +" -m h -o " + output
		)

def base_diffbond_calc(idx_pdb_dict, path):
	folders = os.listdir(path)
	# print(folders)
	index = 0
	for folder in progressbar(folders):
		
		# Use list comprehension to filter for folders starting with "H"
		dir = os.listdir(path + "/" + folder)
		halfs = []
		for d in dir:
			if d[0] == "H" and len(d) >= 4:
				halfs.append(d)
				print(d)
		print(folder)
		print(idx_pdb_dict[folder])

		try:
			calculate_diffbond(str(path + "/" + folder + "/" + halfs[0] + "/hydrogensAdded.pdb"), str(path + "/" + folder + "/" + halfs[1] + "/hydrogensAdded.pdb"), idx_pdb_dict[folder])
			# calculate_diffbond(str(path + "/" + folder + "/" + halfs[0] + "/" +"final_half1.pdb"), str(path + "/" + folder + "/" + halfs[1] + "/" +"final_half2.pdb"), idx_pdb_dict[folder])
		except Exception as error:
			print(error)
			# pass
		
		index = index +1
		if index > 1000:
			break


def wt_diffbond_calc(wt_path):
	folders = os.listdir(wt_path)
	try:
		os.mkdir("Results/wt")
		print("Making new base folder here: Results/wt")
	except OSError as error:
		print(error)

	for folder in folders:
		pdb = folder.split("/")[-1]
		# print(pdb)
		make_output_dir(pdb,"wt")

		halfs = os.listdir(wt_path + "/" + folder)
		half1 = ""
		half2 = ""
		for half in halfs:
			if half[:2] == "H1":
				half1 = half
			if half[:2] == "H2":
				half2 = half
	
		try:
			calculate_diffbond(str(wt_path + "/" + pdb + "/" + half1 + "/hydrogensAdded.pdb"), str(wt_path + "/" + pdb + "/" + half2 + "/hydrogensAdded.pdb"), "wt/" + pdb)
		except Exception as error:
			# print(error)
			with open("error_log.txt", "a") as error_file: 
				error_file.write(error)





def main():
	reader = parse_skempi()
	
	# base0 = "../../SKEMPI_dataset_download/base-0"
	# base_diffbond_calc(reader, base0)
 
	# base1 = "../../SKEMPI_dataset_download/base-1"
	# base_diffbond_calc(reader, base1)
 
	# base2 = "../../SKEMPI_dataset_download/base-2"
	# base_diffbond_calc(reader, base2)

	# base3 = "../../SKEMPI_dataset_download/base-3"
	# base_diffbond_calc(reader, base3)
 
	# base4 = "../../SKEMPI_dataset_download/base-4"
	# base_diffbond_calc(reader, base4)
 
	# base5 = "../../SKEMPI_dataset_download/base-5"
	# base_diffbond_calc(reader, base5)

	wt_path = "../../SKEMPI_dataset_download/wt"
	wt_diffbond_calc(wt_path)




if __name__ == "__main__":
	main()
