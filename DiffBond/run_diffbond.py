import os
from pathlib import Path
import logging
import csv
from tqdm import tqdm
import subprocess
import concurrent.futures
import multiprocessing


def make_output_dir(pdb, path):
    try:
        (Path("Results") / path / pdb).mkdir(exist_ok=True)
    except Exception:
        pass

def parse_skempi():
    # Make a skempi directory for current base
    for i in range(8):
        base_dir = Path("results") / f"base-{i}"
        base_dir.mkdir(exist_ok=True)
        print(f"Making new base folder here: {base_dir}")

    # Open skempi csv file and start reading from base1 and split into folders from PDB and subfolders for index
    idx_pdb_dict = {}
    base_idx = 0
    current_pdb = ""

    with open("../../SKEMPI_dataset_download/skempi_v2.csv", "r") as csv_file:
        for index, row in enumerate(csv.reader(csv_file, delimiter=";")):
            if index == 0: # Skip header
                continue

            pdb = row[0].split("_")[0]
            if pdb != current_pdb:
                current_pdb = pdb
                make_output_dir(pdb, "base-" + str(base_idx))

            if index % 1000 == 0:
                base_idx += 1
                make_output_dir(pdb, "base-" + str(base_idx))

            formatted_index = f"{(index):05d}"
            (Path("results") / f"base-{base_idx}" / pdb / formatted_index).mkdir(exist_ok=True)
            idx_pdb_dict[formatted_index] = f"base-{base_idx}/{pdb}/{formatted_index}"

    return idx_pdb_dict


def calculate_diffbond(file1, file2, output):
    # command = ["python", "DiffBond_v2.py", "-i", file1, file2, "-m", "c", "i", "h" ,"s", "-g", "-o", output]
    command = ["python", "DiffBond_v2.py", "-i", file1, file2, "-m","i" ,"-o", output]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def parallel_process_folder(folder, path, idx_pdb_dict, is_wt=False):
    if os.path.isfile(path + "/" + folder):
        logging.warning(f"{folder} is a file and not a directory.")
        return

    halfs = [d for d in os.listdir(path + "/" + folder) if d[0] == "H" and len(d) >= 4]

    if len(halfs) != 2:
        logging.error(f"{str(folder)} is missing pdb halfs.")
        return

    halfs.sort(reverse=True)

    try:
        if is_wt:
            pdb = folder.split("/")[-1]
            file1 = f"{path}/{folder}/{halfs[0]}/hydrogensAdded.pdb"
            file2 = f"{path}/{folder}/{halfs[1]}/hydrogensAdded.pdb"
            output = f"wt/{pdb}"
        else:
            file1 = f"{path}/{folder}/{halfs[0]}/half1.pdb"
            file2 = f"{path}/{folder}/{halfs[1]}/half2.pdb"
            output = idx_pdb_dict[folder]

        calculate_diffbond(file1, file2, output)
    except Exception as error:
        logging.error("Error in calculating diffbond: %s", error)


def parallel_diffbond_calc(idx_pdb_dict, path, is_wt=False):
    if is_wt:
        try:
            os.makedirs("Results/wt", exist_ok=True)
            print("Making new base folder here: Results/wt")
        except OSError as error:
            print(error)

    folders = os.listdir(path)
    total_folders = len([folder for folder in folders if os.path.isdir(path + "/" + folder)])

    with tqdm(total=total_folders, desc="Processing folders") as progress_bar:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = []
            for index, folder in enumerate(folders):
                futures.append(executor.submit(parallel_process_folder, folder, path, idx_pdb_dict, is_wt))

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    logging.error("Generated an exception: %s", exc)
                finally:
                    progress_bar.update(1)

def main():
    reader = parse_skempi()

    # Process wild type
    wt_path = "datasets/wt"
    # wt_path = "../../MT_Processing_Archive/wt"

    parallel_diffbond_calc(reader, wt_path, is_wt=True)

    # Process base folders
    # for i in range(8):
    #     base_path = f"../../MT_Processing_Archive/base-{i}"
    #     parallel_diffbond_calc(reader, base_path)


if __name__ == "__main__":
    main()