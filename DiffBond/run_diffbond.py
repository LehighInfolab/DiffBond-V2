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
        try:
            (Path("Results") / ("base-" + str(i))).mkdir(exist_ok=True)
            print("Making new base folder here: Results/base-" + str(i))
        except Exception as error:
            print(error)

    # Open skempi csv file and start reading from base1 and split into folders from PDB and subfolders for index
    temp_dir = ""
    idx_pdb_dict = {}
    base_idx = 0
    with open("../../SKEMPI_dataset_download/skempi_v2.csv", "r") as csv_file:
        reader = csv.reader(csv_file, delimiter=";")
        for index, row in enumerate(reader):
            # Skip header
            if index == 0:
                continue

            # Get the pdb from each row - keep track of index as well
            pdb = row[0].split("_")[0]
            if pdb != temp_dir:
                temp_dir = pdb
                make_output_dir(pdb, "base-" + str(base_idx))

            if index % 1000 == 0:
                base_idx += 1
                make_output_dir(pdb, "base-" + str(base_idx))

            # Make a folder for each index in the associated PDB
            try:
                formatted_index = f"{(index):05d}"
                (Path("Results") / ("base-" + str(base_idx)) / pdb / formatted_index).mkdir(exist_ok=True)
            except OSError as error:
                print(error)
                pass

            # Keep track of idx to pdb mapping
            idx_pdb_dict[formatted_index] = "base-" + str(base_idx) + "/" + pdb + "/" + formatted_index

    return idx_pdb_dict


def calculate_diffbond(file1, file2, output):
    # Calculates diffbond using os.system() command
    command = ["python", "DiffBond_v2.py", "-i", file1, file2, "-m", "c", "i", "h" ,"s", "-g", "-o", output]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    # subprocess.run(command, check=True)


def base_diffbond_calc(idx_pdb_dict, path, start_idx):
    folders = os.listdir(path)
    for index, folder in tqdm(enumerate(folders)):
        if index < start_idx:
            continue

        # Use list comprehension to filter for folders starting with "H"
        if os.path.isfile(path + "/" + folder):
            logging.warning(f"{str(folder)} is a file and not a directory.")
            continue

        dir = os.listdir(path + "/" + folder)
        halfs = []
        for d in dir:
            if d[0] == "H" and len(d) >= 4:
                halfs.append(d)

        if len(halfs) != 2:
            logging.error(f"{str(folder)} is missing pdb halfs.")
            continue

        if halfs[0].startswith("H2"):
            temp = halfs[0]
            halfs[0] = halfs[1]
            halfs[1] = temp

        try:
            calculate_diffbond(
                str(path + "/" + folder + "/" + halfs[0] + "/half1.pdb"),
                str(path + "/" + folder + "/" + halfs[1] + "/half2.pdb"),
                idx_pdb_dict[folder],
            )
        except Exception as error:
            logging.error("Error in calculating diffbond:", error)
            continue


def parallel_process_folder(folder, path, idx_pdb_dict):
    if os.path.isfile(path + "/" + folder):
        logging.warning(f"{str(folder)} is a file and not a directory.")
        return

    dir = os.listdir(path + "/" + folder)
    halfs = [d for d in dir if d[0] == "H" and len(d) >= 4]

    if len(halfs) != 2:
        logging.error(f"{str(folder)} is missing pdb halfs.")
        return

    if halfs[0].startswith("H2"):
        halfs[0], halfs[1] = halfs[1], halfs[0]

    try:
        calculate_diffbond(
            str(path + "/" + folder + "/" + halfs[0] + "/half1.pdb"),
            str(path + "/" + folder + "/" + halfs[1] + "/half2.pdb"),
            idx_pdb_dict[folder],
        )
    except Exception as error:
        logging.error("Error in calculating diffbond: %s", error)


def parallel_base_diffbond_calc(idx_pdb_dict, path, start_idx=-1):
    folders = os.listdir(path)
    total_folders = len([folder for folder in folders if os.path.isdir(path + "/" + folder)]) - start_idx

    with tqdm(total=total_folders - start_idx, desc="Processing folders") as progress_bar:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = []
            for index, folder in enumerate(folders):
                # if index > 40:
                #     return
                if index < start_idx:
                    return
                futures.append(executor.submit(parallel_process_folder, folder, path, idx_pdb_dict))

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    logging.error("Generated an exception: %s", exc)
                finally:
                    progress_bar.update(1)


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
        make_output_dir(pdb, "wt")

        halfs = os.listdir(wt_path + "/" + folder)
        half1 = ""
        half2 = ""
        for half in halfs:
            if half[:2] == "H1":
                half1 = half
            if half[:2] == "H2":
                half2 = half

        try:
            calculate_diffbond(
                str(wt_path + "/" + pdb + "/" + half1 + "/hydrogensAdded.pdb"),
                str(wt_path + "/" + pdb + "/" + half2 + "/hydrogensAdded.pdb"),
                "wt/" + pdb,
            )
        except Exception as error:
            logging.error(f"{error}. Failed to calculate diffbond for {str(folder)}.")


def main():
    reader = parse_skempi()

    base0 = "../../MT_Processing_Archive/base-0"
    parallel_base_diffbond_calc(reader, base0)

    base1 = "../../MT_Processing_Archive/base-1"
    parallel_base_diffbond_calc(reader, base1)

    base2 = "../../MT_Processing_Archive/base-2"
    parallel_base_diffbond_calc(reader, base2)

    base3 = "../../MT_Processing_Archive/base-3"
    parallel_base_diffbond_calc(reader, base3)

    base4 = "../../MT_Processing_Archive/base-4"
    parallel_base_diffbond_calc(reader, base4)

    base5 = "../../MT_Processing_Archive/base-5"
    parallel_base_diffbond_calc(reader, base5)

    base6 = "../../MT_Processing_Archive/base-6"
    parallel_base_diffbond_calc(reader, base6)

    base7 = "../../MT_Processing_Archive/base-7"
    parallel_base_diffbond_calc(reader, base7)


if __name__ == "__main__":
    main()
