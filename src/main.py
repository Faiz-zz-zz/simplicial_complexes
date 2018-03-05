"""

Main source file.

Network selection is based on the argument provided. Data downloads are
.gitignore and checked and verified on all initial invocation of a session.
This process is only used to generate graphs
and precompute measures that are algorithmic and computationally heavy.

"""
import os.path
import requests
import sys

from clean_data import clean_data, check_files
import simplicialComplexNetwork as scn


def verify_download():
    location = "../raw_data/"
    urls = {
        "human": ("https://storage.googleapis.com/simplicial-complex-dataset"
                  "/PPI%20Dataset/BIOGRID-Homosapien.csv"),
        "yeast": ("https://storage.googleapis.com/simplicial-complex-dataset"
                  "/PPI%20Dataset/BIOGRID-Saccharomyces-cerevisiae-"
                  "(bakers_yeast).csv")
    }
    for data_set in urls:
        file_name = location + data_set + ".csv"

        if os.path.isfile(file_name):
            print("{}.csv exists, skipping download process.".format(data_set))
            return

        with open(file_name, "wb") as f:
            print("Downloading {}".format(file_name))
            response = requests.get(urls[data_set], stream=True)
            total_length = response.headers.get('content-length')

            if total_length is None:
                f.write(response.content)
            else:
                dl = 0
                total_length = int(total_length)
                for data in response.iter_content(chunk_size=4096):
                    dl += len(data)
                    f.write(data)
                    done = int(50 * dl / total_length)
                    sys.stdout.write(
                        "\r[%s%s]" % ('=' * done, ' ' * (50 - done)))
                    sys.stdout.flush()


def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
        argv = argv[1:]
    return opts


def main(args):
    if not check_files():
        verify_download()
        print("\nDownloads verified.")
        print("Files Downloaded.\nCleaning Datasets...")
        clean_data()
    # human_ppi, yeast_ppi, human_complex, yeast_complex
    data_set = args["-d"].split(",")
    if len(data_set) > 1:
        raise Exception("Only one data set a time allowed")

    methods = args["-m"].split(",")

    if data_set[-3:] == "ppi":
        graph.get_metrics(methods, data_set)
    else:
        scn.generate_metrics(methods, data_set)


if __name__ == '__main__':
    from sys import argv
    args = getopts(argv)
    main(args)
