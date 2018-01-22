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


def verify_download(data_sets):
    location = "../raw_data/"
    urls = {
        "human": ("https://storage.googleapis.com/simplicial-complex-dataset"
                  "/PPI%20Dataset/BIOGRID-Homosapien.csv"),
        "yeast": ("https://storage.googleapis.com/simplicial-complex-dataset"
                  "/PPI%20Dataset/BIOGRID-Saccharomyces-cerevisiae-"
                  "(bakers_yeast).csv")
    }
    for data_set in data_sets:
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
    data_sets = args["-g"].split(",")
    verify_download(data_sets)
    print("\nDownloads verified.")


if __name__ == '__main__':
    from sys import argv
    args = getopts(argv)
    main(args)
