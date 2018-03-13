import csv
import json
import os
from tqdm import tqdm  # progress bar
from filenames import RAW_HUMAN_PPI, RAW_YEAST_PPI, RAW_HUMAN_COMPLEX, \
    RAW_HUMAN_COMPLEX_JSON, RAW_YEAST_COMPLEX, GENE_ID_CONVERSION, \
    CLEANED_LOCATION, LOCATION

header = [
    '#BioGRID Interaction ID',
    'Entrez Gene Interactor A',
    'Entrez Gene Interactor B',
    'BioGRID ID Interactor A',
    'BioGRID ID Interactor B',
    'Systematic Name Interactor A',
    'Systematic Name Interactor B',
    'Official Symbol Interactor A',
    'Official Symbol Interactor B',
    'Synonyms Interactor A',
    'Synonyms Interactor B',
    'Experimental System',
    'Experimental System Type',
    'Author',
    'Pubmed ID',
    'Organism Interactor A',
    'Organism Interactor B',
    'Throughput',
    'Score',
    'Modification',
    'Phenotypes',
    'Qualifications',
    'Tags',
    'Source Database'
]

#  Should check if the experimental systems are the same for human and yeast
expSystem = [
    'Two-hybrid',
    'Affinity Capture-Luminescence',
    'Affinity Capture-MS',
    'Affinity Capture-RNA',
    'Affinity Capture-Western'
]

expSysType = ['physical']


def clean_human_ppi_data():
    organism = ['9606']
    with open(
        LOCATION + RAW_HUMAN_PPI, 'r') as inp, open(
            CLEANED_LOCATION + RAW_HUMAN_COMPLEX, 'w', newline='') as out:
        writer = csv.writer(out)
        reader = csv.reader(inp)
        writer.writerow(header)
        for row in tqdm(reader):
            if (
                row[12] in expSysType
            ) and (
                row[15] in organism
            ) and (
                row[16] in organism
            ) and (
                row[11] in expSystem
            ):
                writer.writerow(row)


def clean_yeast_ppi_data():
    organism = ['559292']
    with open(
        LOCATION + RAW_YEAST_PPI, 'r') as inp, open(
            CLEANED_LOCATION + RAW_YEAST_PPI, 'w', newline='') as out:
        writer = csv.writer(out)
        reader = csv.reader(inp)
        writer.writerow(header)
        for row in tqdm(reader):
            if (
                row[11] in expSystem
            ) and (
                row[12] in expSysType
            ) and (
                row[15] in organism
            ) and (
                row[16] in organism
            ):
                writer.writerow(row)

def clean_human_complex_data():
    data = json.loads(open(CLEANED_LOCATION + RAW_HUMAN_COMPLEX_JSON).read())
    valid_complexes = []
    # filters the mouses (non humans)
    for each in data:
        organisms = set(each["SWISSPROT organism"].split(";"))
        if len(organisms) == 1 and organisms.pop() == "Homo sapiens (Human)":
            valid_complexes.append(each)


    gene_to_entrez = json.loads(open(LOCATION + GENE_ID_CONVERSION).read())
    to_csv = []
    # create edges | [Gene Entrez ID, Complex Name]
    with open(CLEANED_LOCATION + RAW_HUMAN_COMPLEX, "w") as out:
        writer = csv.writer(out)
        writer.writerow(["Name", "Complex"])
        for each in tqdm(valid_complexes):
            genes = each["subunits(Gene name)"].split(";")
            for gene in genes:
                if gene in gene_to_entrez:
                    writer.writerow([gene_to_entrez[gene], each["ComplexName"]])


def check_files():
    for each in [
        RAW_HUMAN_PPI,
        RAW_YEAST_PPI,
        RAW_HUMAN_COMPLEX,
        RAW_YEAST_COMPLEX
    ]:
        if not os.path.isfile(CLEANED_LOCATION + each):
            return False
    return True


def clean_data():
    print("Cleaning BioGrid Human data...")
    clean_human_ppi_data()
    print("Cleaning BioGrid Yeast data...")
    clean_yeast_ppi_data()
    print("Cleaning Corum Human data...")
    clean_human_complex_data()
