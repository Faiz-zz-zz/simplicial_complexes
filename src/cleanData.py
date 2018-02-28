import csv 
from tqdm import tqdm #progress bar

header = ['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'BioGRID ID Interactor A',	'BioGRID ID Interactor B', 'Systematic Name Interactor A', 'Systematic Name Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B', 'Synonyms Interactor A', 'Synonyms Interactor B', 'Experimental System', 'Experimental System Type', 'Author', 'Pubmed ID', 'Organism Interactor A', 'Organism Interactor B', 'Throughput',	'Score', 'Modification', 'Phenotypes', 'Qualifications', 'Tags', 'Source Database']
#Should check if the experimental systems are the same for human and yeast
expSystem = ['Two-hybrid','Affinity Capture-Luminescence','Affinity Capture-MS','Affinity Capture-RNA','Affinity Capture-Western']
expSysType = ['physical']

def cleanHumanData():
	organism = ['9606']

	with open('BIOGRID-Homosapien.csv', 'r') as inp, open('BIOGRID-Homosapien_UPDATED.csv', 'w', newline='') as out:
	    writer = csv.writer(out)
	    reader = csv.reader(inp)
	    writer.writerow(header)
	    for row in tqdm(reader):
	        if row[12] in expSysType and row[15] in organism and row[16] in organism and row[11] in expSystem:
	            writer.writerow(row)

def cleanYeastData():	
	organism = ['559292']

	with open('BIOGRID-Saccharomyces-cerevisiae-(bakers_yeast).csv', 'r') as inp, open('BIOGRID-Saccharomyces-cerevisiae-(bakers_yeast)_UPDATED.csv', 'w', newline='') as out:
	    writer = csv.writer(out)
	    reader = csv.reader(inp)
	    writer.writerow(header)
	    for row in tqdm(reader):
	        if (row[11] in expSystem) and (row[12] in expSysType) and (row[15] in organism) and (row[16] in organism):
	            writer.writerow(row)

cleanHumanData()
cleanYeastData()