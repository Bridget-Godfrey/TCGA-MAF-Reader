import sqlite3
import csv
import re
import os
import pickle
import time
import multiprocessing

con = sqlite3.connect('mutationData.db') #mutatuions3.db
csr = con.cursor()
con.execute("DROP TABLE IF EXISTS mutations")
con.commit()

desiredFields = ["Hugo_Symbol", "UniProt_AApos", "Entrez_Gene_Id", "SwissProt_entry_Id", "i_COSMIC_Gene", "Variant_Classification", "Variant_Type", "Protein_Change", "i_dbNSFP_CADD_phred", "i_dbNSFP_CADD_raw", "i_dbNSFP_CADD_raw_rankscore", "i_dbNSFP_ESP6500_AA_AF", "i_dbNSFP_ESP6500_EA_AF", "i_dbNSFP_Ensembl_geneid", "i_dbNSFP_Ensembl_transcriptid", "i_dbNSFP_FATHMM_pred", "i_dbNSFP_FATHMM_rankscore", "i_dbNSFP_FATHMM_score", "i_dbNSFP_GERPpp_NR", "i_dbNSFP_GERPpp_RS", "i_dbNSFP_GERPpp_RS_rankscore", "i_dbNSFP_Interpro_domain", "i_dbNSFP_LRT_Omega", "i_dbNSFP_LRT_converted_rankscore", "i_dbNSFP_LRT_pred", "i_dbNSFP_LRT_score", "i_dbNSFP_LR_pred", "i_dbNSFP_LR_rankscore", "i_dbNSFP_LR_score", "i_dbNSFP_MutationAssessor_pred", "i_dbNSFP_MutationAssessor_rankscore", "i_dbNSFP_MutationAssessor_score", "i_dbNSFP_MutationTaster_converted_rankscore", "i_dbNSFP_MutationTaster_pred", "i_dbNSFP_MutationTaster_score", "i_dbNSFP_Polyphen2_HDIV_pred", "i_dbNSFP_Polyphen2_HDIV_rankscore", "i_dbNSFP_Polyphen2_HDIV_score", "i_dbNSFP_Polyphen2_HVAR_pred", "i_dbNSFP_Polyphen2_HVAR_rankscore", "i_dbNSFP_Polyphen2_HVAR_score", "i_dbNSFP_RadialSVM_pred", "i_dbNSFP_RadialSVM_rankscore", "i_dbNSFP_RadialSVM_score", "i_dbNSFP_Reliability_index", "i_dbNSFP_SIFT_converted_rankscore", "i_dbNSFP_SIFT_pred", "i_dbNSFP_SIFT_score", "i_dbNSFP_SLR_test_statistic", "i_dbNSFP_SiPhy_29way_logOdds", "i_dbNSFP_SiPhy_29way_logOdds_rankscore", "i_dbNSFP_SiPhy_29way_pi", "i_dbNSFP_UniSNP_ids"]
proteinChangeIndex = 7
casesIndex = len(desiredFields)
ctypeIndex = len(desiredFields) + 1
otherColumns = ["Cases", "Cancer_Types"]

columns = ""
columns2 = ", ".join(desiredFields) + ", "+ ", ".join(otherColumns)
for cn in desiredFields:
	columns += cn + " TEXT, "
for cn in otherColumns:
	columns += cn + " TEXT, "

con.execute("CREATE TABLE mutations ( mutation_id INTEGER UNIQUE, " + columns[:] + " PRIMARY KEY(mutation_id))")
print("\nCREATE TABLE mutations ( mutation_id INTEGER UNIQUE, " + columns[:] + " PRIMARY KEY(mutation_id))")
con.commit()


mutations = {}



def readFile(filename, cancerType, tcgaID):
	global desiredFields
	global con
	global csr
	global mutations
	tsv_file = open(filename, encoding= 'utf-8',errors='ignore' )
	read_tsv = csv.reader(tsv_file, delimiter="\t")
	rowCount = 0
	meaningfulIndexes = []
	for i in range(len(desiredFields)):
		meaningfulIndexes.append(-1)
	for row in read_tsv:
		rowCount += 1
		if rowCount >= 4:
			if rowCount == 4:
				# GET COLUMN NUMBERS HERE
				onIndex = 0

				while onIndex < len(desiredFields): #If so skip it and try again
					for i in range(len(row)):
						if row[i].replace("+", "p") == desiredFields[onIndex]:
							meaningfulIndexes[onIndex] = i
							break
					onIndex += 1

			else:
				# GET ROW DATA HERE
				fields = []
				for n in meaningfulIndexes: 
					val = row[n]
					if n == -1:
						val = "'--"
					fields.append(val)
				
				# if rowCount%100 == 0:
				# 	print(fields, "\n", meaningfulIndexes)
				mutName = fields[0] + " " + fields[proteinChangeIndex]

				try:
					mutations[mutName][casesIndex].append(tcgaID)

					# print(mutName, "added", tcgaID, end= "  ")
					if mutations[mutName][ctypeIndex][len(mutations[mutName][ctypeIndex])-1] != cancerType:
						mutations[mutName][ctypeIndex].append(cancerType)
						# print(cancerType)
					# else:
					# 	print("")

				except:
					fields.append([tcgaID])
					fields.append([cancerType])
					mutations.update({mutName : fields})
					# print("Added  ", mutName)
				
			pass
		pass
	pass
#EOF READFILE






def dumpToDB (startAt, stopAt, mutations, columns2):
	# global con
	# global csr
	# global columns2
	listOfCmds = []
	onKeyNum = 0
	severalCommands = []
	template = 'INSERT INTO mutations ( '+ columns2  +') VALUES ('
	onePercent = len(mutations.keys())//800
	lst = list(mutations.keys())
	for j in range(startAt, stopAt, 1):
		k = lst[j]
		onKeyNum += 1
		temp1 = ", ".join(mutations[k][len(mutations[k])-2])
		temp2 = ", ".join(mutations[k][len(mutations[k])-1])
		temp0 = ", ".join(mutations[k][:-2])

		insertCMD = 'INSERT INTO mutations ( ' + ", ".join(desiredFields) + ', Cases, Cancer_Types ) VALUES ( ' + "?, "*len(desiredFields) + '?, ?)'
		
		mtouple = []
		for i in range(len(desiredFields)):
			mtouple.append(mutations[k][i])
		mtouple.append(temp1)
		mtouple.append(temp2)
		severalCommands.append(tuple(mtouple))

		if onKeyNum%377 == 0:
			con = sqlite3.connect('mutationData.db', timeout = 30.0)
			con.executemany(insertCMD, severalCommands)
			con.commit()
			con.close()
			print("█", end = "")
			severalCommands = []
	con = sqlite3.connect('mutationData.db', timeout = 30.0)
	con.executemany(insertCMD, severalCommands)
	con.commit()
	con.close()
	severalCommands = []
	# print("")


# readFile("MAF/ACC/TCGA-OR-A5J1-01.hg19.oncotator.hugo_entrez_remapped.maf.txt", "ACC", "TCGA-OR-A5J1")

def handleCTs (num):
	cancerTypes = ["ACC", "BRCA", "CHOL", "DLBC", "GBM", "KICH", "KIRP", "LGG", "LUAD", "OV", "PCPG", "READ", "SKCM", "TGCT", "THYM", "UCS", "BLCA", "CESC", "COAD", "ESCA", "HNSC", "KIRC", "LAML", "LIHC", "LUSC", "PAAD", "PRAD", "SARC", "STAD", "THCA", "UCEC", "UVM"]
	types = cancerTypes[:16]
	if num == 1:
		types = cancerTypes[16:]
	for ct in types:
		print("-----------------" + ct + "-----------------------")
		for filename in os.listdir("MAF/" + ct):
			if filename.endswith(".txt") or filename.endswith(".maf"): 
				# print(filename)
				readFile("MAF/" + ct + "/" + filename, ct, filename[0:12])
				pass
			else:
				pass
csr.close()
con.close()


if __name__ == '__main__':
	multiprocessing.freeze_support()
	handleCTs(0)
	handleCTs(1)
		
	for q in range(200):
		print("░", end = "")
	print("")
	qrtr = len(mutations.keys())//4
	print(qrtr)
	t1 = multiprocessing.Process(target=dumpToDB, args=(0, qrtr, mutations, columns2,))
	t1.start()
	t2 = multiprocessing.Process(target=dumpToDB, args=(qrtr, 2*qrtr, mutations, columns2,))
	t2.start()
	t3 = multiprocessing.Process(target=dumpToDB, args=(2*qrtr, 3*qrtr, mutations, columns2,))
	t3.start()
	t4 = multiprocessing.Process(target=dumpToDB, args=(3*qrtr, len(mutations.keys()), mutations, columns2,))
	t4.start()
	t1.join()
	t2.join()
	t3.join()
	t4.join()