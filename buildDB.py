import sqlite3
import csv
import re
import os
con = sqlite3.connect('mutations2.db')
csr = con.cursor()
hugo_Symbol_At = 0
chromosome_At = 0
con.execute("DROP TABLE IF EXISTS mutations")
con.commit()
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
for c in chromosomes:
	con.execute("DROP TABLE IF EXISTS chrom"+ c + "genes")
	con.commit()
	con.execute("DROP TABLE IF EXISTS chrom"+ c + "mutations")
	con.commit()
	con.execute("DROP TABLE IF EXISTS genechroms")
	con.commit()




chromGeneCols = ["hugo_symbol", "UniProt_ID", "Entrez_Gene_Id", "SwissProt_entry_Id", "i_COSMIC_Gene"]
chromGeneCols2 = ["hugo_symbol TEXT", "UniProt_ID TEXT", "Entrez_Gene_Id INT", "SwissProt_entry_Id TEXT", "i_COSMIC_Gene TEXT"]
chromCols = ", ".join(chromGeneCols2)
chromCols2 = ", ".join(chromGeneCols)


isFirst= True
isSecond = True
firstTime = True



ON_CHROMOSOME = {}


def testFunc1 (filename):
	global con
	global csr
	global hugo_Symbol_At
	global chromosome_At
	global chromGeneCols
	global chromGeneCols2
	global chromCols
	global chromCols2
	global isFirst
	global isSecond 
	global firstTime
	tsv_file = open(filename, encoding= 'utf-8',errors='ignore' )

	read_tsv = csv.reader(tsv_file, delimiter="\t")
	columnNames = []
	columnInfo =[]



	startAt = 1
	stopAt = 999
	printPadding = 35

	
	onRow = 0
	onRow2 = 0
	for row in read_tsv:
		onRow +=1
		if len(row) >10 and isFirst:
			isFirst = False
			for i in range(len(row)):
				columnNames.append(row[i])
				
		
		elif len(row) >10 and isSecond:
			isSecond = False
			columnInfo.append([])
			onRow2 = 0
			for i in range(len(row)):
				tmpRowEntry = re.sub('"', "", row[i])
				
				columnInfo[onRow2].append(tmpRowEntry)
				

			
		elif len(row) >10 and onRow < stopAt:
			
			columnInfo.append([])
			for i in range(len(row)):
				tmpRowEntry = re.sub('"', "", row[i])
				print(onRow2, tmpRowEntry )
				columnInfo[onRow2].append(tmpRowEntry)
			onRow2 += 1
	dataType = ["TEXT", "INT", "TEXT", "INT", "TEXT", "INT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "INT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT"]

	columns = ""
	usedNames = []

	helpful = []
	for i in range(len(columnNames)):
		tmpCN = re.sub("\(.*\)", "", columnNames[i])
		tmpCN = re.sub(" +", "_", tmpCN)
		tmpCN = re.sub("\+", "Plus", tmpCN)
		tmpCN = re.sub("-", "Minus", tmpCN)

		if tmpCN.upper() in usedNames:
			print("DUPE AT ", tmpCN)
		else:
			print(i, tmpCN, dataType[i])
			try:
				columns +=  tmpCN + ", "
			except:
				columns +=  tmpCN + ", "

			usedNames.append(tmpCN.upper())
			helpful.append([tmpCN, i, len(helpful)])
	if firstTime:
		con.execute("CREATE TABLE mutations ( mutation_id INTEGER UNIQUE, " + columns[:] + " PRIMARY KEY(mutation_id))")
		print("\nCREATE TABLE mutations ( mutation_id INTEGER UNIQUE, " + columns[:] + " PRIMARY KEY(mutation_id))")
		con.commit()
		con.execute("CREATE TABLE genechroms (" + chromCols + ", chromosome TEXT)")
		print("\n\n\nCREATE TABLE genechroms (" + chromCols + ", chromosome TEXT)")
		con.commit()
		for c in chromosomes:

			con.execute("CREATE TABLE chrom"+ c +"genes (" + chromCols + ")")
			print("\n\n\nCREATE TABLE chrom"+ c +"genes (" + chromCols + ")")
			con.commit()
			con.execute("CREATE TABLE chrom"+ c +"mutations ( mutation_id INTEGER UNIQUE, " + columns[:-2] + " )")
			print("\n\n\nCREATE TABLE chrom"+ c +"mutations ( mutation_id INTEGER UNIQUE, " + columns[:-2] + " )")
			
			con.commit()
		firstTime = False


	onID = startAt

	update_A = "UPDATE mutations SET "
	update_B = " WHERE mutation_id = "
	for i in range( startAt-1, len(columnInfo)):
		
		
		
		vals = ""
		try:
			currentChrom = ON_CHROMOSOME[columnInfo[i][0]][0]
			if currentChrom == None:
				currentChrom = columnInfo[i][4]
				ON_CHROMOSOME[columnInfo[i][0]] = [columnInfo[i][4],  columnInfo[i][194], columnInfo[i][1], columnInfo[i][45], columnInfo[i][117]]
			else:
				pass
		except:
			ON_CHROMOSOME.update({columnInfo[i][0]: [columnInfo[i][4],  columnInfo[i][194], columnInfo[i][1], columnInfo[i][45], columnInfo[i][117]]})
			currentChrom = columnInfo[i][4]
			# print(columnInfo[i][0], "on", columnInfo[i][4])
		
		
		
		# print(columnInfo[i][0])
		for j in range(len(helpful)):
			val = ""
			
			
			try:
				if dataType[helpful[j][1]] != "INT":
					
					val = '"' + columnInfo[i][helpful[j][1]] + '"'
				else:
					
					val = columnInfo[i][helpful[j][1]]
					if val == None  or val == "":
						val = "NULL"
			except:
				try:
					val = '"' + columnInfo[i][helpful[j][1]] + '"'
				except:
					
					val = '""'


			
			vals += val + ", "
			
			
		onID += 1
		
		instCmd = "INSERT INTO mutations (" + columns[:-2] +") VALUES (" + vals[:-2] + ")"
		csr.execute(instCmd)
		lastRow = csr.lastrowid
		instCmd2 = "INSERT INTO chrom" + currentChrom + "mutations (mutation_id, " + columns[:-2] +") VALUES (" + str(lastRow)+  ", " + vals[:-2] + ")"
		
		csr.execute(instCmd2)
		
		
		
		
		
		

	for k in ON_CHROMOSOME.keys():
		# print(k, ON_CHROMOSOME[k])
		instChromCmd = 'INSERT INTO chrom' + ON_CHROMOSOME[k][0] +  'genes ( '+ chromCols2 +' ) VALUES ("' + k + '", "' + ON_CHROMOSOME[k][1] + '", ' + ON_CHROMOSOME[k][2] + ', "' + ON_CHROMOSOME[k][3] + '", "' + ON_CHROMOSOME[k][4] + '" )'
		con.execute(instChromCmd)
		con.commit()
		instChromCmd = 'INSERT INTO genechroms ( '+ chromCols2 +', chromosome ) VALUES ("' + k + '", "' + ON_CHROMOSOME[k][1] + '", ' + ON_CHROMOSOME[k][2] + ', "' + ON_CHROMOSOME[k][3] + '", "' + ON_CHROMOSOME[k][4] + '", "' + ON_CHROMOSOME[k][0] + '" )'
		con.execute(instChromCmd)
		con.commit()


CANCER_TYPE = "ACC"

for filename in os.listdir("MAF/" + CANCER_TYPE):
	if filename.endswith(".txt") or filename.endswith(".maf"): 
		print(filename)
		testFunc1("MAF/" + CANCER_TYPE + "/" + filename)
		pass
	else:
		pass


