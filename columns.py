import csv
tsv_file = open("test.maf.txt", encoding= 'utf-8',errors='ignore' )
#print(csv.list_dialects())
read_tsv = csv.reader(tsv_file, delimiter="\t")# , dialect = "excel-tab", encoding = "utf-8")
columnNames = []
columnInfo =[]
# columnsOfInterest = [4, 7, 0, 5, 6, 9, 34, 8, 39, 40, 41, 10, 11, 12, 13, 47, 48, 49, 50, 194]
columnsOfInterest = [4, 7, 0, 5, 6, 9, 34, 8, 39, 40, 41, 10, 11, 12, 13, 47, 48, 49, 50, 194]
mutIdentifiers = [34, 1, 0, 39, 40, 41, 13, 194, 48]
stopAt = 20
printPadding = 35
#read_tsv = read_tsv.encode('')
isFirst= True
isSecond = True
onRow = 0
for row in read_tsv:
	onRow +=1
	if len(row) >10 and isFirst:
		isFirst = False
		for i in range(len(row)):
			columnNames.append(row[i])
			#print(row[i], "\n")
	
	elif len(row) >10 and isSecond:
		isSecond = False
		for i in range(len(row)):
			columnInfo.append([row[i]])
		#print(columnInfo)
	elif len(row) >10 and onRow < stopAt:
		for i in range(len(row)):
			columnInfo[i].append(row[i])


def printValueAt(c, r):
	if (len(columnInfo[c][r]) >= 1):
		print("    "+str(r) + ".", columnInfo[c][r])
	else:
		print("    "+str(r) + ".", "----")

def printRow(r, colList = None):
	if (colList == None):
		colList = columnsOfInterest

	if (r == -1):
		for i in colList:
			print((columnNames[i] + ":").ljust(printPadding), end = "")
	else:
		for i in colList:
			if (len(columnInfo[i][r]) >= 1):
				print((columnInfo[i][r] + ",").ljust(printPadding), end = "")
			else:
				print("----".ljust(printPadding), end = "")
	print("\n\n")
		


for j in columnsOfInterest:
	print(str(j) + ". ", columnNames[j] + ":")
	printValueAt(j, 0)
	printValueAt(j, 1)
	printValueAt(j, 2)
	printValueAt(j, 3)
	printValueAt(j, 4)
	printValueAt(j, 5)
	print("\n")


# printRow(-1)
# printRow(0)
# printRow(1)
# printRow(2)
# printRow(3)