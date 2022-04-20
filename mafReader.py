# mafReader.py
import sqlite3
import csv
import re
import addClinicalData2

from os import listdir
from os.path import isfile, join


con = sqlite3.connect('mutations.db')
csr = con.cursor()

# CANCER_TYPES = ["ACC"]

dataTypes = {"Hugo_Symbol": "TEXT", "Entrez_Gene_Id": "INT", "Center": "TEXT", "NCBI_Build": "INT", "Chromosome": "TEXT", "Start_position": "INT", "End_position": "INT", "Strand": "TEXT", "Variant_Classification": "TEXT", "Variant_Type": "TEXT", "Reference_Allele": "TEXT", "Tumor_Seq_Allele1": "TEXT", "Tumor_Seq_Allele2": "TEXT", "dbSNP_RS": "TEXT", "dbSNP_Val_Status": "TEXT", "Tumor_Sample_Barcode": "TEXT", "Matched_Norm_Sample_Barcode": "TEXT", "Match_Norm_Seq_Allele1": "TEXT", "Match_Norm_Seq_Allele2": "TEXT", "Tumor_Validation_Allele1": "TEXT", "Tumor_Validation_Allele2": "TEXT", "Match_Norm_Validation_Allele1": "TEXT", "Match_Norm_Validation_Allele2": "TEXT", "Verification_Status": "TEXT", "Validation_Status": "TEXT", "Mutation_Status": "TEXT", "Sequencing_Phase": "TEXT", "Sequence_Source": "TEXT", "Validation_Method": "TEXT", "Score": "TEXT", "BAM_file": "TEXT", "Sequencer": "TEXT", "Tumor_Sample_UUID": "TEXT", "Matched_Norm_Sample_UUID": "TEXT", "Genome_Change": "TEXT", "Annotation_Transcript": "TEXT", "Transcript_Strand": "TEXT", "Transcript_Exon": "TEXT", "Transcript_Position": "TEXT", "cDNA_Change": "TEXT", "Codon_Change": "TEXT", "Protein_Change": "TEXT", "Other_Transcripts": "TEXT", "Refseq_mRNA_Id": "TEXT", "Refseq_prot_Id": "TEXT", "SwissProt_acc_Id": "TEXT", "SwissProt_entry_Id": "TEXT", "Description": "TEXT", "UniProt_AApos": "INT", "UniProt_Region": "TEXT", "UniProt_Site": "TEXT", "UniProt_Natural_Variations": "TEXT", "UniProt_Experimental_Info": "TEXT", "GO_Biological_Process": "TEXT", "GO_Cellular_Component": "TEXT", "GO_Molecular_Function": "TEXT", "COSMIC_overlapping_mutations": "TEXT", "COSMIC_fusion_genes": "TEXT", "COSMIC_tissue_types_affected": "TEXT", "COSMIC_total_alterations_in_gene": "INT", "Tumorscape_Amplification_Peaks": "TEXT", "Tumorscape_Deletion_Peaks": "TEXT", "TCGAscape_Amplification_Peaks": "TEXT", "TCGAscape_Deletion_Peaks": "TEXT", "DrugBank": "TEXT", "ref_context": "TEXT", "gc_content": "TEXT", "CCLE_ONCOMAP_overlapping_mutations": "TEXT", "CCLE_ONCOMAP_total_mutations_in_gene": "TEXT", "CGC_Mutation_Type": "TEXT", "CGC_Translocation_Partner": "TEXT", "CGC_Tumor_Types_Somatic": "TEXT", "CGC_Tumor_Types_Germline": "TEXT", "CGC_Other_Diseases": "TEXT", "DNARepairGenes_Role": "TEXT", "FamilialCancerDatabase_Syndromes": "TEXT", "MUTSIG_Published_Results": "TEXT", "OREGANNO_ID": "TEXT", "OREGANNO_Values": "TEXT", "i_1000gp3_AA": "TEXT", "i_1000gp3_AC": "TEXT", "i_1000gp3_AF": "TEXT", "i_1000gp3_AFR_AF": "TEXT", "i_1000gp3_AMR_AF": "TEXT", "i_1000gp3_AN": "TEXT", "i_1000gp3_CIEND": "TEXT", "i_1000gp3_CIPOS": "TEXT", "i_1000gp3_CS": "TEXT", "i_1000gp3_DP": "TEXT", "i_1000gp3_EAS_AF": "TEXT", "i_1000gp3_END": "TEXT", "i_1000gp3_EUR_AF": "TEXT", "i_1000gp3_IMPRECISE": "TEXT", "i_1000gp3_MC": "TEXT", "i_1000gp3_MEINFO": "TEXT", "i_1000gp3_MEND": "TEXT", "i_1000gp3_MLEN": "TEXT", "i_1000gp3_MSTART": "TEXT", "i_1000gp3_NS": "TEXT", "i_1000gp3_SAS_AF": "TEXT", "i_1000gp3_SVLEN": "TEXT", "i_1000gp3_SVTYPE": "TEXT", "i_1000gp3_TSD": "TEXT", "i_AAChange": "TEXT", "i_ACHILLES_Lineage_Results_Top_Genes": "TEXT", "i_BAM_File": "TEXT", "i_CGC_Cancer_Germline_Mut": "TEXT", "i_CGC_Cancer_Molecular_Genetics": "TEXT", "i_CGC_Cancer_Somatic_Mut": "TEXT", "i_CGC_Cancer_Syndrome": "TEXT", "i_CGC_Chr": "TEXT", "i_CGC_Chr_Band": "TEXT", "i_CGC_GeneID": "TEXT", "i_CGC_Name": "TEXT", "i_CGC_Other_Germline_Mut": "TEXT", "i_CGC_Tissue_Type": "TEXT", "i_COSMIC_Codon": "TEXT", "i_COSMIC_Gene": "TEXT", "i_COSMIC_n_overlapping_mutations": "INT", "i_COSMIC_overlapping_mutation_descriptions": "TEXT", "i_COSMIC_overlapping_primary_sites": "TEXT", "i_ChromChange": "TEXT", "i_ClinVar_ASSEMBLY": "TEXT", "i_ClinVar_HGMD_ID": "TEXT", "i_ClinVar_SYM": "TEXT", "i_ClinVar_TYPE": "TEXT", "i_ClinVar_rs": "TEXT", "i_Drug_Target": "TEXT", "i_ESP_AA": "TEXT", "i_ESP_AAC": "TEXT", "i_ESP_AA_AC": "TEXT", "i_ESP_AA_AGE": "TEXT", "i_ESP_AA_GTC": "TEXT", "i_ESP_AvgAAsampleReadDepth": "TEXT", "i_ESP_AvgEAsampleReadDepth": "TEXT", "i_ESP_AvgSampleReadDepth": "TEXT", "i_ESP_CA": "TEXT", "i_ESP_CDP": "TEXT", "i_ESP_CG": "TEXT", "i_ESP_CP": "TEXT", "i_ESP_Chromosome": "TEXT", "i_ESP_DBSNP": "TEXT", "i_ESP_DP": "TEXT", "i_ESP_EA_AC": "TEXT", "i_ESP_EA_AGE": "TEXT", "i_ESP_EA_GTC": "TEXT", "i_ESP_EXOME_CHIP": "TEXT", "i_ESP_FG": "TEXT", "i_ESP_GL": "TEXT", "i_ESP_GM": "TEXT", "i_ESP_GS": "TEXT", "i_ESP_GTC": "TEXT", "i_ESP_GTS": "TEXT", "i_ESP_GWAS_PUBMED": "TEXT", "i_ESP_MAF": "TEXT", "i_ESP_PH": "TEXT", "i_ESP_PP": "TEXT", "i_ESP_Position": "TEXT", "i_ESP_TAC": "TEXT", "i_ESP_TotalAAsamplesCovered": "INT", "i_ESP_TotalEAsamplesCovered": "INT", "i_ESP_TotalSamplesCovered": "TEXT", "i_Ensembl_so_accession": "TEXT", "i_Ensembl_so_term": "TEXT", "i_Entrez_Gene_Id": "INT", "i_Exon": "TEXT", "i_Familial_Cancer_Genes_Reference": "TEXT", "i_Familial_Cancer_Genes_Synonym": "TEXT", "i_Genome_Plus_Minus_10_Bp": "TEXT", "i_HGNC_Accession_Numbers": "TEXT", "i_HGNC_CCDS_IDs": "TEXT", "i_HGNC_Chromosome": "TEXT", "i_HGNC_Date_Modified": "TEXT", "i_HGNC_Date_Name_Changed": "TEXT", "i_HGNC_Date_Symbol_Changed": "TEXT", "i_HGNC_Ensembl_Gene_ID": "TEXT", "i_HGNC_Ensembl_ID": "TEXT", "i_HGNC_Enzyme_IDs": "TEXT", "i_HGNC_Gene_family_description": "TEXT", "i_HGNC_HGNC_ID": "TEXT", "i_HGNC_Locus_Group": "TEXT", "i_HGNC_Locus_Type": "TEXT", "i_HGNC_Name_Synonyms": "TEXT", "i_HGNC_OMIM_ID": "TEXT", "i_HGNC_Previous_Names": "TEXT", "i_HGNC_Previous_Symbols": "TEXT", "i_HGNC_Primary_IDs": "TEXT", "i_HGNC_Pubmed_IDs": "TEXT", "i_HGNC_Record_Type": "TEXT", "i_HGNC_RefSeq": "TEXT", "i_HGNC_Secondary_IDs": "TEXT", "i_HGNC_Status": "TEXT", "i_HGNC_Synonyms": "TEXT", "i_HGNC_UCSC_ID": "TEXT", "i_HGNC_UniProt_ID": "TEXT", "i_HGNC_VEGA_IDs": "TEXT", "i_HGVS_coding_DNA_change": "TEXT", "i_HGVS_genomic_change": "TEXT", "i_HGVS_protein_change": "TEXT", "i_NTotCov": "TEXT", "i_NVarCov": "TEXT", "i_ORegAnno_bin": "TEXT", "i_TTotCov": "TEXT", "i_TVarCov": "TEXT", "i_Transcript_Id": "TEXT", "i_Trna_alt1": "TEXT", "i_Trna_alt2": "TEXT", "i_Trna_ref": "TEXT", "i_Trna_tot": "TEXT", "i_Trna_var": "TEXT", "i_UniProt_alt_uniprot_accessions": "TEXT", "i_Variant_Classification": "TEXT", "i_Variant_Type": "TEXT", "i_annotation_transcript": "TEXT", "i_build": "TEXT", "i_ccds_id": "TEXT", "i_dbNSFP_1000Gp1_AC": "TEXT", "i_dbNSFP_1000Gp1_AF": "TEXT", "i_dbNSFP_1000Gp1_AFR_AC": "TEXT", "i_dbNSFP_1000Gp1_AFR_AF": "TEXT", "i_dbNSFP_1000Gp1_AMR_AC": "TEXT", "i_dbNSFP_1000Gp1_AMR_AF": "TEXT", "i_dbNSFP_1000Gp1_ASN_AC": "TEXT", "i_dbNSFP_1000Gp1_ASN_AF": "TEXT", "i_dbNSFP_1000Gp1_EUR_AC": "TEXT", "i_dbNSFP_1000Gp1_EUR_AF": "TEXT", "i_dbNSFP_Ancestral_allele": "TEXT", "i_dbNSFP_CADD_phred": "TEXT", "i_dbNSFP_CADD_raw": "TEXT", "i_dbNSFP_CADD_raw_rankscore": "TEXT", "i_dbNSFP_ESP6500_AA_AF": "TEXT", "i_dbNSFP_ESP6500_EA_AF": "TEXT", "i_dbNSFP_Ensembl_geneid": "TEXT", "i_dbNSFP_Ensembl_transcriptid": "TEXT", "i_dbNSFP_FATHMM_pred": "TEXT", "i_dbNSFP_FATHMM_rankscore": "TEXT", "i_dbNSFP_FATHMM_score": "TEXT", "i_dbNSFP_GERPPlusPlus_NR": "TEXT", "i_dbNSFP_GERPPlusPlus_RS": "TEXT", "i_dbNSFP_GERPPlusPlus_RS_rankscore": "TEXT", "i_dbNSFP_Interpro_domain": "TEXT", "i_dbNSFP_LRT_Omega": "TEXT", "i_dbNSFP_LRT_converted_rankscore": "TEXT", "i_dbNSFP_LRT_pred": "TEXT", "i_dbNSFP_LRT_score": "TEXT", "i_dbNSFP_LR_pred": "TEXT", "i_dbNSFP_LR_rankscore": "TEXT", "i_dbNSFP_LR_score": "TEXT", "i_dbNSFP_MutationAssessor_pred": "TEXT", "i_dbNSFP_MutationAssessor_rankscore": "TEXT", "i_dbNSFP_MutationAssessor_score": "TEXT", "i_dbNSFP_MutationTaster_converted_rankscore": "TEXT", "i_dbNSFP_MutationTaster_pred": "TEXT", "i_dbNSFP_MutationTaster_score": "TEXT", "i_dbNSFP_Polyphen2_HDIV_pred": "TEXT", "i_dbNSFP_Polyphen2_HDIV_rankscore": "TEXT", "i_dbNSFP_Polyphen2_HDIV_score": "TEXT", "i_dbNSFP_Polyphen2_HVAR_pred": "TEXT", "i_dbNSFP_Polyphen2_HVAR_rankscore": "TEXT", "i_dbNSFP_Polyphen2_HVAR_score": "TEXT", "i_dbNSFP_RadialSVM_pred": "TEXT", "i_dbNSFP_RadialSVM_rankscore": "TEXT", "i_dbNSFP_RadialSVM_score": "TEXT", "i_dbNSFP_Reliability_index": "TEXT", "i_dbNSFP_SIFT_converted_rankscore": "TEXT", "i_dbNSFP_SIFT_pred": "TEXT", "i_dbNSFP_SIFT_score": "TEXT", "i_dbNSFP_SLR_test_statistic": "TEXT", "i_dbNSFP_SiPhy_29way_logOdds": "TEXT", "i_dbNSFP_SiPhy_29way_logOdds_rankscore": "TEXT", "i_dbNSFP_SiPhy_29way_pi": "TEXT", "i_dbNSFP_UniSNP_ids": "TEXT", "i_dbNSFP_Uniprot_aapos": "TEXT", "i_dbNSFP_Uniprot_acc": "TEXT", "i_dbNSFP_Uniprot_id": "TEXT", "i_dbNSFP_aaalt": "TEXT", "i_dbNSFP_aapos": "TEXT", "i_dbNSFP_aapos_FATHMM": "TEXT", "i_dbNSFP_aapos_SIFT": "TEXT", "i_dbNSFP_aaref": "TEXT", "i_dbNSFP_cds_strand": "TEXT", "i_dbNSFP_codonpos": "TEXT", "i_dbNSFP_foldMinusdegenerate": "TEXT", "i_dbNSFP_genename": "TEXT", "i_dbNSFP_hg18_pos": "TEXT", "i_dbNSFP_phastCons100way_vertebrate": "TEXT", "i_dbNSFP_phastCons100way_vertebrate_rankscore": "TEXT", "i_dbNSFP_phastCons46way_placental": "TEXT", "i_dbNSFP_phastCons46way_placental_rankscore": "TEXT", "i_dbNSFP_phastCons46way_primate": "TEXT", "i_dbNSFP_phastCons46way_primate_rankscore": "TEXT", "i_dbNSFP_phyloP100way_vertebrate": "TEXT", "i_dbNSFP_phyloP100way_vertebrate_rankscore": "TEXT", "i_dbNSFP_phyloP46way_placental": "TEXT", "i_dbNSFP_phyloP46way_placental_rankscore": "TEXT", "i_dbNSFP_phyloP46way_primate": "TEXT", "i_dbNSFP_phyloP46way_primate_rankscore": "TEXT", "i_dbNSFP_refcodon": "TEXT", "i_dbSNPPopFreq": "TEXT", "i_gc_content_full": "TEXT", "i_gencode_transcript_name": "TEXT", "i_gencode_transcript_status": "TEXT", "i_gencode_transcript_tags": "TEXT", "i_gencode_transcript_type": "TEXT", "i_gene_type": "TEXT", "i_havana_transcript": "TEXT", "i_refseq_mrna_id": "TEXT", "i_secondary_variant_classification": "TEXT"}
chromGeneCols = ["hugo_symbol", "UniProt_ID", "Entrez_Gene_Id", "SwissProt_entry_Id", "i_COSMIC_Gene"]
chromGeneCols2 = ["hugo_symbol TEXT", "UniProt_ID TEXT", "Entrez_Gene_Id INT", "SwissProt_entry_Id TEXT", "i_COSMIC_Gene TEXT"]
chromCols = ", ".join(chromGeneCols2)
chromCols2 = ", ".join(chromGeneCols)

startAt = 1
stopAt = 999
columnNames = []
columnInfo =[]

SEEN_MUTATIONS = {}


ROW_NAME_TO_ID = {}
ON_COLUMN = 0
CANCER_TYPE = "ACC"

def readMAF (mafName):
	global columnNames 
	global columnInfo
	global SEEN_MUTATIONS
	global ROW_NAME_TO_ID
	global CANCER_TYPE

	
	hasMuts = []
	ROW_NAME_TO_ID = {}
	onID = 1
	columnNames = []
	columnInfo =[]
	ON_CHROMOSOME = {}
	read_TSV_File(mafName)
	columns = ""
	usedNames = []
	mafColNames2 = []
	helpful = []
	# print(columns[:-2])
	for i in range(len(columnNames)):
		tmpCN = re.sub("\(.*\)", "", columnNames[i])
		tmpCN = re.sub(" +", "_", tmpCN)
		tmpCN = re.sub("\+", "Plus", tmpCN)
		tmpCN = re.sub("-", "Minus", tmpCN)

		if tmpCN.upper() in usedNames:
			# print("DUPE AT ", tmpCN)
			pass
		else:
			ROW_NAME_TO_ID.update({tmpCN : i})
			mafColNames2.append(tmpCN)
			# print(i, tmpCN, dataType(i))
			try:
				columns +=  tmpCN + ", "
			except:
				columns +=  tmpCN + ", "

			usedNames.append(tmpCN.upper())
			helpful.append([tmpCN, i, len(helpful)])

	for i in range( startAt-1, len(columnInfo)):
		vals = ""
		try:
			currentChrom = ON_CHROMOSOME[columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]]][0]
			if currentChrom == None:
				currentChrom = columnInfo[i][ROW_NAME_TO_ID["Chromosome"]]
				ON_CHROMOSOME[columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]]] = [columnInfo[i][ROW_NAME_TO_ID["Chromosome"]],  columnInfo[i][ROW_NAME_TO_ID["i_HGNC_UniProt_ID"]], columnInfo[i][ROW_NAME_TO_ID["Entrez_Gene_Id"]], columnInfo[i][ROW_NAME_TO_ID["SwissProt_acc_Id"]], columnInfo[i][ROW_NAME_TO_ID["i_COSMIC_Gene"]]]
			else:
				pass
		except:
			# print(i)
			tmpLst = []
			try:
				tmpLst.append(columnInfo[i][ROW_NAME_TO_ID["Chromosome"]])
			except:
				tmpLst.append("---")
			try:
				tmpLst.append(columnInfo[i][ROW_NAME_TO_ID["i_HGNC_UniProt_ID"]])
			except:
				tmpLst.append("---")
			try:
				tmpLst.append(columnInfo[i][ROW_NAME_TO_ID["Entrez_Gene_Id"]])
			except:
				tmpLst.append("---")
			try:
				tmpLst.append(columnInfo[i][ROW_NAME_TO_ID["SwissProt_acc_Id"]])
			except:
				tmpLst.append("---")
			try:
				tmpLst.append(columnInfo[i][ROW_NAME_TO_ID["i_COSMIC_Gene"]])
			except:
				tmpLst.append("---")
			ON_CHROMOSOME.update({columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]]: tmpLst})
			currentChrom = columnInfo[i][ROW_NAME_TO_ID["Chromosome"]]
			# print(columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]], "on", columnInfo[i][ROW_NAME_TO_ID["Chromosome"]])
		# instChromCmd = 'INSERT INTO mutations' + columnInfo[i][ROW_NAME_TO_ID["Chromosome"]] +  ' ( Hugo_Symbol ) VALUES ("' + columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]] + '")'
		# con.execute(instChromCmd)
		# con.commit()
		# print(columnInfo[i][ROW_NAME_TO_ID["Hugo_Symbol"]])

		##########################################################################################################################################################################################################################################################################################################
		for clm in mafColNames2:
			val = ""
			try:
				if dataTypes[clm] != "INT":
					val = '"' + columnInfo[i][ROW_NAME_TO_ID[clm]] + '"'
				else:
					val =  columnInfo[i][ROW_NAME_TO_ID[clm]]
					if val == None  or val == "":
						val = "NULL"
			except:
				try:
					val = '"' +  columnInfo[i][ROW_NAME_TO_ID[clm]] + '"'
				except:
					val = '""'
			vals += val + ", "
			#######################################################################################################################################################################################################################################################################################################################
		onID += 1
		seen = False
		try:
			seen = SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]][columnInfo[i][ROW_NAME_TO_ID['Protein_Change']]]
		except:
			try:
				geneSeen = SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]]
				SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]].update({columnInfo[i][ROW_NAME_TO_ID['Protein_Change']] : True})
			except:
				SEEN_MUTATIONS.update({columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]: {columnInfo[i][ROW_NAME_TO_ID['Protein_Change']] : True} })

		if (seen):
			tmpColName = CANCER_TYPE.lower() + "_muts"
			print("WOOHOOO! DUPLICATE MUTATION AT:", columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']], columnInfo[i][ROW_NAME_TO_ID['Protein_Change']])
			sql_query = "UPDATE mutations SET "+ tmpColName +" = "+ tmpColName +" + 1 WHERE mutation_id = " + SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]][columnInfo[i][ROW_NAME_TO_ID['Protein_Change']]]
			csr.execute(sql_query)
			con.commit()
			sql_query = "UPDATE mutations SET total_instances = total_instances + 1 WHERE mutation_id = " + SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]][columnInfo[i][ROW_NAME_TO_ID['Protein_Change']]]

			csr.execute(sql_query)
			con.commit()
			tcganm = columnInfo[i][ROW_NAME_TO_ID['Tumor_Sample_Barcode']][0:12]
			sql_query = "UPDATE mutations SET caseList = caseList" + "|| \", "  + tcganm + "\" WHERE mutation_id = " + SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]][columnInfo[i][ROW_NAME_TO_ID['Protein_Change']]]
			print (sql_query)
			csr.execute(sql_query)
			con.commit()
			lastRow = csr.lastrowid
			hasMuts.append(lastRow)
		else:
			tmpColName = CANCER_TYPE.lower() + "_muts"
			tcganm = columnInfo[i][ROW_NAME_TO_ID['Tumor_Sample_Barcode']][0:12]
			instCmd = "INSERT INTO mutations (" + columns + "caseList, total_instances, " + tmpColName +") VALUES (" + vals + "'" + tcganm + "', 1, 1 )"
			print(instCmd)
			csr.execute(instCmd)
			con.commit()
			# csr.execute(instCmd)
			lastRow = csr.lastrowid
			SEEN_MUTATIONS[columnInfo[i][ROW_NAME_TO_ID['Hugo_Symbol']]][columnInfo[i][ROW_NAME_TO_ID['Protein_Change']]] = str(lastRow)
			hasMuts.append(lastRow)
			instCmd2 = "INSERT INTO chrom" + currentChrom + "mutations (mutation_id, " + columns[:-2]  + ") VALUES (" + str(lastRow)+  ", " + vals[:-2] + ")"
			# print(instCmd2)

			con.execute(instCmd2)
			con.commit()
	tcganm = columnInfo[i][ROW_NAME_TO_ID['Tumor_Sample_Barcode']]
	print(tcganm[0:12], "\n", hasMuts)
	return tcganm[0:12], hasMuts
	# print(con.execute("SELECT TABLE STATUS FROM 'mutations' LIKE AUTO_INCRIMENT"))TCGA-OR-A5L5
	# print(instCmd)
	
	
	
	

	for k in ON_CHROMOSOME.keys():
		# print(k, ON_CHROMOSOME[k])
		instChromCmd = 'INSERT INTO chrom' + ON_CHROMOSOME[k][0] +  'genes ( '+ chromCols2 +' ) VALUES ("' + k + '", "' + ON_CHROMOSOME[k][1] + '", ' + ON_CHROMOSOME[k][2] + ', "' + ON_CHROMOSOME[k][3] + '", "' + ON_CHROMOSOME[k][4] + '" )'
		con.execute(instChromCmd)
		con.commit()
		instChromCmd = 'INSERT INTO genechroms ( '+ chromCols2 +', chromosome ) VALUES ("' + k + '", "' + ON_CHROMOSOME[k][1] + '", ' + ON_CHROMOSOME[k][2] + ', "' + ON_CHROMOSOME[k][3] + '", "' + ON_CHROMOSOME[k][4] + '", "' + ON_CHROMOSOME[k][0] + '" )'
		con.execute(instChromCmd)
		con.commit()


	



def read_TSV_File(tsv_name):
	global columnNames 
	global columnInfo
	isFirst= True
	isSecond = True
	onRow = 0
	onRow2 = 0
	tsv_file = open(tsv_name, encoding= 'utf-8',errors='ignore' )
	read_tsv = csv.reader(tsv_file, delimiter="\t")


	for row in read_tsv:
		onRow +=1
		if len(row) >10 and isFirst:
			isFirst = False
			for i in range(len(row)):
				columnNames.append(row[i])
				#print(row[i], "\n")
		
		elif len(row) >10 and isSecond:
			isSecond = False
			columnInfo.append([])
			onRow2 = 0
			for i in range(len(row)):
				tmpRowEntry = re.sub('"', "", row[i])
				# columnInfo.append([tmpRowEntry])
				columnInfo[onRow2].append(tmpRowEntry)
				# print(row[i])

			#print(columnInfo)
		elif len(row) >10 and onRow < stopAt:
			onRow2 += 1
			columnInfo.append([])
			for i in range(len(row)):
				tmpRowEntry = re.sub('"', "", row[i])
				# columnInfo[i].append(tmpRowEntry)
				columnInfo[onRow2].append(tmpRowEntry)


def dataType(q):
	try:
		ret = dataTypes[str(q)]
		return ret
	except:
		try:
			ret = dataTypes[columnNames[q]]
			return ret
		except:
			return "TEXT"
	pass


def openCancerType(ct):
	global CANCER_TYPE
	newColName = ct.lower() + "_muts"
	CANCER_TYPE = ct
	con.execute("ALTER TABLE mutations ADD "+ newColName+" INT NOT NULL DEFAULT 0")
	con.commit()
	addClinicalData2.readCancerType(ct)
	a = "MAF/" + ct.upper() + "/"
	fileNames = [f for f in listdir(a) if isfile(join(a, f))]
	for patName in fileNames:
		patID, patMutations = readMAF(a + patName)
		for mID in patMutations:
			failed = False
			# try:
			# 	# con.execute("ALTER TABLE cases ADD mut_"+ str(mID) +" INT NOT NULL DEFAULT 0")
			# 	ON_COLUMN += 1
			# 	# con.commit()
			# except  Exception as e:
			# 	# print (e)
			# 	if (str(e) == "too many columns on sqlite_altertab_cases"):
			# 		print("gotcha!")
			# 		break
			# 	# print(type(e))
			# 	# print(e.args) 
			# # print("UPDATE cases SET mut_"+ str(mID) +" = 1 WHERE case_submitter_id = '" + patID + "'")
			try:
				con.execute("UPDATE cases SET mut_"+ str(mID) +" = 1 WHERE case_submitter_id = '" + patID + "'")
				con.commit()
			except:
				pass
		# con.execute() 1997-01-01



if __name__ == "__main__":
	# chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	# con.execute("DROP TABLE IF EXISTS genechroms")
	# con.commit()
	# con.execute("DROP TABLE IF EXISTS mutations")
	# con.commit()
	# con.execute("DROP TABLE IF EXISTS cases")
	# con.commit()
	# for c in chromosomes:
	# 	con.execute("DROP TABLE IF EXISTS chrom"+ c + "genes")
	# 	con.commit()
	# 	con.execute("DROP TABLE IF EXISTS chrom"+ c + "mutations")
	# 	con.commit()
	# 	con.execute("CREATE TABLE chrom" + c + "genes (hugo_symbol TEXT, UniProt_ID TEXT, Entrez_Gene_Id INT, SwissProt_entry_Id TEXT, i_COSMIC_Gene TEXT)")
	# 	con.commit()
	# 	con.execute("CREATE TABLE chrom" + c + "mutations ( mutation_id INTEGER UNIQUE, Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_position, End_position, Strand, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1, Match_Norm_Seq_Allele2, Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2, Verification_Status, Validation_Status, Mutation_Status, Sequencing_Phase, Sequence_Source, Validation_Method, Score, BAM_file, Sequencer, Tumor_Sample_UUID, Matched_Norm_Sample_UUID, Genome_Change, Annotation_Transcript, Transcript_Strand, Transcript_Exon, Transcript_Position, cDNA_Change, Codon_Change, Protein_Change, Other_Transcripts, Refseq_mRNA_Id, Refseq_prot_Id, SwissProt_acc_Id, SwissProt_entry_Id, Description, UniProt_AApos, UniProt_Region, UniProt_Site, UniProt_Natural_Variations, UniProt_Experimental_Info, GO_Biological_Process, GO_Cellular_Component, GO_Molecular_Function, COSMIC_overlapping_mutations, COSMIC_fusion_genes, COSMIC_tissue_types_affected, COSMIC_total_alterations_in_gene, Tumorscape_Amplification_Peaks, Tumorscape_Deletion_Peaks, TCGAscape_Amplification_Peaks, TCGAscape_Deletion_Peaks, DrugBank, ref_context, gc_content, CCLE_ONCOMAP_overlapping_mutations, CCLE_ONCOMAP_total_mutations_in_gene, CGC_Mutation_Type, CGC_Translocation_Partner, CGC_Tumor_Types_Somatic, CGC_Tumor_Types_Germline, CGC_Other_Diseases, DNARepairGenes_Role, FamilialCancerDatabase_Syndromes, MUTSIG_Published_Results, OREGANNO_ID, OREGANNO_Values, i_1000gp3_AA, i_1000gp3_AC, i_1000gp3_AF, i_1000gp3_AFR_AF, i_1000gp3_AMR_AF, i_1000gp3_AN, i_1000gp3_CIEND, i_1000gp3_CIPOS, i_1000gp3_CS, i_1000gp3_DP, i_1000gp3_EAS_AF, i_1000gp3_END, i_1000gp3_EUR_AF, i_1000gp3_IMPRECISE, i_1000gp3_MC, i_1000gp3_MEINFO, i_1000gp3_MEND, i_1000gp3_MLEN, i_1000gp3_MSTART, i_1000gp3_NS, i_1000gp3_SAS_AF, i_1000gp3_SVLEN, i_1000gp3_SVTYPE, i_1000gp3_TSD, i_AAChange, i_ACHILLES_Lineage_Results_Top_Genes, i_BAM_File, i_CGC_Cancer_Germline_Mut, i_CGC_Cancer_Molecular_Genetics, i_CGC_Cancer_Somatic_Mut, i_CGC_Cancer_Syndrome, i_CGC_Chr, i_CGC_Chr_Band, i_CGC_GeneID, i_CGC_Name, i_CGC_Other_Germline_Mut, i_CGC_Tissue_Type, i_COSMIC_Codon, i_COSMIC_Gene, i_COSMIC_n_overlapping_mutations, i_COSMIC_overlapping_mutation_descriptions, i_COSMIC_overlapping_primary_sites, i_ChromChange, i_ClinVar_ASSEMBLY, i_ClinVar_HGMD_ID, i_ClinVar_SYM, i_ClinVar_TYPE, i_ClinVar_rs, i_Drug_Target, i_ESP_AA, i_ESP_AAC, i_ESP_AA_AC, i_ESP_AA_AGE, i_ESP_AA_GTC, i_ESP_AvgAAsampleReadDepth, i_ESP_AvgEAsampleReadDepth, i_ESP_AvgSampleReadDepth, i_ESP_CA, i_ESP_CDP, i_ESP_CG, i_ESP_CP, i_ESP_Chromosome, i_ESP_DBSNP, i_ESP_DP, i_ESP_EA_AC, i_ESP_EA_AGE, i_ESP_EA_GTC, i_ESP_EXOME_CHIP, i_ESP_FG, i_ESP_GL, i_ESP_GM, i_ESP_GS, i_ESP_GTC, i_ESP_GTS, i_ESP_GWAS_PUBMED, i_ESP_MAF, i_ESP_PH, i_ESP_PP, i_ESP_Position, i_ESP_TAC, i_ESP_TotalAAsamplesCovered, i_ESP_TotalEAsamplesCovered, i_ESP_TotalSamplesCovered, i_Ensembl_so_accession, i_Ensembl_so_term, i_Entrez_Gene_Id, i_Exon, i_Familial_Cancer_Genes_Reference, i_Familial_Cancer_Genes_Synonym, i_Genome_Plus_Minus_10_Bp, i_HGNC_Accession_Numbers, i_HGNC_CCDS_IDs, i_HGNC_Chromosome, i_HGNC_Date_Modified, i_HGNC_Date_Name_Changed, i_HGNC_Date_Symbol_Changed, i_HGNC_Ensembl_Gene_ID, i_HGNC_Ensembl_ID, i_HGNC_Enzyme_IDs, i_HGNC_Gene_family_description, i_HGNC_HGNC_ID, i_HGNC_Locus_Group, i_HGNC_Locus_Type, i_HGNC_Name_Synonyms, i_HGNC_OMIM_ID, i_HGNC_Previous_Names, i_HGNC_Previous_Symbols, i_HGNC_Primary_IDs, i_HGNC_Pubmed_IDs, i_HGNC_Record_Type, i_HGNC_RefSeq, i_HGNC_Secondary_IDs, i_HGNC_Status, i_HGNC_Synonyms, i_HGNC_UCSC_ID, i_HGNC_UniProt_ID, i_HGNC_VEGA_IDs, i_HGVS_coding_DNA_change, i_HGVS_genomic_change, i_HGVS_protein_change, i_NTotCov, i_NVarCov, i_ORegAnno_bin, i_TTotCov, i_TVarCov, i_Transcript_Id, i_Trna_alt1, i_Trna_alt2, i_Trna_ref, i_Trna_tot, i_Trna_var, i_UniProt_alt_uniprot_accessions, i_Variant_Classification, i_Variant_Type, i_annotation_transcript, i_build, i_ccds_id, i_dbNSFP_1000Gp1_AC, i_dbNSFP_1000Gp1_AF, i_dbNSFP_1000Gp1_AFR_AC, i_dbNSFP_1000Gp1_AFR_AF, i_dbNSFP_1000Gp1_AMR_AC, i_dbNSFP_1000Gp1_AMR_AF, i_dbNSFP_1000Gp1_ASN_AC, i_dbNSFP_1000Gp1_ASN_AF, i_dbNSFP_1000Gp1_EUR_AC, i_dbNSFP_1000Gp1_EUR_AF, i_dbNSFP_Ancestral_allele, i_dbNSFP_CADD_phred, i_dbNSFP_CADD_raw, i_dbNSFP_CADD_raw_rankscore, i_dbNSFP_ESP6500_AA_AF, i_dbNSFP_ESP6500_EA_AF, i_dbNSFP_Ensembl_geneid, i_dbNSFP_Ensembl_transcriptid, i_dbNSFP_FATHMM_pred, i_dbNSFP_FATHMM_rankscore, i_dbNSFP_FATHMM_score, i_dbNSFP_GERPPlusPlus_NR, i_dbNSFP_GERPPlusPlus_RS, i_dbNSFP_GERPPlusPlus_RS_rankscore, i_dbNSFP_Interpro_domain, i_dbNSFP_LRT_Omega, i_dbNSFP_LRT_converted_rankscore, i_dbNSFP_LRT_pred, i_dbNSFP_LRT_score, i_dbNSFP_LR_pred, i_dbNSFP_LR_rankscore, i_dbNSFP_LR_score, i_dbNSFP_MutationAssessor_pred, i_dbNSFP_MutationAssessor_rankscore, i_dbNSFP_MutationAssessor_score, i_dbNSFP_MutationTaster_converted_rankscore, i_dbNSFP_MutationTaster_pred, i_dbNSFP_MutationTaster_score, i_dbNSFP_Polyphen2_HDIV_pred, i_dbNSFP_Polyphen2_HDIV_rankscore, i_dbNSFP_Polyphen2_HDIV_score, i_dbNSFP_Polyphen2_HVAR_pred, i_dbNSFP_Polyphen2_HVAR_rankscore, i_dbNSFP_Polyphen2_HVAR_score, i_dbNSFP_RadialSVM_pred, i_dbNSFP_RadialSVM_rankscore, i_dbNSFP_RadialSVM_score, i_dbNSFP_Reliability_index, i_dbNSFP_SIFT_converted_rankscore, i_dbNSFP_SIFT_pred, i_dbNSFP_SIFT_score, i_dbNSFP_SLR_test_statistic, i_dbNSFP_SiPhy_29way_logOdds, i_dbNSFP_SiPhy_29way_logOdds_rankscore, i_dbNSFP_SiPhy_29way_pi, i_dbNSFP_UniSNP_ids, i_dbNSFP_Uniprot_aapos, i_dbNSFP_Uniprot_acc, i_dbNSFP_Uniprot_id, i_dbNSFP_aaalt, i_dbNSFP_aapos, i_dbNSFP_aapos_FATHMM, i_dbNSFP_aapos_SIFT, i_dbNSFP_aaref, i_dbNSFP_cds_strand, i_dbNSFP_codonpos, i_dbNSFP_foldMinusdegenerate, i_dbNSFP_genename, i_dbNSFP_hg18_pos, i_dbNSFP_phastCons100way_vertebrate, i_dbNSFP_phastCons100way_vertebrate_rankscore, i_dbNSFP_phastCons46way_placental, i_dbNSFP_phastCons46way_placental_rankscore, i_dbNSFP_phastCons46way_primate, i_dbNSFP_phastCons46way_primate_rankscore, i_dbNSFP_phyloP100way_vertebrate, i_dbNSFP_phyloP100way_vertebrate_rankscore, i_dbNSFP_phyloP46way_placental, i_dbNSFP_phyloP46way_placental_rankscore, i_dbNSFP_phyloP46way_primate, i_dbNSFP_phyloP46way_primate_rankscore, i_dbNSFP_refcodon, i_dbSNPPopFreq, i_gc_content_full, i_gencode_transcript_name, i_gencode_transcript_status, i_gencode_transcript_tags, i_gencode_transcript_type, i_gene_type, i_havana_transcript, i_refseq_mrna_id, i_secondary_variant_classification, total_instances )")
	# 	con.commit()

	# con.execute("CREATE TABLE mutations ( mutation_id INTEGER UNIQUE, Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_position, End_position, Strand, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1, Match_Norm_Seq_Allele2, Tumor_Validation_Allele1, Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2, Verification_Status, Validation_Status, Mutation_Status, Sequencing_Phase, Sequence_Source, Validation_Method, Score, BAM_file, Sequencer, Tumor_Sample_UUID, Matched_Norm_Sample_UUID, Genome_Change, Annotation_Transcript, Transcript_Strand, Transcript_Exon, Transcript_Position, cDNA_Change, Codon_Change, Protein_Change, Other_Transcripts, Refseq_mRNA_Id, Refseq_prot_Id, SwissProt_acc_Id, SwissProt_entry_Id, Description, UniProt_AApos, UniProt_Region, UniProt_Site, UniProt_Natural_Variations, UniProt_Experimental_Info, GO_Biological_Process, GO_Cellular_Component, GO_Molecular_Function, COSMIC_overlapping_mutations, COSMIC_fusion_genes, COSMIC_tissue_types_affected, COSMIC_total_alterations_in_gene, Tumorscape_Amplification_Peaks, Tumorscape_Deletion_Peaks, TCGAscape_Amplification_Peaks, TCGAscape_Deletion_Peaks, DrugBank, ref_context, gc_content, CCLE_ONCOMAP_overlapping_mutations, CCLE_ONCOMAP_total_mutations_in_gene, CGC_Mutation_Type, CGC_Translocation_Partner, CGC_Tumor_Types_Somatic, CGC_Tumor_Types_Germline, CGC_Other_Diseases, DNARepairGenes_Role, FamilialCancerDatabase_Syndromes, MUTSIG_Published_Results, OREGANNO_ID, OREGANNO_Values, i_1000gp3_AA, i_1000gp3_AC, i_1000gp3_AF, i_1000gp3_AFR_AF, i_1000gp3_AMR_AF, i_1000gp3_AN, i_1000gp3_CIEND, i_1000gp3_CIPOS, i_1000gp3_CS, i_1000gp3_DP, i_1000gp3_EAS_AF, i_1000gp3_END, i_1000gp3_EUR_AF, i_1000gp3_IMPRECISE, i_1000gp3_MC, i_1000gp3_MEINFO, i_1000gp3_MEND, i_1000gp3_MLEN, i_1000gp3_MSTART, i_1000gp3_NS, i_1000gp3_SAS_AF, i_1000gp3_SVLEN, i_1000gp3_SVTYPE, i_1000gp3_TSD, i_AAChange, i_ACHILLES_Lineage_Results_Top_Genes, i_BAM_File, i_CGC_Cancer_Germline_Mut, i_CGC_Cancer_Molecular_Genetics, i_CGC_Cancer_Somatic_Mut, i_CGC_Cancer_Syndrome, i_CGC_Chr, i_CGC_Chr_Band, i_CGC_GeneID, i_CGC_Name, i_CGC_Other_Germline_Mut, i_CGC_Tissue_Type, i_COSMIC_Codon, i_COSMIC_Gene, i_COSMIC_n_overlapping_mutations, i_COSMIC_overlapping_mutation_descriptions, i_COSMIC_overlapping_primary_sites, i_ChromChange, i_ClinVar_ASSEMBLY, i_ClinVar_HGMD_ID, i_ClinVar_SYM, i_ClinVar_TYPE, i_ClinVar_rs, i_Drug_Target, i_ESP_AA, i_ESP_AAC, i_ESP_AA_AC, i_ESP_AA_AGE, i_ESP_AA_GTC, i_ESP_AvgAAsampleReadDepth, i_ESP_AvgEAsampleReadDepth, i_ESP_AvgSampleReadDepth, i_ESP_CA, i_ESP_CDP, i_ESP_CG, i_ESP_CP, i_ESP_Chromosome, i_ESP_DBSNP, i_ESP_DP, i_ESP_EA_AC, i_ESP_EA_AGE, i_ESP_EA_GTC, i_ESP_EXOME_CHIP, i_ESP_FG, i_ESP_GL, i_ESP_GM, i_ESP_GS, i_ESP_GTC, i_ESP_GTS, i_ESP_GWAS_PUBMED, i_ESP_MAF, i_ESP_PH, i_ESP_PP, i_ESP_Position, i_ESP_TAC, i_ESP_TotalAAsamplesCovered, i_ESP_TotalEAsamplesCovered, i_ESP_TotalSamplesCovered, i_Ensembl_so_accession, i_Ensembl_so_term, i_Entrez_Gene_Id, i_Exon, i_Familial_Cancer_Genes_Reference, i_Familial_Cancer_Genes_Synonym, i_Genome_Plus_Minus_10_Bp, i_HGNC_Accession_Numbers, i_HGNC_CCDS_IDs, i_HGNC_Chromosome, i_HGNC_Date_Modified, i_HGNC_Date_Name_Changed, i_HGNC_Date_Symbol_Changed, i_HGNC_Ensembl_Gene_ID, i_HGNC_Ensembl_ID, i_HGNC_Enzyme_IDs, i_HGNC_Gene_family_description, i_HGNC_HGNC_ID, i_HGNC_Locus_Group, i_HGNC_Locus_Type, i_HGNC_Name_Synonyms, i_HGNC_OMIM_ID, i_HGNC_Previous_Names, i_HGNC_Previous_Symbols, i_HGNC_Primary_IDs, i_HGNC_Pubmed_IDs, i_HGNC_Record_Type, i_HGNC_RefSeq, i_HGNC_Secondary_IDs, i_HGNC_Status, i_HGNC_Synonyms, i_HGNC_UCSC_ID, i_HGNC_UniProt_ID, i_HGNC_VEGA_IDs, i_HGVS_coding_DNA_change, i_HGVS_genomic_change, i_HGVS_protein_change, i_NTotCov, i_NVarCov, i_ORegAnno_bin, i_TTotCov, i_TVarCov, i_Transcript_Id, i_Trna_alt1, i_Trna_alt2, i_Trna_ref, i_Trna_tot, i_Trna_var, i_UniProt_alt_uniprot_accessions, i_Variant_Classification, i_Variant_Type, i_annotation_transcript, i_build, i_ccds_id, i_dbNSFP_1000Gp1_AC, i_dbNSFP_1000Gp1_AF, i_dbNSFP_1000Gp1_AFR_AC, i_dbNSFP_1000Gp1_AFR_AF, i_dbNSFP_1000Gp1_AMR_AC, i_dbNSFP_1000Gp1_AMR_AF, i_dbNSFP_1000Gp1_ASN_AC, i_dbNSFP_1000Gp1_ASN_AF, i_dbNSFP_1000Gp1_EUR_AC, i_dbNSFP_1000Gp1_EUR_AF, i_dbNSFP_Ancestral_allele, i_dbNSFP_CADD_phred, i_dbNSFP_CADD_raw, i_dbNSFP_CADD_raw_rankscore, i_dbNSFP_ESP6500_AA_AF, i_dbNSFP_ESP6500_EA_AF, i_dbNSFP_Ensembl_geneid, i_dbNSFP_Ensembl_transcriptid, i_dbNSFP_FATHMM_pred, i_dbNSFP_FATHMM_rankscore, i_dbNSFP_FATHMM_score, i_dbNSFP_GERPPlusPlus_NR, i_dbNSFP_GERPPlusPlus_RS, i_dbNSFP_GERPPlusPlus_RS_rankscore, i_dbNSFP_Interpro_domain, i_dbNSFP_LRT_Omega, i_dbNSFP_LRT_converted_rankscore, i_dbNSFP_LRT_pred, i_dbNSFP_LRT_score, i_dbNSFP_LR_pred, i_dbNSFP_LR_rankscore, i_dbNSFP_LR_score, i_dbNSFP_MutationAssessor_pred, i_dbNSFP_MutationAssessor_rankscore, i_dbNSFP_MutationAssessor_score, i_dbNSFP_MutationTaster_converted_rankscore, i_dbNSFP_MutationTaster_pred, i_dbNSFP_MutationTaster_score, i_dbNSFP_Polyphen2_HDIV_pred, i_dbNSFP_Polyphen2_HDIV_rankscore, i_dbNSFP_Polyphen2_HDIV_score, i_dbNSFP_Polyphen2_HVAR_pred, i_dbNSFP_Polyphen2_HVAR_rankscore, i_dbNSFP_Polyphen2_HVAR_score, i_dbNSFP_RadialSVM_pred, i_dbNSFP_RadialSVM_rankscore, i_dbNSFP_RadialSVM_score, i_dbNSFP_Reliability_index, i_dbNSFP_SIFT_converted_rankscore, i_dbNSFP_SIFT_pred, i_dbNSFP_SIFT_score, i_dbNSFP_SLR_test_statistic, i_dbNSFP_SiPhy_29way_logOdds, i_dbNSFP_SiPhy_29way_logOdds_rankscore, i_dbNSFP_SiPhy_29way_pi, i_dbNSFP_UniSNP_ids, i_dbNSFP_Uniprot_aapos, i_dbNSFP_Uniprot_acc, i_dbNSFP_Uniprot_id, i_dbNSFP_aaalt, i_dbNSFP_aapos, i_dbNSFP_aapos_FATHMM, i_dbNSFP_aapos_SIFT, i_dbNSFP_aaref, i_dbNSFP_cds_strand, i_dbNSFP_codonpos, i_dbNSFP_foldMinusdegenerate, i_dbNSFP_genename, i_dbNSFP_hg18_pos, i_dbNSFP_phastCons100way_vertebrate, i_dbNSFP_phastCons100way_vertebrate_rankscore, i_dbNSFP_phastCons46way_placental, i_dbNSFP_phastCons46way_placental_rankscore, i_dbNSFP_phastCons46way_primate, i_dbNSFP_phastCons46way_primate_rankscore, i_dbNSFP_phyloP100way_vertebrate, i_dbNSFP_phyloP100way_vertebrate_rankscore, i_dbNSFP_phyloP46way_placental, i_dbNSFP_phyloP46way_placental_rankscore, i_dbNSFP_phyloP46way_primate, i_dbNSFP_phyloP46way_primate_rankscore, i_dbNSFP_refcodon, i_dbSNPPopFreq, i_gc_content_full, i_gencode_transcript_name, i_gencode_transcript_status, i_gencode_transcript_tags, i_gencode_transcript_type, i_gene_type, i_havana_transcript, i_refseq_mrna_id, i_secondary_variant_classification, total_instances, caseList,  PRIMARY KEY(mutation_id))")
	# con.commit()
	# con.execute("CREATE TABLE genechroms (hugo_symbol TEXT, UniProt_ID TEXT, Entrez_Gene_Id INT, SwissProt_entry_Id TEXT, i_COSMIC_Gene TEXT, chromosome TEXT)")
	# con.commit()


	# CLINICAL_COLUMNS = "internal_case_id INT, case_id TEXT, case_submitter_id TEXT, project_id TEXT, age_at_index TEXT, age_is_obfuscated TEXT, cause_of_death TEXT, cause_of_death_source TEXT, country_of_residence_at_enrollment TEXT, days_to_birth TEXT, days_to_death TEXT, ethnicity TEXT, gender TEXT, occupation_duration_years TEXT, premature_at_birth TEXT, race TEXT, vital_status TEXT, weeks_gestation_at_birth TEXT, year_of_birth TEXT, year_of_death TEXT, age_at_diagnosis TEXT, ajcc_clinical_m TEXT, ajcc_clinical_n TEXT, ajcc_clinical_stage TEXT, ajcc_clinical_t TEXT, ajcc_pathologic_m TEXT, ajcc_pathologic_n TEXT, ajcc_pathologic_stage TEXT, ajcc_pathologic_t TEXT, ajcc_staging_system_edition TEXT, anaplasia_present TEXT, anaplasia_present_type TEXT, ann_arbor_b_symptoms TEXT, ann_arbor_clinical_stage TEXT, ann_arbor_extranodal_involvement TEXT, ann_arbor_pathologic_stage TEXT, best_overall_response TEXT, breslow_thickness TEXT, burkitt_lymphoma_clinical_variant TEXT, child_pugh_classification TEXT, circumferential_resection_margin TEXT, classification_of_tumor TEXT, cog_liver_stage TEXT, cog_neuroblastoma_risk_group TEXT, cog_renal_stage TEXT, cog_rhabdomyosarcoma_risk_group TEXT, days_to_best_overall_response TEXT, days_to_diagnosis TEXT, days_to_last_follow_up TEXT, days_to_last_known_disease_status TEXT, days_to_recurrence TEXT, eln_risk_classification TEXT, enneking_msts_grade TEXT, enneking_msts_metastasis TEXT, enneking_msts_stage TEXT, enneking_msts_tumor_site TEXT, esophageal_columnar_dysplasia_degree TEXT, esophageal_columnar_metaplasia_present TEXT, figo_stage TEXT, figo_staging_edition_year TEXT, first_symptom_prior_to_diagnosis TEXT, gastric_esophageal_junction_involvement TEXT, gleason_grade_group TEXT, gleason_grade_tertiary TEXT, gleason_patterns_percent TEXT, goblet_cells_columnar_mucosa_present TEXT, greatest_tumor_dimension TEXT, gross_tumor_weight TEXT, icd_10_code TEXT, igcccg_stage TEXT, inpc_grade TEXT, inpc_histologic_group TEXT, inrg_stage TEXT, inss_stage TEXT, international_prognostic_index TEXT, irs_group TEXT, irs_stage TEXT, ishak_fibrosis_score TEXT, iss_stage TEXT, largest_extrapelvic_peritoneal_focus TEXT, last_known_disease_status TEXT, laterality TEXT, lymph_node_involved_site TEXT, lymph_nodes_positive TEXT, lymph_nodes_tested TEXT, lymphatic_invasion_present TEXT, margin_distance TEXT, margins_involved_site TEXT, masaoka_stage TEXT, medulloblastoma_molecular_classification TEXT, metastasis_at_diagnosis TEXT, metastasis_at_diagnosis_site TEXT, method_of_diagnosis TEXT, micropapillary_features TEXT, mitosis_karyorrhexis_index TEXT, mitotic_count TEXT, morphology TEXT, non_nodal_regional_disease TEXT, non_nodal_tumor_deposits TEXT, ovarian_specimen_status TEXT, ovarian_surface_involvement TEXT, papillary_renal_cell_type TEXT, percent_tumor_invasion TEXT, perineural_invasion_present TEXT, peripancreatic_lymph_nodes_positive TEXT, peripancreatic_lymph_nodes_tested TEXT, peritoneal_fluid_cytological_status TEXT, pregnant_at_diagnosis TEXT, primary_diagnosis TEXT, primary_gleason_grade TEXT, prior_malignancy TEXT, prior_treatment TEXT, progression_or_recurrence TEXT, residual_disease TEXT, satellite_nodule_present TEXT, secondary_gleason_grade TEXT, site_of_resection_or_biopsy TEXT, sites_of_involvement TEXT, supratentorial_localization TEXT, synchronous_malignancy TEXT, tissue_or_organ_of_origin TEXT, transglottic_extension TEXT, tumor_confined_to_organ_of_origin TEXT, tumor_depth TEXT, tumor_focality TEXT, tumor_grade TEXT, tumor_largest_dimension_diameter TEXT, tumor_regression_grade TEXT, tumor_stage TEXT, vascular_invasion_present TEXT, vascular_invasion_type TEXT, weiss_assessment_score TEXT, who_cns_grade TEXT, who_nte_grade TEXT, wilms_tumor_histologic_subtype TEXT, year_of_diagnosis TEXT, chemo_concurrent_to_radiation TEXT, days_to_treatment_end TEXT, days_to_treatment_start TEXT, initial_disease_status TEXT, number_of_cycles TEXT, reason_treatment_ended TEXT, regimen_or_line_of_therapy TEXT, therapeutic_agents TEXT, treatment_anatomic_site TEXT, treatment_arm TEXT, treatment_dose TEXT, treatment_dose_units TEXT, treatment_effect TEXT, treatment_effect_indicator TEXT, treatment_frequency TEXT, treatment_intent_type TEXT, treatment_or_therapy TEXT, treatment_outcome TEXT, treatment_type TEXT, wilms_tumor_histologic_subtype2 TEXT, year_of_diagnosis2 TEXT, chemo_concurrent_to_radiation2 TEXT, days_to_treatment_end2 TEXT, days_to_treatment_start2 TEXT, initial_disease_status2 TEXT, number_of_cycles2 TEXT, reason_treatment_ended2 TEXT, regimen_or_line_of_therapy2 TEXT, therapeutic_agents2 TEXT, treatment_anatomic_site2 TEXT, treatment_arm2 TEXT, treatment_dose2 TEXT, treatment_dose_units2 TEXT, treatment_effect2 TEXT, treatment_effect_indicator2 TEXT, treatment_frequency2 TEXT, treatment_intent_type2 TEXT, treatment_or_therapy2 TEXT, treatment_outcome2 TEXT, treatment_type2 TEXT, wilms_tumor_histologic_subtype3 TEXT, year_of_diagnosis3 TEXT, chemo_concurrent_to_radiation3 TEXT, days_to_treatment_end3 TEXT, days_to_treatment_start3 TEXT, initial_disease_status3 TEXT, number_of_cycles3 TEXT, reason_treatment_ended3 TEXT, regimen_or_line_of_therapy3 TEXT, therapeutic_agents3 TEXT, treatment_anatomic_site3 TEXT, treatment_arm3 TEXT, treatment_dose3 TEXT, treatment_dose_units3 TEXT, treatment_effect3 TEXT, treatment_effect_indicator3 TEXT, treatment_frequency3 TEXT, treatment_intent_type3 TEXT, treatment_or_therapy3 TEXT, treatment_outcome3 TEXT, treatment_type3 TEXT"
	# CLINICAL_COLUMNS2 = "internal_case_id , case_id, case_submitter_id, project_id, age_at_index , age_is_obfuscated, cause_of_death, cause_of_death_source, country_of_residence_at_enrollment, days_to_birth , days_to_death , ethnicity, gender , occupation_duration_years, premature_at_birth, race, vital_status, weeks_gestation_at_birth , year_of_birth , year_of_death , age_at_diagnosis , ajcc_clinical_m, ajcc_clinical_n, ajcc_clinical_stage, ajcc_clinical_t, ajcc_pathologic_m, ajcc_pathologic_n, ajcc_pathologic_stage, ajcc_pathologic_t, ajcc_staging_system_edition, anaplasia_present, anaplasia_present_type, ann_arbor_b_symptoms, ann_arbor_clinical_stage, ann_arbor_extranodal_involvement, ann_arbor_pathologic_stage, best_overall_response, breslow_thickness, burkitt_lymphoma_clinical_variant, child_pugh_classification, circumferential_resection_margin, classification_of_tumor, cog_liver_stage, cog_neuroblastoma_risk_group, cog_renal_stage, cog_rhabdomyosarcoma_risk_group, days_to_best_overall_response, days_to_diagnosis , days_to_last_follow_up , days_to_last_known_disease_status, days_to_recurrence, eln_risk_classification, enneking_msts_grade, enneking_msts_metastasis, enneking_msts_stage, enneking_msts_tumor_site, esophageal_columnar_dysplasia_degree, esophageal_columnar_metaplasia_present, figo_stage, figo_staging_edition_year, first_symptom_prior_to_diagnosis, gastric_esophageal_junction_involvement, gleason_grade_group, gleason_grade_tertiary, gleason_patterns_percent, goblet_cells_columnar_mucosa_present, greatest_tumor_dimension, gross_tumor_weight, icd_10_code, igcccg_stage, inpc_grade, inpc_histologic_group, inrg_stage, inss_stage, international_prognostic_index, irs_group, irs_stage, ishak_fibrosis_score, iss_stage, largest_extrapelvic_peritoneal_focus, last_known_disease_status, laterality, lymph_node_involved_site, lymph_nodes_positive, lymph_nodes_tested, lymphatic_invasion_present, margin_distance, margins_involved_site, masaoka_stage, medulloblastoma_molecular_classification, metastasis_at_diagnosis, metastasis_at_diagnosis_site, method_of_diagnosis, micropapillary_features, mitosis_karyorrhexis_index, mitotic_count, morphology, non_nodal_regional_disease, non_nodal_tumor_deposits, ovarian_specimen_status, ovarian_surface_involvement, papillary_renal_cell_type, percent_tumor_invasion, perineural_invasion_present, peripancreatic_lymph_nodes_positive, peripancreatic_lymph_nodes_tested, peritoneal_fluid_cytological_status, pregnant_at_diagnosis, primary_diagnosis, primary_gleason_grade, prior_malignancy, prior_treatment, progression_or_recurrence, residual_disease, satellite_nodule_present, secondary_gleason_grade, site_of_resection_or_biopsy, sites_of_involvement, supratentorial_localization, synchronous_malignancy, tissue_or_organ_of_origin, transglottic_extension, tumor_confined_to_organ_of_origin, tumor_depth, tumor_focality, tumor_grade, tumor_largest_dimension_diameter, tumor_regression_grade, tumor_stage, vascular_invasion_present, vascular_invasion_type, weiss_assessment_score, who_cns_grade, who_nte_grade, wilms_tumor_histologic_subtype, year_of_diagnosis , chemo_concurrent_to_radiation, days_to_treatment_end, days_to_treatment_start, initial_disease_status, number_of_cycles, reason_treatment_ended, regimen_or_line_of_therapy, therapeutic_agents, treatment_anatomic_site, treatment_arm, treatment_dose, treatment_dose_units, treatment_effect, treatment_effect_indicator, treatment_frequency, treatment_type, treatment_or_therapy, treatment_outcome, treatment_type, wilms_tumor_histologic_subtype2, year_of_diagnosis2, chemo_concurrent_to_radiation2, days_to_treatment_end2, days_to_treatment_start2, initial_disease_status2, number_of_cycles2, reason_treatment_ended2, regimen_or_line_of_therapy2, therapeutic_agents2, treatment_anatomic_site2, treatment_arm2, treatment_dose2, treatment_dose_units2, treatment_effect2, treatment_effect_indicator2, treatment_frequency2, treatment_type2, treatment_or_therapy2, treatment_outcome2, treatment_type2, wilms_tumor_histologic_subtype3, year_of_diagnosis3, chemo_concurrent_to_radiation3, days_to_treatment_end3, days_to_treatment_start3, initial_disease_status3, number_of_cycles3, reason_treatment_ended3, regimen_or_line_of_therapy3, therapeutic_agents3, treatment_anatomic_site3, treatment_arm3, treatment_dose3, treatment_dose_units3, treatment_effect3, treatment_effect_indicator3, treatment_frequency3, treatment_type3, treatment_or_therapy3, treatment_outcome3, treatment_type3"
	# con.execute("DROP TABLE IF EXISTS cases")
	# con.commit()
	# con.execute("CREATE TABLE cases (" + CLINICAL_COLUMNS + ")")
	# con.commit()



	# openCancerType("ACC")
	con.execute("DROP TABLE IF EXISTS cases")
	con.commit()
	con.close()