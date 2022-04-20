import sqlite3
import csv
import re
CANCER_TYPE = "acc"
NUMREPEAT = 21
num_rows = 154

CASE_NUM = 1
con = sqlite3.connect('clinicalData.db')
csr = con.cursor()
clinical_DataType = []
# clinical_DataType = ["TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "INT", "INT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "INT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT",  "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT"]
for i in range(500):
	clinical_DataType.append("TEXT")

clinical_columnNames = ["case_id", "case_submitter_id", "project_id", "age_at_index", "age_is_obfuscated", "cause_of_death", "cause_of_death_source", "country_of_residence_at_enrollment", "days_to_birth", "days_to_death", "ethnicity", "gender", "occupation_duration_years", "premature_at_birth", "race", "vital_status", "weeks_gestation_at_birth", "year_of_birth", "year_of_death", "age_at_diagnosis", "ajcc_clinical_m", "ajcc_clinical_n", "ajcc_clinical_stage", "ajcc_clinical_t", "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_staging_system_edition", "anaplasia_present", "anaplasia_present_type", "ann_arbor_b_symptoms", "ann_arbor_clinical_stage", "ann_arbor_extranodal_involvement", "ann_arbor_pathologic_stage", "best_overall_response", "breslow_thickness", "burkitt_lymphoma_clinical_variant", "child_pugh_classification", "circumferential_resection_margin", "classification_of_tumor", "cog_liver_stage", "cog_neuroblastoma_risk_group", "cog_renal_stage", "cog_rhabdomyosarcoma_risk_group", "days_to_best_overall_response", "days_to_diagnosis", "days_to_last_follow_up", "days_to_last_known_disease_status", "days_to_recurrence", "eln_risk_classification", "enneking_msts_grade", "enneking_msts_metastasis", "enneking_msts_stage", "enneking_msts_tumor_site", "esophageal_columnar_dysplasia_degree", "esophageal_columnar_metaplasia_present", "figo_stage", "figo_staging_edition_year", "first_symptom_prior_to_diagnosis", "gastric_esophageal_junction_involvement", "gleason_grade_group", "gleason_grade_tertiary", "gleason_patterns_percent", "goblet_cells_columnar_mucosa_present", "greatest_tumor_dimension", "gross_tumor_weight", "icd_10_code", "igcccg_stage", "inpc_grade", "inpc_histologic_group", "inrg_stage", "inss_stage", "international_prognostic_index", "irs_group", "irs_stage", "ishak_fibrosis_score", "iss_stage", "largest_extrapelvic_peritoneal_focus", "last_known_disease_status", "laterality", "lymph_node_involved_site", "lymph_nodes_positive", "lymph_nodes_tested", "lymphatic_invasion_present", "margin_distance", "margins_involved_site", "masaoka_stage", "medulloblastoma_molecular_classification", "metastasis_at_diagnosis", "metastasis_at_diagnosis_site", "method_of_diagnosis", "micropapillary_features", "mitosis_karyorrhexis_index", "mitotic_count", "morphology", "non_nodal_regional_disease", "non_nodal_tumor_deposits", "ovarian_specimen_status", "ovarian_surface_involvement", "papillary_renal_cell_type", "percent_tumor_invasion", "perineural_invasion_present", "peripancreatic_lymph_nodes_positive", "peripancreatic_lymph_nodes_tested", "peritoneal_fluid_cytological_status", "pregnant_at_diagnosis", "primary_diagnosis", "primary_gleason_grade", "prior_malignancy", "prior_treatment", "progression_or_recurrence", "residual_disease", "satellite_nodule_present", "secondary_gleason_grade", "site_of_resection_or_biopsy", "sites_of_involvement", "supratentorial_localization", "synchronous_malignancy", "tissue_or_organ_of_origin", "transglottic_extension", "tumor_confined_to_organ_of_origin", "tumor_depth", "tumor_focality", "tumor_grade", "tumor_largest_dimension_diameter", "tumor_regression_grade", "tumor_stage", "vascular_invasion_present", "vascular_invasion_type", "weiss_assessment_score", "who_cns_grade", "who_nte_grade", "wilms_tumor_histologic_subtype", "year_of_diagnosis", "chemo_concurrent_to_radiation", "days_to_treatment_end", "days_to_treatment_start", "initial_disease_status", "number_of_cycles", "reason_treatment_ended", "regimen_or_line_of_therapy", "therapeutic_agents", "treatment_anatomic_site", "treatment_arm", "treatment_dose", "treatment_dose_units", "treatment_effect", "treatment_effect_indicator", "treatment_frequency", "treatment_intent_type", "treatment_or_therapy", "treatment_outcome", "treatment_type"]
clinical_columnIndex = {}
clinical_dataRows = []
clinical_tsv_file = open("Clinical Data/" +CANCER_TYPE.upper() + "/clinical.tsv", encoding= 'utf-8',errors='ignore' )
clinical_read_tsv = csv.reader(clinical_tsv_file, delimiter="\t")
clinical_onRow = 0
clinical_previousID2 = 0
clinical_previousID = 000


CLINICAL_COLUMNS = "internal_case_id INT, case_id TEXT, case_submitter_id TEXT, project_id TEXT, age_at_index TEXT, age_is_obfuscated TEXT, cause_of_death TEXT, cause_of_death_source TEXT, country_of_residence_at_enrollment TEXT, days_to_birth TEXT, days_to_death TEXT, ethnicity TEXT, gender TEXT, occupation_duration_years TEXT, premature_at_birth TEXT, race TEXT, vital_status TEXT, weeks_gestation_at_birth TEXT, year_of_birth TEXT, year_of_death TEXT, age_at_diagnosis TEXT, ajcc_clinical_m TEXT, ajcc_clinical_n TEXT, ajcc_clinical_stage TEXT, ajcc_clinical_t TEXT, ajcc_pathologic_m TEXT, ajcc_pathologic_n TEXT, ajcc_pathologic_stage TEXT, ajcc_pathologic_t TEXT, ajcc_staging_system_edition TEXT, anaplasia_present TEXT, anaplasia_present_type TEXT, ann_arbor_b_symptoms TEXT, ann_arbor_clinical_stage TEXT, ann_arbor_extranodal_involvement TEXT, ann_arbor_pathologic_stage TEXT, best_overall_response TEXT, breslow_thickness TEXT, burkitt_lymphoma_clinical_variant TEXT, child_pugh_classification TEXT, circumferential_resection_margin TEXT, classification_of_tumor TEXT, cog_liver_stage TEXT, cog_neuroblastoma_risk_group TEXT, cog_renal_stage TEXT, cog_rhabdomyosarcoma_risk_group TEXT, days_to_best_overall_response TEXT, days_to_diagnosis TEXT, days_to_last_follow_up TEXT, days_to_last_known_disease_status TEXT, days_to_recurrence TEXT, eln_risk_classification TEXT, enneking_msts_grade TEXT, enneking_msts_metastasis TEXT, enneking_msts_stage TEXT, enneking_msts_tumor_site TEXT, esophageal_columnar_dysplasia_degree TEXT, esophageal_columnar_metaplasia_present TEXT, figo_stage TEXT, figo_staging_edition_year TEXT, first_symptom_prior_to_diagnosis TEXT, gastric_esophageal_junction_involvement TEXT, gleason_grade_group TEXT, gleason_grade_tertiary TEXT, gleason_patterns_percent TEXT, goblet_cells_columnar_mucosa_present TEXT, greatest_tumor_dimension TEXT, gross_tumor_weight TEXT, icd_10_code TEXT, igcccg_stage TEXT, inpc_grade TEXT, inpc_histologic_group TEXT, inrg_stage TEXT, inss_stage TEXT, international_prognostic_index TEXT, irs_group TEXT, irs_stage TEXT, ishak_fibrosis_score TEXT, iss_stage TEXT, largest_extrapelvic_peritoneal_focus TEXT, last_known_disease_status TEXT, laterality TEXT, lymph_node_involved_site TEXT, lymph_nodes_positive TEXT, lymph_nodes_tested TEXT, lymphatic_invasion_present TEXT, margin_distance TEXT, margins_involved_site TEXT, masaoka_stage TEXT, medulloblastoma_molecular_classification TEXT, metastasis_at_diagnosis TEXT, metastasis_at_diagnosis_site TEXT, method_of_diagnosis TEXT, micropapillary_features TEXT, mitosis_karyorrhexis_index TEXT, mitotic_count TEXT, morphology TEXT, non_nodal_regional_disease TEXT, non_nodal_tumor_deposits TEXT, ovarian_specimen_status TEXT, ovarian_surface_involvement TEXT, papillary_renal_cell_type TEXT, percent_tumor_invasion TEXT, perineural_invasion_present TEXT, peripancreatic_lymph_nodes_positive TEXT, peripancreatic_lymph_nodes_tested TEXT, peritoneal_fluid_cytological_status TEXT, pregnant_at_diagnosis TEXT, primary_diagnosis TEXT, primary_gleason_grade TEXT, prior_malignancy TEXT, prior_treatment TEXT, progression_or_recurrence TEXT, residual_disease TEXT, satellite_nodule_present TEXT, secondary_gleason_grade TEXT, site_of_resection_or_biopsy TEXT, sites_of_involvement TEXT, supratentorial_localization TEXT, synchronous_malignancy TEXT, tissue_or_organ_of_origin TEXT, transglottic_extension TEXT, tumor_confined_to_organ_of_origin TEXT, tumor_depth TEXT, tumor_focality TEXT, tumor_grade TEXT, tumor_largest_dimension_diameter TEXT, tumor_regression_grade TEXT, tumor_stage TEXT, vascular_invasion_present TEXT, vascular_invasion_type TEXT, weiss_assessment_score TEXT, who_cns_grade TEXT, who_nte_grade TEXT, wilms_tumor_histologic_subtype TEXT, year_of_diagnosis TEXT, chemo_concurrent_to_radiation TEXT, days_to_treatment_end TEXT, days_to_treatment_start TEXT, initial_disease_status TEXT, number_of_cycles TEXT, reason_treatment_ended TEXT, regimen_or_line_of_therapy TEXT, therapeutic_agents TEXT, treatment_anatomic_site TEXT, treatment_arm TEXT, treatment_dose TEXT, treatment_dose_units TEXT, treatment_effect TEXT, treatment_effect_indicator TEXT, treatment_frequency TEXT, treatment_intent_type TEXT, treatment_or_therapy TEXT, treatment_outcome TEXT, treatment_type TEXT, wilms_tumor_histologic_subtype2 TEXT, year_of_diagnosis2 TEXT, chemo_concurrent_to_radiation2 TEXT, days_to_treatment_end2 TEXT, days_to_treatment_start2 TEXT, initial_disease_status2 TEXT, number_of_cycles2 TEXT, reason_treatment_ended2 TEXT, regimen_or_line_of_therapy2 TEXT, therapeutic_agents2 TEXT, treatment_anatomic_site2 TEXT, treatment_arm2 TEXT, treatment_dose2 TEXT, treatment_dose_units2 TEXT, treatment_effect2 TEXT, treatment_effect_indicator2 TEXT, treatment_frequency2 TEXT, treatment_intent_type2 TEXT, treatment_or_therapy2 TEXT, treatment_outcome2 TEXT, treatment_type2 TEXT, wilms_tumor_histologic_subtype3 TEXT, year_of_diagnosis3 TEXT, chemo_concurrent_to_radiation3 TEXT, days_to_treatment_end3 TEXT, days_to_treatment_start3 TEXT, initial_disease_status3 TEXT, number_of_cycles3 TEXT, reason_treatment_ended3 TEXT, regimen_or_line_of_therapy3 TEXT, therapeutic_agents3 TEXT, treatment_anatomic_site3 TEXT, treatment_arm3 TEXT, treatment_dose3 TEXT, treatment_dose_units3 TEXT, treatment_effect3 TEXT, treatment_effect_indicator3 TEXT, treatment_frequency3 TEXT, treatment_intent_type3 TEXT, treatment_or_therapy3 TEXT, treatment_outcome3 TEXT, treatment_type3 TEXT"
CLINICAL_COLUMNS2 = "internal_case_id , case_id, case_submitter_id, project_id, age_at_index , age_is_obfuscated, cause_of_death, cause_of_death_source, country_of_residence_at_enrollment, days_to_birth , days_to_death , ethnicity, gender , occupation_duration_years, premature_at_birth, race, vital_status, weeks_gestation_at_birth , year_of_birth , year_of_death , age_at_diagnosis , ajcc_clinical_m, ajcc_clinical_n, ajcc_clinical_stage, ajcc_clinical_t, ajcc_pathologic_m, ajcc_pathologic_n, ajcc_pathologic_stage, ajcc_pathologic_t, ajcc_staging_system_edition, anaplasia_present, anaplasia_present_type, ann_arbor_b_symptoms, ann_arbor_clinical_stage, ann_arbor_extranodal_involvement, ann_arbor_pathologic_stage, best_overall_response, breslow_thickness, burkitt_lymphoma_clinical_variant, child_pugh_classification, circumferential_resection_margin, classification_of_tumor, cog_liver_stage, cog_neuroblastoma_risk_group, cog_renal_stage, cog_rhabdomyosarcoma_risk_group, days_to_best_overall_response, days_to_diagnosis , days_to_last_follow_up , days_to_last_known_disease_status, days_to_recurrence, eln_risk_classification, enneking_msts_grade, enneking_msts_metastasis, enneking_msts_stage, enneking_msts_tumor_site, esophageal_columnar_dysplasia_degree, esophageal_columnar_metaplasia_present, figo_stage, figo_staging_edition_year, first_symptom_prior_to_diagnosis, gastric_esophageal_junction_involvement, gleason_grade_group, gleason_grade_tertiary, gleason_patterns_percent, goblet_cells_columnar_mucosa_present, greatest_tumor_dimension, gross_tumor_weight, icd_10_code, igcccg_stage, inpc_grade, inpc_histologic_group, inrg_stage, inss_stage, international_prognostic_index, irs_group, irs_stage, ishak_fibrosis_score, iss_stage, largest_extrapelvic_peritoneal_focus, last_known_disease_status, laterality, lymph_node_involved_site, lymph_nodes_positive, lymph_nodes_tested, lymphatic_invasion_present, margin_distance, margins_involved_site, masaoka_stage, medulloblastoma_molecular_classification, metastasis_at_diagnosis, metastasis_at_diagnosis_site, method_of_diagnosis, micropapillary_features, mitosis_karyorrhexis_index, mitotic_count, morphology, non_nodal_regional_disease, non_nodal_tumor_deposits, ovarian_specimen_status, ovarian_surface_involvement, papillary_renal_cell_type, percent_tumor_invasion, perineural_invasion_present, peripancreatic_lymph_nodes_positive, peripancreatic_lymph_nodes_tested, peritoneal_fluid_cytological_status, pregnant_at_diagnosis, primary_diagnosis, primary_gleason_grade, prior_malignancy, prior_treatment, progression_or_recurrence, residual_disease, satellite_nodule_present, secondary_gleason_grade, site_of_resection_or_biopsy, sites_of_involvement, supratentorial_localization, synchronous_malignancy, tissue_or_organ_of_origin, transglottic_extension, tumor_confined_to_organ_of_origin, tumor_depth, tumor_focality, tumor_grade, tumor_largest_dimension_diameter, tumor_regression_grade, tumor_stage, vascular_invasion_present, vascular_invasion_type, weiss_assessment_score, who_cns_grade, who_nte_grade, wilms_tumor_histologic_subtype, year_of_diagnosis , chemo_concurrent_to_radiation, days_to_treatment_end, days_to_treatment_start, initial_disease_status, number_of_cycles, reason_treatment_ended, regimen_or_line_of_therapy, therapeutic_agents, treatment_anatomic_site, treatment_arm, treatment_dose, treatment_dose_units, treatment_effect, treatment_effect_indicator, treatment_frequency, treatment_type, treatment_or_therapy, treatment_outcome, treatment_type, wilms_tumor_histologic_subtype2, year_of_diagnosis2, chemo_concurrent_to_radiation2, days_to_treatment_end2, days_to_treatment_start2, initial_disease_status2, number_of_cycles2, reason_treatment_ended2, regimen_or_line_of_therapy2, therapeutic_agents2, treatment_anatomic_site2, treatment_arm2, treatment_dose2, treatment_dose_units2, treatment_effect2, treatment_effect_indicator2, treatment_frequency2, treatment_type2, treatment_or_therapy2, treatment_outcome2, treatment_type2, wilms_tumor_histologic_subtype3, year_of_diagnosis3, chemo_concurrent_to_radiation3, days_to_treatment_end3, days_to_treatment_start3, initial_disease_status3, number_of_cycles3, reason_treatment_ended3, regimen_or_line_of_therapy3, therapeutic_agents3, treatment_anatomic_site3, treatment_arm3, treatment_dose3, treatment_dose_units3, treatment_effect3, treatment_effect_indicator3, treatment_frequency3, treatment_type3, treatment_or_therapy3, treatment_outcome3, treatment_type3"
con.execute("DROP TABLE IF EXISTS cases")
con.commit()
con.execute("CREATE TABLE cases (" + CLINICAL_COLUMNS + ")")
con.commit()
def readCancerType(cncr_type):

	global CANCER_TYPE
	global clinical_columnNames 
	global clinical_columnIndex 
	global clinical_dataRows 
	global clinical_tsv_file 
	global clinical_read_tsv 
	global clinical_onRow 
	global clinical_previousID2 
	global clinical_previousID
	global CASE_NUM
	localColumnNames = []
	localColumnNames = clinical_columnNames[:]
	



	CANCER_TYPE = cncr_type.upper()
	con.execute("DROP TABLE IF EXISTS "+ CANCER_TYPE)
	con.commit()
	con.execute("CREATE TABLE " + CANCER_TYPE +  " (" + CLINICAL_COLUMNS + ")")
	con.commit()
	clinical_columnIndex = {}
	clinical_dataRows = []
	clinical_tsv_file = open("Clinical Data/" +CANCER_TYPE.upper() + "/clinical.tsv", encoding= 'utf-8',errors='ignore' )
	clinical_read_tsv = csv.reader(clinical_tsv_file, delimiter="\t")

	clinical_onRow = 0
	j = 0
	previousID2 = 0
	previousID = 0
	for row in clinical_read_tsv:
		if clinical_onRow == 0:
			for i in range(len(row)):
				clinical_columnIndex.update({row[i]: i})
				j = i
			num_rows = len(row)
			for i in range(len(row)-NUMREPEAT, len(row), 1):
				clinical_columnIndex.update({row[i] + "2": i+NUMREPEAT})
				localColumnNames.append(row[i] + "2")
			for i in range(len(row)-NUMREPEAT, len(row), 1):
				clinical_columnIndex.update({row[i] + "3": i+NUMREPEAT+NUMREPEAT})
				localColumnNames.append(row[i] + "3")
		elif row[0] == previousID:
			if row[0] == previousID2:
				for i in range(num_rows-NUMREPEAT, num_rows, 1):
					clinical_dataRows[len(clinical_dataRows)-1][i+NUMREPEAT+NUMREPEAT] = row[i]
			else:
				for i in range(num_rows-NUMREPEAT, num_rows, 1):
					# print(i, i+NUMREPEAT, num_rows, row[i], columnNames[i+NUMREPEAT])
					clinical_dataRows[len(clinical_dataRows)-1][i+NUMREPEAT] = row[i]
			
			
			# print("treatmentType")

		else:

			previousID2 = previousID
			previousID = row[0]

			tmp = []
			for i in range(len(row)):
				tmp.append(row[i])
			for i in range(num_rows-NUMREPEAT, num_rows, 1):
				tmp.append("'--")
			for i in range(num_rows-NUMREPEAT-NUMREPEAT, num_rows, 1):
				tmp.append("'--")
			
			clinical_dataRows.append(tmp)
		clinical_onRow +=1

		# for i in range(len(localColumnNames)):
		# 	print(localColumnNames[i] + " " + clinical_DataType[i], end = ", ")
	for i in range(len(clinical_dataRows)):
		clinValues = []
		for j in range(len(localColumnNames)):
			if clinical_DataType[j].upper() != "TEXT":
				clinValues.append(clinical_dataRows[i][j])
			else:

				# print("adding text at ", i, j)
				clinValues.append( '"' + clinical_dataRows[i][j] + '"')

		# num = 0		
		# for q in clinValues:
		# 	print(localColumnNames[num], q)
		# 	num += 1
		print(CASE_NUM, clinical_dataRows[i][1])
		con.execute("INSERT INTO cases (" + CLINICAL_COLUMNS2 + ") VALUES (" + str(CASE_NUM) + ", " + ", ".join(clinValues)  + ")")
		# print("(" + str(CASE_NUM) + ", " + ", ".join(clinValues)  + ")")
		con.commit()
		con.execute("INSERT INTO " + CANCER_TYPE + " (" + CLINICAL_COLUMNS2 + ") VALUES (" + str(CASE_NUM) + ", " + ", ".join(clinValues)  + ")")
		con.commit()
		CASE_NUM += 1
	# con.close()




# readCancerType("acc")
# readCancerType("blca")
# readCancerType("blca")


readCancerType('acc')
readCancerType('blca')
readCancerType('brca')
readCancerType('cesc')
readCancerType('chol')
readCancerType('coad')
readCancerType('dlbc')
readCancerType('esca')
readCancerType('gbm')
readCancerType('hnsc')
readCancerType('kich')
readCancerType('kirc')
readCancerType('kirp')
readCancerType('laml')
readCancerType('lgg')
readCancerType('lihc')
readCancerType('luad')
readCancerType('lusc')
readCancerType('ov')
readCancerType('paad')
readCancerType('pcpg')
readCancerType('prad')
readCancerType('read')
readCancerType('sarc')
readCancerType('skcm')
readCancerType('stad')
readCancerType('tgct')
readCancerType('thca')
readCancerType('thym')
readCancerType('ucec')
readCancerType('ucs')
readCancerType('uvm')