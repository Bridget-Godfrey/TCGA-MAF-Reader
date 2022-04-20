import sqlite3
import csv
import re
# con = sqlite3.connect('mutations.db')
# cur = con.cursor()

CANCER_TYPE = "acc"
NUMREPEAT = 21
num_rows = 154
# con.execute("DROP TABLE IF EXISTS " + CANCER_TYPE)
# con.execute("DROP TABLE IF EXISTS allTypes")

dataType = ["text", "text", "text", "int", "TEXT", "text", "text", "text", "int", "int", "text", "int", "text", "text", "text", "text", "int", "int", "int", "int", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "int", "int", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "int", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text",  "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text"]
for i in range(50):
	dataType.append("text")
columnNames = ["case_id", "case_submitter_id", "project_id", "age_at_index", "age_is_obfuscated", "cause_of_death", "cause_of_death_source", "country_of_residence_at_enrollment", "days_to_birth", "days_to_death", "ethnicity", "gender", "occupation_duration_years", "premature_at_birth", "race", "vital_status", "weeks_gestation_at_birth", "year_of_birth", "year_of_death", "age_at_diagnosis", "ajcc_clinical_m", "ajcc_clinical_n", "ajcc_clinical_stage", "ajcc_clinical_t", "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_staging_system_edition", "anaplasia_present", "anaplasia_present_type", "ann_arbor_b_symptoms", "ann_arbor_clinical_stage", "ann_arbor_extranodal_involvement", "ann_arbor_pathologic_stage", "best_overall_response", "breslow_thickness", "burkitt_lymphoma_clinical_variant", "child_pugh_classification", "circumferential_resection_margin", "classification_of_tumor", "cog_liver_stage", "cog_neuroblastoma_risk_group", "cog_renal_stage", "cog_rhabdomyosarcoma_risk_group", "days_to_best_overall_response", "days_to_diagnosis", "days_to_last_follow_up", "days_to_last_known_disease_status", "days_to_recurrence", "eln_risk_classification", "enneking_msts_grade", "enneking_msts_metastasis", "enneking_msts_stage", "enneking_msts_tumor_site", "esophageal_columnar_dysplasia_degree", "esophageal_columnar_metaplasia_present", "figo_stage", "figo_staging_edition_year", "first_symptom_prior_to_diagnosis", "gastric_esophageal_junction_involvement", "gleason_grade_group", "gleason_grade_tertiary", "gleason_patterns_percent", "goblet_cells_columnar_mucosa_present", "greatest_tumor_dimension", "gross_tumor_weight", "icd_10_code", "igcccg_stage", "inpc_grade", "inpc_histologic_group", "inrg_stage", "inss_stage", "international_prognostic_index", "irs_group", "irs_stage", "ishak_fibrosis_score", "iss_stage", "largest_extrapelvic_peritoneal_focus", "last_known_disease_status", "laterality", "lymph_node_involved_site", "lymph_nodes_positive", "lymph_nodes_tested", "lymphatic_invasion_present", "margin_distance", "margins_involved_site", "masaoka_stage", "medulloblastoma_molecular_classification", "metastasis_at_diagnosis", "metastasis_at_diagnosis_site", "method_of_diagnosis", "micropapillary_features", "mitosis_karyorrhexis_index", "mitotic_count", "morphology", "non_nodal_regional_disease", "non_nodal_tumor_deposits", "ovarian_specimen_status", "ovarian_surface_involvement", "papillary_renal_cell_type", "percent_tumor_invasion", "perineural_invasion_present", "peripancreatic_lymph_nodes_positive", "peripancreatic_lymph_nodes_tested", "peritoneal_fluid_cytological_status", "pregnant_at_diagnosis", "primary_diagnosis", "primary_gleason_grade", "prior_malignancy", "prior_treatment", "progression_or_recurrence", "residual_disease", "satellite_nodule_present", "secondary_gleason_grade", "site_of_resection_or_biopsy", "sites_of_involvement", "supratentorial_localization", "synchronous_malignancy", "tissue_or_organ_of_origin", "transglottic_extension", "tumor_confined_to_organ_of_origin", "tumor_depth", "tumor_focality", "tumor_grade", "tumor_largest_dimension_diameter", "tumor_regression_grade", "tumor_stage", "vascular_invasion_present", "vascular_invasion_type", "weiss_assessment_score", "who_cns_grade", "who_nte_grade", "wilms_tumor_histologic_subtype", "year_of_diagnosis", "chemo_concurrent_to_radiation", "days_to_treatment_end", "days_to_treatment_start", "initial_disease_status", "number_of_cycles", "reason_treatment_ended", "regimen_or_line_of_therapy", "therapeutic_agents", "treatment_anatomic_site", "treatment_arm", "treatment_dose", "treatment_dose_units", "treatment_effect", "treatment_effect_indicator", "treatment_frequency", "treatment_intent_type", "treatment_or_therapy", "treatment_outcome", "treatment_type"]
columnIndex = {}
dataRows = []
tsv_file = open("Clinical Data/" +CANCER_TYPE.upper() + "/clinical.tsv", encoding= 'utf-8',errors='ignore' )
read_tsv = csv.reader(tsv_file, delimiter="\t")
onRow = 0
j = 0
previousID2 = 0
previousID = 000
for row in read_tsv:
	if onRow == 0:
		for i in range(len(row)):
			columnIndex.update({row[i]: i})
			j = i
		num_rows = len(row)
		for i in range(len(row)-NUMREPEAT, len(row), 1):
			columnIndex.update({row[i] + "2": i+NUMREPEAT})
			columnNames.append(row[i] + "2")
		for i in range(len(row)-NUMREPEAT, len(row), 1):
			columnIndex.update({row[i] + "3": i+NUMREPEAT+NUMREPEAT})
			columnNames.append(row[i] + "3")
	elif row[0] == previousID:
		if row[0] == previousID2:
			for i in range(num_rows-NUMREPEAT, num_rows, 1):
				dataRows[len(dataRows)-1][i+NUMREPEAT+NUMREPEAT] = row[i]
		else:
			for i in range(num_rows-NUMREPEAT, num_rows, 1):
				# print(i, i+NUMREPEAT, num_rows, row[i], columnNames[i+NUMREPEAT])
				dataRows[len(dataRows)-1][i+NUMREPEAT] = row[i]
		
		
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
		
		dataRows.append(tmp)
	onRow +=1





for col in range(num_rows+NUMREPEAT + NUMREPEAT):
	print(str(col) + ". " + columnNames[col] + "(" + dataType[col]+ ") :")
	for i in range(10):
		print(dataRows[i][col])
	print("\n\n")


