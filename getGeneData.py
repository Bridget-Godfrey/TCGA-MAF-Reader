# getGeneData.py
import sqlite3
import csv
import re
import os
import pickle
import time
import multiprocessing

seenIDs = {}


# SELECT * FROM "mutations" WHERE "Hugo_Symbol" LIKE 'PIK3CA'
forwardTable = "internal_case_id, case_id, case_submitter_id, project_id, age_at_index , age_is_obfuscated, cause_of_death, cause_of_death_source, country_of_residence_at_enrollment, days_to_birth , days_to_death , ethnicity, gender , occupation_duration_years, premature_at_birth, race, vital_status, weeks_gestation_at_birth , year_of_birth , year_of_death , age_at_diagnosis , ajcc_clinical_m, ajcc_clinical_n, ajcc_clinical_stage, ajcc_clinical_t, ajcc_pathologic_m, ajcc_pathologic_n, ajcc_pathologic_stage, ajcc_pathologic_t, ajcc_staging_system_edition, anaplasia_present, anaplasia_present_type, ann_arbor_b_symptoms, ann_arbor_clinical_stage, ann_arbor_extranodal_involvement, ann_arbor_pathologic_stage, best_overall_response, breslow_thickness, burkitt_lymphoma_clinical_variant, child_pugh_classification, circumferential_resection_margin, classification_of_tumor, cog_liver_stage, cog_neuroblastoma_risk_group, cog_renal_stage, cog_rhabdomyosarcoma_risk_group, days_to_best_overall_response, days_to_diagnosis , days_to_last_follow_up , days_to_last_known_disease_status, days_to_recurrence, eln_risk_classification, enneking_msts_grade, enneking_msts_metastasis, enneking_msts_stage, enneking_msts_tumor_site, esophageal_columnar_dysplasia_degree, esophageal_columnar_metaplasia_present, figo_stage, figo_staging_edition_year, first_symptom_prior_to_diagnosis, gastric_esophageal_junction_involvement, gleason_grade_group, gleason_grade_tertiary, gleason_patterns_percent, goblet_cells_columnar_mucosa_present, greatest_tumor_dimension, gross_tumor_weight, icd_10_code, igcccg_stage, inpc_grade, inpc_histologic_group, inrg_stage, inss_stage, international_prognostic_index, irs_group, irs_stage, ishak_fibrosis_score, iss_stage, largest_extrapelvic_peritoneal_focus, last_known_disease_status, laterality, lymph_node_involved_site, lymph_nodes_positive, lymph_nodes_tested, lymphatic_invasion_present, margin_distance, margins_involved_site, masaoka_stage, medulloblastoma_molecular_classification, metastasis_at_diagnosis, metastasis_at_diagnosis_site, method_of_diagnosis, micropapillary_features, mitosis_karyorrhexis_index, mitotic_count, morphology, non_nodal_regional_disease, non_nodal_tumor_deposits, ovarian_specimen_status, ovarian_surface_involvement, papillary_renal_cell_type, percent_tumor_invasion, perineural_invasion_present, peripancreatic_lymph_nodes_positive, peripancreatic_lymph_nodes_tested, peritoneal_fluid_cytological_status, pregnant_at_diagnosis, primary_diagnosis, primary_gleason_grade, prior_malignancy, prior_treatment, progression_or_recurrence, residual_disease, satellite_nodule_present, secondary_gleason_grade, site_of_resection_or_biopsy, sites_of_involvement, supratentorial_localization, synchronous_malignancy, tissue_or_organ_of_origin, transglottic_extension, tumor_confined_to_organ_of_origin, tumor_depth, tumor_focality, tumor_grade, tumor_largest_dimension_diameter, tumor_regression_grade, tumor_stage, vascular_invasion_present, vascular_invasion_type, weiss_assessment_score, who_cns_grade, who_nte_grade, wilms_tumor_histologic_subtype, year_of_diagnosis , chemo_concurrent_to_radiation, days_to_treatment_end, days_to_treatment_start, initial_disease_status, number_of_cycles, reason_treatment_ended, regimen_or_line_of_therapy, therapeutic_agents, treatment_anatomic_site, treatment_arm, treatment_dose, treatment_dose_units, treatment_effect, treatment_effect_indicator, treatment_frequency, treatment_type, treatment_or_therapy, treatment_outcome, treatment_type, wilms_tumor_histologic_subtype2, year_of_diagnosis2, chemo_concurrent_to_radiation2, days_to_treatment_end2, days_to_treatment_start2, initial_disease_status2, number_of_cycles2, reason_treatment_ended2, regimen_or_line_of_therapy2, therapeutic_agents2, treatment_anatomic_site2, treatment_arm2, treatment_dose2, treatment_dose_units2, treatment_effect2, treatment_effect_indicator2, treatment_frequency2, treatment_type2, treatment_or_therapy2, treatment_outcome2, treatment_type2, wilms_tumor_histologic_subtype3, year_of_diagnosis3, chemo_concurrent_to_radiation3, days_to_treatment_end3, days_to_treatment_start3, initial_disease_status3, number_of_cycles3, reason_treatment_ended3, regimen_or_line_of_therapy3, therapeutic_agents3, treatment_anatomic_site3, treatment_arm3, treatment_dose3, treatment_dose_units3, treatment_effect3, treatment_effect_indicator3, treatment_frequency3, treatment_type3, treatment_or_therapy3, treatment_outcome3, treatment_type3".split(", ")
desiredFields = ["Hugo_Symbol", "UniProt_AApos", "Entrez_Gene_Id", "SwissProt_entry_Id", "i_COSMIC_Gene", "Variant_Classification", "Variant_Type", "Protein_Change", "i_dbNSFP_CADD_phred", "i_dbNSFP_CADD_raw", "i_dbNSFP_CADD_raw_rankscore", "i_dbNSFP_ESP6500_AA_AF", "i_dbNSFP_ESP6500_EA_AF", "i_dbNSFP_Ensembl_geneid", "i_dbNSFP_Ensembl_transcriptid", "i_dbNSFP_FATHMM_pred", "i_dbNSFP_FATHMM_rankscore", "i_dbNSFP_FATHMM_score", "i_dbNSFP_GERPpp_NR", "i_dbNSFP_GERPpp_RS", "i_dbNSFP_GERPpp_RS_rankscore", "i_dbNSFP_Interpro_domain", "i_dbNSFP_LRT_Omega", "i_dbNSFP_LRT_converted_rankscore", "i_dbNSFP_LRT_pred", "i_dbNSFP_LRT_score", "i_dbNSFP_LR_pred", "i_dbNSFP_LR_rankscore", "i_dbNSFP_LR_score", "i_dbNSFP_MutationAssessor_pred", "i_dbNSFP_MutationAssessor_rankscore", "i_dbNSFP_MutationAssessor_score", "i_dbNSFP_MutationTaster_converted_rankscore", "i_dbNSFP_MutationTaster_pred", "i_dbNSFP_MutationTaster_score", "i_dbNSFP_Polyphen2_HDIV_pred", "i_dbNSFP_Polyphen2_HDIV_rankscore", "i_dbNSFP_Polyphen2_HDIV_score", "i_dbNSFP_Polyphen2_HVAR_pred", "i_dbNSFP_Polyphen2_HVAR_rankscore", "i_dbNSFP_Polyphen2_HVAR_score", "i_dbNSFP_RadialSVM_pred", "i_dbNSFP_RadialSVM_rankscore", "i_dbNSFP_RadialSVM_score", "i_dbNSFP_Reliability_index", "i_dbNSFP_SIFT_converted_rankscore", "i_dbNSFP_SIFT_pred", "i_dbNSFP_SIFT_score", "i_dbNSFP_SLR_test_statistic", "i_dbNSFP_SiPhy_29way_logOdds", "i_dbNSFP_SiPhy_29way_logOdds_rankscore", "i_dbNSFP_SiPhy_29way_pi", "i_dbNSFP_UniSNP_ids"]

def getGeneData (Hugo_Symbol, AA_change):
	seenIDs = {}

	mutationDB = sqlite3.connect('filledDB/mutationData.db')
	mcsr = mutationDB.cursor()
	clinicalDB = sqlite3.connect('filledDB/clinical_data')
	ccsr = clinicalDB.cursor()

	mcsr.execute('SELECT * FROM "mutations" WHERE "Hugo_Symbol" LIKE \'' + Hugo_Symbol + '\' AND "Protein_Change" LIKE \'' + AA_change + '\'')
	output = mcsr.fetchall()
	
	for row in output:

		tmp = row[len(desiredFields)+1].split(", ")
		# print (tmp)
		for patName in tmp:
			try:
				seenIDs[patName].append(row)
			except:

				output2 = ccsr.execute('SELECT * FROM "cases" WHERE "case_submitter_id" LIKE \'' + patName + '\'')
				# output2 = ccsr.fetchall()
				for caseData in output2:
					seenIDs.update({patName: [caseData, row]})
					# print("here")
					# break

			
	for k in seenIDs.keys():
		# print(k)
		caseData = seenIDs[k][0]
		print("\t".join(str(d) for d in caseData), end = "\t")
		for i in range(1, len(seenIDs[k])):
			tmpRow = seenIDs[k][i]
			print("\t".join(str(d) for d in tmpRow[1:]), end = "\t")
	print("")




def getGeneData2 (Hugo_Symbol, AA_change):
	global seenIDs
	seenIDs = {}

	mutationDB = sqlite3.connect('filledDB/mutationData.db')
	mcsr = mutationDB.cursor()
	clinicalDB = sqlite3.connect('filledDB/clinical_data')
	ccsr = clinicalDB.cursor()

	mcsr.execute('SELECT * FROM "mutations" WHERE "Hugo_Symbol" LIKE \'' + Hugo_Symbol + '\' AND "Protein_Change" LIKE \'' + AA_change + '\'')
	output = mcsr.fetchall()
	
	for row in output:

		tmp = row[len(desiredFields)+1].split(", ")
		# print (tmp)
		for patName in tmp:
			try:
				seenIDs[patName].append([caseData, row])
			except:

				output2 = ccsr.execute('SELECT * FROM "cases" WHERE "case_submitter_id" LIKE \'' + patName + '\'')
				# output2 = ccsr.fetchall()
				for caseData in output2:
					seenIDs.update({patName: [[caseData, row]]})
					# print("here")
					# break
	print(AA_change)
	for k in seenIDs.keys():
		caseData = seenIDs[k]
		for row in caseData:
			print("\t".join(str(d) for d in row[0]), end = "\t")
			print("\t".join(str(d) for d in row[1]))
			# print("")
		# print("\t".join(str(d) for d in caseData), end = "\t")
		# for i in range(1, len(seenIDs[k])):
		# 	tmpRow = seenIDs[k][i]
		# 	print("\t".join(str(d) for d in tmpRow[1:]), end = "\t")
	print("")


print("\t".join(forwardTable) + "\t" + (("\t".join(desiredFields))))
aa_vars = ["p.E545K", "p.H1047R", "p.E542K", "p.R88Q", "p.H1047L", "p.N345K", "p.E726K", "p.G118D", "p.C420R", "p.Q546R", "p.E453K", "p.E545A", "p.R108H", "p.Q546K", "p.R93Q", "p.H1047Y", "p.M1043V", "p.E81K", "p.Q546P", "p.Y1021C", "p.G1049R", "p.K111E", "p.R38C", "p.R38H", "p.G106V", "p.T1025A", "p.E545Q", "p.E365K", "p.M1043I", "p.R93W", "p.E545G", "p.M1043I", "p.G1007R", "p.V344M", "p.R357Q", "p.K111N", "p.D350N", "p.V344G", "p.E453Q", "p.E970K", "p.C901F", "p.R115L", "p.E39K", "p.E600K", "p.D350G", "p.E418K", "p.E542A", "p.Q546E", "p.C604R", "p.N345I", "p.P539R", "p.P471L", "p.D939G", "p.E542V", "p.G364R", "p.E545D", "p.K111N", "p.F83S", "p.F667L", "p.G106R", "p.Y1021H", "p.N1044K", "p.E542G", "p.N1044K", "p.D725G", "p.N345T", "p.C378F", "p.M1043I", "p.R115P", "p.P104L", "p.E542Q", "p.Q75E", "p.M1004V", "p.R852Q", "p.M1004I", "p.R992P", "p.C378R", "p.H1047Q", "p.P366R", "p.G106S", "p.P104T", "p.N1044S", "p.H1047Q", "p.M1004I", "p.H1048R", "p.R38S", "p.C378Y", "p.R401Q", "p.G106D", "p.Q546H", "p.M282V", "p.M1043L", "p.T1052K", "p.V151M", "p.L339I", "p.P449S", "p.E545D", "p.G914R", "p.N107S", "p.L569I", "p.R818C", "p.E365V", "p.N1044Y", "p.E103G", "p.R617Q", "p.R979G", "p.G865D", "p.R951C", "p.S379T", "p.L956F", "p.G903E", "p.Y392H", "p.P953S", "p.P609H", "p.E722K", "p.Y182H", "p.E39G", "p.M123I", "p.T957P", "p.D520V", "p.D603H", "p.P266T", "p.D1017N", "p.L445I", "p.N380S", "p.S773F", "p.E737K", "p.W11S", "p.D1029N", "p.Q643H", "p.W446C", "p.L929M", "p.M732I", "p.C971R", "p.L866F", "p.E474A", "p.E978K", "p.P471A", "p.G12D", "p.E726G", "p.D84H", "p.Q546H", "p.S292I", "p.P449L", "p.T322A", "p.F744I", "p.M1040I", "p.P124A", "p.I1058M", "p.R777M", "p.E710Q", "p.I841V", "p.D258N", "p.F1002L", "p.L156H", "p.R519G", "p.K337N", "p.D454Y", "p.E600V", "p.C90Y", "p.V636L", "p.F684L", "p.E547K", "p.E849K", "p.E545V", "p.H495Y", "p.E547D", "p.E417Q", "p.Q958R", "p.L339V", "p.R398H", "p.K711N", "p.G363A", "p.I351T", "p.L239R", "p.G359R", "p.R38L", "p.E80K", "p.R818H", "p.N467K", "p.K886E", "p.D746Y", "p.G1009A", "p.Q981R", "p.M1043T", "p.F945C", "p.L658F", "p.F930V", "p.M1005V", "p.Y432C", "p.D300V", "p.D725N", "p.Y606C", "p.R274K", "p.Q296E", "p.C90G", "p.W328S", "p.N345Y", "p.N345S", "p.L1006R", "p.S66C", "p.N107T", "p.E474K", "p.S405F", "p.T1025S", "p.R412Q", "p.I112F", "p.E522A", "p.R108L", "p.I13T", "p.P27T", "p.K111R", "p.R310C", "p.M1040V", "p.V71I", "p.G1049S", "p.E791Q", "p.D589N", "p.F909C", "p.D1029H", "p.L766F", "p.L989V", "p.C905S", "p.G451V", "p.E65K", "p.G451R", "p.L671V", "p.D1045V", "p.A224S", "p.Y165H", "p.S499F", "p.R335G", "p.S576Y", "p.A1020T", "p.L752V", "p.R19I", "p.N170S", "p.V101A", "p.H495L", "p.W11L", "p.R617W", "p.L628R", "p.M811I", "p.C90R", "p.C36Y", "p.R93P", "p.E85G", "p.S629C", "p.E1012Q", "p.S673T", "p.R770Q", "p.R777K", "p.F1016C", "p.F614I", "p.T342S", "p.L866W", "p.R140Q", "p.G914R", "p.M1004R", "p.V344A", "p.T86S", "p.E1037K", "p.W783L", "p.S1015Y", "p.N345H", "p.L531V", "p.I406V", "p.Q879R", "p.L279I", "p.P104R", "p.E417K", "p.R693H", "p.P57L"]

for aa in aa_vars:
	# print(aa, end = "	")
	# print(aa)
	getGeneData2("PIK3CA", aa)
	# print("")
