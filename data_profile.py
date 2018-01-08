# -*- coding: utf-8 -*-
# Instructions for execution: 3 arguments to be supplied at run time
# 1.driver name (e.g. Fabry etc.)
# 2.driver type (diagnosis/procedure/drugs)
# 3.driver codes file (text file)

import psycopg2
import pandas as pd
import time, datetime
import sqlalchemy
import sys
from NDC_Label_And_Code_Mapping_v2 import ndc_code_to_label_mapping

#%%
def execute_sql(sql):
    conn = psycopg2.connect(database="mydb", user="postgres", password="satyam", host="127.0.0.1", port="5432")
    cur = conn.cursor()
    cur.execute(sql)
    rows = cur.fetchall()
    # Extract the column names
    col_names = []
    for elt in cur.description:
        col_names.append(elt[0])
    result = pd.DataFrame(rows, columns=col_names)
    return(result)

#%%
def write_to_table(dataframe, tablename):
    '''Returns a connection and a metadata object'''
    user="postgres"; password="satyam"; host="127.0.0.1"; port="5432"; db="mydb"
    url  = 'postgresql://{}:{}@{}:{}/{}'
    url = url.format(user, password, host, port, db)
    # The return value of create_engine() is our connection object
    con = sqlalchemy.create_engine(url, client_encoding='utf8')
    # We then bind the connection to MetaData()
    meta = sqlalchemy.MetaData(bind=con, reflect=True)
    meta = sqlalchemy.MetaData(con, schema='public')
    meta.reflect()
    pdsql = pd.io.sql.SQLDatabase(con, meta=meta)
    pdsql.to_sql(dataframe, tablename, if_exists='append', index=False)
    return ("DONE")

#%%
def update_sql(sql):
    conn = psycopg2.connect(database="mydb", user="postgres", password="satyam", host="127.0.0.1", port="5432")
    updated_rows = 0
    try:
        cur = conn.cursor()
        # execute the UPDATE  statement
        cur.execute(sql)
        # get the number of updated rows
        updated_rows = cur.rowcount
        # Commit the changes to the database
        conn.commit()
        # Close communication with the PostgreSQL database
        cur.close()
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
    finally:
        if conn is not None:
            conn.close()
    return (updated_rows)

#%%
def mimic_diagnosis_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
    code_list = str(code_list).strip('[]')
    # fetch all data required for given codes
    mimic_diag_codes_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_DESC FROM mimic_diagnoses A, mimic_d_icd_diagnoses B WHERE A.ICD9_CODE=B.ICD9_CODE AND A.ICD9_CODE IN (" + code_list + ")"
    mimic_diag_codes = execute_sql(mimic_diag_codes_sql)
    mimic_diag = mimic_diag_codes[['subject_id', 'icd9_code', 'icd9_desc']]
    mimic_diag = mimic_diag.drop_duplicates()
    # patient count per icd9_code
    mimic_diag_count = mimic_diag.groupby(['icd9_code', 'icd9_desc'],as_index = False).size().reset_index()
    mimic_diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
    # patient list
    mimic_diag_subject = list(mimic_diag.subject_id.unique())
    # extract all data for the patients fetched above
    mimic_diag_all_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_DESC FROM mimic_diagnoses A, mimic_d_icd_diagnoses B WHERE A.ICD9_CODE=B.ICD9_CODE AND A.SUBJECT_ID IN (" + str(mimic_diag_subject).strip('[]') + ")"
    mimic_diag_all = execute_sql(mimic_diag_all_sql)
    return(mimic_diag_all, mimic_diag_count, mimic_diag_subject)

#%%
def mimic_diagnosis_from_subject_id(subject_ids):
    mimic_subject = []
    for c in subject_ids:
        mimic_subject.append(str(c))
    mimic_subject = str(mimic_subject).strip('[]')
    mimic_diag_all_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_DESC FROM mimic_diagnoses A, mimic_d_icd_diagnoses B WHERE A.ICD9_CODE=B.ICD9_CODE AND A.SUBJECT_ID IN (" + mimic_subject + ")"
    # fetch all data for required subjects
    mimic_diag_all = execute_sql(mimic_diag_all_sql)
    mimic_diag = mimic_diag_all[['subject_id', 'icd9_code', 'icd9_desc']]
    mimic_diag = mimic_diag.drop_duplicates()
    # patient count per icd9_code
    mimic_diag_count = mimic_diag.groupby(['icd9_code', 'icd9_desc'],as_index = False).size().reset_index()
    mimic_diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
    return(mimic_diag_all, mimic_diag_count)

#%%
def cms_diagnosis_summary(cmsin_diag_all, cmsout_diag_all, codes=None):
    # append the two dataframes
    cms_diag_all = cmsin_diag_all.append(pd.DataFrame(data = cmsout_diag_all), ignore_index=True)
    # to be executed if codes are supplied
    if codes is not None:
        cmsin_diag_all = cmsin_diag_all.loc[cmsin_diag_all['icd9_code'].isin(codes)]
        cmsout_diag_all = cmsout_diag_all.loc[cmsout_diag_all['icd9_code'].isin(codes)]
    # combine the two data frames
    cms_diag = cmsin_diag_all.append(pd.DataFrame(data = cmsout_diag_all), ignore_index=True)
    cms_diag = cms_diag[['desynpuf_id', 'icd9_code', 'icd9_desc']]
    # remove duplicates
    cms_combined = cms_diag.drop_duplicates()
    # patient count per icd9_code from both inpatient & outpatient
    cms_diag_count = cms_combined.groupby(['icd9_code', 'icd9_desc'],as_index = False).size().reset_index()
    cms_diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
    # patient list
    cms_diag_desynpuf = list(cms_combined.desynpuf_id.unique())
    return(cms_diag_all, cms_diag_count, cms_diag_desynpuf)

#%%
def cms_diagnosis_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
    code_list = str(code_list).strip('[]')
    cms_diag_list = ["icd9_dgns_cd_1", "icd9_dgns_cd_2", "icd9_dgns_cd_3", "icd9_dgns_cd_4", "icd9_dgns_cd_5", "icd9_dgns_cd_6", "icd9_dgns_cd_7", "icd9_dgns_cd_8", "icd9_dgns_cd_9", "icd9_dgns_cd_10"]
    # fetch all data for required codes
    cms_diag_sql = []
    for cms in cms_diag_list:
        cms_diag_sql_temp = cms + " in (" + code_list + ")"
        cms_diag_sql.append(cms_diag_sql_temp)
    cms_diag_sql = ' or '.join(cms_diag_sql)
    # inpatient
    cmsin_diag_sql = "select a.desynpuf_id, a.colnames, 'in' as io_flag, a.icd9_code, b.long_title as icd9_desc from(SELECT desynpuf_id, unnest(array['icd9_dgns_cd_1', 'icd9_dgns_cd_2', 'icd9_dgns_cd_3', 'icd9_dgns_cd_4', 'icd9_dgns_cd_5', 'icd9_dgns_cd_6', 'icd9_dgns_cd_7', 'icd9_dgns_cd_8', 'icd9_dgns_cd_9', 'icd9_dgns_cd_10']) AS colnames,   unnest(array[icd9_dgns_cd_1, icd9_dgns_cd_2, icd9_dgns_cd_3, icd9_dgns_cd_4, icd9_dgns_cd_5, icd9_dgns_cd_6, icd9_dgns_cd_7, icd9_dgns_cd_8, icd9_dgns_cd_9, icd9_dgns_cd_10]) AS icd9_code FROM cms_ipclaims_combined where " + cms_diag_sql + ")a left outer join mimic_d_icd_diagnoses b on a.icd9_code = b.icd9_code" 
    cmsin_diag_all = execute_sql(cmsin_diag_sql)
    cmsin_diag_all = cmsin_diag_all.dropna(subset=['icd9_code'])
    # outpatient
    cmsout_diag_sql = "select a.desynpuf_id, 'out' as io_flag, a.colnames, a.icd9_code, b.long_title as icd9_desc from(SELECT desynpuf_id, unnest(array['icd9_dgns_cd_1', 'icd9_dgns_cd_2', 'icd9_dgns_cd_3', 'icd9_dgns_cd_4', 'icd9_dgns_cd_5', 'icd9_dgns_cd_6', 'icd9_dgns_cd_7', 'icd9_dgns_cd_8', 'icd9_dgns_cd_9', 'icd9_dgns_cd_10']) AS colnames,   unnest(array[icd9_dgns_cd_1, icd9_dgns_cd_2, icd9_dgns_cd_3, icd9_dgns_cd_4, icd9_dgns_cd_5, icd9_dgns_cd_6, icd9_dgns_cd_7, icd9_dgns_cd_8, icd9_dgns_cd_9, icd9_dgns_cd_10]) AS icd9_code FROM cms_opclaims_combined where " + cms_diag_sql + ")a left outer join mimic_d_icd_diagnoses b on a.icd9_code = b.icd9_code" 
    cmsout_diag_all = execute_sql(cmsout_diag_sql)
    cmsout_diag_all = cmsout_diag_all.dropna(subset=['icd9_code'])
    # call the summary function
    cms_diag_all, cms_diag_count, cms_diag_desynpuf = cms_diagnosis_summary(cmsin_diag_all, cmsout_diag_all, codes)
    return(cms_diag_all, cms_diag_count, cms_diag_desynpuf)

#%%
def cms_diagnosis_from_desynpuf_id(desynpuf_ids):
    cms_desynpuf = []
    for c in desynpuf_ids:
        cms_desynpuf.append(str(c))
    cms_desynpuf = str(cms_desynpuf).strip('[]')
    # inpatient
    cmsin_diag_sql = "select a.desynpuf_id, 'in' as io_flag, a.colnames, a.icd9_code, b.long_title as icd9_desc from(SELECT desynpuf_id, unnest(array['icd9_dgns_cd_1', 'icd9_dgns_cd_2', 'icd9_dgns_cd_3', 'icd9_dgns_cd_4', 'icd9_dgns_cd_5', 'icd9_dgns_cd_6', 'icd9_dgns_cd_7', 'icd9_dgns_cd_8', 'icd9_dgns_cd_9', 'icd9_dgns_cd_10']) AS colnames, unnest(array[icd9_dgns_cd_1, icd9_dgns_cd_2, icd9_dgns_cd_3, icd9_dgns_cd_4, icd9_dgns_cd_5, icd9_dgns_cd_6, icd9_dgns_cd_7, icd9_dgns_cd_8, icd9_dgns_cd_9, icd9_dgns_cd_10]) AS icd9_code FROM cms_ipclaims_combined where desynpuf_id in (" + cms_desynpuf + "))a left outer join mimic_d_icd_diagnoses b on a.icd9_code = b.icd9_code"
    cmsin_diag_all = execute_sql(cmsin_diag_sql)
    cmsin_diag_all = cmsin_diag_all.dropna(subset=['icd9_code'])

    # outpatient
    cmsout_diag_sql = "select a.desynpuf_id, 'out' as io_flag, a.colnames, a.icd9_code, b.long_title as icd9_desc from(SELECT desynpuf_id, unnest(array['icd9_dgns_cd_1', 'icd9_dgns_cd_2', 'icd9_dgns_cd_3', 'icd9_dgns_cd_4', 'icd9_dgns_cd_5', 'icd9_dgns_cd_6', 'icd9_dgns_cd_7', 'icd9_dgns_cd_8', 'icd9_dgns_cd_9', 'icd9_dgns_cd_10']) AS colnames, unnest(array[icd9_dgns_cd_1, icd9_dgns_cd_2, icd9_dgns_cd_3, icd9_dgns_cd_4, icd9_dgns_cd_5, icd9_dgns_cd_6, icd9_dgns_cd_7, icd9_dgns_cd_8, icd9_dgns_cd_9, icd9_dgns_cd_10]) AS icd9_code FROM cms_opclaims_combined where desynpuf_id in (" + cms_desynpuf + "))a left outer join mimic_d_icd_diagnoses b on a.icd9_code = b.icd9_code"
    cmsout_diag_all = execute_sql(cmsout_diag_sql)
    cmsout_diag_all = cmsout_diag_all.dropna(subset=['icd9_code'])
    cms_diag_all, cms_diag_count, cms_diag_desynpuf = cms_diagnosis_summary(cmsin_diag_all, cmsout_diag_all)
    return(cms_diag_all, cms_diag_count)

#%%
def mimic_procedures_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
    code_list = str(code_list).strip('[]')
    # fetch all data required for given codes
    mimic_proc_codes_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_PROC_DESC FROM mimic_procedures A, mimic_d_icd_procedures B WHERE A.ICD9_PROCEDURE=B.ICD9_PROCEDURE AND A.ICD9_PROCEDURE IN (" + code_list + ")"
    mimic_proc_codes = execute_sql(mimic_proc_codes_sql)
    mimic_proc = mimic_proc_codes[['subject_id', 'icd9_procedure', 'icd9_proc_desc']]
    mimic_proc = mimic_proc.drop_duplicates()
    # patient count per icd9_code
    mimic_proc_count = mimic_proc.groupby(['icd9_procedure', 'icd9_proc_desc'], as_index = False).size().reset_index()
    mimic_proc_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
    # patient list
    mimic_proc_subject = list(mimic_proc.subject_id.unique())
	# extract all data for the patients fetched above
    mimic_proc_all_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_PROC_DESC FROM mimic_procedures A, mimic_d_icd_procedures B WHERE A.ICD9_PROCEDURE=B.ICD9_PROCEDURE AND A.SUBJECT_ID IN (" + mimic_proc_subject + ")"
    mimic_proc_all = execute_sql(mimic_proc_all_sql)
    return(mimic_proc_all, mimic_proc_count, mimic_proc_subject)

#%%
def mimic_procedures_from_subject_id(subject_ids):
    mimic_subject = []
    for c in subject_ids:
        mimic_subject.append(str(c))
    mimic_subject = str(mimic_subject).strip('[]')
    mimic_proc_all_sql = "SELECT A.*, B.LONG_TITLE AS ICD9_PROC_DESC FROM mimic_procedures A, mimic_d_icd_procedures B WHERE A.ICD9_PROCEDURE=B.ICD9_PROCEDURE AND A.SUBJECT_ID IN (" + mimic_subject + ")"
    # fetch all data for required subjects
    mimic_proc_all = execute_sql(mimic_proc_all_sql)
    mimic_proc = mimic_proc_all[['subject_id', 'icd9_procedure', 'icd9_proc_desc']]
    mimic_proc = mimic_proc.drop_duplicates()
    # patient count per procedure code
    mimic_proc_count = mimic_proc.groupby(['icd9_procedure', 'icd9_proc_desc'],as_index = False).size().reset_index()
    mimic_proc_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
    return(mimic_proc_all, mimic_proc_count)

#%%
def cms_procedures_summary(cmsin_proc_all, cmsout_proc_all, codes=None):
    # append the two dataframes
    cms_proc_all = cmsin_proc_all.append(pd.DataFrame(data = cmsout_proc_all), ignore_index=True)
    # to be executed if codes are supplied
    if codes is not None:
        cmsin_proc_all = cmsin_proc_all.loc[cmsin_proc_all['value'].isin(codes)]
        cmsout_proc_all = cmsout_proc_all.loc[cmsout_proc_all['value'].isin(codes)]
    # combine the two data frames
    cms_proc = cmsin_proc_all.append(pd.DataFrame(data = cmsout_proc_all), ignore_index=True)
    cms_proc = cms_proc[['desynpuf_id', 'icd9_procedure', 'icd9_proc_desc']]
    # remove duplicates
    cms_proc_combined = cms_proc.drop_duplicates()
    # patient count per icd9_code from both inpatient & outpatient
    cms_proc_count = cms_proc_combined.groupby(['icd9_procedure', 'icd9_proc_desc'],as_index = False).size().reset_index()
    cms_proc_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
    # patient list
    cms_proc_desynpuf = list(cms_proc_combined.desynpuf_id.unique())
    return(cms_proc_all, cms_proc_count, cms_proc_desynpuf)

#%%
def cms_procedures_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
        code_list = str(code_list).strip('[]')
    cms_proc_list = ["icd9_prcdr_cd_1", "icd9_prcdr_cd_2", "icd9_prcdr_cd_3", "icd9_prcdr_cd_4", "icd9_prcdr_cd_5", "icd9_prcdr_cd_6"]
    # fetch all data for required codes
    cms_proc_sql = []
    for cms in cms_proc_list:
        cms_proc_sql_temp = cms + " in (" + code_list + ")"
        cms_proc_sql.append(cms_proc_sql_temp)
    cms_proc_sql = ' or '.join(cms_proc_sql)
    # inpatient
    cmsin_proc_sql = "select a.desynpuf_id, 'in' as io_flag, a.colnames, a.icd9_procedure, b.long_title as icd9_proc_desc from(SELECT desynpuf_id, unnest(array['icd9_prcdr_cd_1', 'icd9_prcdr_cd_2', 'icd9_prcdr_cd_3', 'icd9_prcdr_cd_4', 'icd9_prcdr_cd_5', 'icd9_prcdr_cd_6']) AS colnames, unnest(array[icd9_prcdr_cd_1, icd9_prcdr_cd_2, icd9_prcdr_cd_3, icd9_prcdr_cd_4, icd9_prcdr_cd_5, icd9_prcdr_cd_6]) AS icd9_procedure FROM cms_ipclaims_combined where " + cms_proc_sql + ")a left outer join mimic_d_icd_procedures b on a.icd9_procedure = b.icd9_procedure"
    cmsin_proc_all = execute_sql(cmsin_proc_sql)
    cmsin_proc_all = cmsin_proc_all.dropna(subset=['icd9_procedure'])

    # outpatient
    cmsout_proc_sql = "select a.desynpuf_id, 'out' as io_flag, a.colnames, a.icd9_procedure, b.long_title as icd9_proc_desc from(SELECT desynpuf_id, unnest(array['icd9_prcdr_cd_1', 'icd9_prcdr_cd_2', 'icd9_prcdr_cd_3', 'icd9_prcdr_cd_4', 'icd9_prcdr_cd_5', 'icd9_prcdr_cd_6']) AS colnames, unnest(array[icd9_prcdr_cd_1, icd9_prcdr_cd_2, icd9_prcdr_cd_3, icd9_prcdr_cd_4, icd9_prcdr_cd_5, icd9_prcdr_cd_6]) AS icd9_procedure FROM cms_opclaims_combined where " + cms_proc_sql + ")a left outer join mimic_d_icd_procedures b on a.icd9_procedure = b.icd9_procedure"    
    cmsout_proc_all = execute_sql(cmsout_proc_sql)
    cmsout_proc_all = cmsout_proc_all.dropna(subset=['icd9_procedure'])
    # call the summary function
    cms_proc_all, cms_proc_count, cms_proc_desynpuf = cms_procedures_summary(cmsin_proc_all, cmsout_proc_all, codes)
    return(cms_proc_all, cms_proc_count, cms_proc_desynpuf)

#%%
def cms_procedures_from_desynpuf_id(desynpuf_ids):
    cms_desynpuf = []
    for c in desynpuf_ids:
        cms_desynpuf.append(str(c))
    cms_desynpuf = str(cms_desynpuf).strip('[]')
    # inpatient
    cmsin_proc_sql = "select a.desynpuf_id, 'in' as io_flag, a.colnames, a.icd9_procedure, b.long_title as icd9_proc_desc from(SELECT desynpuf_id, unnest(array['icd9_prcdr_cd_1', 'icd9_prcdr_cd_2', 'icd9_prcdr_cd_3', 'icd9_prcdr_cd_4', 'icd9_prcdr_cd_5', 'icd9_prcdr_cd_6']) AS colnames, unnest(array[icd9_prcdr_cd_1, icd9_prcdr_cd_2, icd9_prcdr_cd_3, icd9_prcdr_cd_4, icd9_prcdr_cd_5, icd9_prcdr_cd_6]) AS icd9_procedure FROM cms_ipclaims_combined where desynpuf_id in (" + cms_desynpuf + "))a left outer join mimic_d_icd_procedures b on a.icd9_procedure = b.icd9_procedure"
    cmsin_proc_all = execute_sql(cmsin_proc_sql)
    cmsin_proc_all = cmsin_proc_all.dropna(subset=['icd9_procedure'])

    # outpatient
    cmsout_proc_sql = "select a.desynpuf_id, 'out' as io_flag, a.colnames, a.icd9_procedure, b.long_title as icd9_proc_desc from(SELECT desynpuf_id, unnest(array['icd9_prcdr_cd_1', 'icd9_prcdr_cd_2', 'icd9_prcdr_cd_3', 'icd9_prcdr_cd_4', 'icd9_prcdr_cd_5', 'icd9_prcdr_cd_6']) AS colnames, unnest(array[icd9_prcdr_cd_1, icd9_prcdr_cd_2, icd9_prcdr_cd_3, icd9_prcdr_cd_4, icd9_prcdr_cd_5, icd9_prcdr_cd_6]) AS icd9_procedure FROM cms_opclaims_combined where desynpuf_id in (" + cms_desynpuf + "))a left outer join mimic_d_icd_procedures b on a.icd9_procedure = b.icd9_procedure"
    cmsout_proc_all = execute_sql(cmsout_proc_sql)
    cmsout_proc_all = cmsout_proc_all.dropna(subset=['icd9_procedure'])
    cms_proc_all, cms_proc_count, cms_proc_desynpuf = cms_procedures_summary(cmsin_proc_all, cmsout_proc_all)
    return(cms_proc_all, cms_proc_count)

#%%
def mimic_drugs_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
    code_list = str(code_list).strip('[]')
    # fetch all data for required subjects
    mimic_drug_all_sql = "SELECT * FROM mimic_prescriptions WHERE NDC IN (%s) and drug_type = 'MAIN'" %(code_list)
    mimic_drug_all = execute_sql(mimic_drug_all_sql)
    mimic_drug = mimic_drug_all[['subject_id', 'ndc', 'drug']]
    mimic_drug = mimic_drug.drop_duplicates()
    mimic_drug = mimic_drug.loc[mimic_drug['ndc'] != '0']
    # patient count per drug code
    mimic_drug_count = mimic_drug.groupby(['ndc', 'drug'],as_index = False).size().reset_index()
    mimic_drug_count.columns = ['ndc', 'drug', 'NoOfPatients']
	# patient list
    mimic_drug_subject = list(mimic_drug.subject_id.unique())
    return(mimic_drug_all, mimic_drug_count, mimic_drug_subject)

#%%
def mimic_drugs_from_subject_id(subject_ids):
    mimic_subject = []
    for c in subject_ids:
        mimic_subject.append(str(c))
    mimic_subject = str(mimic_subject).strip('[]')
    mimic_drug_all_sql = "SELECT * FROM mimic_prescriptions WHERE SUBJECT_ID IN (%s)" %(mimic_subject)
    # fetch all data for required subjects
    mimic_drug_all = execute_sql(mimic_drug_all_sql)
    mimic_drug = mimic_drug_all[['subject_id', 'ndc', 'drug']]
    mimic_drug = mimic_drug.drop_duplicates()
    mimic_drug = mimic_drug.loc[mimic_drug['ndc'] != '0']
    # patient count per drug code
    mimic_drug_count = mimic_drug.groupby(['ndc', 'drug'],as_index = False).size().reset_index()
    mimic_drug_count.columns = ['ndc', 'drug', 'NoOfPatients']
    return(mimic_drug_all, mimic_drug_count)

#%%
def cms_drugs_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
        code_list = str(code_list).strip('[]')
    cms_drug_all_sql = "select * from cms_prescriptions_combined where prod_srvc_id in (" + code_list + ")"
    # fetch all data for required subjects
    cms_drug_all = execute_sql(cms_drug_all_sql)
    cms_drug = cms_drug_all[['desynpuf_id', 'prod_srvc_id', 'drug']]
    cms_drug = cms_drug.drop_duplicates()
    cms_drug = cms_drug.loc[cms_drug['prod_srvc_id'] != '0']
    cms_drug_ref = cms_drug[['prod_srvc_id', 'drug']].drop_duplicates()
    cms_drug_ref.columns = ['ndc', 'drug']
    # patient count per drug code
    cms_drug_count = cms_drug.groupby('prod_srvc_id', as_index = False).size().reset_index()
    cms_drug_count.columns = ['ndc', 'NoOfPatients']
    # for missing drug labels
    print("\nFetching labels for missing drug labels...")
    for i, row in cms_drug_count.iterrows():
        if row.drug == '' or pd.isnull(row.drug):
            ndc, ndc_lbl = ndc_code_to_label_mapping(row.ndc)
            if ndc==ndc_lbl or ndc_lbl=="NOTFOUND":
                ndc_12d, ndc_lbl_12d = ndc_code_to_label_mapping(row.ndc.zfill(12))
                if ndc_12d==ndc_lbl_12d or ndc_lbl_12d=="NOTFOUND":
                    pass
                else:
                    ndc = ndc_12d
                    ndc_lbl = ndc_lbl_12d
            cms_drug_count.loc[i,'ndc'], cms_drug_count.loc[i,'drug']  = ndc, ndc_lbl

    # adding labels to count frame
    cms_drug_count = pd.merge(cms_drug_count, cms_drug_ref)
    # adding labels to raw frame
    cms_drug_all = cms_drug_all.drop('drug',1)
    cms_drug_all = pd.merge(cms_drug_all, cms_drug_ref, left_on='prod_srvc_id', right_on='ndc')
	# patient list
    cms_drug_desynpuf = list(cms_drug.desynpuf_id.unique())
    return(cms_drug_all, cms_drug_count, cms_drug_desynpuf)

#%%
def cms_drugs_from_desynpuf_id(desynpuf_ids):
    cms_desynpuf = []
    for c in desynpuf_ids:
        cms_desynpuf.append(str(c))
    cms_desynpuf = str(cms_desynpuf).strip('[]')
    cms_drug_all_sql = "select * from cms_prescriptions_combined where desynpuf_id in (" + cms_desynpuf + ")"
    # fetch all data for required subjects
    cms_drug_all = execute_sql(cms_drug_all_sql)
    cms_drug = cms_drug_all[['desynpuf_id', 'prod_srvc_id', 'drug']]
    cms_drug = cms_drug.drop_duplicates()
    cms_drug = cms_drug.loc[cms_drug['prod_srvc_id'] != '0']
    cms_drug_ref = cms_drug[['prod_srvc_id', 'drug']].drop_duplicates()
    cms_drug_ref.columns = ['ndc', 'drug']
    # patient count per drug code
    cms_drug_count = cms_drug.groupby('prod_srvc_id', as_index = False).size().reset_index()
    cms_drug_count.columns = ['ndc', 'NoOfPatients']
    # for missing drug labels
    print("\nFetching labels for missing drug labels...")
    for i, row in cms_drug_ref.iterrows():
        if row.drug == '' or pd.isnull(row.drug):
            ndc, ndc_lbl = ndc_code_to_label_mapping(row.ndc)
            if ndc==ndc_lbl or ndc_lbl=="NOTFOUND":
                ndc_12d, ndc_lbl_12d = ndc_code_to_label_mapping(row.ndc.zfill(12))
                if ndc_12d==ndc_lbl_12d or ndc_lbl_12d=="NOTFOUND":
                    pass
                else:
                    ndc = ndc_12d
                    ndc_lbl = ndc_lbl_12d
            cms_drug_ref.loc[i,'ndc'], cms_drug_ref.loc[i,'drug']  = ndc, ndc_lbl
    # adding labels to count frame
    cms_drug_count = pd.merge(cms_drug_count, cms_drug_ref)
    # adding labels to raw frame
    cms_drug_all = cms_drug_all.drop('drug',1)
    cms_drug_all = pd.merge(cms_drug_all, cms_drug_ref, left_on='prod_srvc_id', right_on='ndc')
    return(cms_drug_all, cms_drug_count)
                
#%%
def mimic_labs_from_codes(codes):
    code_list = []
    for c in codes:
        code_list.append(str(c))
    code_list = str(code_list).strip('[]')
    # fetch all data for required subjects
    mimic_labs_all_sql = "SELECT * FROM mimic_labevents A, mimic_d_labitems B WHERE A.ITEMID = B.ITEMID AND B.LOINC_CODE IN (" + code_list + ")"
    mimic_labs_all = execute_sql(mimic_labs_all_sql)
    mimic_labs = mimic_labs_all[['subject_id', 'loinc_code', 'label']]
    mimic_labs = mimic_labs.drop_duplicates()
    # patient count per lab code
    mimic_labs_count = mimic_labs.groupby(['loinc_code', 'label'],as_index = False).size().reset_index()
    mimic_labs_count.columns = ['loinc_code', 'label', 'NoOfPatients']
	# patient list
    mimic_labs_subject = list(mimic_labs.subject_id.unique())
    return(mimic_labs_all, mimic_labs_count, mimic_labs_subject)
	
#%%
def mimic_labs_from_subject_id(subject_ids):
    mimic_subject = []
    for c in subject_ids:
        mimic_subject.append(str(c))
    mimic_subject = str(mimic_subject).strip('[]')
    mimic_labs_all_sql = "SELECT A.*, B.LABEL, B.LOINC_CODE FROM mimic_labevents A, mimic_d_labitems B WHERE A.ITEMID = B.ITEMID AND A.SUBJECT_ID IN (" + mimic_subject + ")"
    # fetch all data for required subjects
    mimic_labs_all = execute_sql(mimic_labs_all_sql)
    mimic_labs = mimic_labs_all[['subject_id', 'loinc_code', 'label']]
    mimic_labs = mimic_labs.drop_duplicates()
    # patient count per lab code
    mimic_labs_count = mimic_labs.groupby(['loinc_code', 'label'],as_index = False).size().reset_index()
    mimic_labs_count.columns = ['loinc_code', 'label', 'NoOfPatients']
    return(mimic_labs_all, mimic_labs_count)

#%%
def create_profile(uniq_count, diag_count, mimic_diag_all, cms_diag_all, proc_count, mimic_proc_all, cms_proc_all, drug_count, mimic_drug_all, cms_drug_all, labs_count, mimic_labs_all):
    # a dataframe to add blank columns between adjacent dataframes
    blnk_df = pd.DataFrame({'' : []})
    # create the count profile
    profile_df = pd.concat([uniq_count, blnk_df, diag_count, blnk_df, proc_count, blnk_df, drug_count, blnk_df, labs_count], axis=1)
    # create the combined dataframe
    # Diagnosis
    mimic_diag_df = mimic_diag_all[['subject_id', 'icd9_code', 'icd9_desc']]
    mimic_diag_df.columns = ['patient_id', 'code', 'code_desc']
    cms_diag_df = cms_diag_all[['desynpuf_id', 'icd9_code', 'icd9_desc']]
    cms_diag_df.columns = ['patient_id', 'code', 'code_desc']
    diag_df = mimic_diag_df.append(pd.DataFrame(data = cms_diag_df), ignore_index=True)
    diag_df['code_type'] = 'diagnosis'
    # Procedures
    mimic_proc_df = mimic_proc_all[['subject_id', 'icd9_procedure', 'icd9_proc_desc']]
    mimic_proc_df.columns = ['patient_id', 'code', 'code_desc']
    cms_proc_df = cms_proc_all[['desynpuf_id', 'icd9_procedure', 'icd9_proc_desc']]
    cms_proc_df.columns = ['patient_id', 'code', 'code_desc']
    proc_df = mimic_proc_df.append(pd.DataFrame(data = cms_proc_df), ignore_index=True)
    proc_df['code_type'] = 'procedure'
    # Drugs
    mimic_drug_df = mimic_drug_all[['subject_id', 'ndc', 'drug']]
    mimic_drug_df.columns = ['patient_id', 'code', 'code_desc']
    cms_drug_df = cms_drug_all[['desynpuf_id', 'prod_srvc_id', 'drug']]
    cms_drug_df.columns = ['patient_id', 'code', 'code_desc']
    drug_df = mimic_drug_df.append(pd.DataFrame(data = cms_drug_df), ignore_index=True)
    drug_df['code_type'] = 'drug'
    # Labs
    labs_df = mimic_labs_all[['subject_id', 'loinc_code', 'label']]
    labs_df.columns = ['patient_id', 'code', 'code_desc']
    labs_df['code_type'] = 'labs'
    # combine all the data frames
    combined_df = diag_df.append(pd.DataFrame(data = proc_df), ignore_index=True)
    combined_df = combined_df.append(pd.DataFrame(data = drug_df), ignore_index=True)
    combined_df = combined_df.append(pd.DataFrame(data = labs_df), ignore_index=True)
    combined_df = combined_df.drop_duplicates()
    return(profile_df, combined_df)
    
#%%
def update_patient_combinations(combined_df, profile_df):
    patient_list = combined_df.patient_id.unique()
    combn_summ = []
    for patient_id in patient_list:
        # diagnosis codes
        diag_codes = combined_df[(combined_df.patient_id==patient_id) & (combined_df.code_type=='diagnosis')].code
        diag = []
        for d in diag_codes:
            diag.append(str(d))
        diag = str(diag).strip('[]')
        diag_codes_sql = "select icd9_std from mimic_d_icd_diagnoses where icd9_code in (%s)" %(diag)
        diag_codes = execute_sql(diag_codes_sql)
        diag_codes = "IDX_" + diag_codes.iloc[:,0]
        # drug codes
        drug_codes = combined_df[(combined_df.patient_id==patient_id) & (combined_df.code_type=='drug')].code
        drug_codes = "N_" + drug_codes[drug_codes != '0']
        # lab codes
        lab_codes = combined_df[(combined_df.patient_id==patient_id) & (combined_df.code_type=='labs')].code
        lab_codes = "L_" + lab_codes
        # diag-drug pairs
        diag_drug_list = []
        for diag in diag_codes:
            for drug in drug_codes:
                diag_drug_list.append([diag,drug])
        # diag-lab pairs
        diag_lab_list = []
        for diag in diag_codes:
            for lab in lab_codes:
                diag_lab_list.append([diag,lab])
        # append the two lists
        pt_combn = diag_drug_list + diag_lab_list
        lbl_pt_combn = ['pair_l1','pair_l2']
        pt_combn = pd.DataFrame(pt_combn, columns=lbl_pt_combn)
        pt_combn['patient_id'] = patient_id
        pt_combn['combination_flag'] = 'N'
        pt_combn = pt_combn[['patient_id','pair_l1','pair_l2','combination_flag']]
        # generate combination summary
        combn_summ.append([patient_id, len(diag_codes), len(drug_codes), len(lab_codes), len(diag_drug_list), len(diag_lab_list), len(pt_combn)])
        # store combination pairs in patient_combination table
        print("Writing combinations to table for: %s" %(patient_id))
        write_to_table(pt_combn, 'patient_combinations')
    # combination summary for all patients
    combn_summ = pd.DataFrame(combn_summ)
    combn_summ.columns = ['patient_id','diag_cnt','drug_cnt','lab_cnt','diag_drug','diag_lab','combn_pland']
    # identify the actual combinations to be run
    pat_list = []
    for c in patient_list:
        pat_list.append(str(c))
    pat_list = str(pat_list).strip('[]')
    # update patient_combination table and mark combinations done as 'Y' which are already present
    print("Updating the patient_combination table for existing combinations...")
    update_pat_comb_sql = "update patient_combinations b set combination_flag = 'Y' from (select pair_l1, pair_l2 from combinations_scores group by pair_l1, pair_l2) a where b.pair_l1 = a.pair_l1 and b.pair_l2 = a.pair_l2 and b.patient_id in(%s);" %(pat_list)
    print("Rows updated = ",update_sql(update_pat_comb_sql))
    # to be executed counts
    comb_to_be_exec_sql = "select patient_id, count(*) as comb_exec from patient_combinations where combination_flag='N' and patient_id in (%s) group by patient_id" %(pat_list)
    comb_to_be_exec = execute_sql(comb_to_be_exec_sql)
    combn_summ = pd.merge(combn_summ,comb_to_be_exec,how='inner',on='patient_id')
    # a dataframe to add blank columns between adjacent dataframes
    blnk_df = pd.DataFrame({'' : []})
    # update the profile
    profile_df = pd.concat([profile_df, blnk_df, combn_summ], axis=1)
    # update patient_combination table and mark combinations done as 'Y'
    return(profile_df)

#%%
def out_to_excel(driver_name, profile_df, mimic_diag_all, cms_diag_all, mimic_proc_all, cms_proc_all, mimic_drug_all, cms_drug_all, mimic_labs_all, combined_df):
    # generate file name
    st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%b%d_%H%M')
    fname = "Data_Profiling_" + driver_name + "_" + st + ".xlsx"
    # write output to file
    writer = pd.ExcelWriter(fname, engine='xlsxwriter')
    print(fname)
    # profile summary
    profile_df.to_excel(writer, sheet_name='Data_Profile')
    # diagnosis data
    mimic_diag_all.to_excel(writer, sheet_name='mimic_diagnosis')
    cms_diag_all.to_excel(writer, sheet_name='cms_diagnosis')
    # procedure data
    mimic_proc_all.to_excel(writer, sheet_name='mimic_procedure')
    cms_proc_all.to_excel(writer, sheet_name='cms_procedure')
    # drugs data
    mimic_drug_all.to_excel(writer, sheet_name='mimic_drugs')
    cms_drug_all.to_excel(writer, sheet_name='cms_drugs')
    # labs data
    mimic_labs_all.to_excel(writer, sheet_name='mimic_labs')
    # combined dataframe
    combined_df.to_excel(writer, sheet_name='patients_all_data')
    writer.save()
    # write combined frame to database
    write_to_table(combined_df, 'combinations_input_codes')
    return(writer)
    
#%%
# main function
if __name__ == "__main__":
    path = 'C:/Users/310223340/OneDrive - Philips/OtherProject/Personal_Project/DrugNeighborhood/'
    name = sys.argv[1]
    driver = sys.argv[2]
    print(driver, " = ", name)
    with open(path+sys.argv[3], 'r') as temp_file:
        contents = [line.rstrip('\n') for line in temp_file]
    print("codes = ",contents)
    #### start point: DIAGNOSIS ####
    if driver == 'diagnosis':
        # DIAGNOSIS
        print("Executing mimic diagnosis from codes...\n")
        mimic_diag_all, mimic_diag_count, mimic_diag_subject = mimic_diagnosis_from_codes(contents)
        print("Executing cms diagnosis from codes...\n")
        cms_diag_all, cms_diag_count, cms_diag_desynpuf = cms_diagnosis_from_codes(contents)
        # append frames to get the final count per diagnosis code
        diag_count = mimic_diag_count.append(cms_diag_count, ignore_index = True)
        diag_count = diag_count.groupby(['icd9_code', 'icd9_desc']).sum().reset_index()
        diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
        # number of unique patients
        uniq_pat = len(mimic_diag_subject) + len(cms_diag_desynpuf)
        uniq_count = pd.DataFrame({'NoOfUniqPatients':uniq_pat}, index=[0])
        # PROCEDURES
        print("Executing mimic procedures using subject_ids(diagnosis)...\n")
        mimic_proc_all, mimic_proc_count = mimic_procedures_from_subject_id(mimic_diag_subject)
        print("Executing cms procedures using desynpuf_ids(diagnosis)...\n")
        cms_proc_all, cms_proc_count = cms_procedures_from_desynpuf_id(cms_diag_desynpuf)
        # append frames to get the final count per procedure code
        proc_count = mimic_proc_count.append(cms_proc_count, ignore_index = True)
        proc_count = proc_count.groupby(['icd9_procedure', 'icd9_proc_desc']).sum().reset_index()
        proc_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
        # DRUGS
        print("Executing mimic drugs using subject_ids(diagnosis)...\n")
        mimic_drug_all, mimic_drug_count = mimic_drugs_from_subject_id(mimic_diag_subject)
        print("Executing cms drugs using desynpuf_ids(diagnosis)...\n")
        cms_drug_all, cms_drug_count = cms_drugs_from_desynpuf_id(cms_diag_desynpuf)
        # append frames to get the final count per procedure code
        drug_count = mimic_drug_count.append(cms_drug_count, ignore_index = True)
        drug_count = drug_count.groupby(['ndc', 'drug']).sum().reset_index()
        drug_count.columns = ['ndc', 'drug', 'NoOfPatients']
        # LABS
        print("Executing mimic labs using subject_ids(diagnosis)...\n")
        mimic_labs_all, mimic_labs_count = mimic_labs_from_subject_id(mimic_diag_subject)
        labs_count = mimic_labs_count
        # combined dataframe consisting of all
        print("Creating the profile & combining dataframes...")
        profile_df, combined_df = create_profile(uniq_count, diag_count, mimic_diag_all, cms_diag_all, proc_count, mimic_proc_all, cms_proc_all, drug_count, mimic_drug_all, cms_drug_all, labs_count, mimic_labs_all)
        # update profile to create combination summary and write combinations to table
        print("Creating combination summary and writing combinations to table...")
        profile_df = update_patient_combinations(combined_df, profile_df)
        # passing arguments for table to be written to excel
        print("Writing to spreadsheet & table...")
        out_to_excel(name, profile_df, mimic_diag_all, cms_diag_all, mimic_proc_all, cms_proc_all, mimic_drug_all, cms_drug_all, mimic_labs_all, combined_df)
                
    #### start point: PROCEDURES ####
    elif driver == 'procedure':
        # PROCEDURES
        print("Executing mimic procedures using codes...\n")
        mimic_proc_all, mimic_proc_count, mimic_proc_subject = mimic_procedures_from_codes(contents)
        print("Executing cms procedures using codes...\n")
        cms_proc_all, cms_proc_count, cms_proc_desynpuf = cms_procedures_from_codes(contents)
        # append frames to get the final count per procedure code
        proc_count = mimic_proc_count.append(cms_proc_count, ignore_index = True)
        proc_count = proc_count.groupby(['icd9_procedure', 'icd9_proc_desc']).sum().reset_index()
        diag_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
        # number of unique patients
        uniq_pat = len(mimic_proc_subject) + len(cms_proc_desynpuf)
        uniq_count = pd.DataFrame({'NoOfUniqPatients':uniq_pat}, index=[0])
        # DIAGNOSIS
        print("Executing mimic diagnosis using subject_ids(procedures)...\n")
        mimic_diag_all, mimic_diag_count = mimic_diagnosis_from_subject_id(mimic_proc_subject)
        print("Executing cms diagnosis using desynpuf_ids(procedures)...\n")
        cms_diag_all, cms_diag_count = cms_diagnosis_from_desynpuf_id(cms_proc_desynpuf)
        # append frames to get the final count per diagnosis code
        diag_count = mimic_diag_count.append(cms_diag_count, ignore_index = True)
        diag_count = diag_count.groupby(['icd9_code', 'icd9_desc']).sum().reset_index()
        diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
        # DRUGS
        print("Executing mimic drugs using subject_ids(procedures)...\n")
        mimic_drug_all, mimic_drug_count = mimic_drugs_from_subject_id(mimic_proc_subject)
        print("Executing cms drugs using desynpuf_ids(procedures)...\n")
        cms_drug_all, cms_drug_count = cms_drugs_from_desynpuf_id(cms_proc_desynpuf)
        # append frames to get the final count per procedure code
        drug_count = mimic_drug_count.append(cms_drug_count, ignore_index = True)
        drug_count = drug_count.groupby(['ndc', 'drug']).sum().reset_index()
        drug_count.columns = ['ndc', 'drug', 'NoOfPatients']
        # LABS
        print("Executing mimic labs using subject_ids(procedures)...\n")
        mimic_labs_all, mimic_labs_count = mimic_labs_from_subject_id(mimic_proc_subject)
        labs_count = mimic_labs_count
        # combined dataframe consisting of all
        print("Creating the profile & combining dataframes...")
        profile_df, combined_df = create_profile(uniq_count, diag_count, mimic_diag_all, cms_diag_all, proc_count, mimic_proc_all, cms_proc_all, drug_count, mimic_drug_all, cms_drug_all, labs_count, mimic_labs_all)
        # update profile to create combination summary and write combinations to table
        print("Creating combination summary and writing combinations to table...")
        profile_df = update_patient_combinations(combined_df, profile_df)
        # passing arguments for table to be written to excel
        print("Writing to spreadsheet & table...")
        out_to_excel(name, profile_df, mimic_diag_all, cms_diag_all, mimic_proc_all, cms_proc_all, mimic_drug_all, cms_drug_all, mimic_labs_all, combined_df)
        
    #### start point: DRUGS ####
    elif driver == 'drugs':
        # DRUGS
        print("Executing mimic drugs using codes...\n")
        mimic_drug_all, mimic_drug_count, mimic_drug_subject = mimic_drugs_from_codes(contents)
        print("Executing cms drugs using codes...\n")
        cms_drug_all, cms_drug_count, cms_drug_desynpuf = cms_drugs_from_codes(contents)
        # append frames to get the final count per procedure code
        drug_count = mimic_drug_count.append(cms_drug_count, ignore_index = True)
        drug_count = drug_count.groupby(['ndc', 'drug']).sum().reset_index()
        drug_count.columns = ['ndc', 'drug', 'NoOfPatients']
        # number of unique patients
        uniq_pat = len(mimic_drug_subject) + len(cms_drug_desynpuf)
        uniq_count = pd.DataFrame({'NoOfUniqPatients':uniq_pat}, index=[0])
        # PROCEDURES
        print("Executing mimic procedures using subject_ids(drugs)...\n")
        mimic_proc_all, mimic_proc_count = mimic_procedures_from_subject_id(mimic_drug_subject)
        print("Executing cms procedures desynpuf_ids(drugs)...\n")
        cms_proc_all, cms_proc_count = cms_procedures_from_desynpuf_id(cms_drug_desynpuf)
        # append frames to get the final count per procedure code
        proc_count = mimic_proc_count.append(cms_proc_count, ignore_index = True)
        proc_count = proc_count.groupby(['icd9_procedure', 'icd9_proc_desc']).sum().reset_index()
        diag_count.columns = ['icd9_procedure', 'icd9_proc_desc', 'NoOfPatients']
        # DIAGNOSIS
        print("Executing mimic diagnosis using subject_ids(drugs)...\n")
        mimic_diag_all, mimic_diag_count = mimic_diagnosis_from_subject_id(mimic_drug_subject)
        print("Executing cms diagnosis using desynpuf_ids(drugs)...\n")
        cms_diag_all, cms_diag_count = cms_procedures_from_desynpuf_id(cms_drug_desynpuf)
        # append frames to get the final count per diagnosis code
        diag_count = mimic_diag_count.append(cms_diag_count, ignore_index = True)
        diag_count = diag_count.groupby(['icd9_code', 'icd9_desc']).sum().reset_index()
        diag_count.columns = ['icd9_code', 'icd9_desc', 'NoOfPatients']
        # LABS
        print("Executing mimic labs using subject_ids(drugs)...\n")
        mimic_labs_all, mimic_labs_count = mimic_labs_from_subject_id(mimic_drug_subject)
        labs_count = mimic_labs_count
        # combined dataframe consisting of all
        print("Creating the profile & combining dataframes...")
        profile_df, combined_df = create_profile(uniq_count, diag_count, mimic_diag_all, cms_diag_all, proc_count, mimic_proc_all, cms_proc_all, drug_count, mimic_drug_all, cms_drug_all, labs_count, mimic_labs_all)
        # update profile to create combination summary and write combinations to table
        print("Creating combination summary and writing combinations to table...")
        profile_df = update_patient_combinations(combined_df, profile_df)
        # passing arguments for table to be written to excel
        print("Writing to spreadsheet & table...")
        out_to_excel(name, profile_df, mimic_diag_all, cms_diag_all, mimic_proc_all, cms_proc_all, mimic_drug_all, cms_drug_all, mimic_labs_all, combined_df)
