# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 19:24:21 2017

@author: 310223340
"""
import pandas as pd
import numpy as np
from data_profile import execute_sql

#%% create master diagnosis frame
def master_icd9():
    print("\nCreating MASTER_icd9...")
    ALZHDMTA_icd9 = pd.DataFrame(['3310','3311','33111','33119','3312','3317','2900','2901','29010','29011','29012','29013','29020','29021','2903','29040','29041','29042','29043','2940','2941','29410','29411','2948','797'],columns=['ICD9'])
    ALZHDMTA_icd9['Diagnosis'] = 'ALZHDMTA'
    
    DIABETES_icd9 = pd.DataFrame(['24900', '24901', '24910', '24911', '24920', '24921', '24930', '24931', '24940', '24941', '24950', '24951', '24960', '24961', '24970', '24971', '24980', '24981', '24990', '24991', '25000', '25001', '25002', '25003', '25010', '25011', '25012', '25013', '25020', '25021', '25022', '25023', '25030', '25031', '25032', '25033', '25040', '25041', '25042', '25043', '25050', '25051', '25052', '25053', '25060', '25061', '25062', '25063', '25070', '25071', '25072', '25073', '25080', '25081', '25082', '25083', '25090', '25091', '25092', '25093', '3572', '36201', '36202', '36641'],columns=['ICD9'])
    DIABETES_icd9['Diagnosis'] = 'DIABETES'
    
    CHF_icd9 = pd.DataFrame(['39891','40201','40211','40291','40401','40411','40491','40403','40413','40493','4280','4281','42820','42821','42822','42823','42830','42831','42832','42833','42840','42841','42842','42843','4289'],columns=['ICD9'])
    CHF_icd9['Diagnosis'] = 'CHF'
    
    CHRNKIDN_icd9 = pd.DataFrame(['01600','01601','01602','01603','01604','01605','01606','0954','1890','1899','2230','23691','24940','24941','25040','25041','25042','25043','2714','2741','27410','28311','40301','40311','40391','40402','40403','40412','40413','40492','40493','4401','4421','5724','5800','5804','58081','58089','5809','5810','5811','5812','5813','58181','58189','5819','5820','5821','5822','5824','58281','58289','5829','5830','5831','5832','5834','5836','5837','58381','58389','5839','5845','5846','5847','5848','5849','585','5851','5852','5853','5854','5855','5856','5859','586','587','5880','5881','58881','58889','5889','591','75312','75313','75314','75315','75316','75317','75319','75320','75321','75322','75323','75329','7944'],columns=['ICD9'])
    CHRNKIDN_icd9['Diagnosis'] = 'CHRNKIDN'
    
    COPD_icd9 = pd.DataFrame(['4910', '4911', '49120', '49121', '49122', '4918', '4919', '4920', '4928', '4940', '4941', '496'],columns=['ICD9'])
    COPD_icd9['Diagnosis'] = 'COPD'
    
    DEPRESSN_icd9 = pd.DataFrame(['29620','29621','29622','29623','29624','29625','29626','29630','29631','29632','29633','29634','29635','29636','29650','29651','29652','29653','29654','29655','29656','29660','29661','29662','29663','29664','29665','29666','29689','2980','3004','3091','311'],columns=['ICD9'])
    DEPRESSN_icd9['Diagnosis'] = 'DEPRESSN'
    
    ISCHMCHT_icd9 = pd.DataFrame(['41000','41001','41002','41010','41011','41012','41020','41021','41022','41030','41031','41032','41040','41041','41042','41050','41051','41052','41060','41061','41062','41070','41071','41072','41080','41081','41082','41090','41091','41092','4110','4111','41181','41189','412','4130','4131','4139','41400','41401','41402','41403','41404','41405','41406','41407','41410','41411','41412','41419','4142','4143','4148','4149','Proc','0066','3601','3602','3603','3604','3605','3606','3607','3609','3610','3611','3612','3613','3614','3615','3616','3617','3619','362','3631','3632','HCPCS','33510','33511','33512','33513','33514','33515','33516','33517','33518','33519','33521','33522','33523','33533','33534','33535','33536','33542','33545','33548','92975','92977','92980','92982','92995','33140','33141'],columns=['ICD9'])
    ISCHMCHT_icd9['Diagnosis'] = 'ISCHMCHT'
    
    OSTEOPRS_icd9 = pd.DataFrame(['73300','73301','73302','73303','73309'],columns=['ICD9'])
    OSTEOPRS_icd9['Diagnosis'] = 'OSTEOPRS'
    
    RA_OA_icd9 = pd.DataFrame(['7140','7141','7142','71430','71431','71432','71433','71500','71504','71509','71510','71511','71512','71513','71514','71515','71516','71517','71518','71520','71521','71522','71523','71524','71525','71526','71527','71528','71530','71531','71532','71533','71534','71535','71536','71537','71538','71580','71589','71590','71598'],columns=['ICD9'])
    RA_OA_icd9['Diagnosis'] = 'RA_OA'
    
    STRKETIA_icd9 = pd.DataFrame(['430','431','43400','43401','43410','43411','43490','43491','4350','4351','4353','4358','4359','436','99702'],columns=['ICD9'])
    STRKETIA_icd9['Diagnosis'] = 'STRKETIA'
    
    CNCR_BRST_icd9 = pd.DataFrame(['1740','1741','1742','1743','1744','1745','1746','1748','1749','2330'],columns=['ICD9'])
    CNCR_BRST_icd9['Diagnosis'] = 'CNCR_BRST'
    
    CNCR_COLR_icd9 = pd.DataFrame(['1540','1541','1530','1531','1532','1533','1534','1535','1536','1537','1538','1539','2303','2304'],columns=['ICD9'])
    CNCR_COLR_icd9['Diagnosis'] = 'CNCR_COLR'
    
    CNCR_PROS_icd9 = pd.DataFrame(['185','2334'],columns=['ICD9'])
    CNCR_PROS_icd9['Diagnosis'] = 'CNCR_PROS'
    
    CNCR_LUNG_icd9 = pd.DataFrame(['1620','1622','1623','1624','1625','1628','1629','2312'],columns=['ICD9'])
    CNCR_LUNG_icd9['Diagnosis'] = 'CNCR_LUNG'
    
    # combine all dataframes into one
    MASTER_icd9 = pd.concat([ALZHDMTA_icd9,CHF_icd9,CHRNKIDN_icd9,COPD_icd9,DEPRESSN_icd9,DIABETES_icd9,ISCHMCHT_icd9,OSTEOPRS_icd9,RA_OA_icd9,STRKETIA_icd9,CNCR_BRST_icd9,CNCR_COLR_icd9,CNCR_PROS_icd9,CNCR_LUNG_icd9], ignore_index=True)
    
    MASTER_icd9 = MASTER_icd9[['Diagnosis', 'ICD9']]
    return (MASTER_icd9)

#%% Function to create diagnosis summaries
def diagnosis_summary(diag_data, comorbidity, MASTER_icd9):
    #The diagnosis related dataset    
    print("\n\tDiagnosis summary:",comorbidity)
    codes = list(MASTER_icd9.loc[MASTER_icd9.Diagnosis==comorbidity,'ICD9'])
    inpat_diag_only = diag_data.loc[diag_data.icd9_code.isin(codes),]
    inpat_diag_only = pd.merge(inpat_diag_only, admission[['hadm_id','visit_year','admittime','dischtime']])
    inpat_diag_only = inpat_diag_only.loc[inpat_diag_only.visit_year.isin([2008,2009,2010]),]
    inpat_diag_only = inpat_diag_only.sort_values(['subject_id','admittime'])
    
    #Date of first diagnosis#
    inpat_first_date = inpat_diag_only.drop_duplicates(subset='subject_id',keep='first')
    inpat_first_date = inpat_first_date[['subject_id', 'admittime']]
    inpat_first_date.columns = ['subject_id', comorbidity+'_First_Diag']
  
    #Diagnosis in 2008#
    inpat_diag_08 = inpat_diag_only.loc[inpat_diag_only.visit_year==2008,]
    inpat_diag_visit_08 = inpat_diag_08[['subject_id','hadm_id']].drop_duplicates()
    inpat_diag_visit_08 = inpat_diag_visit_08.groupby('subject_id').agg({'hadm_id':'count'}).reset_index()
    inpat_diag_visit_08.columns = ['subject_id', comorbidity+'_Visit_2008']
    inpat_diag_cnt_08 = inpat_diag_08[['subject_id','icd9_code']]
    inpat_diag_cnt_08 = inpat_diag_cnt_08.groupby('subject_id').agg({'icd9_code':'count'}).reset_index()
    inpat_diag_cnt_08.columns = ['subject_id', comorbidity+'_Count_2008']
    inpat_diag_08 = pd.merge(inpat_diag_visit_08,inpat_diag_cnt_08)
    
    #Diagnosis in 2009#
    inpat_diag_09 = inpat_diag_only.loc[inpat_diag_only.visit_year==2009,]
    inpat_diag_visit_09 = inpat_diag_09[['subject_id','hadm_id']].drop_duplicates()
    inpat_diag_visit_09 = inpat_diag_visit_09.groupby('subject_id').agg({'hadm_id':'count'}).reset_index()
    inpat_diag_visit_09.columns = ['subject_id', comorbidity+'_Visit_2009']
    inpat_diag_cnt_09 = inpat_diag_09[['subject_id','icd9_code']]
    inpat_diag_cnt_09 = inpat_diag_cnt_09.groupby('subject_id').agg({'icd9_code':'count'}).reset_index()
    inpat_diag_cnt_09.columns = ['subject_id', comorbidity+'_Count_2009']
    inpat_diag_09 = pd.merge(inpat_diag_visit_09, inpat_diag_cnt_09)
    
    #Diagnosis in 2010#
    inpat_diag_10 = inpat_diag_only.loc[inpat_diag_only.visit_year==2010,]
    inpat_diag_visit_10 = inpat_diag_10[['subject_id','hadm_id']].drop_duplicates()
    inpat_diag_visit_10 = inpat_diag_visit_10.groupby('subject_id').agg({'hadm_id':'count'}).reset_index()
    inpat_diag_visit_10.columns = ['subject_id', comorbidity+'_Visit_2010']
    inpat_diag_cnt_10 = inpat_diag_10[['subject_id','icd9_code']]
    inpat_diag_cnt_10 = inpat_diag_cnt_10.groupby('subject_id').agg({'icd9_code':'count'}).reset_index()
    inpat_diag_cnt_10.columns = ['subject_id', comorbidity+'_Count_2010']
    inpat_diag_10 = pd.merge(inpat_diag_visit_10,inpat_diag_cnt_10)
  
    #Condition Era within 30,60,90,180 days ##
    inpat_cera = inpat_diag_only[['subject_id','admittime','dischtime','visit_year']]
    def cal_days_btwn_cond(group):
                delta = abs(group.admittime - group.dischtime.shift(1))
                return(pd.DataFrame(delta.dt.days))
    inpat_cera['days_btwn_cond'] = inpat_cera.groupby(['subject_id','visit_year']).apply(cal_days_btwn_cond)
    #dataframe to be used further for calculating mtbe
    inpat_mtbe = inpat_cera
    #condition era with persistence window of 30 days
    inpat_cera['cera_30'] = inpat_cera['days_btwn_cond'].apply(lambda x: 0 if x <= 30 else 1)
    inpat_cera['cera_60'] = inpat_cera['days_btwn_cond'].apply(lambda x: 0 if x <= 60 else 1)
    inpat_cera['cera_90'] = inpat_cera['days_btwn_cond'].apply(lambda x: 0 if x <= 90 else 1)
    inpat_cera['cera_180'] = inpat_cera['days_btwn_cond'].apply(lambda x: 0 if x <= 180 else 1)
    #condition era based on years
    inpat_cera = inpat_cera.groupby(['subject_id','visit_year']).agg({'cera_30':'sum', 'cera_60':'sum','cera_90':'sum', 'cera_180':'sum'}).reset_index()
    inpat_cera.visit_year = comorbidity + '_' + inpat_cera.visit_year.astype(str)
    inpat_cera = pd.pivot_table(inpat_cera, values=['cera_30','cera_60','cera_90','cera_180'], index='subject_id', columns='visit_year',aggfunc=np.sum,fill_value=0).reset_index()
    colnames = ["_".join((j,i)) for i,j in inpat_cera.columns[1:]]
    inpat_cera.columns = ['subject_id'] + colnames
    
    #MTBE
    inpat_mtbe = inpat_mtbe.groupby('subject_id').agg({'days_btwn_cond':'mean'}).reset_index()
    inpat_mtbe.columns = ['subject_id', comorbidity+'_MTBE']
    
    #Merging all computed datasets
    inpat_diag_summary = inpat_first_date.merge(inpat_diag_08,on='subject_id',how='outer').merge(inpat_diag_09,on='subject_id',how='outer').merge(inpat_diag_10,on='subject_id',how='outer').merge(inpat_cera,on='subject_id',how='outer').merge(inpat_mtbe,on='subject_id',how='outer')
    return(inpat_diag_summary)
    
#%% Function to create drug summaries
def drug_summary(pres_data, pres_name, adm_data):
    print("\nDrug summary:",pres_name)
    drug_era = pres_data[['subject_id', 'hadm_id', 'starttime', 'endtime', 'ndc', 'drug', 'dose_val_rx', 'dose_unit_rx', 'drug_type']]
    drug_era.columns = ['subject_id','hadm_id','drug_strt_orig','drug_end_orig','ndc','drug','dose_val_rx', 'dose_unit_rx', 'drug_type']
    drug_era = drug_era.dropna(subset=['drug_strt_orig','drug_end_orig'])
    drug_era.drug_strt_orig = pd.to_datetime(drug_era.drug_strt_orig, format='%Y-%m-%d')
    drug_era.drug_end_orig = pd.to_datetime(drug_era.drug_end_orig, format='%Y-%m-%d')
    drug_era = drug_era.sort_values(['drug_strt_orig', 'drug_end_orig']).reset_index(drop=True)
    drug_era = pd.merge(drug_era, adm_data, on=['subject_id','hadm_id'])
    # check for reversed start & end date, and make changes if required
    for i, row in drug_era.iterrows():
        if row.drug_strt_orig > row.drug_end_orig:
            drug_era.loc[i,'drug_strt_trns'] = row.drug_end_orig
            drug_era.loc[i,'drug_end_trns'] = row.drug_strt_orig
        else:
            drug_era.loc[i,'drug_strt_trns'] = row.drug_strt_orig
            drug_era.loc[i,'drug_end_trns'] = row.drug_end_orig
    drug_era['drug_days_supply'] = (drug_era.drug_end_trns - drug_era.drug_strt_trns).dt.days + 1
    drug_era['days_drug_strt_orig_admittime_orig'] = (drug_era.drug_strt_trns - drug_era.admittime_orig).dt.days
    drug_era['shftd_drug_days_supply'] = drug_era.drug_days_supply * drug_era.orig_adm_compression.astype(float)
    drug_era['shftd_days_drug_strt_orig_admittime_orig'] = drug_era.days_drug_strt_orig_admittime_orig * drug_era.orig_adm_compression.astype(float)
    drug_era['drug_strt_shftd'] = drug_era.admittime + drug_era.shftd_days_drug_strt_orig_admittime_orig.apply(lambda x: pd.Timedelta(days=round(x,0)))
    drug_era['drug_end_shftd'] = drug_era.drug_strt_shftd + drug_era.shftd_drug_days_supply.apply(lambda x: pd.Timedelta(days=round(x,0)))
    # count of prescription
    inpat_drug_count = drug_era.groupby(['subject_id','visit_year']).agg({'drug':'count'}).reset_index()
    inpat_drug_count.visit_year = pres_name + '_' + inpat_drug_count.visit_year.astype(str)
    inpat_drug_count = pd.pivot_table(inpat_drug_count, values='drug', index='subject_id', columns='visit_year',aggfunc=np.sum,fill_value=0).reset_index()
    # era calculation
    inpat_dera = drug_era[['subject_id','drug_strt_shftd','drug_end_shftd','visit_year']]
    inpat_dera = inpat_dera.sort_values(['subject_id','drug_strt_shftd','drug_end_shftd']).reset_index( drop=True)
    #calculate days between prescriptions
    def cal_days_btwn_pres(group):
        delta = group.drug_strt_shftd - group.drug_end_shftd.shift(1)
        return(pd.DataFrame(delta.dt.days))
    inpat_dera['days_btwn_pres'] = inpat_dera.groupby(['subject_id','visit_year']).apply(cal_days_btwn_pres)
    #drug era with perseistence window of 30 days
    inpat_dera['dera_30'] = inpat_dera.days_btwn_pres.apply(lambda x: 0 if x <= 30 else 1)
    inpat_dera['dera_60'] = inpat_dera.days_btwn_pres.apply(lambda x: 0 if x <= 60 else 1)
    inpat_dera['dera_90'] = inpat_dera.days_btwn_pres.apply(lambda x: 0 if x <= 90 else 1)
    #drug era based on years
    inpat_dera = inpat_dera.groupby(['subject_id','visit_year']).agg({'dera_30':'sum', 'dera_60':'sum','dera_90':'sum'}).reset_index()
    inpat_dera.visit_year = pres_name + '_' + inpat_dera.visit_year.astype(str)
    inpat_dera = pd.pivot_table(inpat_dera, values=['dera_30','dera_60','dera_90'], index='subject_id', columns='visit_year',aggfunc=np.sum,fill_value=0).reset_index()
    #flatten multiindex column names
    colnames = ["_".join((j,i)) for i,j in inpat_dera.columns[1:]]
    inpat_dera.columns = ['subject_id'] + colnames
    inpat_dera = pd.merge(inpat_drug_count, inpat_dera, on='subject_id')
    return(inpat_dera)

#%% feature creation from raw data
def feature_creation(admission, diagnosis, procedures, MASTER_icd9):
    # Total number of comorbidities 
    print("\nFeature Creation: Total Comorbidities")
    inpat_totcom = diagnosis.loc[diagnosis.hadm_id.isin(admission.hadm_id),['subject_id','hadm_id','seq_num']]
    inpat_totcom = inpat_totcom.groupby(['subject_id','hadm_id']).agg({'seq_num':'max'}).reset_index()
    inpat_totcom = inpat_totcom.groupby('subject_id').agg({'seq_num':'sum'}).reset_index().rename(columns={'seq_num':'TotalComorb'})
    
    # For first visit ##
    firstadm = admission.sort_values(['subject_id','dischtime']).reset_index(drop = True)
    firstadm = firstadm.drop_duplicates(subset='subject_id',keep='first')
    
    # Number of comorbidities in 1st visit
    print("\nFeature Creation: Total Comorbidities in 1st visit")
    inpat_com = diagnosis.loc[diagnosis.hadm_id.isin(firstadm.hadm_id),]
    inpat_com = inpat_com.groupby('subject_id').agg({'seq_num':'max'}).reset_index().rename(columns={'seq_num':'ComorFirstVisit'})
    
    # Total number of comorbidities in each year
    print("\nFeature Creation: Total Comorbidities in each year")
    inpat_totcom_yr = pd.merge(admission[['hadm_id','visit_year']],diagnosis,on='hadm_id')
    inpat_totcom_yr = inpat_totcom_yr.groupby(['subject_id','hadm_id','visit_year']).agg({'seq_num':'max'}).reset_index()
    inpat_totcom_yr = inpat_totcom_yr.groupby(['subject_id','visit_year']).agg({'seq_num':'sum'}).reset_index().rename(columns={'seq_num':'TotalComorb'})
    inpat_totcom_yr = pd.pivot_table(inpat_totcom_yr, values='TotalComorb', index='subject_id', columns='visit_year',aggfunc=np.sum,fill_value=0).reset_index()
    colnames = list('comorb_'+ inpat_totcom_yr.columns[1:].astype(str))
    inpat_totcom_yr.columns = inpat_totcom_yr.columns[:1].tolist() + colnames
    
    # Number of procedures in 1st visit
    print("\nFeature Creation: Total Procedures in 1st visit")
    inpat_proc = procedures.loc[procedures.hadm_id.isin(firstadm.hadm_id),]
    inpat_proc = inpat_proc.groupby('subject_id').agg({'seq_num': 'max'}).reset_index().rename(columns={'seq_num': 'ProcsFirstVisit'})
    
    # Total number of procedures 
    print("\nFeature Creation: Total Procedures")
    inpat_totproc = procedures.loc[procedures.hadm_id.isin(admission.hadm_id),['subject_id','hadm_id','seq_num']]
    inpat_totproc = inpat_totproc.groupby(['subject_id','hadm_id']).agg({'seq_num':'max'}).reset_index()
    inpat_totproc = inpat_totproc.groupby('subject_id').agg({'seq_num':'sum'}).reset_index().rename(columns={'seq_num':'TotalProcs'})
    
    # Total number of admissions in each year
    print("\nFeature Creation: Total admissions in each year")
    inpat_visit = admission.groupby(['subject_id','visit_year']).agg({'hadm_id':'count'}).reset_index()
    inpat_visit = pd.pivot_table(inpat_visit, values='hadm_id', index='subject_id', columns='visit_year',aggfunc=np.sum,fill_value=0).reset_index()
    colnames = list('TotalVisits_'+ inpat_visit.columns[1:].astype(str))
    inpat_visit.columns = inpat_visit.columns[:1].tolist() + colnames

    ##compute the diagnosis summaries
    print("\nFeature Creation: Diagosis summary per diagnosis")
    inpat_diagnosis = pd.DataFrame(columns=['subject_id'])
    for com in MASTER_icd9.Diagnosis.unique():
        inpat_diagnosis_temp = diagnosis_summary(diagnosis,com,MASTER_icd9)
        inpat_diagnosis = inpat_diagnosis.merge(inpat_diagnosis_temp,on='subject_id',how='outer')
        
    ## Length of stay
    print("\nFeature Creation: Length of stay in each year")
    admission['length_of_stay'] = (admission.dischtime - admission.admittime).dt.days + 1
    inpat_los = admission.groupby(['subject_id','visit_year']).agg({'length_of_stay':'sum'}).reset_index()
    inpat_los = pd.pivot_table(inpat_los, values='length_of_stay', index='subject_id', columns='visit_year', aggfunc=np.sum,fill_value=0).reset_index()
    colnames = list('los_'+ inpat_los.columns[1:].astype(str))
    inpat_los.columns = inpat_los.columns[:1].tolist() + colnames
    
    ## If there's a new ICD9 code in the last visit ##
    print("\nFeature Creation: New diagnosis in last visit")
    #keep patients with more than 1 visit
    inpat_last = admission[['subject_id','hadm_id','admittime']].groupby('subject_id').filter(lambda x: len(x) > 1)
    ## For last visit ##
    inpat_last = inpat_last.sort_values('admittime',ascending=False).reset_index(drop=True)
    #last admission
    inpat_last_adm = inpat_last.drop_duplicates(subset='subject_id',keep='first')
    #all but last admission
    inpat_all_but_last_adm = inpat_last.loc[~inpat_last.hadm_id.isin(inpat_last_adm.hadm_id),]
    #last admission diagnosis
    inpat_last_diag = diagnosis.loc[diagnosis.hadm_id.isin(inpat_last_adm.hadm_id),['subject_id','icd9_code']]
    #all but last admission diagnosis
    inpat_all_but_last_diag = diagnosis.loc[diagnosis.hadm_id.isin(inpat_all_but_last_adm.hadm_id),['subject_id','icd9_code']]
    #define function to find new diagnosis in last visit
    def new_diag_last_adm(subject_id):
        last_diag = inpat_last_diag.loc[inpat_last_diag.subject_id==subject_id,]
        all_but_last_diag = inpat_all_but_last_diag.loc[inpat_all_but_last_diag.subject_id==subject_id,]
        last_new = pd.concat([last_diag,all_but_last_diag,all_but_last_diag]).drop_duplicates(keep=False)
        return(last_new)
    #call function for each patient
    inpat_last_new = pd.DataFrame()    
    for subj in inpat_last_diag.subject_id.unique():
        inpat_last_new_temp = new_diag_last_adm(subj)
        inpat_last_new = pd.concat([inpat_last_new, inpat_last_new_temp])
    #match ICD9 codes with the names, for unatched fill as OTHERS
    inpat_last_new = pd.merge(inpat_last_new, MASTER_icd9, left_on='icd9_code', right_on='ICD9', how='left').fillna('OTHERS')
    inpat_last_new = inpat_last_new.drop('ICD9',1)
    #rename new diagnosis in last visit columns
    inpat_last_new.Diagnosis = 'NewICD9LastVisit_' + inpat_last_new.Diagnosis
    #pivot data in the required format
    inpat_last_new = pd.pivot_table(inpat_last_new, values='icd9_code', index='subject_id', columns='Diagnosis', aggfunc=np.size,fill_value=0).reset_index()

    ### merging all the features created
    inpat_final = inpat_totcom.merge(inpat_totcom, on='subject_id', how='outer').merge(inpat_totcom_yr, on='subject_id', how='outer').merge(inpat_proc, on='subject_id', how='outer').merge(inpat_totproc, on='subject_id', how='outer').merge(inpat_visit, on='subject_id', how='outer').merge(inpat_diagnosis, on='subject_id', how='outer').merge(inpat_los,on='subject_id',how='outer').merge(inpat_last_new,on='subject_id',how='outer')
    
    ##add inpatient flag
    inpat_final['io_flag'] = 'in'

    # Procedures processing
    print("\nFeature Creation: Procedure related features")
    #merge with admissions for relevant hadm_id
    proc = pd.merge(procedures, admission[['hadm_id','admittime','dischtime','visit_year']])
    proc['visit_year'] = proc.visit_year.astype(str)
    #consider procedures carried out in 2008, 2009 & 2010 only
    proc = proc.loc[proc.visit_year.isin(['2008','2009','2010']),]
    #add the reporting summary & reporting label details
    proc_report_details = execute_sql('select * from mimic_d_icd9_proc_reference')
    proc = pd.merge(proc, proc_report_details[['icd9_proc', 'icd9_proc_report_summ', 'icd9_proc_report_label']], left_on='icd9_procedure', right_on='icd9_proc')
    #pivoting dataset in the required format
    proc_final = pd.pivot_table(proc, values='icd9_procedure', index='subject_id', columns=['visit_year', 'icd9_proc_report_summ', 'icd9_proc_report_label'], aggfunc=np.size,fill_value=0)
    #flatten multiindex column names
    proc_final.columns = ['_'.join(col).strip() for col in proc_final.columns]
    proc_final = proc_final.reset_index()
    #### meging with inpatient set above
    inpat_final = pd.merge(inpat_final,proc_final, on='subject_id', how='left')
    
    # Drugs processing (added on 01/08/2018)
    pres = execute_sql("select * from mimic_prescriptions where drug_type='MAIN'")
    pres_unique = pres[['drug']].drop_duplicates()
    inpat_pres = pd.DataFrame(columns=['subject_id'])
    for i, row in pres_unique.iterrows():
        pres_data = pres.loc[pres.drug == row.drug,]
        adm_data = admission.loc[admission.hadm_id.isin(pres_data.hadm_id),]
        inpat_pres_temp = drug_summary(pres_data, row.drug, adm_data)
        inpat_pres = inpat_pres.merge(inpat_pres_temp,on='subject_id',how='outer')
    #### meging with inpatient set above
    inpat_final = pd.merge(inpat_final, inpat_pres, on='subject_id', how='left')
    return(inpat_final)

#%% derived features based on features summary
def derived_summary_features(inpat_final):
    print("\nFeature Creation: Summarised features")
    inpat_final['NoOfCom_14_2008'] = inpat_final[["ALZHDMTA_Count_2008","CHF_Count_2008","CHRNKIDN_Count_2008","COPD_Count_2008","DEPRESSN_Count_2008","DIABETES_Count_2008","ISCHMCHT_Count_2008","OSTEOPRS_Count_2008","RA_OA_Count_2008","STRKETIA_Count_2008","CNCR_BRST_Count_2008","CNCR_COLR_Count_2008","CNCR_PROS_Count_2008","CNCR_LUNG_Count_2008"]].sum(axis=1)
    inpat_final['NoOfCom_14_2009'] = inpat_final[["ALZHDMTA_Count_2009","CHF_Count_2009","CHRNKIDN_Count_2009","COPD_Count_2009","DEPRESSN_Count_2009","DIABETES_Count_2009","ISCHMCHT_Count_2009","OSTEOPRS_Count_2009","RA_OA_Count_2009","STRKETIA_Count_2009","CNCR_BRST_Count_2009","CNCR_COLR_Count_2009","CNCR_PROS_Count_2009","CNCR_LUNG_Count_2009"]].sum(axis=1)
    inpat_final['NoOfCom_14_2010'] = inpat_final[["ALZHDMTA_Count_2010","CHF_Count_2010","CHRNKIDN_Count_2010","COPD_Count_2010","DEPRESSN_Count_2010","DIABETES_Count_2010","ISCHMCHT_Count_2010","OSTEOPRS_Count_2010","RA_OA_Count_2010","STRKETIA_Count_2010","CNCR_BRST_Count_2010","CNCR_COLR_Count_2010","CNCR_PROS_Count_2010","CNCR_LUNG_Count_2010"]].sum(axis=1)
    inpat_final['Total_Comorb_14'] = inpat_final[["NoOfCom_14_2008", "NoOfCom_14_2009", "NoOfCom_14_2010"]].sum(axis=1)
    inpat_final['Total_Comorb_14_0809'] = inpat_final[["NoOfCom_14_2008","NoOfCom_14_2009"]].sum(axis=1)
    
    ### ALZHDMTA ###
    inpat_final['ALZHDMTA_Prc_Dstrb_2008'] = inpat_final.ALZHDMTA_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['ALZHDMTA_Prc_Dstrb_2009'] = inpat_final.ALZHDMTA_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['ALZHDMTA_Prc_Dstrb_2010'] = inpat_final.ALZHDMTA_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['ALZHDMTA_Count_2008_2009'] = inpat_final[["ALZHDMTA_Count_2008","ALZHDMTA_Count_2009"]].sum(axis=1)
    inpat_final['ALZHDMTA_Prc_2008_2009'] = inpat_final.ALZHDMTA_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['ALZHDMTA_Count_2008_2010'] = inpat_final[["ALZHDMTA_Count_2008","ALZHDMTA_Count_2009","ALZHDMTA_Count_2010"]].sum(axis=1)
    inpat_final['ALZHDMTA_Prc_2008_2010'] = inpat_final.ALZHDMTA_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CHF ###
    inpat_final['CHF_Prc_Dstrb_2008'] = inpat_final.CHF_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CHF_Prc_Dstrb_2009'] = inpat_final.CHF_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CHF_Prc_Dstrb_2010'] = inpat_final.CHF_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CHF_Count_2008_2009'] = inpat_final[["CHF_Count_2008","CHF_Count_2009"]].sum(axis=1)
    inpat_final['CHF_Prc_2008_2009'] = inpat_final.CHF_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CHF_Count_2008_2010'] = inpat_final[["CHF_Count_2008","CHF_Count_2009","CHF_Count_2010"]].sum(axis=1)
    inpat_final['CHF_Prc_2008_2010'] = inpat_final.CHF_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CHRNKIDN ###
    inpat_final['CHRNKIDN_Prc_Dstrb_2008'] = inpat_final.CHRNKIDN_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CHRNKIDN_Prc_Dstrb_2009'] = inpat_final.CHRNKIDN_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CHRNKIDN_Prc_Dstrb_2010'] = inpat_final.CHRNKIDN_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CHRNKIDN_Count_2008_2009'] = inpat_final[["CHRNKIDN_Count_2008","CHRNKIDN_Count_2009"]].sum(axis=1)
    inpat_final['CHRNKIDN_Prc_2008_2009'] = inpat_final.CHRNKIDN_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CHRNKIDN_Count_2008_2010'] = inpat_final[["CHRNKIDN_Count_2008","CHRNKIDN_Count_2009","CHRNKIDN_Count_2010"]].sum(axis=1)
    inpat_final['CHRNKIDN_Prc_2008_2010'] = inpat_final.CHRNKIDN_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### COPD ###
    inpat_final['COPD_Prc_Dstrb_2008'] = inpat_final.COPD_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['COPD_Prc_Dstrb_2009'] = inpat_final.COPD_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['COPD_Prc_Dstrb_2010'] = inpat_final.COPD_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['COPD_Count_2008_2009'] = inpat_final[["COPD_Count_2008","COPD_Count_2009"]].sum(axis=1)
    inpat_final['COPD_Prc_2008_2009'] = inpat_final.COPD_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['COPD_Count_2008_2010'] = inpat_final[["COPD_Count_2008","COPD_Count_2009","COPD_Count_2010"]].sum(axis=1)
    inpat_final['COPD_Prc_2008_2010'] = inpat_final.COPD_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### DEPRESSN ###
    inpat_final['DEPRESSN_Prc_Dstrb_2008'] = inpat_final.DEPRESSN_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['DEPRESSN_Prc_Dstrb_2009'] = inpat_final.DEPRESSN_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['DEPRESSN_Prc_Dstrb_2010'] = inpat_final.DEPRESSN_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['DEPRESSN_Count_2008_2009'] = inpat_final[["DEPRESSN_Count_2008","DEPRESSN_Count_2009"]].sum(axis=1)
    inpat_final['DEPRESSN_Prc_2008_2009'] = inpat_final.DEPRESSN_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['DEPRESSN_Count_2008_2010'] = inpat_final[["DEPRESSN_Count_2008","DEPRESSN_Count_2009","DEPRESSN_Count_2010"]].sum(axis=1)
    inpat_final['DEPRESSN_Prc_2008_2010'] = inpat_final.DEPRESSN_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### DIABETES ###
    inpat_final['DIABETES_Prc_Dstrb_2008'] = inpat_final.DIABETES_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['DIABETES_Prc_Dstrb_2009'] = inpat_final.DIABETES_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['DIABETES_Prc_Dstrb_2010'] = inpat_final.DIABETES_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['DIABETES_Count_2008_2009'] = inpat_final[["DIABETES_Count_2008","DIABETES_Count_2009"]].sum(axis=1)
    inpat_final['DIABETES_Prc_2008_2009'] = inpat_final.DIABETES_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['DIABETES_Count_2008_2010'] = inpat_final[["DIABETES_Count_2008","DIABETES_Count_2009","DIABETES_Count_2010"]].sum(axis=1)
    inpat_final['DIABETES_Prc_2008_2010'] = inpat_final.DIABETES_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### ISCHMCHT ###
    inpat_final['ISCHMCHT_Prc_Dstrb_2008'] = inpat_final.ISCHMCHT_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['ISCHMCHT_Prc_Dstrb_2009'] = inpat_final.ISCHMCHT_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['ISCHMCHT_Prc_Dstrb_2010'] = inpat_final.ISCHMCHT_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['ISCHMCHT_Count_2008_2009'] = inpat_final[["ISCHMCHT_Count_2008","ISCHMCHT_Count_2009"]].sum(axis=1)
    inpat_final['ISCHMCHT_Prc_2008_2009'] = inpat_final.ISCHMCHT_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['ISCHMCHT_Count_2008_2010'] = inpat_final[["ISCHMCHT_Count_2008","ISCHMCHT_Count_2009","ISCHMCHT_Count_2010"]].sum(axis=1)
    inpat_final['ISCHMCHT_Prc_2008_2010'] = inpat_final.ISCHMCHT_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### OSTEOPRS ###
    inpat_final['OSTEOPRS_Prc_Dstrb_2008'] = inpat_final.OSTEOPRS_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['OSTEOPRS_Prc_Dstrb_2009'] = inpat_final.OSTEOPRS_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['OSTEOPRS_Prc_Dstrb_2010'] = inpat_final.OSTEOPRS_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['OSTEOPRS_Count_2008_2009'] = inpat_final[["OSTEOPRS_Count_2008","OSTEOPRS_Count_2009"]].sum(axis=1)
    inpat_final['OSTEOPRS_Prc_2008_2009'] = inpat_final.OSTEOPRS_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['OSTEOPRS_Count_2008_2010'] = inpat_final[["OSTEOPRS_Count_2008","OSTEOPRS_Count_2009","OSTEOPRS_Count_2010"]].sum(axis=1)
    inpat_final['OSTEOPRS_Prc_2008_2010'] = inpat_final.OSTEOPRS_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### RA_OA ###
    inpat_final['RA_OA_Prc_Dstrb_2008'] = inpat_final.RA_OA_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['RA_OA_Prc_Dstrb_2009'] = inpat_final.RA_OA_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['RA_OA_Prc_Dstrb_2010'] = inpat_final.RA_OA_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['RA_OA_Count_2008_2009'] = inpat_final[["RA_OA_Count_2008","RA_OA_Count_2009"]].sum(axis=1)
    inpat_final['RA_OA_Prc_2008_2009'] = inpat_final.RA_OA_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['RA_OA_Count_2008_2010'] = inpat_final[["RA_OA_Count_2008","RA_OA_Count_2009","RA_OA_Count_2010"]].sum(axis=1)
    inpat_final['RA_OA_Prc_2008_2010'] = inpat_final.RA_OA_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### STRKETIA ###
    inpat_final['STRKETIA_Prc_Dstrb_2008'] = inpat_final.STRKETIA_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['STRKETIA_Prc_Dstrb_2009'] = inpat_final.STRKETIA_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['STRKETIA_Prc_Dstrb_2010'] = inpat_final.STRKETIA_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['STRKETIA_Count_2008_2009'] = inpat_final[["STRKETIA_Count_2008","STRKETIA_Count_2009"]].sum(axis=1)
    inpat_final['STRKETIA_Prc_2008_2009'] = inpat_final.STRKETIA_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['STRKETIA_Count_2008_2010'] = inpat_final[["STRKETIA_Count_2008","STRKETIA_Count_2009","STRKETIA_Count_2010"]].sum(axis=1)
    inpat_final['STRKETIA_Prc_2008_2010'] = inpat_final.STRKETIA_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CNCR_BRST ###
    inpat_final['CNCR_BRST_Prc_Dstrb_2008'] = inpat_final.CNCR_BRST_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CNCR_BRST_Prc_Dstrb_2009'] = inpat_final.CNCR_BRST_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CNCR_BRST_Prc_Dstrb_2010'] = inpat_final.CNCR_BRST_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CNCR_BRST_Count_2008_2009'] = inpat_final[["CNCR_BRST_Count_2008","CNCR_BRST_Count_2009"]].sum(axis=1)
    inpat_final['CNCR_BRST_Prc_2008_2009'] = inpat_final.CNCR_BRST_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CNCR_BRST_Count_2008_2010'] = inpat_final[["CNCR_BRST_Count_2008","CNCR_BRST_Count_2009","CNCR_BRST_Count_2010"]].sum(axis=1)
    inpat_final['CNCR_BRST_Prc_2008_2010'] = inpat_final.CNCR_BRST_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CNCR_COLR ###
    inpat_final['CNCR_COLR_Prc_Dstrb_2008'] = inpat_final.CNCR_COLR_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CNCR_COLR_Prc_Dstrb_2009'] = inpat_final.CNCR_COLR_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CNCR_COLR_Prc_Dstrb_2010'] = inpat_final.CNCR_COLR_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CNCR_COLR_Count_2008_2009'] = inpat_final[["CNCR_COLR_Count_2008","CNCR_COLR_Count_2009"]].sum(axis=1)
    inpat_final['CNCR_COLR_Prc_2008_2009'] = inpat_final.CNCR_COLR_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CNCR_COLR_Count_2008_2010'] = inpat_final[["CNCR_COLR_Count_2008","CNCR_COLR_Count_2009","CNCR_COLR_Count_2010"]].sum(axis=1)
    inpat_final['CNCR_COLR_Prc_2008_2010'] = inpat_final.CNCR_COLR_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CNCR_PROS ###
    inpat_final['CNCR_PROS_Prc_Dstrb_2008'] = inpat_final.CNCR_PROS_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CNCR_PROS_Prc_Dstrb_2009'] = inpat_final.CNCR_PROS_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CNCR_PROS_Prc_Dstrb_2010'] = inpat_final.CNCR_PROS_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CNCR_PROS_Count_2008_2009'] = inpat_final[["CNCR_PROS_Count_2008","CNCR_PROS_Count_2009"]].sum(axis=1)
    inpat_final['CNCR_PROS_Prc_2008_2009'] = inpat_final.CNCR_PROS_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CNCR_PROS_Count_2008_2010'] = inpat_final[["CNCR_PROS_Count_2008","CNCR_PROS_Count_2009","CNCR_PROS_Count_2010"]].sum(axis=1)
    inpat_final['CNCR_PROS_Prc_2008_2010'] = inpat_final.CNCR_PROS_Count_2008_2010/inpat_final.Total_Comorb_14
    
    ### CNCR_LUNG ###
    inpat_final['CNCR_LUNG_Prc_Dstrb_2008'] = inpat_final.CNCR_LUNG_Count_2008/inpat_final.NoOfCom_14_2008
    inpat_final['CNCR_LUNG_Prc_Dstrb_2009'] = inpat_final.CNCR_LUNG_Count_2009/inpat_final.NoOfCom_14_2009
    inpat_final['CNCR_LUNG_Prc_Dstrb_2010'] = inpat_final.CNCR_LUNG_Count_2010/inpat_final.NoOfCom_14_2010
    
    inpat_final['CNCR_LUNG_Count_2008_2009'] = inpat_final[["CNCR_LUNG_Count_2008","CNCR_LUNG_Count_2009"]].sum(axis=1)
    inpat_final['CNCR_LUNG_Prc_2008_2009'] = inpat_final.CNCR_LUNG_Count_2008_2009/inpat_final.Total_Comorb_14_0809
    
    inpat_final['CNCR_LUNG_Count_2008_2010'] = inpat_final[["CNCR_LUNG_Count_2008","CNCR_LUNG_Count_2009","CNCR_LUNG_Count_2010"]].sum(axis=1)
    inpat_final['CNCR_LUNG_Prc_2008_2010'] = inpat_final.CNCR_LUNG_Count_2008_2010/inpat_final.Total_Comorb_14
    return(inpat_final)

#%%# main function
if __name__ == "__main__":
    print("\nExtracting the admissions, diagnoses & procedures dataset...")
    # admission dataset
    admission = execute_sql("select * from mimic_admissions_shftd")
    # remove Records with NA in dischtime
    admission = admission.dropna(subset=['dischtime'])
    # add visit_year column
    admission['visit_year'] = admission['admittime'].dt.year
    # diagnosis dataset
    diagnosis = execute_sql("select a.*, b.long_title, b.icd9_std from mimic_diagnoses a, mimic_d_icd_diagnoses b where a.icd9_code = b.icd9_code")
    # procedures dataset
    procedures = execute_sql('select a.*, b.long_title from mimic_procedures a, mimic_d_icd_procedures b where a.icd9_procedure = b.icd9_procedure')
    # master_icd9
    MASTER_icd9 = master_icd9()
    #feature creation
    inpat_features = feature_creation(admission, diagnosis, procedures, MASTER_icd9)
    #calculating & adding summarised features
    mimic_summary = derived_summary_features(inpat_features)
    #write output to a spreadsheet
    print("\nGenerating the spreadsheet - mimic_all_patients.xlsx")
    writer = pd.ExcelWriter("mimic_all_patients.xlsx", engine='xlsxwriter')
    mimic_summary.to_excel(writer, sheet_name='mimic_all_patients')
    writer.save()
