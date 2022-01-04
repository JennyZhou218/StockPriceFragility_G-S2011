# -*- coding: utf-8 -*-
"""
Calculate stock price fragility (Greenwood Thesmar 2011)

Input data:
    1. CRSP MSF: 19902020_crsp_msf.dta
    2. CRSP Mutual fund: 19902021_crsp_fund.dta
    3. Thonsom s12type3: s12type3.sas7bdat
    4. MFLink1: MFLink1.sas7bdat
    5. MFLink2: MFLink2.sas7bdat
Output data: fragility.dta at stock-quarter level

Created on Thu Oct 21 15:49:46 2021
@author: JennyZhou
"""

import numpy as np 
import pandas as pd
import os
import gc
gc.collect()
os.chdir("yourdirectory")

#############################################################################
# Step 1                                                                    #
# Get the crsp stock price and shrout from  CRSP MSF                        #
#############################################################################
crsp_m = pd.read_stata("data/19902020_crsp_msf.dta")
crsp_m = crsp_m[(crsp_m['SHRCD']==10)|(crsp_m['SHRCD']==11)]
# get month and quarter-end dates
crsp_m['mdate']=crsp_m['date']+pd.offsets.MonthEnd(0)
crsp_m['qdate']=crsp_m['date']+pd.offsets.QuarterEnd(0)

# calculate adjusted price, total shares and market cap
crsp_m['p']=crsp_m['PRC'].abs()/crsp_m['CFACPR'] # price adjusted
crsp_m['tso']=crsp_m['SHROUT']*crsp_m['CFACSHR']*1e3 # total shares out adjusted
crsp_m['me'] = crsp_m['p']*crsp_m['tso']/1e6 # market cap in $mil

# keep only relevant columns
crsp_m = crsp_m[['PERMNO','mdate','qdate','date','CFACSHR', 'p', 'tso','me','CUSIP','NCUSIP']]

# For each stock(permno), each quarter (qdate), find the last monthly date (mdate)
qend = crsp_m[['PERMNO','mdate','qdate']].groupby(['PERMNO','qdate'])['mdate'].max().reset_index()

# Merge back to keep last monthly observation for each quarter
crsp_qend = pd.merge(crsp_m, qend, how='inner', on=['PERMNO','qdate','mdate'])

del [[crsp_m,qend]]
crsp_qend = crsp_qend.replace([np.inf, -np.inf], np.nan)  
crsp_qend.to_stata("data/intermediate/crsp_qend.dta")
gc.collect()

#############################################################################
# Step 2                                                                    #
# Get the Thomson s12type3 holdings data, merge with crsp stock information #
# Use s12type3                                                              #
#############################################################################
s12 = pd.read_sas("data/s12type3.sas7bdat")
crsp_qend = pd.read_stata("data/intermediate/crsp_qend.dta")

s12['CUSIP1'] = s12['CUSIP'].str.decode('utf-8')
# merge with crsp by CUSIP first, then by NCUSIP
vint = pd.merge(s12,crsp_qend,how='left',left_on=['CUSIP1','FDATE'],right_on=['CUSIP','mdate'])
vint = pd.merge(vint,crsp_qend,how='left',left_on=['CUSIP1','FDATE'],right_on=['NCUSIP','mdate'])
del [[s12,crsp_qend]]

vint['p'] = vint['p_x'].mask(pd.isnull,vint['p_y'])
vint['cfacshr'] = vint['CFACSHR_x'].mask(pd.isnull,vint['CFACSHR_y'])
vint['tso'] = vint['tso_x'].mask(pd.isnull,vint['tso_y'])
vint['me'] = vint['me_x'].mask(pd.isnull,vint['me_y'])
vint['CUSIP'] = vint['NCUSIP_x'].mask(pd.isnull,vint['NCUSIP_y'])
vint['PERMNO'] = vint['PERMNO_x'].mask(pd.isnull,vint['PERMNO_y'])

# select the columns needed for calculation
vint = vint.loc[:,['FDATE','CUSIP','FUNDNO','SHARES','PERMNO','cfacshr','p','tso','me']]
vint.to_stata("data/intermediate/s12_crsp.dta")

#############################################################################
# Step 3                                                                    #
# USE CRSP MF and MFLINK1 TO GET RAW DATA FOR CALCULATE OMEGA               #
#############################################################################
crspmf =  pd.read_stata("data/19902021_crsp_fund.dta")
mflink1 = pd.read_sas("data/s12/data2019/MFLink1.sas7bdat")
mflink1 = mflink1[['crsp_fundno','wficn']]
crspmf = pd.merge(crspmf,mflink1,how='inner',on='crsp_fundno') 
# drop those with missing tna and missing return
crspmf = crspmf.dropna()  

# portfolio level tna & return
# use last tna (total net asset) obeservation as portw (portfolio weight)
crspmf['portw'] = crspmf.sort_values('caldt').groupby(['crsp_fundno'])['mtna'].shift(1)
crspmf_tna = crspmf.groupby(['wficn','caldt'])['mtna'].sum().reset_index().rename(columns={'mtna':'tna'})
crspmf_portw = crspmf.groupby(['wficn','caldt'])['portw'].sum().reset_index().rename(columns={'portw':'sum_portw'})
crspmf = pd.merge(crspmf,crspmf_tna,how = 'left', on=['wficn','caldt'])
crspmf = pd.merge(crspmf,crspmf_portw,how = 'left', on=['wficn','caldt'])
crspmf['ret'] = crspmf['mret']*crspmf['portw']/crspmf['sum_portw']

crspmf1 =crspmf.drop_duplicates(subset=['wficn','caldt']) 
crspmf1 = crspmf1[['wficn','caldt','tna','ret']]
del [[crspmf,crspmf_portw,crspmf_tna]]

# convert monthly fund return into quarterly fund return
crspmf1['qdate']=crspmf1['caldt']+pd.offsets.QuarterEnd(0)
crspmf1['lnret'] = np.log(crspmf1['ret']+1)
crspmf_qtr = crspmf1.groupby(['qdate','wficn'])['lnret'].sum().reset_index()
crspmf_qtr['ret'] = np.exp(crspmf_qtr['lnret'])-1

crspmf1 = crspmf1.sort_values('caldt')
#sanity check, keep the last observation of tna as quarter end tna
crspmf1 = crspmf1.drop_duplicates(subset=['wficn','qdate']) 
crspmf1 = crspmf1[['qdate','wficn','tna']]
crspmf_qtr = pd.merge(crspmf1,crspmf_qtr,on=['wficn','qdate'],how='inner')  
del [crspmf1]

# quarterly fund flow and tna
crspmf_qtr = crspmf_qtr.sort_values(['wficn','qdate'], ascending=[True, True])
crspmf_qtr['lag_tna'] = crspmf_qtr.groupby(['wficn'])['tna'].shift(1)
crspmf_qtr['f'] = (crspmf_qtr['tna'] - crspmf_qtr['lag_tna']*(crspmf_qtr['ret']+1))/crspmf_qtr['lag_tna']
crspmf_qtr = crspmf_qtr[['qdate','wficn','f','tna']] 


crspmf_qtr['f'] = crspmf_qtr['f'].replace(np.inf, np.nan)
crspmf_qtr['f'] = crspmf_qtr['f'].replace(-np.inf, np.nan)

crspmf_qtr.to_stata("data/intermediate/211101_omega_raw.dta")

#############################################################################
# Step 4                                                                    #
# Calculate Omega_hat t for each quarter from 2001Q1 to 2019Q4              #
# Use CRSP Mutual fund database                                             #
#############################################################################
# generate a list of quarter end values 
date =  pd.date_range(start='2000-01-01', end='2021-01-01', freq='Q') # you can adjust the time period to fit your purpose
crspmf_qtr = pd.read_stata("data/intermediate/211101_omega_raw.dta")
crspmf_qtr = crspmf_qtr.dropna()   
# trim fund flow before regression to avoid extreme values: I capped the fund flow to -90% to 1000%
crspmf_qtr = crspmf_qtr[crspmf_qtr['f']>-0.9]
crspmf_qtr = crspmf_qtr[crspmf_qtr['f']<10]
# delete if total net asset <= 0
crspmf_qtr = crspmf_qtr[crspmf_qtr['tna']>0] 

#  calculate the omega matrix for each quarter
for t in range(len(date)):
    # calculate pct omega
    collist = crspmf_qtr[crspmf_qtr['qdate']==date[t]]['wficn']
    temp = crspmf_qtr[crspmf_qtr['qdate']<=date[t]]
    temp = temp[['qdate','wficn','f']]
    temp = pd.pivot_table(temp,index='qdate',columns='wficn')
    temp.columns = temp.columns.droplevel(0)
    temp = temp.drop(temp.columns.difference(collist),axis=1)
    omega_pct =  temp.cov()
    
    # as long as there is a missing value in covariance matrix, the whole outcome matrix will be NaNs
    # delete the missing rows with the ones that miss the most amount of covariance first
    find_col = pd.DataFrame(omega_pct.count())
    minimum = find_col.min()[0]
    while(len(find_col)>minimum):
        drop_col = find_col.loc[find_col[0]==minimum]
        omega_pct = omega_pct.drop(drop_col.index,axis=1)     
        omega_pct = omega_pct.drop(drop_col.index,axis=0)     
        find_col = pd.DataFrame(omega_pct.count())
        minimum = find_col.min()[0]
    
    # tna_diag
    tna = crspmf_qtr[crspmf_qtr['qdate']==date[t]]
    tna.index = tna['wficn']
    tna = tna.reindex(index=omega_pct.index)    
    tna = tna['tna']
    tna_diag =  np.diag(tna)
    
    # omega
    omega = tna_diag.dot(omega_pct).dot(tna_diag)
    omega = pd.DataFrame(omega,columns=omega_pct.columns,index=omega_pct.index)
    omega.to_stata("data/intermediate/omega/omega{index}.dta".format(index=t))

#############################################################################
# Step 5                                                                    #
# Calculate the stock level fragility                                       #                 
#############################################################################
# merge crspmf_qtr with mflink2, and then merge with s12_crsp
mflink2 = pd.read_sas("data/s12/data2019/MFLink2.sas7bdat")
crspmf_qtr = pd.read_stata("data/intermediate/211101_omega_raw.dta")
mf = pd.merge(mflink2,crspmf_qtr,how = 'inner', left_on=['wficn','fdate'],right_on=['wficn','qdate'])
mf = mf[['fundno','wficn','fdate','tna']] 
mf.to_stata("data/intermediate/mf.dta")

del [[mflink2,crspmf_qtr]]
gc.collect()

vint = pd.read_stata("data/intermediate/s12_crsp2.dta")
mf = pd.read_stata("data/intermediate/mf.dta")
vint = vint.rename(columns={'FUNDNO':'fundno'})
vint = vint.rename(columns={'FDATE':'fdate'})
vintmf = pd.merge(vint,mf,how='inner',on=['fdate','fundno']) 
vintmf['w'] = vintmf['p']*vintmf['SHARES']*vintmf['cfacshr']/(vintmf['tna']*1e6) # tna is in millions of dollars
vintmf = vintmf.sort_values(['PERMNO','fdate','wficn','tna'],ascending=[True,True,True, False])
vintmf = vintmf.drop_duplicates(subset=['PERMNO','fdate','wficn']) #sanity check
# vintmf is fund-stock-quarter level
del [[vint,mf]]
gc.collect()

# for each quarter, create a W matrix
# W is a n*m matrix, n is number of stocks and m is number of funds
# W*Omega*W is a n*n fragility-cofragility matrix, stock price fragility is the diagnal value

date =  pd.date_range(start='2001-01-01', end='2021-01-01', freq='Q')
for t in range(len(date)):
    omega = pd.read_stata("data/intermediate/omega/omega{index}.dta".format(index=t))
    omega = omega.set_index('wficn')
    omega.columns = omega.index
    
    # w is at fund-stock level for each quarter
    w = vintmf[vintmf['fdate']==date[t]] 
    w = w[['PERMNO','wficn','w']]
    # turn long w to wide w
    w = pd.pivot_table(w,index='PERMNO',columns='wficn')
    w.columns = w.columns.droplevel(0)
    # change na obs to 0 - if a fund doesn't hold stock, it will show 0 share holdings
    w = w.replace(np.nan,0)
    # get the common funds between W and omega
    fund_index = np.intersect1d(w.columns,omega.index)
    w = w.reindex(columns=fund_index)
    omega = omega.reindex(columns=fund_index)
    omega = omega.reindex(index=fund_index)
    
    # get the fragility-cofragility matrix
    fragility = w.dot(omega).dot(np.transpose(w))
    # the diagnal value is the stock price fragility we want
    frag = pd.Series(np.diag(fragility), index=[fragility.index])
    frag = pd.DataFrame(frag)
    frag.columns = ['pre_frag']
    frag['fdate'] = date[t]
    frag = frag.reset_index()
    frag.to_stata("data/intermediate/frag/frag{index}.dta".format(index=t))

# construct a panel of stock price fragility (pre_frag) before adjusting for its market cap
frag = pd.read_stata("data/intermediate/frag/frag{index}.dta".format(index=0)) 
for t in range(1,len(date)):
    frag = frag.append(pd.read_stata("data/intermediate/frag/frag{index}.dta".format(index=t)))
del [[t,date]]

# merge fragility data with stock characteristics to get the final product
crsp = pd.read_stata("data/intermediate/crsp_qend.dta")
crsp = crsp.rename(columns={'qdate':'fdate'})  
crsp = pd.merge(crsp,frag,on=['PERMNO','fdate'],how='inner')
# sanity check
crsp = crsp.sort_values(['PERMNO','fdate','me'],ascending=[True,True,False])
crsp = crsp.drop_duplicates(subset=['PERMNO','fdate'])
# calculate the stock price fragility using pre_frag/(market cap)^2
crsp['frag'] = crsp['pre_frag']/(crsp['me']*crsp['me']) 
# calculate the square root of stock price fragility
crsp['sqrt_frag'] = crsp['frag']**(1/2)
crsp.to_stata("data/fragilit.dta")








