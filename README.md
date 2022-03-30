# StockPriceFragility_GreenwoodThesmar2011
   
The chain event caused by the collapse of Archegos Capital Management raised investors' attention on the stock price fragility. The stock price fragility measures the likelyhood of stock prices to fluctuate due to non-fundamental demand shocks from correlated fund flow.   
This Python code calculates the stock price fragility following Greenwood and Thesmar (2011).      

The input data are accessible from wrds database. I used CRSP MSF, CRSP Mutual fund dataset, Thomson s12 type3, Thomson MFLink1, Thomson MFLink2. The output file is fragility.dta, which is a stock-quarter level dataset.   
          
## Methodology 
### Step 1
Step 1 calculates the quaterly stock return, shares outstanding, and market capitalization using CRSP MSF.   
### Step 2
Step 2 merge the stock-quarter level information with Thomson s12 type3, which is a stock-fund-quarter level data that shows fund holdings.
### Step 3
Step 3 calculates the percentage fund flow using CRSP Mututal fund dataset and Thomson MFLink1. Thomson MFLink1 generates a fund-level identifier for the crsp portfolio, which can be merged with MFLink2 in order to merge with Thomson s12 type3.    
The percentage fund flow is calculated as follows:   
![img](http://latex.codecogs.com/svg.latex?f%5E%5C%25_%7Bi%2Ct%7D%3D%5Cfrac%7BTNA_%7Bi%2Ct%7D-TNA_%7Bi%2Ct-1%7D%281%2BR_%7Bi%2Ct-1%7D%29%7D%7BTNA_%7Bi%2Ct-1%7D%7D)    
for fund i in quarter t. Where TNA is total net assets and R is the fund return.   
### Step 4
Step 4 calculats omega_hat for each quarter. Omega_hat is the covariance matrix transformed from percentage fund flows, calculated as follows:  
![img](http://latex.codecogs.com/svg.latex?%5Chat%7B%5COmega%7D_%7Bt%7D%3Ddiag%28TNA_%7Bt%7D%29%5COmega%5E%7B%5C%25%7D_%7Bt%7Ddiag%28TNA_%7Bt%7D%29)
### Step 5   
Finally, stock-level fragility is estimated by taking the diagonal value of Gt:   
![img](http://latex.codecogs.com/svg.latex?G_t%3DW%27_t%5Chat%7B%5COmega%7D_%7Bt%7DW_t)    
where Wt is a stock-fund matrix that measures the stock holdings for each fund at quarter t. Last, I scale the value by squared market capitalization value to get the stock price fragility.

