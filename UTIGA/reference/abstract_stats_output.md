
# SECTION 1: Data Integrity & Linkage


## 1.1 File Timestamps



|File                            |mtime               |   size|
|:-------------------------------|:-------------------|------:|
|results/clinical/status_map.csv |2025-12-02 12:38:46 |  29520|
|results/mlst/mlst_with_meta.csv |2025-12-02 12:43:56 | 306884|

## 1.2 Row Counts Loaded

Status Map: 274 rowsMLST File:  382 rows
## 1.3 Linkage Yield

Total Episodes (status_map): 274ASB/UTI Episodes: 233Episodes with ST data linked: 183 (66.8%)ASB/UTI Episodes with ST linked: 152 (65.2%)Duplicate isolates collapsed: 199
# SECTION 2: Cohort Size by Follow-up Depth



|Recall_Depth   | Participants| ASB_UTI_Episodes| Linked_Episodes|TP_Range |Epi_Per_Ppt_Median_Range |
|:--------------|------------:|----------------:|---------------:|:--------|:------------------------|
|>=2 Timepoints |           92|              230|             150|2-6      |2 (0-6)                  |
|>=3 Timepoints |           54|              176|             104|3-6      |3 (1-6)                  |
|>=4 Timepoints |           24|               97|              54|4-6      |4 (1-6)                  |

# SECTION 3: ST Distribution (ASB/UTI Episodes)


## K>=2 Subset



|ST | Episodes| Participants| Pct_Episodes| Pct_Participants|
|:--|--------:|------------:|------------:|----------------:|
|-  |       24|           13|         16.0|             14.1|
|43 |       20|           10|         13.3|             10.9|
|4  |       12|            7|          8.0|              7.6|
|6  |       11|            4|          7.3|              4.3|
|1  |        6|            4|          4.0|              4.3|
|22 |        6|            3|          4.0|              3.3|
|2  |        5|            3|          3.3|              3.3|
|87 |        5|            4|          3.3|              4.3|
|10 |        4|            2|          2.7|              2.2|
|33 |        4|            2|          2.7|              2.2|
Distinct STs observed: 34Top 5 STs comprise 48.7% of all linked episodes.
## K>=3 Subset



|ST | Episodes| Participants| Pct_Episodes| Pct_Participants|
|:--|--------:|------------:|------------:|----------------:|
|-  |       17|            8|         16.3|             14.8|
|6  |       11|            4|         10.6|              7.4|
|43 |       10|            5|          9.6|              9.3|
|1  |        6|            4|          5.8|              7.4|
|4  |        6|            3|          5.8|              5.6|
|22 |        4|            2|          3.8|              3.7|
|33 |        4|            2|          3.8|              3.7|
|73 |        4|            2|          3.8|              3.7|
|87 |        4|            3|          3.8|              5.6|
|2  |        3|            2|          2.9|              3.7|
Distinct STs observed: 27Top 5 STs comprise 48.1% of all linked episodes.
## K>=4 Subset



|ST  | Episodes| Participants| Pct_Episodes| Pct_Participants|
|:---|--------:|------------:|------------:|----------------:|
|6   |       11|            4|         20.4|             16.7|
|-   |        6|            2|         11.1|              8.3|
|33  |        4|            2|          7.4|              8.3|
|43  |        4|            2|          7.4|              8.3|
|36  |        3|            2|          5.6|              8.3|
|506 |        3|            2|          5.6|              8.3|
|722 |        3|            1|          5.6|              4.2|
|87  |        3|            2|          5.6|              8.3|
|1   |        2|            1|          3.7|              4.2|
|197 |        2|            1|          3.7|              4.2|
Distinct STs observed: 17Top 5 STs comprise 51.9% of all linked episodes.
# SECTION 4: UTI Proportion by ST


## K>=2 Subset



|ST | n_total| n_UTI| n_ASB| UTI_Prop|
|:--|-------:|-----:|-----:|--------:|
|87 |       5|     1|     4|     20.0|
|6  |      11|     2|     9|     18.2|
|-  |      24|     3|    21|     12.5|
|43 |      20|     1|    19|      5.0|
|1  |       6|     0|     6|      0.0|
|2  |       5|     0|     5|      0.0|
|22 |       6|     0|     6|      0.0|
|4  |      12|     0|    12|      0.0|
> Among STs with >=5 episodes, ST87 had the highest UTI proportion (1/5 = 20%), compared with... 
## K>=3 Subset



|ST | n_total| n_UTI| n_ASB| UTI_Prop|
|:--|-------:|-----:|-----:|--------:|
|6  |      11|     2|     9|     18.2|
|-  |      17|     2|    15|     11.8|
|43 |      10|     1|     9|     10.0|
|1  |       6|     0|     6|      0.0|
|4  |       6|     0|     6|      0.0|
> Among STs with >=5 episodes, ST6 had the highest UTI proportion (2/11 = 18.2%), compared with... 
## K>=4 Subset



|ST | n_total| n_UTI| n_ASB| UTI_Prop|
|:--|-------:|-----:|-----:|--------:|
|6  |      11|     2|     9|     18.2|
|-  |       6|     1|     5|     16.7|
> Among STs with >=5 episodes, ST6 had the highest UTI proportion (2/11 = 18.2%), compared with... 
# SECTION 5: Within-Participant Stability


## K>=2 Subset

Participants with >=2 linked episodes: 65Stable ST: 56 (86.2%)At least one switch: 9 (13.8%)Median switches (IQR): 0 (0)

Table: Top ST Switches

|Prev_ST |ST   |  n|
|:-------|:----|--:|
|-       |1106 |  1|
|-       |87   |  1|
|1112    |4    |  1|
|2       |83   |  1|
|3       |658  |  1|
|36      |6    |  1|
|567     |506  |  1|
|6       |36   |  1|
|6       |87   |  1|
|938     |882  |  1|

## K>=3 Subset

Participants with >=2 linked episodes: 45Stable ST: 39 (86.7%)At least one switch: 6 (13.3%)Median switches (IQR): 0 (0)

Table: Top ST Switches

|Prev_ST |ST   |  n|
|:-------|:----|--:|
|-       |1106 |  1|
|2       |83   |  1|
|36      |6    |  1|
|567     |506  |  1|
|6       |36   |  1|
|6       |87   |  1|
|938     |882  |  1|

## K>=4 Subset

Participants with >=2 linked episodes: 23Stable ST: 20 (87%)At least one switch: 3 (13%)Median switches (IQR): 0 (0)

Table: Top ST Switches

|Prev_ST |ST  |  n|
|:-------|:---|--:|
|36      |6   |  1|
|567     |506 |  1|
|6       |36  |  1|
|6       |87  |  1|

# SECTION 6: VF Burden Analysis

Loaded VF data. Rows:  32182 


|ST | N_Episodes| Median_Burden| IQR_Burden| UTI_Burden_Med| ASB_Burden_Med|
|:--|----------:|-------------:|----------:|--------------:|--------------:|
|-  |         25|          69.0|      24.00|           61.0|           69.0|
|1  |          6|          81.0|       2.00|             NA|           81.0|
|10 |          4|         100.5|       1.00|          100.5|          100.5|
|2  |          5|          53.0|      15.00|             NA|           53.0|
|22 |          6|         100.0|       2.25|             NA|          100.0|
|33 |          4|          89.5|      19.00|             NA|           89.5|
|4  |         12|          91.0|      18.25|             NA|           91.0|
|43 |         21|          83.0|       8.00|           80.0|           83.0|
|6  |         11|          77.0|       6.50|           80.0|           77.0|
|87 |          5|          68.0|      18.00|           86.0|           59.0|

# FINAL DELIVERABLE: Abstract Numbers

Run the script to see the generated lists above. Copy specific rows as needed.
Methods One-Liner:
'We analyzed E. coli ST dynamics in 24 participants with >=X follow-up episodes. STs were determined from whole-genome assemblies and linked to clinical episodes (UTI/ASB) via timepoint-matched sequence typing.'