
# Tiger snake venom transcriptomics project

2019-05-21

Jenna Crowe-Riddell
Stephen Pederson
Vick Thompson

## To do list
- [X] Upload data
- [X] Trim reads
- [X] Align transcripts
- [ ] Assemble transcripts:  Stringtie using ggf file for TS genome*
- [ ] DE analysis

*This accounts for transcripts that don't 'fit' & will output a gtf file
of each sample then create a combined gtf with location of those transcripts
We will have to manually inspect & adjust parameters 

## Workflow 
1. Assemble transcripts for each aligned reads & merge into singe gtf file (StringTie)
2. Manually inspect & adjust parameters/transcripts, consider populations/individuals, multiple layers of merging
3. Estimate transcript abundances (Feature counts)
4. Visualise samples - PCA plot, are samples grouping by population / prey? Differences in triplicate samples? If no, then merge
5. Perform differential expression analysis for mainland versus island populations (EdgeR) 

## Samples
- TASMANIA = 5
- MAINLAND = 5
- REEVESBY = 5
- FRANKLIN = 5

Data located in Jenna's FAST on Phoenix:
/fast/users/a1662801/20190520_TigerSnake_RNASeq/0_rawData/fastq
