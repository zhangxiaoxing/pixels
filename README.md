# FNC-lab codes for Neuropixels data analysis

---
This repository hosts necessary codes for the statistics and visualization in FNC lab neuropixels studies. 

### Dependencies
- MATLAB (>R2020a)
	- [FieldTrip Toolbox](https://www.fieldtriptoolbox.org)
	- [npy-matlab](https://github.com/kwikteam/npy-matlab)
- Python
	- [Elephant - Electrophysiology Analysis Toolkit](https://github.com/INM-6/elephant)
- Java SE
	- [Gephi and GephiToolkit](https://gephi.org)  
  
  
---

### Bleeding edge data pipeline

- `sync.py`  
For extract behavior events from neuropixels binary data


- `jpsth/extract_waveform.m`  
Extract waveform from neuropixels binary data  

- `jpsth/+ephys/pixFlatFR.m`  
Generate *FR_All_{bin_width}.hdf5* firing rate file.  
Include tagged behavioral trials, no criteria applied yet  

- `per_sec/per_sec_stats.gen_align_files()`  
Cluster_id to brain region localization
- `per_sec/per_sec_stats.gen_selectivity_stats`  
Generate *transient_{delay}.hdf5* dataset  
Well-trained performace criteria applied  

- `+bz/x_corr_bz.m`  
Function coupling based on ccg
- `+bz/sums_conn.m`  
(map->)reduce functional coupling data  

-  `+bz/reg_conn_bz.m`  
Assign selectivity & region for function coupling