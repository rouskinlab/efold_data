# eFold training data

Here's a breakdown of the data we used to train eFold.

| Training stage   | Name on HuggingFace   | Source                     | Method                          | Number of sequences we used   | Families                      | L ∈ [10, 199]   | L ∈ [200, 499]   | L ∈ [500, 999]   | L ∈ [1000, 1999]   |   L ∈ [2000, inf] |
|:-----------------------|:----------------------|:---------------------------|:--------------------------------|:----------------------|:------------------------------|:----------------|:-----------------|:-----------------|:-------------------|------------------:|
| Pre-training     | rnacentral_synthetic  | Sequences from RNA central | RNAstructure                    | 226'729               | All known families            | 176'486         | 49'463           | 780              | 0                  |                 0 |
| Pre-training     | ribo500-blast         | Ribonanza Competition      | RNAstructure + DMS and/or SHAPE | 46'060                | Unlabelled                    | 46'049          | 11               | 0                | 0                  |                 0 |
| Pre-training     | bpRNA-1m              | bpRNA-1m                   | Covariance analysis             | 66'715                | Unlabelled, sRNA, tRNA        | 48'090          | 6'167            | 2'829            | 9'260              |               369 |
| Fine-tuning      | pri_miRNA             | This work                  | RNAstructure + DMS              | 1'098                 | pri-miRNA                     | 0               | 1'098            | 0                | 0                  |                 0 |
| Fine-tuning      | human_mRNA            | This work                  | RNAstructure + DMS              | 1'456                 | mRNA                          | 0               | 493              | 882              | 81                 |                 0 |
| Testing          | PDB                   | PDB                        | NMR, crystallography            | 356                   | Short non-coding RNA          | 343             | 6                | 6                | 1                  |                 0 |
| Testing          | viral_fragments       | Peer-reviewed literature   | RNAstructure + DMS              | 40                    | Viral RNA                     | 12              | 17               | 11               | 0                  |                 0 |
| Testing          | lncRNA                | Bugnon and al, 2022        | RNAstructure + DMS              | 10                    | Long non-coding RNA           | 0               | 2                | 1                | 7                  |                 0 |
| Testing          | archiveII_blast       | Archive II                 | Covariance analysis             | 355                   | rRNA, tRNA, tmRNA, unlabelled | 242             | 65               | 43               | 5                  |                 0 |

## Notes
- all of the raw data in on [this google drive](https://drive.google.com/drive/folders/1pKUBGlvcft4WsUSztaUCOXcyGi9a8NUy)
- the normalisation, filtering and replicates matching of SEISMIC data is made using `seismic__postprocessing.py`
