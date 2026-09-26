# Same-accession size H1 weakens under stricter cell-count checks

The check was locked in `results/preregistration_block_eligibility_sensitivity.md` at b9978eb before computation. Data are the original OrgaSegment DIS accession, not independent replication. The only changed parameter is the number of measured organoids required in **each** of four drug/control x small/large cells within a donor/plate/forskolin-dose block. The originally registered primary test used 3, gave 31 blocks, 12 donors, 9 positive donor medians, exact one-sided p=.072998 and **FAIL**.

| Minimum per cell | Blocks | Eligible donors | Positive donors | One-sided sign p | Frozen >=8-donor verdict |
|---:|---:|---:|---:|---:|---|
| 3 | 31 | 12 | 9 | .072998 | FAIL |
| 5 | 28 | 12 | 8 | .193848 | FAIL |
| 10 | 15 | 6 | 6 | .015625 | UNINTERPRETABLE |
| 20 | 9 | 3 | 3 | .125000 | UNINTERPRETABLE |

The nominal p=.015625 at 10 per cell is **not** a discovery: only 6 donors remain, below the preregistered minimum 8, and the threshold was examined as sensitivity after the failed primary result. At 5 per cell, four of 12 donor medians are nonpositive and the sign test weakens. Cell eligibility is consequential and selective, especially at stricter cutoffs. The analysis does not identify whether the missing donors differ biologically or due to imaging; it adds no mechanism or clinical use. Full donor and threshold outputs: `results/block_eligibility_sensitivity.json`.
