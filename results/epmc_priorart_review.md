# Europe PMC prior-art review (fixed search commit cdb1051)

The three pre-registered title/abstract searches returned 4, 3 and 1 records, six distinct articles/preprints. I read their returned abstracts in results/epmc_priorart.json.

- https://europepmc.org/article/MED/37633792 : pancreatic duct organoids from CF pigs were small with no identifiable lumen. Related biological prior art for size/lumen, not human intestinal modulator-by-size FIS.
- https://europepmc.org/article/MED/32485957 : nasal human organoid lumen formation predicts baseline CFTR function. Related and importantly weakens any claim that lumen or size effects are new across all CF organoids.
- https://europepmc.org/article/MED/35086832 : 173-person longitudinal clinical association of *well-level* FIS and progression; not per-organoid starting-size moderation.
- https://europepmc.org/article/MED/32377875 : human epididymis FIS, not intestinal per-organoid modulator-by-size.
- https://europepmc.org/article/MED/40700044 : review discussing SLA, ROMA and FIS for CFSPID/CRMS but not the specific within-organoid drug-by-starting-size test in its abstract.
- https://europepmc.org/article/PPR/PPR1218915 : canine intestinal organoids, not human CF donor modulator-by-size.

H1 novelty-falsification condition: no returned title/abstract explicitly reports the exact result. **Novelty remains unresolved, not established.** Search terms were narrow and titles/abstracts omit details; full-text and wider systematic review would be required to claim novelty. The biologically related nasal and pancreatic findings should remain visible in the paper.

## Wider full-text prior-art correction, 25 Sep 2026
The six-record preregistered title/abstract search is NOT a complete novelty search. Calucho et al. 2021, https://www.nature.com/articles/s41598-021-94798-x , already tested initial nasospheroid size versus individual FSK shrinking and found no significant association among 138 WT responders; its supplementary slope result likewise failed (r=.04124, p=.6298). It also studied modulators in seven CF subjects. Kim et al. 2020, https://www.slas-discovery.org/article/S2472-5552(22)06605-9/fulltext , assessed initial size and growth response in colorectal tumor organoids and warned that small-object segmentation may inflate variation. These do not duplicate the proposed donor-matched intestinal *modulator-by-size* interaction, but they bar a broad claim that organoid starting-size effects or single-object size analysis are new. The stronger matched-block intestinal test still fails, and no established novel discovery is claimed.

## Source-level airway FIS prior art, 26 Sep 2026
Berical et al. 2022 (https://www.nature.com/articles/s41467-022-31854-8) tracked individual iPSC-derived airway spheroid responses to FSK/vehicle and studied modulator FIS in multiple patient lines. The publicly released source workbook (https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-31854-8/MediaObjects/41467_2022_31854_MOESM7_ESM.xlsx) has Figure 2 experiment-level response summaries, Figure 3B separate baseline sphere sizes, and Figure 3D patient-line/condition experiment values, not joinable object tracks by size and treatment. The workbook was inspected after preregistration at commit 096acf5; detailed audit: `results/berical_source_suitability.md`. It limits broad novelty, not an external test of the intestinal donor-matched size interaction.
