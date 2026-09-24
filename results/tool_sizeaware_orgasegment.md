# vorganoid sizeaware on OrgaSegment DIS (default settings)
Command: PYTHONPATH=src python -m vorganoid.cli sizeaware data/raw/orgasegment/dis_merged_A0.csv --drug VX445_VX661_VX770 --forskolin-col forskolin_concentration_µM
Default settings: 4 quartile bins on pooled drug+DMSO A0, at least 5 organoids per donor x condition x bin, 1000 donor-bootstrap resamples.
Result: 14 evaluable donors (3 lack enough small organoids), mean attenuation (top-quartile minus bottom-quartile effect) 0.249, 95% CI [0.165, 0.326], positive in 12/14.
This agrees with the octile analysis in finding_size_dependent_swelling.md (top minus bottom 0.251 [0.184, 0.323]).
