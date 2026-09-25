# Pre-registration: organoid deficit genes vs Open Targets organ-disease genes (committed before any Open Targets download)
Input: results/strict_organ_deficits.csv (top-100 deficit genes per organ, 5 organs).
Data: Open Targets Platform GraphQL, top-500 targets by overall association score for liver disorder MONDO_0005154, kidney disorder MONDO_0005240, lung disorder MONDO_0005275, brain disorder MONDO_0005560, colonic disorder MONDO_0003409.
H1 (disease relevance, organ-specific): for each organ, its deficit list overlaps its own organ-disease set more than the other four organs' deficit lists do (one-sided Fisher, own list vs pooled other lists). Primary statistic: diagonal-dominance D = mean own-overlap minus mean off-diagonal overlap in the 5x5 overlap matrix, p by 10,000 permutations of organ labels across deficit genes.
Pass: permutation p < 0.05. Otherwise H1 fails and is reported as a negative.
