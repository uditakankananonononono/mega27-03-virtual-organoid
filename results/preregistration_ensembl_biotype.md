# Pre-registration: coding-composition check of missing organoid genes

Registered before querying Ensembl REST. Source: Ensembl REST POST /lookup/id, https://rest.ensembl.org/documentation/info/lookup_post . Held-fixed gene sets: top-100 Ensembl IDs per organ in results/strict_organ_deficits.csv. Map each ID through the live lookup endpoint, count `biotype == protein_coding`, exclude missing IDs.

Primary H1: liver top-100 deficit genes have a higher protein-coding fraction than brain top-100 deficit genes, one-sided Fisher exact p < 0.05. This tests whether the metabolic-parenchyma description may reflect composition differences in coding vs non-coding top deficits. This is not proof of a biological mechanism, because gene biotype is a coarse annotation and the organ contrast was selected after seeing deficits. Report all five coding fractions as descriptive. If either mapping <80/100, verdict UNINTERPRETABLE. Preserve response provenance and endpoint/version metadata when provided.
