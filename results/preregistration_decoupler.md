# Pre-registration: TF/pathway activity of organoid deficits (decoupler 2.1.4, CollecTRI + PROGENy via OmniPath); committed before running
Input: per-organ consensus deficit vector delta_g = z(GTEx organ) - z(organoid), strict series (as src/strict_gsea.py).
Method: decoupler univariate linear model (ulm) per organ.
Positive-control hypothesis H: the organ's lineage master TFs rank in the top 10% of CollecTRI TFs by deficit activity:
liver HNF4A, HNF1A, NR1H4, CEBPA; kidney HNF1B, HNF4A, PAX2, PAX8; lung NKX2-1, FOXA2; colon CDX2, CDX1, HNF4A; brain NEUROD2, NEUROD6, TBR1.
Score = fraction of listed TFs (present in CollecTRI) in top 10%; binomial one-sided test vs 0.10 pooled over organs. PASS if p < 0.05.
PROGENy pathway activities are exploratory (no hypothesis).
