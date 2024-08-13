import allel
import numpy as np
import pandas as pd

# import vcf
callset = allel.read_vcf('/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/vcfs/subsetch15.mediterranean.vcf.gz',  fields=['calldata/GT', 'samples'])

n_components = 5

g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT'])) # genotype array
ac = g.count_alleles() # allele count array
samples = callset['samples']

# filter singletons and multi-allelic sites
mask = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
masked_g = g.compress(mask, axis=0)

# create the PCA input (i.e. number of alternative alleles per individual at each site)
gn = masked_g.to_n_alt()

# run the PCA
coords, model = allel.pca(gn, n_components=n_components, copy=True, scaler='patterson', ploidy=2)

# Coordinates of PC1 and PC2 are recorded in coords[:, 0] and coords[:, 1]. The variance explained by each component is recorded in model.explained_variance_ratio_[n-1], where n is the number of the component.
coords_df = pd.DataFrame(coords, columns=["PC" + str(i + 1) for i in range(n_components)])
samples_df = pd.DataFrame(samples, columns=["id"])
pca = pd.concat([samples_df, coords_df], axis=1)


# output to file
pca.to_csv("/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/pca.eu40.vs.dchry2.3.Inversion.csv", index=False)

with open('/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/pca.eu40.vs.dchry2.3.Inversion.explained_variance.txt', 'w') as f:
    f.write("component\tvariance_explained\n")
    for idx, line in enumerate(model.explained_variance_ratio_.tolist()):
        f.write(f"{idx + 1}\t{line}\n")
