from argparse import ArgumentParser

from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.decomposition import PCA as PCA
from pandas import DataFrame, concat
from seaborn import violinplot, catplot, distplot
from numpy import log, array, cov
from numpy.linalg import LinAlgError, cholesky


import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/3810865/matplotlib-unknown-projection-3d-error?noredirect=1
from mpl_toolkits.mplot3d import Axes3D

#from cavefinomics import AstyanaxMe, MaMammalianMe
from cavefinomics import AstyanaxMe
from itertools import repeat, chain
from scipy.spatial.distance import mahalanobis
from scipy.linalg import inv
from pandas import Series
from scipy.stats import chi
from seaborn import swarmplot, violinplot

parser = ArgumentParser(description="Mahalanobis outlier.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--out-plot", type=str, help="Outplot.")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

#https://franciscomorales.org/2017/03/09/how-to-python-calculate-mahalanobis-distance/

# more rigorous method:
#https://towardsdatascience.com/multivariate-outlier-detection-in-high-dimensional-spectral-data-45878fd0ccb8

# get data and compute common metabolites
astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
# remove "pools" columns
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]

astyanax_data.columns = [f'{c} {n}' for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,7),3*3*3)))]
#print(astyanax_data)
#print(astyanax_data.columns)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['30d Starved', '4d Starved', 'Refed']

# mean within groups
maha_local_means = []
maha_local_means_labels = []
for tissue in tissues:
    tissue_subset = astyanax_data.loc[:,astyanax_data.columns.str.contains(tissue)]
    v = inv(cov(tissue_subset))
    print(tissue)
    #print(tissue_subset.transpose().values)
    try:
        x = cholesky(v)
        print('pos def')
    except LinAlgError:
        print('Not pos def')
    for pop in pops:
        for condition in conditions:
            subset = tissue_subset.loc[:,tissue_subset.columns.str.contains(pop) & tissue_subset.columns.str.contains(condition)]
            #print(subset)
            mean = subset.mean(axis=1)
            #print(mean)
            print(len(mean))
            print(v.shape)
            #stop
            #print(len(subset[x]))
            maha_local_means += [float(mahalanobis(subset[x].values,mean,v)) for x in subset.columns]
            print(maha_local_means[-1])
            maha_local_means_labels += list(subset.columns)

maha_local_means = Series(maha_local_means,index=maha_local_means_labels).sort_values(ascending=False)
print(maha_local_means)
#maha_local_means.name = 'Mahalanobis'
maha_local_means.to_csv('/tmp/maha_local_means.csv')

maha_local_means_flat = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            subset = maha_local_means.loc[maha_local_means.index.str.contains(pop) & maha_local_means.index.str.contains(tissue) & maha_local_means.index.str.contains(condition)]
            for x in subset:
                maha_local_means_flat.append({'Population':pop,'Tissue':tissue,'Condition':condition,'MD (within group)':x})
maha_local_means_flat = DataFrame(maha_local_means_flat)
#print(maha_local_means_flat)

# chi square cutoff
#maha_local_means = maha_local_means.apply(lambda x: x*x)
cutoff = chi.ppf(0.975,len(astyanax_data.index)/3)
#https://seaborn.pydata.org/examples/scatterplot_categorical.html
#https://seaborn.pydata.org/examples/simple_violinplots.html
fig,ax = plt.subplots()
swarmplot(x='Tissue',y='MD (within group)',hue='Population',data=maha_local_means_flat,palette=['firebrick','goldenrod','dodgerblue'])
plt.gca().axhline(cutoff,0.,1.,linestyle='--',color='k',alpha=0.5)
plt.gca().set_xlabel('')

plt.gca().patch.set_visible(False)
for s in ["top", "bottom", "left", "right"]:
    plt.gca().spines[s].set_visible(False)

plt.gca().annotate(
    f'Tinaja Liver\nRefed 6',
    xy=(0.825, 0.9),
    xycoords="axes fraction",
    xytext=(0.65, 0.97),
    textcoords="axes fraction",
    arrowprops=dict(facecolor="black", shrink=0.02, linewidth=0.0, alpha=0.75),
    horizontalalignment="center",
    verticalalignment="top",
    bbox=dict(fc="white", alpha=0.5, boxstyle="round"),
)

plt.gca().annotate(
    f'Pachon Muscle Refed 5',
    xy=(0.5, 0.955),
    xycoords="axes fraction",
    xytext=(0.2, 0.975),
    textcoords="axes fraction",
    arrowprops=dict(facecolor="black", shrink=0.05, linewidth=0.0, alpha=0.75),
    horizontalalignment="center",
    verticalalignment="top",
    bbox=dict(fc="white", alpha=0.5, boxstyle="round"),
)

plt.gca().annotate(
    f'Pachon Liver\n30d Fasted 3',
    xy=(0.825, 0.87),
    xycoords="axes fraction",
    xytext=(0.68, 0.75),
    textcoords="axes fraction",
    arrowprops=dict(facecolor="black", shrink=0.00, linewidth=0.0, alpha=0.75),
    horizontalalignment="center",
    verticalalignment="top",
    bbox=dict(fc="white", alpha=0.5, boxstyle="round"),
)

plt.savefig(args.out_plot)
