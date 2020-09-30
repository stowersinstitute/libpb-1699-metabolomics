from argparse import ArgumentParser

from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.decomposition import PCA as PCA
from pandas import DataFrame, concat
from seaborn import violinplot, catplot, distplot
from numpy import log, array
from numpy import log10

import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/3810865/matplotlib-unknown-projection-3d-error?noredirect=1
from mpl_toolkits.mplot3d import Axes3D
import os
#https://www.delftstack.com/howto/matplotlib/how-to-make-different-subplot-sizes-in-matplotlib/
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap

from cavefinomics import AstyanaxMe
from cavefinomics import AstyanaxLi

parser = ArgumentParser(description="PCA.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--out-dir", type=str, help="Out dir.")
args = parser.parse_args()

ame = AstyanaxMe(args.astyanax, kegg_compounds_file=args.compounds, sample_sheet_path=args.sample_sheet, hmdb_file=args.hmdb)

ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
# remove "pools" columns
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data = astyanax_data.apply(log10)
astyanax_data.columns = (f'{c} {n}' for c,n in zip(astyanax_data.columns,list(range(1,7))*27))
astyanax_data = astyanax_data.reindex(sorted(astyanax_data.columns), axis=1)
#print(astyanax_data)

# get data and compute common metabolites
astyanax_lipids = ali.lmdata.iloc[:,ali.non_numeric_cols-1:]
astyanax_lipids = astyanax_lipids.apply(log10)
astyanax_lipids.columns = (c[0] + ' ' + c[-2] + ' ' + ' '.join(c[1:-2]) + ' ' + c[-1] for c in (cc.split() for cc in astyanax_lipids.columns))
astyanax_lipids = astyanax_lipids.reindex(sorted(astyanax_lipids.columns), axis=1)
#astyanax_lipids.columns = astyanax_data.columns
#print(astyanax_lipids)

#print(astyanax_data.columns,astyanax_lipids.columns)
astyanax_data = concat((astyanax_data,astyanax_lipids),axis=0)
#print(astyanax_data)

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

if args.exclude_outlier:
    for o in outliers:
        astyanax_data = astyanax_data.loc[:,~astyanax_data.columns.str.contains(o)]

#astyanax_data = DataFrame(StandardScaler().fit_transform(astyanax_data.values), index=astyanax_data.index, columns=astyanax_data.columns)

#gridspec_kw = {"height_ratios":[1.,1.], "width_ratios" : [1.]}
#fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8, 12), gridspec_kw=gridspec_kw, projection='3d')
#fig = plt.figure()
#fig = plt.figure(figsize=plt.figaspect(0.5))

first_cmap = LinearSegmentedColormap.from_list('First', ["#224b67", "#3498db", "#98cbee"])
second_cmap = LinearSegmentedColormap.from_list('Second', ["#493054", "#9b59b6", "#e2bdf1"])
third_cmap = LinearSegmentedColormap.from_list('Third', ["#1f4736", "#57997d", "#8eddbb"])

def plot_pca(astyanax_data,selectors,subset_category,ax):
    global handles, labels
    joined_dataset = astyanax_data

    pca_data = StandardScaler().fit_transform(joined_dataset.transpose())

    # https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
    pca = PCA(n_components=3)
    components = DataFrame(pca.fit_transform(pca_data), columns=['Comp1', 'Comp2', 'Comp3'])


    import matplotlib
    for color_theme,tissue in zip([first_cmap, second_cmap, third_cmap], ['Brain', 'Muscle', 'Liver']):
        cmap = matplotlib.cm.get_cmap(color_theme)
        for color_pos,var2 in selectors:
            p = components.iloc[joined_dataset.columns.isin(astyanax_data.columns) & joined_dataset.columns.str.contains(tissue) & joined_dataset.columns.str.contains(var2)]
            label = ', '.join((tissue, var2))
            #print(p)
            #https://www.thetopsites.net/article/52439590.shtml
            ax.scatter(p['Comp1'],p['Comp2'],p['Comp3'], label=label, depthshade=False, color=cmap(color_pos))
            handles, labels = ax.get_legend_handles_labels()

    ax.set_xlabel('Comp1')
    ax.set_ylabel('Comp2')
    ax.set_zlabel('Comp3')

    ax.set_title(f'Lipid & Primary Metabolite PCA\n by Tissue & {subset_category}')

    #ax.legend()

    #plt.show()
    #plt.cla()



for subset_category,selectors,filename in zip(
    ['Condition','Population'],
    [
      list(zip([0.0, 0.5, 1.0], ['4d Starved', '30d Starved', 'Refed'])),
      list(zip([0.0, 0.5, 1.0], ['Pachon', 'Tinaja', 'Surface']))
    ],
    ['pca-global-condition.pdf','pca-global-pop.pdf']):
    fig = plt.figure(figsize=(9, 7))
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                         width_ratios=[2, 1])
    ax = fig.add_subplot(spec[0],projection='3d')
    handles, labels = None, None
    plot_pca(astyanax_data,
            selectors,
            subset_category,
            ax)
    ax2 = fig.add_subplot(spec[1])
    ax2.legend(handles, labels, loc='center right')
    ax2.set_yticks([], minor=[])
    ax2.set_xticks([], minor=[])
    ax2.patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax2.spines[s].set_visible(False)
    plt.savefig(os.path.join(args.out_dir,filename),bbox_inches='tight',transparent=True,pad_inches=0)
