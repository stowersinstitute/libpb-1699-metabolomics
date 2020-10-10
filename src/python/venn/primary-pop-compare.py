from argparse import ArgumentParser
import os
from pandas import read_csv, concat, DataFrame
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
import matplotlib.patches as patches

parser = ArgumentParser(description="Venn diagrams.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--output-dir", type=str, help="Out dir")
parser.add_argument("--level", type=float, help="p value cutoff.")
args = parser.parse_args()

if args.exclude_outlier:
    outlier = 'kein-Ausreißern'
else:
    outlier = 'mit-Ausreißern'

#colors = ['firebrick','goldenrod','mediumseagreen','rosybrown']
#https://colorhunt.co/palette/184190
#colors = ['#0e9aa7','#3da4ab','#f6cd61','#fe8a71']
#https://colorhunt.co/palette/189889
colors = ['#d92027','#ff9234','#ffcd3c','#35d0ba']
#https://colorhunt.co/palette/183823
#colors = ['#ffc4a3','#ff9a76','#f96d80','#bb596b']


conditions = {'4d':'4d Starved', '30d':'30d Starved', 'Ref':'Refed'}
pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
comparisons = {'PvS':('Pachon','Surface'),'TvS':('Tinaja','Surface'),'PvT':('Pachon','Tinaja')}
primary_categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}

metabolites = {}
for i,cond in zip(range(3),conditions):
    for j,tissue in zip(range(1,4),tissues):
        for comp,groups in comparisons.items():
            sig = []
            for cat,category in primary_categories.items():
                data = read_csv(f"out/work/primary/glm/singlefactor/{outlier}/{cat}/{tissue}/{cond}/{comp}.csv").rename({'Pr(>|z|)':'p'},axis=1)
                data = data.loc[data['p']<args.level]
                cols = list(data.columns)
                cols[0] = 'Name'
                data.columns = cols
                data = data.set_index('Name')
                sig.append(data)
            sig = concat(sig)
            metabolites[cond,tissue,comp] = sig

lipids = {}
for i,cond in zip(range(3),conditions):
    for j,tissue in zip(range(1,4),tissues):
        for comp,groups in comparisons.items():
            data = read_csv(f"out/work/lipids/glm/singlefactor/{outlier}/{tissue}/{cond}/{comp}.csv").rename({'Pr(>|z|)':'p'},axis=1)
            data = data.loc[data['p']<args.level]
            cols = list(data.columns)
            cols[0] = 'InChIKey'
            data.columns = cols
            data = data.set_index('InChIKey')
            lipids[cond,tissue,comp] = data

def make_merged_table(cond,tissue):
    c = list(comparisons.keys())
    merged = metabolites[cond,tissue,c[0]].merge(metabolites[cond,tissue,c[1]],on='Name',suffixes=c[:2])
    d = DataFrame(metabolites[cond,tissue,c[2]]).add_suffix(c[2])
    d = d.rename({f'Name{c[2]}':'Name'},axis=1)
    merged = merged.merge(d,on='Name')
    merged = merged.dropna()
    return merged

height_ratios = [1]*len(conditions)
gridspec_kw = {"height_ratios":height_ratios, "width_ratios" : [1.,]+[3.,3.,3.]*2}
fig,ax = plt.subplots(nrows=len(conditions),ncols=len(tissues)*2+1,figsize=(26, 12), sharex=False, sharey=False, gridspec_kw=gridspec_kw)
for i,cond in zip(range(3),conditions):
    ax[i,0].text(0.5,0.5,conditions[cond],size=16,rotation=90.,ha='center',va='center',transform=ax[i,0].transAxes)
    # hide graphics
    ax[i,0].set_yticks([], minor=[])
    ax[i,0].set_xticks([], minor=[])
    ax[i,0].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,0].spines[s].set_visible(False)
    for j,tissue in zip(range(1,7,2),tissues):
        def get_sig_metabs2(cond1,cond2):
            table = make_merged_table(cond,tissue)
            return set(table.loc[(table[f'Estimate{cond1}'] > 0.) == (table[f'Estimate{cond2}'] > 0.)].index)

        def get_sig_metabs3():
            table = make_merged_table(cond,tissue)
            return set(table.loc[
              ((table[f'Estimate{cond1}'] > 0.) == (table[f'Estimate{cond2}'] > 0.)) & ((table[f'Estimate{cond2}'] > 0.) == (table[f'Estimate{cond3}'] > 0.))].index)

        # directed version
        sets = tuple(set(metabolites[cond,tissue,c].apply(lambda u: f'{u.name}u' if u['Estimate'] > 0. else f'{u.name}d', axis=1)) for c in comparisons)
        # directionless version
        #sets = tuple(set(metabolites[cond,tissue,c].apply(lambda u: f'{u.name}' if u['Estimate'] > 0. else f'{u.name}', axis=1)) for c in comparisons)

        #https://stackoverflow.com/questions/35320437/drawing-brackets-over-plot

        v = venn3(subsets = sets,
                      set_labels = [' vs \n'.join(c) for c in comparisons.values()], alpha=0.5, ax=ax[i,j])
        if v.patches[0] is not None:
            v.patches[0].set_color(colors[0])
        if v.patches[1] is not None:
            v.patches[1].set_color(colors[1])
        if v.patches[2] is not None:
            v.patches[2].set_color(colors[2])
        if v.patches[3] is not None:
            v.patches[3].set_color(colors[3])
        if v.patches[4] is not None:
            v.patches[4].set_color('grey')
        if v.patches[5] is not None:
            v.patches[5].set_color('grey')
        if v.patches[6] is not None:
            v.patches[6].set_color('silver') # dodgerblue

        if i == 0:
            ax[i,j].set_title(f'$\\bf{{Primary}}$')

    for j,tissue in zip(range(2,8,2),tissues):
        # directed version
        sets = tuple(set(lipids[cond,tissue,c].apply(lambda u: f'{u.name}u' if u['Estimate'] > 0. else f'{u.name}d', axis=1)) for c in comparisons)
        # directionless version
        #sets = tuple(set(lipids[cond,tissue,c].apply(lambda u: f'{u.name}' if u['Estimate'] > 0. else f'{u.name}', axis=1)) for c in comparisons)

        v = venn3(subsets = sets,
                      set_labels = [' vs \n'.join(c) for c in comparisons.values()], alpha=0.5, ax=ax[i,j])
        if v.patches[0] is not None:
            v.patches[0].set_color(colors[0])
        if v.patches[1] is not None:
            v.patches[1].set_color(colors[1])
        if v.patches[2] is not None:
            v.patches[2].set_color(colors[2])
        if v.patches[3] is not None:
            v.patches[3].set_color(colors[3])
        if v.patches[4] is not None:
            v.patches[4].set_color('grey')
        if v.patches[5] is not None:
            v.patches[5].set_color('grey')
        if v.patches[6] is not None:
            v.patches[6].set_color('silver') # dodgerblue

        if i == 0:
            ax[i,j].set_title(f'$\\bf{{Lipids}}$')

ax2 = fig.add_subplot(111)
rect = patches.FancyBboxPatch((0.29,0.065),0.125,0.85,linewidth=0,boxstyle='round,pad=0.0,rounding_size=0.01',edgecolor=None,facecolor='gainsboro',mutation_aspect=1./2.,transform=fig.transFigure)
ax2.add_patch(rect)
rect = patches.FancyBboxPatch((0.29+0.25,0.065),0.125,0.85,linewidth=0,boxstyle='round,pad=0.0,rounding_size=0.01',edgecolor=None,facecolor='gainsboro',mutation_aspect=1./2.,transform=fig.transFigure)
ax2.add_patch(rect)
rect = patches.FancyBboxPatch((0.29+0.25*2,0.065),0.125,0.85,linewidth=0,boxstyle='round,pad=0.0,rounding_size=0.01',edgecolor=None,facecolor='gainsboro',mutation_aspect=1./2.,transform=fig.transFigure)
ax2.add_patch(rect)

ax2.set_yticks([], minor=[])
ax2.set_xticks([], minor=[])
ax2.patch.set_visible(False)
ax2.set_zorder(-1)
for s in ["top", "bottom", "left", "right"]:
    ax2.spines[s].set_visible(False)
plt.axis('off')
plt.margins(0,0)
ax2.set_frame_on(False)
for o in fig.findobj():
    o.set_clip_on(False)

#https://stackoverflow.com/questions/35320437/drawing-brackets-over-plot
plt.annotate('Brain', xy=(0.3-0.005, 0.95), xytext=(0.3-0.005, 0.97), xycoords='figure fraction', textcoords='figure fraction',
            fontsize=14.0*1.5, ha='center', va='bottom',
            bbox=None,
            arrowprops=dict(arrowstyle='-[, widthB=10.0, lengthB=1.5, angleB=0.0', lw=2.0))
plt.annotate('Muscle', xy=(0.3+0.245, 0.95), xytext=(0.3+0.245, 0.97), xycoords='figure fraction', textcoords='figure fraction',
            fontsize=14.0*1.5, ha='center', va='bottom',
            bbox=None,
            arrowprops=dict(arrowstyle='-[, widthB=10.0, lengthB=1.5, angleB=0.0', lw=2.0))
plt.annotate('Liver', xy=(0.3+0.245*2+0.005, 0.95), xytext=(0.3+0.245*2+0.005, 0.97), xycoords='figure fraction', textcoords='figure fraction',
            fontsize=14.0*1.5, ha='center', va='bottom',
            bbox=None,
            arrowprops=dict(arrowstyle='-[, widthB=10.0, lengthB=1.5, angleB=0.0', lw=2.0))


filename=os.path.join(args.output_dir,outlier,f'venn-diagrams-cross-pop.pdf')
plt.savefig(filename)
