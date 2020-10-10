from cavefinomics import AstyanaxMe
from argparse import ArgumentParser
from numpy import log10, array
from itertools import repeat, chain
from seaborn import (catplot, FacetGrid, palplot,
                     color_palette, set_palette, heatmap)
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib.ticker as plticker
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from pandas import concat

SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

# palettes
#https://colorhunt.co/palette/207147
bloodmeridian = ['#931a25','#e97171','#ffcb8e','#f5efef'][::-1]
#https://colorhunt.co/palette/207237
barney = ['#d789d7','#9d65c9','#5d54a4','#2a3d66'][::-1]
#https://colorhunt.co/palette/202154
fireinthesky = ['#ed6663','#b52b65','#3c2c3e','#59405c'][::-1]
#https://colorhunt.co/palette/195946
oceansunset = ['#848ccf','#93b5e1','#ffe4e4','#be5683']
#https://colorhunt.co/palette/192115
oceansunset2 = ['#ddf3f5','#a6dcef','#f2aaaa','#e36387']
#https://colorhunt.co/palette/175167
desertsunset = ['#d63447','#f57b51','#f6eedf','#d1cebd']
#https://colorhunt.co/palette/192164
watermelon = ['#e7305b','#e2979c','#f7f5dd','#9bdeac']
#https://colorhunt.co/palette/13607
watermelon2 = ['#9bb899','#fcceaa','#f5827d','#ea4961']
#https://colorhunt.co/palette/191947
dusk = ['#111d5e','#c70039','#f37121','#ffbd69']
#https://colorhunt.co/palette/178936
dusk2 = ['#202040','#543864','#ff6363','#ffbd69']
#https://colorhunt.co/palette/175358
pinkdew = ['#ffb2a7','#e6739f','#cc0e74','#790c5a']
#https://colorhunt.co/palette/172856
twilight = ['#4f3961','#ea728c','#fc9d9d','#f3d4d4']
#https://colorhunt.co/palette/168310
pastel = ['#ffb6b9','#fae3d9','#bbded6','#8ac6d1'][::-1]
#https://colorhunt.co/palette/167893
bigblue = ['#1b262c','#0f4c75','#3282b8','#bbe1fa']
#https://colorhunt.co/palette/159679
capecanaveral = ['#f0134d','#ff6f5e','#f5f0e3','#40bfc1'][::-1]
#https://colorhunt.co/palette/159617
pastelblue = ['#ecfcff','#b2fcff','#5edfff','#3e64ff']
#https://colorhunt.co/palette/137194
thunderbirds = ['#c82121','#dee1ec','#becbff','#0d0cb5'][::-1]
#https://colorhunt.co/palette/3460
thunderbirds2 = ['#eb586f','#d8e9f0','#4aa0d5','#454553'][::-1]
#https://colorhunt.co/palette/80280
thunderbirds3 = ['#e2eff1','#65799b','#555273','#e23e57']
#https://colorhunt.co/palette/116523
sandysunset = ['#84b9ef','#fbe4c9','#ff5d5d','#952e4b']
#https://colorhunt.co/palette/14114
#brightsandy
#https://colorhunt.co/palette/13313
reddwarf = ['#2a363b','#e84a5f','#ff847c','#fecea8']
#https://colorhunt.co/palette/131292
savannahsunrise = ['#a26ea1','#f18a9b','#ffb480','#ffff9d']
#https://colorhunt.co/palette/125241
savannahsunrise = ['#ff9900','#ca431d','#8b104e','#520556']
#https://colorhunt.co/palette/44449
savannahsunrise = ['#2e94b9','#fffdc1','#f0b775','#fd5959']
#https://colorhunt.co/palette/134419
highnoon = ['#feffdf','#ffe79a','#ffa952','#ef5a5a']

#dusk2
mycmap = LinearSegmentedColormap.from_list('thecmap', thunderbirds2)

#https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
#colors = ['firebrick','goldenrod','mediumseagreen','rosybrown']
#https://colorhunt.co/palette/184190
#colors = ['#0e9aa7','#3da4ab','#f6cd61','#fe8a71']
#https://colorhunt.co/palette/189889
colors = ['#d92027','#ff9234','#ffcd3c','#35d0ba']
#https://colorhunt.co/palette/183823
#colors = ['#ffc4a3','#ff9a76','#f96d80','#bb596b']

#palplot(color_palette('Set1'))
set_palette('Set2')

parser = ArgumentParser(description="Sugar phosphate heatmap.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--level", type=float, help="P value cutoff.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--output", type=str, help="Output.")
parser.add_argument("--output-cbar", type=str, help="Output cbar.")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

compounds = [
    'glucose-1-phosphate',
    'glucose-6-phosphate',
    'galactose-6-phosphate',
    'fructose-1-phosphate',
    'fructose-6-phosphate',
    'ribose-5-phosphate',
    'ribulose-5-phosphate',
    'glucose',
    'fructose',
    'glucuronic acid',
    'gluconic acid',
  ]

sugar_phosphates = [
    'glucose-1-phosphate',
    'glucose-6-phosphate',
    'galactose-6-phosphate',
    'fructose-1-phosphate',
    'fructose-6-phosphate',
    'ribose-5-phosphate',
    'ribulose-5-phosphate',
    ]

sugars = [
    'glucose',
    'fructose',
    ]

uronic_acids = [
    'glucuronic acid',
    'gluconic acid',
    ]

compounds_by_class = [sugar_phosphates, sugars, uronic_acids]

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,6+1),9*3))))
astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)
astyanax_data = astyanax_data.loc[compounds]

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    return subset

astyanax_data = process_outlier(args.exclude_outlier,astyanax_data)

pops = ['Pachon', 'Tinaja', 'Surface']
pops2 = ['Pachón', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
#conditions = ['30d Starved', '4d Starved', 'Refed']
conditions = {'4d':'4d Starved', '30d':'30d Starved', 'Ref':'Refed'}
condmap = {v:k for k,v in {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}.items()}
comparisons = ['PvS','TvS','PvT']
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}

data = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions.values():
            for compound in compounds:
                subset = astyanax_data.loc[compound,
                  astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
                for val in subset:
                    data.append({'Population':pop if pop != 'Pachon' else 'Pachón','Tissue':tissue,'Condition':condmap[condition],'Compound':compound,'Value':val})
data = DataFrame(data)
# normalize to condition
#cond_maxes = data[['Condition','Compound','Value']].groupby(['Condition','Compound']).max()
#data['Value'] = data.apply(lambda u: u['Value']/cond_maxes.loc[u['Condition'],u['Compound']],axis=1)
averaged_data = data.pivot_table(index=['Condition','Compound'],columns=['Population','Tissue'])
averaged_data = averaged_data.max(axis=1)
print(averaged_data)
#stop
data['Value'] = data.apply(lambda u: u['Value']/averaged_data.loc[u['Condition'],u['Compound']],axis=1)
print(data)


# significance data
datasets = []
sig = {}
up = {}
if args.exclude_outlier:
    outlier_text = 'kein-Ausreißern'
else:
    outlier_text = 'mit-Ausreißern'
for cat in categories:
    for tissue in tissues:
        for cond in conditions:
            for comp in comparisons:
                d = read_csv(f'out/work/primary/glm/singlefactor/{outlier_text}/{cat}/{tissue}/{cond}/{comp}.csv',index_col=0)
                d['Category'] = cat
                d['Tissue'] = tissue
                d['Condition'] = cond
                d['Comparison'] = comp
                for m in d.index:
                    if d.loc[m,'Pr(>|z|)'] < args.level:
                        sig[m,tissue,cond,comp] = True
                    else:
                        sig[m,tissue,cond,comp] = False
                    if d.loc[m,'Estimate'] > 0.:
                        up[m,tissue,cond,comp] = True
                    else:
                        up[m,tissue,cond,comp] = False
                datasets.append(d)
sig_data = concat(datasets,axis=0).dropna()
#print(sig_data)


#fig,ax = plt.subplots(nrows=len(set(data['Compound'])),ncols=len(tissues)*len(conditions),figsize=(12.,12.))
gridspec_kw = {"height_ratios":[3.,1.,1.], "width_ratios" : [3.,3.,3.,2.]*2+[3.,3.,3.]}
fig,ax = plt.subplots(nrows=3,ncols=11,figsize=(12.,10.),gridspec_kw=gridspec_kw)

#print(data)
#j = 0
for kc,(cond_label,condition) in enumerate(conditions.items()):
    cond_data = data.loc[(data['Condition'] == cond_label)]
    #print(cond_data)
    #stop
    #cond_data.iloc[2:] = cond_data.div(cond_data.iloc[2:].max(axis=1),axis=0)
    for kt,tissue in enumerate(tissues):
        for i,cpds in enumerate(compounds_by_class):
            data2d = cond_data.loc[(cond_data['Tissue'] == tissue)]
            data2d = data2d.drop('Tissue',axis=1).drop('Condition',axis=1)
            data2d = data2d.pivot_table(index='Compound',columns='Population')
            #https://stackoverflow.com/questions/35678874/normalize-rows-of-pandas-data-frame-by-their-sums/35679163
            #https://stackoverflow.com/questions/39273441/flatten-pandas-pivot-table
            data2d.columns = data2d.columns.to_series().str.join('_').str.replace('Value_','')
            #print(list(data2d.columns))
            data2d = data2d[pops2]
            data2d = data2d.loc[cpds,:]
            #print(data2d)
            #print(len(data2d.index))

            #print(list(sig_data.columns))
            #print(sig_data['Comparison'] == 'PvS')
            sig2d = concat(
              [sig_data.loc[(sig_data['Comparison'] == comp) & (sig_data['Condition'] == cond_label) & (sig_data['Tissue'] == tissue)]
              .loc[data2d.index,'Pr(>|z|)']
              .loc[cpds] for comp in ['PvS','TvS']],
              axis=1)
            sig2d.insert(2,'',0.)
            #print(sig2d)
            #print(data2d)
            #stop
            j = kt*4+kc
            data2d.index.name = ''
            #heatmap(data2d, vmin=0., vmax=1., annot=sig2d, fmt = '.2f', ax=ax[j], cbar=False, xticklabels=j==0, yticklabels=True, cmap='coolwarm')
            heatmap(data2d, vmin=0., vmax=1., annot=None, fmt = '.2f', ax=ax[i,j], cbar=False, xticklabels=i==2, yticklabels=j==0, cmap=mycmap)
            #https://stackoverflow.com/questions/26540035/rotate-label-text-in-seaborn-factorplot/34722235
            ax[i,j].set_yticklabels(ax[i,j].get_yticklabels(),fontsize='x-large',rotation=0)
            ax[i,j].set_xticklabels(ax[i,j].get_xticklabels(),fontsize='large')
            if i==0:
                ax[i,j].set_title(f'{cond_label}',fontsize='x-large')
            #break

for i in range(3):
    for j in [3,7]:
        ax[i,j].set_yticks([], minor=[])
        ax[i,j].set_xticks([], minor=[])
        ax[i,j].patch.set_visible(False)
        for s in ["top", "bottom", "left", "right"]:
            ax[i,j].spines[s].set_visible(False)

#https://stackoverflow.com/questions/27037241/changing-the-rotation-of-tick-labels-in-seaborn-heatmap
#plt.yticks(rotation=0)
for x,tissue in zip([0.2325,0.515,0.7925],tissues):
    fig.text(x,0.95,tissue,size=20,ha='center',va='center',fontweight='bold',transform=fig.transFigure)
for y,cls in zip([0.675,0.3575,0.1715],['Sugar\nphos-\nphates','Sugars','Uronic\nacids']):
    fig.text(0.915,y,cls,size=20,ha='left',va='center',rotation=0.,fontweight='bold',transform=fig.transFigure)

ax2 = fig.add_subplot(111)
rect = patches.FancyBboxPatch((0.38,0.0),0.275,0.985,linewidth=0,boxstyle='round,pad=0.0,rounding_size=0.01',edgecolor=None,facecolor='gainsboro',mutation_aspect=1./2.,transform=fig.transFigure)
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

plt.savefig(args.output,bbox_inches='tight',transparent=True,pad_inches=0)



fig,ax = plt.subplots(figsize=(2, 10))
# draw c bar
norm = matplotlib.colors.Normalize(vmin=-1., vmax=1.)
sm = matplotlib.cm.ScalarMappable(cmap=mycmap, norm=norm)
sm.set_array([])
cax = fig.add_axes([0.1,0.1,0.25,0.8])
fig.colorbar(sm, cax=cax)
cax.yaxis.set_major_locator(plt.FixedLocator([-1., 1.]))
cax.set_yticklabels(['0','Max/\ngroup'],fontsize='xx-large',fontweight='bold')
#cax.text(0.5,1.125,'Max of mTIC\npeak avg.\nwithin group\n(6 reps.)',size=20,ha='center',va='center',fontweight='bold',transform=ax.transAxes)
#cax.text(0.5,-0.1,'Down in \ncave+starved',size=20,ha='center',va='center',fontweight='bold',transform=ax.transAxes)
# clear extra axes
# hide graphics
ax.set_yticks([], minor=[])
ax.set_xticks([], minor=[])
ax.patch.set_visible(False)
for s in ["top", "bottom", "left", "right"]:
    ax.spines[s].set_visible(False)
filename=f'{args.output_cbar}'
plt.savefig(filename,bbox_inches='tight',transparent=True,pad_inches=0)
