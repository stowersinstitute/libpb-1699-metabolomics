from pandas import read_csv, DataFrame, concat, Series
import attr
import os
from json import load
from numpy import log10, isfinite, any, all, quantile, nan, isnan
from numbers import Number

tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']

@attr.s(auto_attribs=True)
class AstyanaxLi:
    lipids_normalized: str
    lipidmaps_js: str
    lipids_unnormalized: str = None

    def __attrs_post_init__(self):
        if self.lipids_unnormalized is not None:
            # unnormalized
            self.unnormalized = {}
            self.unnormalized_polarity = {}
            for tissue in tissues:
                filename = os.path.join(self.lipids_unnormalized,f'{tissue.lower()}.csv')
                column_table = (
                    read_csv(filename, skiprows=0, nrows=6, index_col=6, header=0)
                    .iloc[:, 6:]
                    .transpose()
                    )
                #print(column_table)

                row_table = read_csv(filename, skiprows=6, usecols=list(range(6)), index_col=0, header=0)
                #print(row_table)

                data = read_csv(filename, skiprows=7, header=None).iloc[:, 7:]
                data.insert(0,'InChIKey',list(row_table['InChI Key']))
                data = data.set_index('InChIKey')
                data.columns = list(column_table['Treatment'])
                print(f'data before {tissue} {data.shape}')
                if tissue != 'Muscle':
                    t = (row_table.loc[list(data.index.notna())]['ESI mode'] == 'ESI (+)').shape
                    #print(f't {t}')
                    self.unnormalized_polarity[tissue] = (row_table.loc[list(data.index.notna())]['ESI mode'] == 'ESI (+)')
                    #print(f'the {self.unnormalized_polarity[tissue].shape}')
                else:
                    self.unnormalized_polarity[tissue] = (row_table.loc[list(data.index.notna())]['ESI (+)'] == 'ESI (+)')
                print(f'un table {tissue} {self.unnormalized_polarity[tissue].shape}')
                data = data.loc[data.index.notna()]
                #print(f'data 1 {tissue} {data.shape}')
                data = data.loc[:,data.columns.notna()]
                self.unnormalized[tissue] = data
                #print(f'unnormalized data {tissue} {self.unnormalized[tissue].shape}')

        # normalized
        self.normalized = {}
        for tissue in tissues:
            for polarity in polarities:
                filename = os.path.join(self.lipids_normalized,polarity,f'{tissue.lower()}.csv')
                column_table = (
                    read_csv(filename, skiprows=0, nrows=6, index_col=5, header=0)
                    .iloc[:, 5:]
                    .transpose()
                    )

                row_table = read_csv(filename, skiprows=6, usecols=list(range(5)), index_col=0, header=0)
                #if tissue == 'Brain' and polarity == 'negative':
                    #print(print(row_table.loc[row_table['InChI Key'] == 'JBDGKEXQKCCQFK-JWQIMADESA-N',:]))
                    #stop

                data = read_csv(filename, skiprows=7, header=None).iloc[:, 6:]
                data.insert(0,'InChIKey',list(row_table['InChI Key']))
                data = data.set_index('InChIKey')
                data.columns = list(column_table['Treatment'])
                data = data.loc[data.index.notna()]
                data = data.loc[:,data.columns.notna()]
                #if tissue == 'Brain' and polarity == 'negative':
                    #print(print(data.loc['JBDGKEXQKCCQFK-JWQIMADESA-N',:]))
                    #stop
                #data = data.drop_duplicates()
                #https://stackoverflow.com/questions/13035764/remove-rows-with-duplicate-indices-pandas-dataframe-and-timeseries
                #data = data.loc[~data.index.duplicated(keep='first')]
                data = data.groupby(data.index).sum()
                self.normalized[tissue,polarity] = data

        with open(self.lipidmaps_js) as f:
            self.lipidmaps = load(f)

        self.lipidmaps_inchikey = {v['INCHI_KEY']:u for u,v in self.lipidmaps.items()}
        self.lipidmaps_inchikey2 = {'-'.join(v['INCHI_KEY'].split('-')[:2]):u for u,v in self.lipidmaps.items()}
        self.lipidmaps_inchikey1 = {v['INCHI_KEY'].split('-')[0]:u for u,v in self.lipidmaps.items()}

        def get_lipidmap_id(inchikey):
            if inchikey in self.lipidmaps_inchikey:
                return self.lipidmaps_inchikey[inchikey]
            elif '-'.join(inchikey.split('-')[:2]) in self.lipidmaps_inchikey2:
                return self.lipidmaps_inchikey2['-'.join(inchikey.split('-')[:2])]
            elif inchikey.split('-')[0] in self.lipidmaps_inchikey1:
                return self.lipidmaps_inchikey1[inchikey.split('-')[0]]
            return nan

        def fix_category(c):
            if '[' in c:
                return c.split('[')[0].strip()
            else:
                return c

        tissuedata = {}
        for tissue in tissues:
            d = []
            for polarity in polarities:
                lmd = DataFrame(self.normalized[tissue,polarity].apply(lambda u: get_lipidmap_id(u.name),axis=1), columns=['LMID'])
                lmd['Category'] = lmd.apply(lambda u: fix_category(self.lipidmaps[u[0]]['CATEGORY']) if isinstance(u[0],str) else nan,axis=1)
                lmd['Class'] = lmd.apply(lambda u: fix_category(self.lipidmaps[u[0]]['MAIN_CLASS']) if isinstance(u[0],str) else nan,axis=1)
                lmd['Tissue'] = tissue
                lmd['Polarity'] = polarity
                lmd['InChIKey'] = self.normalized[tissue,polarity].index
                lmd = lmd.set_index("InChIKey")
                d.append(lmd)

            #lmids = (set(d[0]['LMID']) | set(d[1]['LMID'])) - set([nan])
            inchikeys = set(set(d[0].index) | set(d[1].index)) - set([nan])
            lmdata = []
            from itertools import repeat, chain
            n = 0
            for inchikey in inchikeys:
                try:
                    #print(self.normalized[tissue,polarities[0]].columns)
                    if inchikey in d[0].index:
                        k = 0
                    elif inchikey in d[1].index:
                        k = 1
                    else:
                        assert False
                    r = {'LMID':d[k]['LMID'].loc[inchikey], 'Category': d[k]['Category'].loc[inchikey], 'Class': d[k]['Class'].loc[inchikey], 'Tissue': d[k]['Tissue'].loc[inchikey], 'InChIKey': inchikey}
                    if isinstance(r['LMID'],Number) and isnan(r['LMID']):
                        continue
                    n += 1
                    self.non_numeric_cols= len(r)
                    r.update({f'{c} {tissue} {k}': v1+v2 for k,c,v1,v2 in zip(
                      chain.from_iterable(repeat(range(1,1+6),9)),
                      self.normalized[tissue,polarities[0]].columns,
                      self.normalized[tissue,polarities[0]].loc[inchikey,:] if inchikey in self.normalized[tissue,polarities[0]].index else [0.]*6*9,
                      self.normalized[tissue,polarities[1]].loc[inchikey,:] if inchikey in self.normalized[tissue,polarities[1]].index else [0.]*6*9,
                      )})
                    lmdata.append(r)
                    #print(n,len(lmdata))
                    #print(n,lmdata)
                    assert n == len(lmdata)
                except KeyError:
                    pass
            #print(f'l {len(lmdata)}')
            #print(f'n {n}')
            lmdata = DataFrame(lmdata)
            #lmdata = lmdata.set_index('LMID')
            lmdata = lmdata.set_index('InChIKey')
            #lmdata = lmdata.groupby('LMID').sum()
            tissuedata[tissue] = lmdata
            #print(lmdata)
            #print(lmdata.iloc[:,self.non_numeric_cols-1:])
            #print(len(set(lmdata.index)))
        lmdata = concat(tissuedata.values(), axis=1)
        #lmdata = lmdata.drop_duplicates(axis=1)
        lmdata = lmdata.loc[:,~lmdata.columns.duplicated()]
        lmdata = lmdata.dropna()
        #lmdata = lmdata.fillna(0.)
        self.lmdata = lmdata
        #print(lmdata)
        #print(len(set(tissuedata['Liver'].index) & set(tissuedata['Brain'].index)))
        #print(lmdata.loc['LMFA01030132'])
