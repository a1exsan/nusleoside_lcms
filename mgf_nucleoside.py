
import pandas as pd
import msvis
import numpy as np
from tqdm import tqdm
import pickle
import pyopenms as poms


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)


def getAverageMolecularMass(formula):
    return poms.EmpiricalFormula(formula).getAverageWeight()

class MSspecData():
    def __init__(self):
        super().__init__()
        self.specTab = []
        self.index = -1
        self.tol = 0.1

    def __call__(self, index):
        return self.specTab[index]

    def __iter__(self):
        return self

    def __next__(self):
        self.index += 1
        if self.index < len(self.specTab):
            return self.specTab[self.index]
        raise StopIteration

    def add_spec(self, spec):
        for i in spec:
            self.specTab.append(i)

    def spec_uniq_union(self, ms2index, tol=0.1):
        ret = {}
        data = []
        for t in self.specTab:
            if t['ms2index'] == ms2index:
                data.append(t)

        for k in list(data[0].keys()):
            ret[k] = data[0][k]

        initMZ = [i for i in data[0]['ms2']]
        initInt = [i for i in data[0]['intens']]
        for t in data:
            for j in range(len(t['ms2'])):
                ctrl = True
                for i in range(len(initMZ)):
                    d = abs(initMZ[i] - t['ms2'][j])
                    if d <= tol:
                        ctrl = False
                        if t['intens'][j] > initInt[i]:
                            initInt[i] = t['intens'][j]
                if ctrl:
                    initMZ.append(t['ms2'][j])
                    initInt.append(t['intens'][j])

        ret['ms2'] = [i for i in initMZ]
        ret['intens'] = [i for i in initInt]

        return ret

    def spec_uniq_union_data(self, data, tol=0.1):
        ret = {}

        for k in list(data[0].keys()):
            ret[k] = data[0][k]

        initMZ = [i for i in data[0]['ms2']]
        initInt = [i for i in data[0]['intens']]
        rt = []
        for t in data:
            rt.append(t['rt'])
            for j in range(len(t['ms2'])):
                ctrl = True
                for i in range(len(initMZ)):
                    d = abs(initMZ[i] - t['ms2'][j])
                    if d <= tol:
                        ctrl = False
                        if t['intens'][j] > initInt[i]:
                            initInt[i] = t['intens'][j]
                if ctrl:
                    initMZ.append(t['ms2'][j])
                    initInt.append(t['intens'][j])

        ret['ms2'] = [i for i in initMZ]
        ret['intens'] = [i for i in initInt]
        ret['rt-min'] = min(rt)
        ret['rt-max'] = max(rt)

        return ret

    def spec_uniq_union_data_v2(self, data, tol=0.1):
        ret = {}

        for k in list(data[0].keys()):
            ret[k] = data[0][k]

        mz_list, int_list, rt = [], [], []
        for d in data:
            # print(d)
            rt.append(d['rt'])
            for i in range(len(d['ms2'])):
                mz_list.append(round(d['ms2'][i], 1))
                int_list.append(d['intens'][i])

        df = pd.DataFrame({'index': mz_list, 'ms2': mz_list, 'intens': int_list})
        df = df.groupby('index').agg(lambda x: max(list(x)))

        ret['ms2'] = list(df['ms2'])
        ret['intens'] = list(df['intens'])
        ret['rt-min'] = min(rt)
        ret['rt-max'] = max(rt)

        return ret

    def group_by_mz_charge(self, intens_treshold=100, tol=0.1):

        def treshold_by_intens(ms2, intens, treshold):
            rms2, rint = [], []
            for i in range(len(ms2)):
                if intens[i] >= treshold:
                    rint.append(intens[i])
                    rms2.append(ms2[i])
            return rms2, rint

        df = pd.DataFrame(self.specTab)
        df = df.sort_values(by='rt', ascending=True)

        df = df.groupby('mzcharge').agg(lambda x: list(x))

        df = list(df.T.to_dict().values())

        # print(df[211])

        ret = []
        index = 0
        for t in df:
            # print(index, t['charge'])
            data = []
            for i in range(len(t['mz'])):
                d = {}
                d['mz'] = t['mz'][i]
                d['charge'] = t['charge'][i]
                d['ms2'], d['intens'] = treshold_by_intens(t['ms2'][i], t['intens'][i], intens_treshold)
                d['rt'] = t['rt'][i]
                d['index'] = index
                data.append(d)
            index += 1
            group = self.spec_uniq_union_data_v2(data, tol=tol)
            # print(group)
            if len(group['ms2']) > 0:
                ret.append(group)

        # print(ret[0])
        return ret

    def loadMGF(self, fn, seqCtrl=False):

        f = open(fn, 'r')
        ln = f.readlines()
        f.close()

        v = {'ms2': []}
        k = 0
        ctrl = False
        ms2, intens = [], []
        for i in tqdm(ln, desc='loading ms/ms data: '):
            #print(i)
            if i.find('TITLE') > -1:
                v['title'] = i[6:i.find('\n')]
                # print(v['title'])
                p1 = v['title'].find('at ') + 2
                p2 = v['title'].find('mins')
                if (p1 > -1) and (p2 > -1):
                    v['rt'] = float(v['title'][p1:p2])
                else:
                    v['rt'] = 0

                if seqCtrl:
                    v['seq'] = i[6:i.find('\n')]
            if i.find('SEQ') > -1:
                v['seq'] = i[4:i.find('\n')]
            if i.find('PEPMASS') > -1:
                v['mz'] = i[8:i.find('\n')]
            if i.find('CHARGE') > -1:
                v['charge'] = i[7:8]
            else:
                v['charge'] = 1

                try:
                    v['mzcharge'] = str(v['mz']) + '_' + str(v['charge'])
                except:
                    v['mzcharge'] = ''

            if i.find('END IONS') > -1:
                ctrl = False
                if len(ms2) > 0:
                    v['ms2'] = ms2
                    v['intens'] = intens
                    ms2, intens = [], []
                    self.specTab.append(v)
                    v = {'ms2': []}

            if i.find('RTINSECONDS') > -1:
                ctrl = True
            elif ctrl:
                pos = 0
                if i.find(' ') > -1:
                    pos = i.find(' ')
                if i.find('\t') > -1:
                    pos = i.find('\t')
                mz = i[0:pos]
                ii = i[pos + 1:i.find('\n')]
                ms2.append(float(mz))
                intens.append(float(ii))

            #if (k + 1) == len(ln):
            #    self.specTab.append(v)
            k += 1

    def write_mgf(self, fn, tab):
        f = open(fn, 'w')
        ln = []
        ln.append('COM=convert from MSP NIST ' + fn + '\n')
        k = 0
        for t in tab:
            ln.append('BEGIN IONS\n')
            ln.append('TITLE ' + t['title'] + '\n')
            # ln.append('SEQ='+t['title']+'\n')

            pm = 'PEPMASS=' + t['mz'] + '\n'
            ln.append(pm)

            ch = 'CHARGE=' + t['charge'] + '+\n'
            ln.append(ch)

            ln.append('INSTRUMENT=Default\n')
            rt = 'RTINSECONDS=' + str(t['rt']) + '\n'
            ln.append(rt)

            for j in range(len(t['ms2'])):
                ions = str(t['ms2'][j]) + ' ' + str(t['intens'][j]) + '\n'
                ln.append(ions)

            ln.append('END IONS\n\n')
        f.writelines(ln)

    def PMZ_filter(self, tab, key_seq='seq'):
        ret = []

        for item in tab:
            peptide = poms.AASequence.fromString(item[key_seq])
            mz = peptide.getMonoWeight(poms.Residue.ResidueType.Full, int(item['charge'])) / int(item['charge'])
            if abs(mz - float(item['mz'])) < self.tol:
                ret.append(item)

        return ret

class mgf2DMapConverter():
    def __init__(self, ms2Tab):
        self.ms2Tab = ms2Tab
        self.mzData = None
        self.max_mz = 350
        self.min_mz = 50
        self.ms2Data = None

    def convert(self):
        self.mzData = []
        for i in range(len(self.ms2Tab)):
            d = {}
            d['rt'] = float(self.ms2Tab[i]['rt'])
            d['mz'] = round(float(self.ms2Tab[i]['mz']), 1)
            d['intens'] = max(self.ms2Tab[i]['intens'])
            d['uid'] = i
            self.mzData.append(d)
        self.mzData = pd.DataFrame(self.mzData)
        return self.mzData

    def mzFiltrate(self, min_mz, max_mz):
        self.mzData = self.mzData[self.mzData['mz'] >= min_mz]
        self.mzData = self.mzData[self.mzData['mz'] <= max_mz]

    def rtFiltrate(self, min_rt, max_rt):
        self.mzData = self.mzData[self.mzData['rt'] >= min_rt]
        self.mzData = self.mzData[self.mzData['rt'] <= max_rt]

    def intFiltrate(self, min_int, max_int):
        if max_int == -1:
            max_int = self.mzData['intens'].max()
        self.mzData = self.mzData[self.mzData['intens'] >= min_int]
        self.mzData = self.mzData[self.mzData['intens'] <= max_int]

    def group_by_mz_rt(self):
        mz_list = list(set(self.mzData['mz']))
        self.mzData['class'] = [str(0) for i in range(self.mzData.shape[0])]

        index = 1
        for mz in mz_list:
            ind_array = self.mzData[self.mzData['mz'] == mz].index
            self.mzData.loc[ind_array[0], 'class'] = str(index)
            for i in range(1, ind_array.shape[0]):
                if abs(self.mzData['rt'].loc[ind_array[i]] - self.mzData['rt'].loc[ind_array[i - 1]]) < 0.1:
                    self.mzData.loc[ind_array[i], 'class'] = str(index)
                else:
                    index += 1
                    self.mzData.loc[ind_array[i], 'class'] = str(index)
            index += 1

        self.mzData['name'] = self.mzData['class']

    def __ms2ToVec(self, ms2, intens, min_mz=50, max_mz=300):
        ret = [0. for mz in range(max_mz + 1)]
        for i, mz in zip(intens, ms2):
            if mz <= max_mz and mz >= min_mz:
                ret[int(round(mz, 0))] = i
        return np.array(ret)

    def average_ms2(self):
        classes = list(set(self.mzData['class']))
        self.ms2Data = {'mz': [], 'rt': [], 'name': [],
                        'ms2': [], 'intens': [], 'class': [],
                        'max_int': [], 'point_count': []}
        for c in tqdm(classes, desc='average ms2:'):
            spec = self.__ms2ToVec([], [], self.min_mz, self.max_mz)
            count = 0
            for uid in self.mzData[self.mzData['class'] == c]['uid']:
                spec += self.__ms2ToVec(self.ms2Tab[uid]['ms2'], self.ms2Tab[uid]['intens'], self.min_mz, self.max_mz)
                count += 1
            #spec = spec / count
            self.ms2Data['ms2'].append([i for i in range(spec.shape[0]) if spec[i] != 0])
            self.ms2Data['intens'].append([spec[i] for i in range(spec.shape[0]) if spec[i] != 0])
            self.ms2Data['class'].append(c)
            self.ms2Data['max_int'].append(max(self.ms2Data['intens'][-1]))
            self.ms2Data['point_count'].append(count)
            self.ms2Data['mz'].append(self.mzData[self.mzData['class'] == c]['mz'].max())
            self.ms2Data['rt'].append(self.mzData[self.mzData['class'] == c]['rt'].mean())
            self.ms2Data['name'].append(c)

        self.ms2Data = pd.DataFrame(self.ms2Data)

    def reset_intensity(self):
        for i, c in zip(self.ms2Data['max_int'], self.ms2Data['class']):
            self.mzData.loc[self.mzData['class'] == c, 'intens'] = i

    def identification_by_modomix(self, base):

        for intens, c, ms2, mz in zip(self.ms2Data['intens'],
                                        self.ms2Data['class'],
                                        self.ms2Data['ms2'],
                                        self.ms2Data['mz']):
            #max_index = intens.index(max(intens))
            #prod_ion = ms2[max_index]

            df = pd.DataFrame({'ms2': ms2, 'intens': intens})
            df = df.sort_values(by='intens', ascending=False)
            df = df[:3]
            prod_ion = list(df['ms2'])

            result = base.search_simple(mz, prod_ion, 0.5, 1)

            names = ''
            if result.shape[0] > 0:
                for name in result['short_name']:
                    names += name + '/'
            names = names + c

            if names != '':
                self.ms2Data.loc[self.ms2Data['class'] == c, 'name'] = names
                self.mzData.loc[self.mzData['class'] == c, 'name'] = names




class modomix_database():
    def __init__(self, db_filename='data/nucleoside_modif_base_dict.pkl'):
        self.db_filename = db_filename

        with open(self.db_filename, 'rb') as f:
            self.base = pickle.load(f)

        self.base = pd.DataFrame(self.base)
        self.base = self.base[self.base['mass_avg'] != 'null']
        self.base = self.base[self.base['mass_avg'] != '0.0']
        self.base = self.base.reset_index()
        self.base = self.base.drop(['index'], axis=1)
        self.base['mass_avg'] = [float(m) for m in self.base['mass_avg']]
        self.base['mass_prot'] = self.base['mass_avg'] + 1

    def search_simple(self, prot_mass, product_ion, treshold1=0.5, treshold2=1):

        df = self.base[self.base['mass_prot'] >= prot_mass - treshold1]
        df = df[df['mass_prot'] <= prot_mass + treshold1]

        ret = []
        for index in df.index:
            if df['product_ions'].loc[index].find('null') == -1:
                ions = [float(i) for i in df['product_ions'].loc[index].split("/")]
                ctrl = False
                for i in ions:
                    for p in product_ion:
                        if abs(i - p) <= treshold2:
                            ctrl = True
                            break
                if ctrl:
                    ret.append(dict(df.loc[index]))

        return pd.DataFrame(ret)





def main():
    path = r'C:\Users\Alex\Documents\LCMS\nucleoside\eGFP_JB.mgf'
    path = r'C:\Users\Alex\Documents\LCMS\nucleoside\xRNA\xRNA_20cid_1_1.mgf'
    path = r'C:\Users\Alex\Documents\LCMS\nucleoside\xRNA\xRNA_10cid.mgf'

    path = r'C:\Users\Alex\Documents\LCMS\nucleoside\trna_ecoli.mgf'
    #path = r'C:\Users\Alex\Documents\LCMS\nucleoside\trna_yeast.mgf'

    #path = r'C:\Users\Alex\Documents\LCMS\nucleoside\rna\v7_6ul.mgf'
    path = r'C:\Users\Alex\Documents\LCMS\nucleoside\rna\yRNA_221221.mgf'

    mgf = MSspecData()
    mgf.loadMGF(path)

    conv = mgf2DMapConverter(mgf.specTab)
    conv.convert()
    conv.mzFiltrate(200, 400)
    conv.intFiltrate(100, -1)
    conv.group_by_mz_rt()
    conv.average_ms2()
    mapData = conv.mzData

    #index = conv.ms2Data[conv.ms2Data['class'] == 'sel_class'].index
    #d = dict(conv.ms2Data[conv.ms2Data['class'] == sel_class].loc[index[0]])
    #ms2, intens = d['ms2'], d['intens']

    conv.identification_by_modomix(modomix_database())
    print(conv.ms2Data['name'])

    #spec_vis = msvis.massSpec(ms2, intens)
    #spec_vis.draw_spec()

    print(mapData)

    #viewer = msvis.simple_ms_map(mapData['rt'], mapData['mz'], mapData['intens'])
    viewer = msvis.clusters_ms_map(mapData['rt'], mapData['mz'], mapData['intens'], mapData['class'])
    viewer.transperancy = 0.1
    viewer.draw_map()

def main1():
    base = modomix_database()
    print(base.base)
    print(base.search_simple(284.1, 152, 0.5))

if __name__ == '__main__':
    main()