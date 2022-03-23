import pandas as pd
import streamlit as st

import msvis
import mgf_nucleoside as mn

st.set_page_config(layout="wide")

st.markdown("""
<style>
.big-font {
    font-size: 20px !important;
    font-weight: bold;
}
</style>
""", unsafe_allow_html=True)

@st.cache
def upload_mgf(path):
    mgf = mn.MSspecData()
    mgf.loadMGF(path)
    return mgf

@st.cache
def load_modomix_base():
    return mn.modomix_database()

@st.cache
def load_pubchem_base():
    return pd.read_csv('data/nuclesides_pubchem_group_inchi_mass.csv')

def find_mass(dataset, msData, sel_class, treshold=1):
    index = msData.ms2Data[msData.ms2Data['class'] == sel_class].index
    d = dict(msData.ms2Data[msData.ms2Data['class'] == sel_class].loc[index[0]])
    massH = d['mz']

    df = dataset[dataset['M+H'] <= massH + treshold]
    df = df[df['M+H'] >= massH - treshold]
    base = 'https://pubchem.ncbi.nlm.nih.gov/#query='

    inchi = []
    for i in df['INCHI']:
        if i.find('InChI', 2) != -1:
            inchi.append(i[0: i.find('InChI', 2)])
        else:
            inchi.append(i)

    df['link'] = [f"{base}{i}" for i in inchi]
    out = pd.DataFrame({'chemical': df['chemical'], 'M+H': df['M+H'], 'link': df['link']})
    return out

def get_msData(mgf):
    conv = mn.mgf2DMapConverter(mgf.specTab)
    conv.convert()
    return conv

def draw_msMap(msData):
    viewer = msvis.bokeh_ms_map_class(msData.mzData['rt'],
                                      msData.mzData['mz'],
                                      msData.mzData['intens'],
                                      msData.mzData['class'],
                                      msData.mzData['name'],
                                    rt_position=0, title='MS:Time Map')
    viewer.xaxis_label = 'Retention time, min'
    viewer.transperancy = 0.5
    viewer.draw_map(is_show=False)
    st.bokeh_chart(viewer.plot, use_container_width=True)

def draw_ms_spectra(msData, sel_class):

    index = msData.ms2Data[msData.ms2Data['class'] == sel_class].index
    d = dict(msData.ms2Data[msData.ms2Data['class'] == sel_class].loc[index[0]])
    ms2, intens = d['ms2'], d['intens']

    viewer = msvis.bokeh_ms_spectra_simple(ms2, intens, title='MS Spectrum')
    viewer.draw_map(is_show=False)
    st.bokeh_chart(viewer.plot, use_container_width=True)

def find_classes(msData):
    msData.group_by_mz_rt()
    msData.average_ms2()
    #msData.reset_intensity()
    return msData

def filtrate_data_1(msData, min_mz, max_mz, int_treshold, min_rt, max_rt):
    msData.rtFiltrate(min_rt, max_rt)
    msData.mzFiltrate(min_mz, max_mz)
    msData.intFiltrate(int_treshold, -1)
    return msData

def get_class_info(msData, selected_class):
    index = msData.ms2Data[msData.ms2Data['class'] == selected_class].index
    results = dict(msData.ms2Data[msData.ms2Data['class'] == selected_class].loc[index[0]])

    df = pd.DataFrame({'ms2': results['ms2'], 'intens': results['intens']})
    df = df.sort_values(by='intens', ascending=False)

    results['top5'] = df[:5]
    return results

def filtrate_data_2(msData, low_points, low_intens):

    msData.ms2Data = msData.ms2Data[msData.ms2Data['point_count'] > low_points]
    msData.ms2Data = msData.ms2Data[msData.ms2Data['max_int'] > low_intens]

    in_list = list(msData.ms2Data['class'])

    msData.mzData = msData.mzData[msData.mzData['class'].isin(in_list)]

    return msData

class openfile():
    name = None


st.sidebar.markdown('<p class="big-font">Nucleoside modification analysis tool</p>', unsafe_allow_html=True)
st.sidebar.write("Beta version")
st.sidebar.markdown('<p class="big-font">Upload mass spectrometry data:</p>', unsafe_allow_html=True)


msData = None
selected_class = None

uploaded_file = st.sidebar.file_uploader("Choose a file (*.mgf)")

demo_data = st.sidebar.radio('Demo data', ("E.coli tRNA data", "Yeast tRNA data"))

if demo_data == "E.coli tRNA data":
    mgfData = upload_mgf(f'data/trna_ecoli.mgf')
    msData = get_msData(mgfData)

elif demo_data == "Yeast tRNA data":
    mgfData = upload_mgf(f'data/trna_yeast.mgf')
    msData = get_msData(mgfData)


pubchem_db = load_pubchem_base()
modomix_db = load_modomix_base()

if uploaded_file is not None:
    with open(f'data/temp/{uploaded_file.name}', 'wb') as f:
        f.write(uploaded_file.getvalue())
    mgfData = upload_mgf(f'data/temp/{uploaded_file.name}')

    msData = get_msData(mgfData)


col1, col2 = st.columns([2, 1])

with col2:
    if msData != None:

        min_rt, max_rt = st.select_slider(
            'Select a range of mass / charge',
            options=range(0, 21, 1),
            value=(0, 15))
        #st.write('Range of time', min_rt, '--', max_rt)

        min_mz, max_mz = st.select_slider(
        'Select a range of mass / charge',
        options=range(50, 1610, 10),
        value=(200, 400))
        #st.write('Range of mass / charge', min_mz, '--', max_mz)

        intens_treshold = 500

        class_points_tr = st.slider(
            'Select a class points number low treshold:',
            0, 50, 0)

        class_intens_tr = st.slider(
            'Select a class Intensity low treshold:',
            0, 10000, 100)

        msData = filtrate_data_1(msData, min_mz, max_mz, intens_treshold, min_rt, max_rt)

        msData = find_classes(msData)

        msData = filtrate_data_2(msData, class_points_tr, class_intens_tr)

        st.markdown('<p class="big-font">Select one of the class: </p>', unsafe_allow_html=True)
        selected_class = st.selectbox(
            '',
            tuple(msData.ms2Data['class']))

        if st.checkbox('modomix identification'):
            msData.identification_by_modomix(modomix_db)

        url = "https://iimcb.genesilico.pl/modomics/"
        st.write("Modomix [link](%s)" % url)

        class_info = get_class_info(msData, selected_class)

        st.markdown('<p class="big-font">Class information: </p>', unsafe_allow_html=True)

        st.write("Points number: ", class_info['point_count'],
                 " mass/charge ", class_info['mz'],
                 " Retentin time: ", round(class_info['rt'], 2))
        st.write("Intensity: ", round(class_info['max_int'], 0),
                 "Class name: ", class_info['name'])
        st.write("Top 5 product Ions:", class_info['top5'])

        mFormula = st.text_input("Enter Molecular Formula:", 'C5N5H6')
        try:
            st.write('Molecular weight: ', round(mn.getAverageMolecularMass(mFormula), 2))
        except:
            st.write('failed formula')



with col1:
    if msData != None:
        draw_msMap(msData)

    if selected_class != None:
        draw_ms_spectra(msData, selected_class)

    if msData != None:
        st.markdown('<p class="big-font">Possible candidates (PubChem):</p>', unsafe_allow_html=True)
        pubchem_df = find_mass(pubchem_db, msData, selected_class, treshold=0.5)

        if pubchem_df.shape[0] > 0:
            index = 1
            for n, m, url in zip(pubchem_df['chemical'], pubchem_df['M+H'], pubchem_df['link']):
                st.markdown(f"{index}: {n} [M+H] = {m}; [link](%s)" % url)
                index += 1
