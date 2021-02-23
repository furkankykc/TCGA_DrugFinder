from Bio import Entrez
import re
import pickle
import numpy as np
from matplotlib import pyplot as plt
import subprocess
import pandas as pd


class PubChemCompound:

    def __init__(self, mail="furkanfbr@gmail.com", save_file='save.p'):
        Entrez.email = mail
        self.save_file = save_file
        self.known_drugs = self.load()

    def fetch_data(self, compound, max_drug: int = 10):
        idList = self.search_cid_from_compound(compound, max_drug)
        meshIds = self.get_mesh_term_id(idList)
        results = self.get_mesh_term_from_id(meshIds)
        self.known_drugs['c_id'] = dict(self.known_drugs['c_id'].items() | dict(zip(compound, idList)).items())
        self.known_drugs['term'] = dict(self.known_drugs['term'].items() | dict(zip(compound, results)).items())

    def load(self):
        import os
        if not os.path.exists(self.save_file + 'd'):
            self.known_drugs = {'c_id': {}, 'term': {}}
            self.save()
        return {'c_id': pickle.load(open(self.save_file + 'c', "rb")),
                'term': pickle.load(open(self.save_file + 'd', "rb"))}

    def save(self):
        # self.known_drugs['c_id'] = list(map(lambda x: x.strip(), self.known_drugs['c_id']))
        # self.known_drugs['term'] = list(map(lambda x: x.strip(), self.known_drugs['term']))
        pickle.dump(self.known_drugs['c_id'], open(self.save_file + 'c', "wb"),
                    protocol=pickle.HIGHEST_PROTOCOL)  # save it into a file named save.p
        pickle.dump(self.known_drugs['term'], open(self.save_file + 'd', "wb"),
                    protocol=pickle.HIGHEST_PROTOCOL)  # save it into a file named save.p

    def get_pick(self, compound: list):
        compound = list(map(lambda x: x.lower(), compound))
        temp = [i for i in set(compound) if not i in self.known_drugs['c_id'].keys()]
        print(temp)
        self.fetch_data(temp, 1)
        self.save()
        return [self.known_drugs['c_id'][i] if i in self.known_drugs['c_id'] else np.nan for i in compound], [
            self.known_drugs['term'][i] if i in self.known_drugs['term'] else np.nan for i in compound]

    def extract_cid_n_term(self, data: dict):
        return list(data.keys()), list(map(lambda x: x.strip(), list(list(data.values()).keys())))

    def search_cid_from_compound(self, compound, max_drug: int = 10):

        idList = []
        for c in compound:
            handle = Entrez.esearch(db='pccompound', term=f'"{c}"[Synonym]', rettype='xml', retMax=max_drug,
                                    sort='relevance')
            res = Entrez.read(handle)
            idList.extend(map(str, res['IdList']))
        return idList

    def get_mesh_term_id(self, idList):
        handle = Entrez.elink(db='mesh', dbfrom='pccompound', id=idList, linkname='pccompound_mesh')
        result = Entrez.read(handle)
        results = []
        for id, i in zip(idList, result):
            links = []
            for j in i["LinkSetDb"]:
                for link in j["Link"]:
                    links.append(link['Id'])
            if len(links) > 0:
                results.append(links)
        return results

    def extract_important_section(self, sample: str):
        match = re.search(r'\bEntry Terms:\n', sample)
        literal_name = sample.split('\n')[1]
        lname_macth = re.search(r':', literal_name)
        if lname_macth:
            literal_name = literal_name[lname_macth.span()[1]:]
        if match:
            return literal_name, list(
                map(lambda x: x.strip(), sample[match.span()[1]:].split('\n\n')[0].strip().split('\n')))
        # else:
        #     print('Notfound', sample)
        return literal_name.strip(), []

    def get_mesh_term_from_id(self, results) -> list:
        res = []
        for v in results:
            handle = Entrez.efetch(db='mesh', id=list(v), rettype='xml')
            temp = handle.read().split('\n\n\n')[:-1]
            result = dict(map(self.extract_important_section, temp)).keys()
            res.append(list(result)[0].strip())
        return res

    def pprint(self, result: dict):
        if len(result) > 0:
            for k, v in result.items():
                for s, d in v.items():
                    print(f'Mesh Entry={s} for cID={k}')
                    # print('Mesh Terms:')
                    # [print(f'\t {r}') for r in d]
        else:
            print('No results found')

    def cprint(self, result: dict):
        if len(result) > 0:
            for k, v in result.items():
                for s, d in v.items():
                    self.title_print(f'Mesh Entry={s} for cID={k}')
                    self.text_print('Entry Terms:')
                    [self.ctext_print(f'\t {r}') for r in d]
        else:
            print('No results found')

    def get_color(self):
        c = self.color_num
        self.color_num += 1
        return c

    def title_print(self, str):
        print(f'<h5>{str}</h5')

    def text_print(self, str):
        print(f'<a>{str}</a')

    def ctext_print(self, str):
        bcolor = ['ff9292', 'ffb4b4', 'ffdcdc', 'ffe8e8']
        print(f'<a style="background-color:{bcolor[self.get_color() % len(bcolor)]}" >{str}</a')

    def cprint(self, result: dict):
        if len(result) > 0:
            for k, v in result.items():
                for s, d in v.items():
                    print('<div class="jumbotron"><div class="container">')

                    self.title_print(f'Mesh Terms for {s} (cID={k})')
                    # self.text_print('Mesh Terms:')
                    [self.ctext_print(f'{r}') for r in d]
                    print('</div></div>')
        else:
            print('No results found')


def fetch_projects():
    import requests
    import json
    url = "https://api.gdc.cancer.gov/projects"
    # The 'fields' parameter is passed as a comma-separated string of single names.
    fields = [
        "project.project_id"
    ]

    fields = ','.join(fields)

    params = {
        # "fields": fields,
        # "format": "TSV",
        "from": "0",
        "size": "1000"
    }

    response = requests.get(url, params=params)
    data = response.json()  # Check the JSON Response Content documentation below
    # print(json.dumps(data, indent=4, sort_keys=True))

    project_ids = [i['project_id'] for i in data['data']['hits'] if i['project_id'].startswith("TCGA")]
    return project_ids


def statisctics(data, filter_vals=None):
    pdData = pd.DataFrame.from_dict(data).dropna(subset=['term'])
    # filter_vals = ['Carboplatin', 'transplatin [Supplementary Concept]']
    if filter_vals != None:
        pdData = pdData[pdData['term'].isin(filter_vals)]

    # %%
    def survived(x):
        return x.isnull().sum()

    def mean_survived(x):
        return np.mean(x)

    def smoking(x):
        return x[x == 'Smoking'].count()

    # %%
    ex = pdData[['therapy_types', 'days_to_death', 'tobacco_smoking_history', 'term']].groupby("term")
    # %%
    aggs = {'days_to_death': survived, 'tobacco_smoking_history': smoking}
    ex2 = ex.agg(aggs)
    ex2.columns = ['survived', 'smoking']
    print(ex2)
    # %%
    ex2.plot(kind='bar', figsize=(10, 19))
    plt.show()
    # %%
    ex3 = pdData[['therapy_types', 'days_to_death', 'term']].groupby(['term', 'therapy_types']).agg(
        mean_survived).sort_values('days_to_death', ascending=False)
    print(ex3)
    ex3.plot(kind='bar', figsize=(10, 19))
    plt.show()


if __name__ == '__main__':
    import sys
    import os

    cancer_type = "TCGA-LUAD"
    project_file = f"{cancer_type}.csv"
    if len(sys.argv) > 1:
        cancer_type = sys.argv[1]
        project_file = f"{cancer_type}.csv"
        if not os.path.exists(project_file):
            command = f"./tcga.R \"{cancer_type}\""
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            process.wait()

        data = pd.read_csv(project_file).dropna(subset=['drug_name']).sort_values(by=['days_to_death'],
                                                                                  ascending=True).to_dict('list')
        pc = PubChemCompound()
        data['c_id'], data['term'] = pc.get_pick(data['drug_name'])
        statisctics(data)
