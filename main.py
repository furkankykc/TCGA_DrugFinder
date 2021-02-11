from Bio import Entrez
import re


class PubChemCompound:

    def __init__(self, mail="furkanfbr@gmail.com"):
        Entrez.email = mail
        self.color_num = 0

    def fetch_data(self, compound, max_drug: int = 10, max_comp: int = 10):
        idList = self.search_cid_from_compound(compound, max_drug)
        meshIds = self.get_mesh_term_id(idList)
        results = self.get_mesh_term_from_id(meshIds, max_comp)
        return results

    def extract_cid_n_term(self, data: dict):
        return list(data.keys()), list(map(lambda x: x.strip(), list(list(data.values()).keys())))

    def search_cid_from_compound(self, compound, max_drug: int = 10):

        idList = []
        for c in compound:
            handle = Entrez.esearch(db='pccompound', term=f'"{c}"[Synonym]', rettype='xml', retMax=max_drug,
                                    sort='relevance')
            res = Entrez.read(handle)
            idList.extend(res['IdList'])
        return idList

    def get_mesh_term_id(self, idList):
        handle = Entrez.elink(db='mesh', dbfrom='pccompound', id=idList, linkname='pccompound_mesh')
        result = Entrez.read(handle)
        results = {}
        for id, i in zip(idList, result):
            links = {}
            for j in i["LinkSetDb"]:
                for link in j["Link"]:
                    links[link['Id']] = ''
            if len(links) > 0:
                results[id] = links
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
        return literal_name, []

    def get_mesh_term_from_id(self, results, max_comp: int = 10) -> dict:
        for k, v in results.items():
            handle = Entrez.efetch(db='mesh', id=list(v.keys())[:max_comp], rettype='xml')
            temp = handle.read().split('\n\n\n')[:-1]
            result = dict(map(self.extract_important_section, temp))
            results[k] = result
        return results

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


if __name__ == '__main__':
    import sys
    import pandas as pd

    # if len(sys.argv) > 1:
    data = pd.read_csv('temp.csv').sort_values(by=['days_to_death'], ascending=True).to_dict('list')
    pc = PubChemCompound()
    print(pc.fetch_data(data['drug_name'][:10], 1, 1))

    # result = pc.fetch_data(compound=str(sys.argv[1]), max_drug=10, max_comp=1)
    # pc.pprint(result)
