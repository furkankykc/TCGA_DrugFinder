#!/usr/bin python
import cgi
import cgitb

cgitb.enable()

print("Content-Type: text/html")
print()

from Bio import Entrez


class PubChemCompound:
    def __init__(self, mail="furkanfbr@gmail.com"):
        Entrez.email = mail
        # handle = Entrez.einfo()
        # record = Entrez.read(handle)
        # print(record["DbList"])

    def fetch_data(self, compound):
        idList = self.__search_cid_from_compound(compound)
        meshIds = self.__get_mesh_term_id(idList)
        results = self.__get_mesh_term_from_id(meshIds)
        return results

    def __search_cid_from_compound(self, compound):

        idList = []
        handle = Entrez.esearch(db='pccompound', term=f'"{compound}"[Synonym]', rettype='xml', retMax=10,
                                sort='relevance')
        res = Entrez.read(handle)
        idList.extend(res['IdList'])
        return idList

    def __get_mesh_term_id(self, idList):
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

    def extract_important_section(self, sample):
        import re
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

    def __get_mesh_term_from_id(self, results) -> dict:
        for k, v in results.items():
            handle = Entrez.efetch(db='mesh', id=v.keys(), rettype='xml')
            temp = handle.read().split('\n\n\n')[:-1]
            result = dict(map(self.extract_important_section, temp))
            results[k] = result
        return results

    def pprint(self, result: dict):
        for k, v in result.items():
            for s, d in v.items():
                print(f'Mesh Entry={s} for cID={k}')
                print('Mesh Terms:')
                [print(f'\t {r}') for r in d]

        # [print(r) for r in result['2244']['68001241'].split('\n')]


form = cgi.FieldStorage()
compound = form.getvalue('compound')
rettype = form.getvalue('rettype')

pc = PubChemCompound('huseyinfurkankiyikci@posta.mu.edu.tr')
print("<pre>")
pc.pprint(pc.fetch_data(compound=compound))
print("</pre>")
