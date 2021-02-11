import cgi
import cgitb
from Bio import Entrez

cgitb.enable()

print("Content-Type: text/html")
print()
print('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
print('''<!-- Latest compiled and minified CSS -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">

<!-- Optional theme -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">

<!-- Latest compiled and minified JavaScript -->
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>''')
print('<body>')


class PubChemCompound:
    def __init__(self, mail="furkanfbr@gmail.com"):
        Entrez.email = mail
        # handle = Entrez.einfo()
        # record = Entrez.read(handle)
        self.color_num = 0

    def fetch_data(self, compound, max_drugs: int = 10, max_compounds: int = 10):
        idList = self.__search_cid_from_compound(compound, max_drugs)
        meshIds = self.__get_mesh_term_id(idList)
        results = self.__get_mesh_term_from_id(meshIds, max_compounds)
        return results

    def __search_cid_from_compound(self, compound, max_drugs: int = 10):

        idList = []
        handle = Entrez.esearch(db='pccompound', term=f'"{compound}"[Synonym]', rettype='xml', retMax=max_drugs,
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

    def __get_mesh_term_from_id(self, results, max_compounds: int = 10) -> dict:
        for k, v in results.items():
            handle = Entrez.efetch(db='mesh', id=list(v.keys())[:max_compounds], rettype='xml')
            temp = handle.read().split('\n\n\n')[:-1]
            result = dict(map(self.extract_important_section, temp))
            results[k] = result
        return results

    def pprint(self, result: dict):
        if len(result) > 0:
            for k, v in result.items():
                for s, d in v.items():
                    print(f'Mesh Entry={s} for cID={k}')
                    print('Mesh Terms:')
                    [print(f'\t {r}') for r in d]
        else:
            print('No results found')

    def cprint(self, result: dict):
        if len(result) > 0:
            for k, v in result.items():
                for s, d in v.items():
                    print('<div class="jumbotron">')
                    print('<div class="container">')
                    self.title_print(f'Mesh Term for cID={k} is {s} ')
                    # self.text_print('Mesh Terms:')
                    # [self.ctext_print(f'{r}') for r in d]
                    print('</div></div>')
        else:
            print('No results found')

    def get_color(self):
        c = self.color_num
        self.color_num += 1
        return c

    def title_print(self, str):
        print(f'<h2 style="color:red">{str}</h2>')

    def text_print(self, str):
        print(f'<a>{str}</a>')

    def ctext_print(self, str):
        bcolor = ['ff9292', 'ffb4b4', 'ffdcdc', 'ffe8e8', 'ffdcdc', 'ffb4b4']
        print(f'<p style="margin-left:50;display:list-item;" >{str}</p>')


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        pc = PubChemCompound()
        result = pc.fetch_data(compound=str(sys.argv[1]))
        pc.pprint(result)

print("""<form class="form-inline mt-2 mt-md-0" action = \"main.py\" method = get> <input style="width:100%" placeholder="Enter a drug name" class="form-control mr-sm-2" type = text name = compound />  <input style=" width: 100%;background-color:green;color:white"type = submit class="btn btn-outline-success my-2 my-sm-0" value = Search /> </form>
""")
import sys

if len(sys.argv) > 1:
    pc = PubChemCompound()
    result = pc.fetch_data(compound=str(sys.argv[1]), max_drugs=10, max_compounds=1)
    pc.pprint(result)
else:
    form = cgi.FieldStorage()
    compound = form.getvalue('compound')
    rettype = form.getvalue('rettype')
    if compound is not None:
        print(f"<h3>Showing results for {compound.title()}</h3>")
        pc = PubChemCompound('huseyinfurkankiyikci@posta.mu.edu.tr')
        pc.cprint(pc.fetch_data(compound=compound, max_drugs=10, max_compounds=1))

print('</body>')
