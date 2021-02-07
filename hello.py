#!/usr/bin python
import cgi
import cgitb

cgitb.enable()
from main import PubChemCompound

print("Content-Type: text/html")
print()
form = cgi.FieldStorage()
compound = form.getvalue('compound')
rettype = form.getvalue('rettype')
print("""<form action = \"~hf.kiyikci19/cgi-bin/main.py\" method = get> Enter a drug name: <input type = text name = db /><br /> <input type = submit value = Search /> </form>"
""")
pc = PubChemCompound('huseyinfurkankiyikci@posta.mu.edu.tr')
print("<pre>")
pc.pprint(pc.fetch_data(compound=compound))
print("</pre>")
