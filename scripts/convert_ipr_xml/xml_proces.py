#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xml.etree.ElementTree as ET
import re
import sys
import os
inputfile = sys.argv[1]
outdir = sys.argv[2]

import os
if not os.path.exists(outdir):
    os.makedirs(outdir)


ET.register_namespace('', "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5")
tree = ET.parse(inputfile)
root = tree.getroot()
for child in root:
    filename = re.sub("_.....$","", child[1].attrib["id"])
    string=ET.tostring(child, encoding='utf8').decode('utf8')
    f=open(outdir+"/"+filename+".xml", "a+")
    f.write(string.replace("<?xml version=\'1.0\' encoding=\'utf8\'?>\n<protein xmlns=\"http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5\">","    <protein>"))


#root[3][1].attrib["id"]
#re.sub("_.....$","",a)
#t=root[0]
#string=ET.tostring(t, encoding='utf8').decode('utf8')
#print(string.replace("<?xml version=\'1.0\' encoding=\'utf8\'?>\n<protein xmlns=\"http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5\">","    <protein>"))
#print(string)
#ET.dump(t)
