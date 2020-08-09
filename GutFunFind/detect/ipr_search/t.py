import csv
from collections import defaultdict
from Bio import SearchIO
from typing import Dict, IO, List, Set, OrderedDict, Callable
from Bio.SearchIO._model.query import QueryResult


qresults = SearchIO.parse("MGYG-HGUT-02420.xml", "interproscan-xml")
ortho_pair_file="domain_precision.txt"

from collections import defaultdict
OrthScore_dict = defaultdict(dict)
#############################################################
#  User can change the last column to indicate specificity  #
#############################################################
with open(ortho_pair_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    header = next(csv_reader)
    col_num = len(header)
    if col_num == 2:
        OrthScore_dict[header[1]] = {"orthoID":header[0],"precision": 1}
        for row in csv_reader:
            OrthScore_dict[row[1]] = {"orthoID":row[0], "precision":1}
    else:
        OrthScore_dict[header[1]] = {"orthoID":header[0],"precision": float(header[2])}
        for row in csv_reader:
            OrthScore_dict[row[1]] = {"orthoID":row[0],"precision": float(row[2])}


#OrthScore_dict.keys()
#def hsp_filter_func(hsp):
#    return(hsp.hit_id in OrthScore_dict.keys())
#q_list = [i for _, i in SearchIO.to_dict(qresults).items() if len(i) > 0]
#filter_res = [blast_filter(config=filter_cf, qres=i) for i in q_list]
#qres.hsp_filter(hsp_filter_func)
#q_list = [ i.hsp_filter(hsp_filter_func) for _, i in SearchIO.to_dict(qresults).items() if len(i.hsp_filter(hsp_filter_func)) >0 ]
#for i in q_list:
#    dtlist = [OrthScore_dict[x.id] for x in i]
#    print(i.id)
#    print(sorted(dtlist, key = lambda i: i['precision'],reverse=True)[0])

q_list=[]
for qres in qresults:

    # remove QueryResult that doest not hit any domain in function-related domain list
    i = qres.hsp_filter(lambda hsp: hsp.hit_id in OrthScore_dict.keys())


    # for those without any hit match to the domain
    if len(i) > 0 :
        
        # sort the hits based on precision
        i.sort(key = lambda hit:OrthScore_dict[hit.id]["precision"],reverse=True)

        max_dict = OrthScore_dict[i.hits[0].id]

        # set the QueryResult attribution
        setattr(i, "orthoID", max_dict["orthoID"])
        setattr(i, "orthoID_weight", max_dict["precision"])
        q_list.append(i)

