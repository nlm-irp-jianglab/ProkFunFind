import operator
import csv
from collections import defaultdict

from Bio.File import as_handle
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from ProkFunFind.toolkit.utility import check_path_existence


class baktaTabParser:
    """Class for a Parser of the bakta table format.

       Attributes:
          handle: str
              open file object
          line: str
              next line in the file

       Methods:
          __iter__: Function to iterate over file and return QueryResults
          _parse_row: Function to parse and split a row in the file
          _parse_qresult: Function to take row data and make QueryResult

    """

    def __init__(self, handle):
        """Initialize the class.

           Arguments:
               handle: open bakta tabular output file.
        """
        self.handle = handle
        self.line = self.handle.readline().rstrip("\n")

    def __iter__(self):
        """Iterate over baktaTabParser, yields query results."""
        header_mark = "#"
        # read through the header if it exists
        while self.line.startswith(header_mark):
            self.line = self.handle.readline()
        # if we have result rows, parse it
        if self.line:
            yield from self._parse_qresult()

    def _parse_row(self):
        """Return a dictionary of parsed row values (PRIVATE)."""
        cols = self.line.strip("\n").split("\t")

        # Number of columns should be 7 in the standard bakta table
        if len(cols) != 9:
            raise ValueError("Less columns than expected, only %i" % len(cols))

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[5]  # gene name
        qresult['program'] = "bakta"

        hit = {}
        if cols[8] == '':
            xrefs = ['']
        else:
            xrefs_full = cols[8].split(", ")
            xrefs = [i.split(":")[1] for i in xrefs_full]
        hit['id'] = xrefs  # DB ID
        hit['description'] = cols[7]  # description of target
        hit['query_id'] = cols[5]  # query name
        hsp = {}
        # evalue or score should be float but sometimes not
        hsp['evalue'] = None

        frag = {}

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp, 'frag': frag}

    def _parse_qresult(self):
        """Parse query results (PRIVATE)."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3

        # initial value dummies
        qres_state = None
        file_state = None

        cur_qid = None
        prev_qid = None

        cur, prev = None, None
        hsp_dict = defaultdict(list)
        hit_list = []

        while True:
            # store previous line's parsed values for all lines after the first
            if cur is not None:
                prev = cur
                prev_qid = cur_qid

            if self.line and not self.line.startswith("#"):
                cur = self._parse_row()
                cur_qid = cur['qresult']['id']
            else:
                file_state = state_EOF
                # mock values for cur_qid since the line is empty
                cur_qid = None
            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME

            if prev is not None:
                for cog in prev['hit']['id']:
                    prev_hid = cog

                    frag = HSPFragment(prev_hid, prev_qid)

                    for attr, value in prev['frag'].items():
                        setattr(frag, attr, value)
                    hsp = HSP([frag])

                    for attr, value in prev['hsp'].items():
                        setattr(hsp, attr, value)
                    hsp_dict[prev_hid].append(hsp)

                    hit = Hit()
                    for attr, value in prev['hit'].items():
                        if attr == "id":
                            setattr(hit, "id", cog)
                        else:
                            setattr(hit, attr, value)
                    if hit.id not in [i.id for i in hit_list]:
                        hit_list.append(hit)

                # create qresult and yield if we're at a new qresult or at EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:

                    for hit in hit_list:
                        for hsp in hsp_dict[hit.id]:
                            hit.hsps.append(hsp)
                    qresult = QueryResult(hit_list, prev_qid)
                    for attr, value in prev['qresult'].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    # if we're at EOF, break
                    if file_state == state_EOF:
                        break
                    hsp_dict = defaultdict(list)
                    hit_list = []

            self.line = self.handle.readline()


def bakta_tab_parse(handle, **kwargs):
    """Parse bakta table and yield results

       Arguments:
           handle: open bakta tabular file
           kwargs:

       Yields:
           yields from iterator of the source file.
    """
    # get the iterator object and do error checking
    mod = __import__("ProkFunFind.detect.bakta_search.bakta_parse",
                     fromlist=[""])
    obj_name = "baktaTabParser"
    iterator = getattr(mod, obj_name)

    # and start iterating
    with as_handle(handle) as source_file:
        generator = iterator(source_file, **kwargs)
        yield from generator


