# Copyright 2020 by Xiaofang Jiang. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SearchIO parser for InterProScan tab output formats."""
# for more info: https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats

from collections import defaultdict

from Bio.File import as_handle
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
import operator

class InterproscanTabParser:
    """Parser for the InterProScan table format.

       Attributes:
           handle: str
               open file object
           line: str
               next line in the file

       Methods:
           __iter__: Function to iterate over file and return QueryResults
           _parse_row: Function to parse and split row in file
           _parse_qresult: Function to take row data and make QueryResult
    """

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = self.handle.readline().rstrip("\n")

    def __iter__(self):
        """Iterate over InterproscanTabParser, yields query results."""
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
        if len(cols) <5:
            raise ValueError("Less columns than expected, only %i" % len(cols))

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[0]  # query name
        qresult['seq_len'] = int(cols[2])  # query length
        qresult['program'] = "InterProScan"  #

        hit = {}
        hit['id'] = cols[4]  # target name
        hit['description'] = cols[5]  # description of target
        hit['query_id'] = cols[0]  # query name
        hit['attributes'] = {'Target': str(cols[4])}
        if len(cols) == 15:
            xrefs = ["IPR:" + cols[11]] if cols[11] else []
            xrefs += cols[13].split("|") if cols[13] else []
            xrefs += cols[14].replace(" ", "").split("|") if cols[14] else []
            hit['dbxrefs'] = xrefs
        hit['evalue'] = float(cols[8]) if cols[8] != "-" else None

        hsp = {}
        # evalue or score should be float but sometime not
        hsp['evalue'] = float(cols[8]) if cols[8] != "-" else None

        frag = {}
        frag['query_start'] = int(cols[6]) - 1  # query start, zero-based
        frag['query_end'] = int(cols[7])  # query end

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
                # since domain tab formats only have 1 Hit per line
                # we always create HSPFragment, HSP, and Hit per line
                prev_hid = prev['hit']['id']

                frag = HSPFragment(prev_hid, prev_qid)

                for attr, value in prev['frag'].items():
                    setattr(frag, attr, value)
                hsp = HSP([frag])

                for attr, value in prev['hsp'].items():
                    setattr(hsp, attr, value)
                hsp_dict[prev_hid].append(hsp)

                hit = Hit()
                for attr, value in prev['hit'].items():
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


def ipr_tab_parse(handle, **kwargs):
    """Fucntion to handle parsing of IPR file

       Arguments:
           handle: IPR file object.
    """
    # get the iterator object and do error checking
    mod = __import__("ProkFunFind.detect.ipr_search.interproscan_tab",
                     fromlist=[""])
    obj_name = "InterproscanTabParser"
    iterator = getattr(mod, obj_name)

    # and start iterating
    with as_handle(handle) as source_file:
        generator = iterator(source_file, **kwargs)
        yield from generator


def ipr_filter(config: dict, qres: QueryResult, basedir: str, filter_dict: dict) -> QueryResult:
    global_evalue = config['interproscan'].get('evalue', 0.01)
    ops = {
        '<=': operator.le,
        '>=': operator.ge,
        '>': operator.gt,
        '<': operator.lt,
        '==': operator.eq,
        '!=': operator.ne
    }

    def hsp_filter_func(hsp):
        status = True
        if hsp.hit_id in filter_dict:
            for one in filter_dict[hsp.hit_id]:
                if one['cpfun'](getattr(hsp, one['attr']), one['value']):
                    pass
                else:
                    status = False
                    break
        else:
            if hsp.evalue is None:
                status = False
            elif hsp.evalue > float(global_evalue):
                status = False
        return status
    return qres.hsp_filter(hsp_filter_func)
