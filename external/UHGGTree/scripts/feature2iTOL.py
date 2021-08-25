#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ete3 import Tree
import math

################################################################################
from color import ggsci_db

default_col = ggsci_db["default"]

def dt2str(dt):
    return "|".join([i + "-" + str(dt[i]) for i in dt])


def str2dt(st):
    return {
        string.split("-")[0]: int(string.split("-")[1])
        for string in st.split("|")
    }


####################################################
#  return if the rank is lower than the tip_rank   #
####################################################


def rank_lower(n, tip_rank):
    order_rank = [
        "root", "domain", "phylum", "class", "order", "family", "genus",
        "species", "strain"
    ]
    return order_rank.index(tip_rank) < order_rank.index(n.rank)


#############################################################
#  return if all the children are of the same state or not  #
#############################################################


def same_state(state):
    state_dt = str2dt(state)
    val = set(state_dt.values())
    if val == {0}:
        return True
    elif len(val) == 2 and (0 in val):
        return True
    else:
        return False


################################################################################
################################################################################


def main(args):
    infile = args.infile
    attfile = args.attfile          # provide gene names if annotfile=extshape (@dombraccia)
    outfile = args.outfile
    min_leaves = args.min_leaves
    tip_rank = args.tip_rank
    collapse = args.collapse
    annotfile = args.annotfile      # specify piechart, externalshape or barchart (@dombraccia)

    tree = Tree(infile, format=3)

    if (annotfile == 'piechart'):
        feature_dt = {
            line.rstrip('\n').split("\t")[0]: line.rstrip('\n').split("\t")[1]
            for line in open(attfile)
        }
        state_set = list(set(feature_dt.values()))
    elif (annotfile in ['extshape', 'barchart']):
        feature_dt = {
            line.rstrip('\n').split("\t")[0]: line.rstrip('\n').split("\t")[1].split(',')
            for line in open(attfile)
        }
        state_set = list(set(sum(feature_dt.values(), [])))
    
    leaves_set = feature_dt.keys()

    #############################################################
    #  only report the leaves that only present in the attfile  #
    #############################################################
    #
    #node_set = []
    #for leaf_name in leaves_set:
    #    leaf_node = tree.get_leaves_by_name(leaf_name)[0]
    #    node_set += leaf_node.get_ancestors()
    #    node_set.append(leaf_node)
    #
    #tree.prune(node_set,preserve_branch_length = True)
    #
    ##################################################################
    #  read the state from the attfile and assign the state to tree  #
    ##################################################################

    if (annotfile == 'piechart'):
        for node in tree.traverse("postorder"):
            state_dt = {i: 0 for i in state_set}
            if node.is_leaf():
                if node.name in feature_dt:
                    state_dt[feature_dt[node.name]] = 1
            else:
                for child in node.get_children():
                    child_dt = str2dt(child.state)
                    for key in state_dt:
                        state_dt[key] += child_dt[key]
            node.add_features(leaves=len(
                node.get_leaves()))  # add the number of leaves for internal node
            node.add_features(state=dt2str(state_dt))

    ##### ========= modified for externalshape (@dombraccia) ========= #####
    
    if (annotfile in ['extshape', 'barchart']):
        for node in tree.traverse("postorder"):
            state_dt = {i: 0 for i in state_set}
            if node.is_leaf():
                if node.name in feature_dt:
                    for gene in feature_dt[node.name]:
                        state_dt[gene] = 1
            else:
                for child in node.get_children():
                    child_dt = str2dt(child.state)
                    for key in state_dt:
                        state_dt[key] += child_dt[key]
            node.add_features(leaves=len(
                node.get_leaves()))  # add the number of leaves for internal node
            node.add_features(state=dt2str(state_dt))

    ####################################################
    #  process the tree and select the node to report  #
    ####################################################

    for node in tree.traverse("preorder"):
        if node.is_root(): continue
        if int(node.leaves) < int(min_leaves):
            #print("remove due to less than "+ str(min_leaves) + ":" + node.name)
            node.detach()
        elif rank_lower(node, tip_rank):
            #print("remove due to less than "+ tip_rank + ":"  +node.name+"\t"+node.rank)
            node.detach()
        if collapse and same_state(node.state):
            #print("remove due to same state "+ node.name +":"  +"\t".join([i.name for i in node.get_children()]))
            [i.detach() for i in node.get_children()]

    #######################################
    #  generate annotation file for iTOL  #
    #######################################

    # check if PIECHART annotation file specified
    if (annotfile == 'piechart'):

        col_len = len(state_set)

        f = open(outfile + ".piechart.txt","w")
        f.write("DATASET_PIECHART\nSEPARATOR TAB\nDATASET_LABEL\tpiechart\nCOLOR\t#ff0000\n")
        f.write("FIELD_COLORS\t")
        if col_len <= 12:
            f.write("\t".join([ i + "80" for i in default_col[0:col_len] ])) 
        else: 
            f.write("# TODO: more than 12 type; Please add color palette here; separated by tab")
        f.write("\nFIELD_LABELS\t" + "\t".join(state_set) + "\nDATA\n")

        for node in tree.traverse("preorder"):
            if node.is_leaf():
                state_dt = str2dt(node.state)
                state_list = [str(state_dt[i]) for i in state_set]
                f.write(node.name + "\t1\t" + str(math.sqrt(node.leaves)) +
                    "\t" + "\t".join(state_list)+"\n")

        f.close()

    #####################################################
    #  generate EXTERNALSHAPE annotation file for iTOL  #
    #  author: domenick j. braccia
    #####################################################

    # check if EXTERNALSHAPE annotation file specified
    if (annotfile == 'extshape'):

        col_len = len(state_set)

        f = open(outfile + ".extshape.txt","w")
        f.write("DATASET_EXTERNALSHAPE\nSEPARATOR TAB\nDATASET_LABEL\textshape\nCOLOR\t#ff0000\n")
        f.write("FIELD_COLORS\t")
        if col_len <= 12:
            f.write("\t".join([ i + "80" for i in default_col[0:col_len] ])) 
        else:
            f.write("# TODO: more than 12 type; Please add color palette here; separated by tab")
        f.write("\nFIELD_LABELS\t" + "\t".join(state_set) + "\nDATA\n")

        for node in tree.traverse("preorder"):
            if node.is_leaf():
                state_dt = str2dt(node.state)
                state_list = [str(round(math.sqrt(state_dt[i]/node.leaves), 3)) for i in state_set]
                f.write(node.name + "\t" +
                    "\t".join(state_list)+"\n")

        f.close()

    #####################################################
    #  generate BARCHART annotation file for iTOL  #
    #  author: domenick j. braccia
    #####################################################

    # check if EXTERNALSHAPE annotation file specified
    if (annotfile == 'barchart'):

        col_len = len(state_set)

        f = open(outfile + ".barchart.txt","w")
        f.write("DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tbarchart\nCOLOR\t#ff0000\n")
        f.write("FIELD_COLORS\t")
        if col_len <= 12:
            f.write("\t".join([ i + "80" for i in default_col[0:col_len] ])) 
        else:
            f.write("# TODO: more than 12 type; Please add color palette here; separated by tab")
        f.write("\nFIELD_LABELS\t" + "\t".join(state_set) + "\nDATA\n")

        for node in tree.traverse("preorder"):
            if node.is_leaf():
                f.write(node.name + "\t" + str(node.leaves) + "\n")

        f.close()

    ############################
    #  generate tree for iTOL  #
    ############################

    tree.write(format=3,
               outfile=outfile+".tree",
               features=["state", "leaves", "rank"],
               format_root_node=True)


if __name__ == "__main__":
    from argparse import ArgumentParser, FileType

    parser = ArgumentParser(
        description='Generate tree as well as PIECHART annotation file for iTOL'
    )
    parser.add_argument('-t',
                        help='input tree path',
                        #type=FileType('r'),
                        required=True,
                        dest='infile',
                        metavar='tree')
    parser.add_argument('-f',
                        help='path to feature file',
                        #type=FileType('r'),
                        required=True,
                        dest='attfile',
                        metavar='table')
    parser.add_argument('-p',
                        help='the prefix of output tree and annotation',
                        type=str,
                        required=True,
                        dest='outfile',
                        metavar='prefix')
    parser.add_argument('-m',
                        help='report nodes with minimal number of leaves',
                        type=int,
                        default=1,
                        required=False,
                        dest='min_leaves',
                        metavar='num')
    parser.add_argument('-r',
                        help='remove leaves lower than the rank',
                        type=str,
                        default="species",
                        required=False,
                        dest='tip_rank',
                        choices=[
                            "root", "domain", "phylum", "class", "order",
                            "family", "genus", "species", "strain"
                        ],
                        metavar='rank')
    parser.add_argument('-c',
                        help='collapse the clades whoes children are of same state',
                        action='store_true',
                        required=False,
                        dest='collapse')
    parser.add_argument('-a',
                        help='type of output dataset desired (piechart or barchart)',
                        default='piechart',
                        required=False,
                        dest='annotfile',
                        choices=['piechart', 'extshape', 'barchart'])
    args = parser.parse_args()

    main(args)
