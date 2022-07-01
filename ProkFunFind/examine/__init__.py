import json
import logging


def read_json(system_file):
    with open(system_file) as f:
        system_dict = json.load(f)
    return system_dict


def find_gene_with_orthoID(orthoID, genomeObject, detect_tools):
    """TODO: find the genes in the genomeObject with the orthoID
    :returns: A list of genes fit the above requirement
    """
    gene_list = []
    for detect_tool in detect_tools:
        gene_list += [
            gene.id for gene in genomeObject.genes.values()
            if hasattr(gene, detect_tool)
            and getattr(getattr(gene, detect_tool), "orthoID") == orthoID
        ]
    return gene_list


def check_gene_in_subsystem(system_dict,
                            genomeObject,
                            detect_tools,
                            fun_string=[]):
    """
    Assign genes assigned with orthoID to the system_dict
    """

    completeness = {
        'essential': 0,
        'nonessential': 0,
        'essential_presence': 0,
        'nonessential_presence': 0
    }

    status = 0
    fun_list = fun_string.copy()

    if "analogs" in system_dict:  # if the system has analogs_status
        _, analogs_status = check_gene_in_subsystem(system_dict["analogs"],
                                                    genomeObject,
                                                    detect_tools,
                                                    fun_string=fun_list)
        if analogs_status:
            status = 1

    if "name" in system_dict:
        fun_list.append(system_dict['name'])
        for sub_component in system_dict['components']:
            _, sub_status = check_gene_in_subsystem(sub_component,
                                                    genomeObject,
                                                    detect_tools,
                                                    fun_string=fun_list)
            completeness[sub_component['presence']] += 1
            completeness[sub_component['presence'] + "_presence"] += sub_status
        system_dict['completeness'] = completeness
        if completeness['essential'] == completeness[
                'essential_presence'] and completeness['essential'] > 0:
            status = 1
        if system_dict['completeness']['essential'] == 0:
            logging.warning("Component {} was defined without any essential \
                            subcomponents. Please define at least one \
                            subcomponent as essential.".format(
                            system_dict['name']))
            quit()

    elif "orthoID" in system_dict:
        fun_list.append(system_dict['orthoID'])
        gene_list = find_gene_with_orthoID(system_dict['orthoID'],
                                           genomeObject, detect_tools)
        if gene_list:
            status = 1
            system_dict['genes'] = gene_list
            for gene in gene_list:
                if hasattr(genomeObject.genes[gene], "Functions"):
                    genomeObject.genes[gene].Functions.append(
                        "/".join(fun_list))
                else:
                    setattr(genomeObject.genes[gene],
                            "Functions", ["/".join(fun_list)])

    return (system_dict, status)


def pipeline(system_file, genome_object, detect_tools):
    system_dict = read_json(system_file)
    system_dict, status = check_gene_in_subsystem(system_dict,
                                                  genomeObject=genome_object,
                                                  detect_tools=detect_tools)
    return (system_dict, status, genome_object)
