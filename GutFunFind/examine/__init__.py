import json


def read_json(system_file):
    with open(system_file) as f:
        system_dict = json.load(f)
    return system_dict


def find_gene_with_orthoID(orthoID, genomeObject, detect_tool):
    """TODO: find the genes in the genomeObject with the orthoID assigned by detect_tool
    :returns: A list of genes fit the above requirement
    """
    gene_list = [
        gene.id for gene in genomeObject.genes.values() if hasattr(
            gene,
            detect_tool) and getattr(
            getattr(
                gene,
                detect_tool),
            "orthoID") == orthoID]
    return(gene_list)


def check_gene_in_subsystem(
        system_dict,
        genomeObject,
        detect_tool,
        fun_string=[]):
    """
    Assign genes assigned with orthoID to the system_dict
    """
    fun_list = fun_string.copy()
    completeness = {
        "essential": 0,
        "nonessential": 0,
        "essential_presence": 0,
        "nonessential_presence": 0}
    if "analogs" in system_dict:
        _, analogs_status = check_gene_in_subsystem(
            system_dict["analogs"], genomeObject, detect_tool, fun_string=fun_list)

    if "name" in system_dict:
        fun_list.append(system_dict["name"])

        for sub_component in system_dict["components"]:
            _, status = check_gene_in_subsystem(
                sub_component, genomeObject, detect_tool, fun_string=fun_list)

            completeness[sub_component["presence"]] += 1
            completeness[sub_component["presence"] + "_presence"] += status

    elif "orthoID" in system_dict:
        fun_list.append(system_dict["orthoID"])
        gene_list = find_gene_with_orthoID(
            system_dict["orthoID"], genomeObject, detect_tool)
        status = 0
        if gene_list:
            system_dict["genes"] = gene_list
            status = 1
            for gene in gene_list:
                if hasattr(genomeObject.genes[gene], "Functions"):
                    genomeObject.genes[gene]["Functions"].append(
                        "/".join(fun_list))
                else:
                    setattr(
                        genomeObject.genes[gene], "Functions", [
                            "/".join(fun_list)])
        return (system_dict, status)
        fun_list.pop()

    status = 0
    if completeness["essential"] == completeness["essential_presence"] and completeness["essential"] > 0:
        status = 1
    elif analogs_status:
        status = 1

    if "name" in system_dict:
        system_dict["completeness"] = completeness

    return (system_dict, status)


def pipeline(system_file, genome_object, detect_tool):
    system_dict = read_json(system_file)
    system_dict, status = check_gene_in_subsystem(
        system_dict, genomeObject=genome_object, detect_tool=detect_tool)
    return (system_dict, status, genome_object)
