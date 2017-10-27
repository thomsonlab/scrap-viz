import urllib3 as urllib
import pandas

gene_summary_marker = "Entrez Gene Summary for"
gene_name_marker = "aliasMainName"


def get_gene_summaries(genes):

    gene_name_descriptions = []

    http = urllib.PoolManager()

    for gene in genes:

        url = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" % gene
        response = http.request("GET", url)

        if response.status != 200:
            gene_name_descriptions.append((gene, "No response"))
            continue

        string = response.data.decode("UTF-8")

        if string.find("was not found") != -1:
            gene_name_descriptions.append((gene, "Gene not found"))
            continue

        gene_name_marker_location = string.find(gene_name_marker)

        if gene_name_marker_location == -1:
            gene_name_descriptions.append((gene, "Alias not found"))
            continue

        name_start_index = string.find(">", gene_name_marker_location)
        name_end_index = string.find("<", name_start_index)

        gene_name = string[name_start_index+1:name_end_index]

        gene_summary_marker_location = string.find(gene_summary_marker)

        if gene_summary_marker_location == -1:
            gene_name_descriptions.append((gene, "Entrez summary not found"))
            continue

        summary_start_index = string.find("<p>", gene_summary_marker_location)

        summary_end_index = string.find("</p>", summary_start_index)

        gene_description = string[summary_start_index+3:summary_end_index]

        gene_name_descriptions.append((gene_name, gene_description))


    return gene_name_descriptions
