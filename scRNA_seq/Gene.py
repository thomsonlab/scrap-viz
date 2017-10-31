import urllib3 as urllib

Entrez_gene_summary_marker = " Entrez Gene Summary for "
gene_summary_marker = " Summary for "
gene_name_marker = "aliasMainName"
gene_cache = {}


def get_gene_summaries(genes):

    global gene_cache

    gene_name_descriptions = []

    http = urllib.PoolManager()

    for gene in genes:

        if gene in gene_cache and gene_cache[gene] != (gene, "No response"):
            gene_name_descriptions.append(gene_cache[gene])
            continue

        url = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" % gene
        response = http.request("GET", url)

        if response.status != 200:
            gene_name_descriptions.append((gene, "No response"))
            continue

        string = response.data.decode("UTF-8")

        if string.find("was not found") != -1:
            gene_cache[gene] = (gene, "Gene not found")
            gene_name_descriptions.append(gene_cache[gene])
            continue

        gene_name_marker_location = string.find(gene_name_marker)

        if gene_name_marker_location == -1:
            gene_name = gene
        else:
            name_start_index = string.find(">", gene_name_marker_location)
            name_end_index = string.find("<", name_start_index)

            gene_name = string[name_start_index+1:name_end_index]

        Entrez_gene_summary_location = string.find(Entrez_gene_summary_marker)

        if Entrez_gene_summary_location == -1:
            gene_summary_marker_location = string.find(gene_summary_marker)

            if gene_summary_marker_location == -1:
                gene_cache[gene] = (gene_name, "Gene summary not found")
                gene_name_descriptions.append(gene_cache[gene])
                continue
        else:
            gene_summary_marker_location = Entrez_gene_summary_location

        summary_start_index = string.find("<p>", gene_summary_marker_location)

        summary_end_index = string.find("</p>", summary_start_index)

        gene_description = string[summary_start_index+3:summary_end_index]
        gene_description = " ".join(gene_description.split())

        gene_cache[gene] = (gene_name, gene_description)
        gene_name_descriptions.append(gene_cache[gene])

    return gene_name_descriptions
