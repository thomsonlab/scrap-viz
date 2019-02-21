import os
import urllib3 as urllib
import json

ENTREZ_GENE_SUMMARY_MARKER = " Entrez Gene Summary for "
GENE_SUMMARY_MARKER = " Summary for "
GENE_NAME_MARKER = "aliasMainName"


class Gene_Metadata:

    def __init__(self, dir_path):

        self._gene_cache = {}
        self._gene_dir_path = dir_path

        self.load()

    def load(self):

        gene_cache_path = os.path.join(self._gene_dir_path, "gene_cache.json")

        if os.path.isfile(gene_cache_path):
            with open(gene_cache_path, "r") as gene_cache_file:
                self._gene_cache = json.load(gene_cache_file)

    def save(self):

        gene_cache_path = os.path.join(self._gene_dir_path, "gene_cache.json")

        with open(gene_cache_path, "w") as gene_cache_file:
            json.dump(self._gene_cache, gene_cache_file, indent=4)

    def get_gene_summaries(self, genes):

        gene_name_descriptions = []

        for _ in genes:
            gene_name_descriptions.append("No response")

        return gene_name_descriptions

        http = urllib.PoolManager()

        added_gene = False

        for gene in genes:

            if gene in self._gene_cache and self._gene_cache[gene] != \
                    (gene, "No response"):
                gene_name_descriptions.append(self._gene_cache[gene])
                continue

            url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" % gene
            response = http.request("GET", url)

            if response.status != 200:
                gene_name_descriptions.append((gene, "No response"))
                continue

            added_gene = True

            string = response.data.decode("UTF-8")

            if string.find("was not found") != -1:
                self._gene_cache[gene] = (gene, "Gene not found")
                gene_name_descriptions.append(self._gene_cache[gene])
                continue

            gene_name_marker_location = string.find(GENE_NAME_MARKER)

            if gene_name_marker_location == -1:
                gene_name = gene
            else:
                name_start_index = string.find(">", gene_name_marker_location)
                name_end_index = string.find("<", name_start_index)

                gene_name = string[name_start_index + 1:name_end_index]

            Entrez_gene_summary_location = string.find(
                ENTREZ_GENE_SUMMARY_MARKER)

            if Entrez_gene_summary_location == -1:
                gene_summary_marker_location = string.find(GENE_SUMMARY_MARKER)

                if gene_summary_marker_location == -1:
                    self._gene_cache[gene] = \
                        (gene_name, "Gene summary not found")
                    gene_name_descriptions.append(self._gene_cache[gene])
                    continue
            else:
                gene_summary_marker_location = Entrez_gene_summary_location

            summary_start_index = string.find("<p>",
                                              gene_summary_marker_location)

            summary_end_index = string.find("</p>", summary_start_index)

            gene_description = string[summary_start_index + 3:summary_end_index]
            gene_description = " ".join(gene_description.split())

            self._gene_cache[gene] = (gene_name, gene_description)
            gene_name_descriptions.append(self._gene_cache[gene])

        if added_gene:
            self.save()

        return gene_name_descriptions
