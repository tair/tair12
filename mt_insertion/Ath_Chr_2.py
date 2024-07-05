#!/usr/bin/python3
"""
This program annotates MtDNA in Chromosome 2
Author: Asher Pasha
Date: April 2024
"""
import os
import subprocess
from BCBio import GFF
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

BLAST_PATH = "/home/asher/BAR/Blast/ncbi-blast-2.15.0+/bin/"
NUMT_START = 8654381
NUMT_END = 9294939


def download_sequences():
    Entrez.email = "asher.pasha@utoronto.ca"

    if not os.path.exists("input"):
        os.mkdir("input")

    if not os.path.exists("input/CP116281.2.fasta"):
        # Download Chr 2
        with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="CP116281.2") as handle:
            seq_record = SeqIO.read(handle, "fasta")

        with open("input/CP116281.2.fasta", "w") as handle:
            SeqIO.write(seq_record, handle, "fasta")

    if not os.path.exists("input/NC_037304.1.gbk"):
        # Download MtDNA
        with Entrez.efetch(db="nucleotide", rettype="gb", id="NC_037304.1") as handle:
            seq_record = SeqIO.read(handle, "gb")

        with open("input/NC_037304.1.gbk", "w") as handle:
            SeqIO.write(seq_record, handle, "gb")


def extract_sequences():
    mt_genes = []

    # There is only one record in this file
    # So, we can use read instead of parse
    mt_dna_record = SeqIO.read("input/NC_037304.1.gbk", "genbank")
    mt_dna_seq = mt_dna_record.seq

    for feature in mt_dna_record.features:
        # Extracting gene, tRNA, rRNA only for MtDNA
        # Skip nad1, nad2, nad5
        if feature.type == "gene" and feature.qualifiers["gene"][0] not in ["nad1", "nad2", "nad5"]:
            gene_id = feature.qualifiers["db_xref"][0].split(":")[1]
            header = "gi|" + gene_id + "|NC_037304.1|"
            # Adding GI id again show it shows up in GFF files for now
            description = feature.type + ": " + feature.qualifiers["gene"][0]

            # Note: This requires at least BioPython 1.78
            mt_record = SeqRecord(Seq(feature.extract(mt_dna_seq)), id=header, description=description)
            mt_genes.append(mt_record)

    if not os.path.exists("output"):
        os.mkdir("output")

    SeqIO.write(mt_genes, open("output/mt_genes.fasta", "w"), "fasta")

    # Copy Chr2 to output
    subprocess.call("cp input/CP116281.2.fasta output/CP116281.2.fasta", shell=True)


def extract_nad_genes():
    # Extract nad1, nad2, nad5

    mt_genes = []

    # Get the DNA sequence
    mt_dna_record = SeqIO.read("input/NC_037304.1.gbk", "genbank")
    mt_dna_seq = mt_dna_record.seq

    # nad1_1
    # complement(82937..83321)
    mt_gene = mt_dna_seq[82937 - 1 : 83321]
    mt_gene = mt_gene.reverse_complement()
    header = "gi|36335702|NC_037304.1|"
    description = "gene: nad1_1"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad1_2
    # 58315..59481
    mt_gene = mt_dna_seq[58315 - 1 : 59481]
    header = "gi|36335702|NC_037304.1|"
    description = "gene: nad1_2"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad1_3
    # 230304..234132
    mt_gene = mt_dna_seq[230304 - 1 : 234132]
    header = "gi|36335702|NC_037304.1|"
    description = "gene: nad1_3"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad2_1
    # 161832..163357
    mt_gene = mt_dna_seq[161832 - 1 : 163357]
    header = "gi|36335684|NC_037304.1|"
    description = "gene: nad2_1"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad2_2
    # complement(92819..98042)
    mt_gene = mt_dna_seq[92819 - 1 : 98042]
    mt_gene = mt_gene.reverse_complement()
    header = "gi|36335684|NC_037304.1|"
    description = "gene: nad2_2"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad5_1
    # 234352..236628
    mt_gene = mt_dna_seq[234352 - 1 : 236628]
    header = "gi|36335705|NC_037304.1|"
    description = "gene: nad5_1"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad5_2
    # 360107..360128
    mt_gene = mt_dna_seq[360107 - 1 : 360128]
    header = "gi|36335705|NC_037304.1|"
    description = "gene: nad5_2"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # nad5_3
    # 26817..28336
    mt_gene = mt_dna_seq[26817 - 1 : 28336]
    header = "gi|36335705|NC_037304.1|"
    description = "gene: nad5_3"
    mt_record = SeqRecord(Seq(mt_gene), id=header, description=description)
    mt_genes.append(mt_record)

    # Append
    SeqIO.write(mt_genes, open("output/mt_genes.fasta", "a"), "fasta")


def make_blast_db():
    make_blast_db_app = BLAST_PATH + "makeblastdb"
    cmd = make_blast_db_app + " -dbtype nucl -in CP116281.2.fasta -parse_seqids"

    os.chdir("output")
    subprocess.run(cmd, shell=True)
    os.chdir("..")


def run_blast():
    blastn = BLAST_PATH + "blastn"
    cmd = blastn + " -query mt_genes.fasta -db CP116281.2.fasta -out blastn_results.xml -evalue 0.001 -outfmt 5"

    os.chdir("output")
    subprocess.run(cmd, shell=True)
    os.chdir("..")


def parse_blast_output():
    # Here is what can be read for the BLAST XML
    # https://github.com/biopython/biopython/blob/master/Bio/SearchIO/BlastIO/__init__.py
    # Create Sequence record for GFF3
    # Documentation: https://biopython.org/wiki/GFF_Parsing

    # Only needed for GFF writer
    # I am open to suggestions.
    chr2 = SeqIO.read("input/CP116281.2.fasta", "fasta")
    seq = chr2.seq
    rec = SeqRecord(seq, id="Chr2")
    rec.features = []

    # Now parse Blast results
    blast_qresults = SearchIO.parse("output/blastn_results.xml", "blast-xml")

    for blast_qresult in blast_qresults:
        # Take the top hit of each query for now
        for hsp in blast_qresult.hsps:
            # Only process hsp in the numt region
            if hsp.hit_start > NUMT_START and hsp.hit_end < NUMT_END:
                if ((hsp.aln_span / blast_qresult.seq_len) > 0.95) and ((hsp.ident_num / hsp.aln_span) > 0.95):
                    # Now filter by Seq coverage and identity.
                    # Right now both are 95% coverage
                    top_feature = make_feature(blast_qresult, hsp)

                    # Check if a gene already exists at that location. If not add.
                    if any(hsp.hit_start == feature.location.start for feature in rec.features) and any(
                        hsp.hit_end == feature.location.end for feature in rec.features
                    ):
                        pass
                    else:
                        rec.features.append(top_feature)

    # Sort GFF based on start location
    rec.features.sort(key=lambda r: r.location.start)

    # Write the GFF file
    with open("output/MtDNA_in_Chr2.gff", "w") as handle:
        GFF.write([rec], handle)


def make_feature(blast_qresult, hsp):
    """
    Make GFF feature
    :param blast_qresult:
    :param hsp:
    :return:
    """
    description = blast_qresult.description
    description = description.replace("(", " ")
    description = description.replace(")", "")

    qualifiers = {"source": "Arabidopsis12", "ID": "AT2G00000", "other": [description]}

    feature = SeqFeature(
        FeatureLocation(hsp.hit_start, hsp.hit_end, hsp.hit_frame),
        type="nuclear_mt_pseudogene",
        qualifiers=qualifiers,
    )

    return feature


def main():
    # download_sequences()
    extract_sequences()
    extract_nad_genes()
    make_blast_db()
    run_blast()
    parse_blast_output()


if __name__ == "__main__":
    main()
