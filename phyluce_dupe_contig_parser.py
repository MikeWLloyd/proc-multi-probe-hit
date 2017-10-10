from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import align_info
#customized to output consensus with IUPAC ambigious bases. 
import sys
import re
import os
from collections import defaultdict
import argparse
from shutil import copyfile

def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Find, align, generate consensus, and add the 'multiple-contigs-to-probe' contigs to the master fasta file."""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory in which to store the output files"""
    )
    parser.add_argument(
        "--contigs",
        required=True,
        type=is_dir,
        action=FullPaths,
        default=None,
        help="""Directory containing assembled contigs"""
    )
    parser.add_argument(
        "--dupes",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The path and name of the '--keep-duplicates' file output during the 'match_to_probe' phyluce script step."""
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The path and name of the FASTA file output during the 'get_fasta' phyluce script step."""
    )
    parser.add_argument(
        "--conf",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The path and name of the CONF file output during the 'get_match_counts' phyluce script step."""
    )
    parser.add_argument(
        "--incomplete",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The path and name of the INCOMPLETE file output during the 'get_match_counts' phyluce script step."""
    )
    return parser.parse_args()

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  
        # False when `fasta` is empty, i.e. wasn't a FASTA file
#####

def printoutincomplete(taxa, loci, d):
#taxa : unique list of taxa
#loci : unique list of loci
#d : "taxon - locus" dictionary.

# Print incomplete conf file and incomplete file for use in phyluce pipeline.
# Taxa, loci are unique lists, d is "taxon - locus" dictionary.
    taxa = sorted(taxa)
    loci = sorted(loci)
    incomplete_conf = open(args.output + '/incomplete_matrix.conf', 'w')

    incomplete_conf.write("[Organisms]\n")
    incomplete_conf.write("\n".join(taxa))
    incomplete_conf.write("\n[Loci]\n")
    incomplete_conf.write("\n".join(loci))

    incomplete_matrix = open(args.output + '/incomplete_matrix.incomplete', 'w')

    for t in d.items():
        incomplete_matrix.write("[" + str(t[0]) + "]\n")
        incomplete_matrix.write("\n".join(d[t[0]]))
        incomplete_matrix.write("\n")
    incomplete_matrix.close()
    incomplete_conf.close()
######


def generateconseq():

    dupes = open(args.dupes,"r")

    if is_fasta(args.fasta):
        pass
    else:
        sys.exit("Found " + args.fasta + " but it does not appear to be a FASTA file. You should check that one.")

    copyfile(args.fasta, args.output + '/amended_incomplete_matrix.fasta')    
    master_fasta = open(args.output + '/amended_incomplete_matrix.fasta', 'a')

    multi_probe = False
    added_loci = []
    a = defaultdict(list)
    for aline in dupes:
        values = aline.rstrip('\n').split(":")
        if (re.search("contigs hitting multiple probes", values[0])):
            multi_probe = False
        if (re.search("probes hitting multiple contigs", values[0])):
            multi_probe = True
            taxon = values[0].split(" ")
            taxon = taxon[0].replace("[","")
            print("Working on {}".format(taxon))
            contig_file = args.contigs + '/' + taxon + ".contigs.fasta"
            if is_fasta(args.contigs + '/' + taxon + ".contigs.fasta"):
                pass
            else:
                sys.exit("Found " + args.contigs + '/' + taxon + ".contigs.fasta" + " but it does not appear to be a FASTA file. You should check that one.")

        if multi_probe and (values[0] != '') and not re.search("probes hitting multiple contigs", values[0]):
            #values[0] = locus

            added_loci.append(values[0])
            a[taxon].append(values[0])
            #build list of 'new' loci, and associate those loci with correct taxon. 
            #these are used for re-building the .conf and .incomplete files downstream. 

            contigs = set(values[1].split(", "))
            records = (r for r in SeqIO.parse(contig_file, "fasta") if r.id in contigs)
            count = SeqIO.write(records, "temp.file", "fasta")
            mafft_cline = MafftCommandline(input="temp.file")
            stdout, stderr = mafft_cline()
            outfile = args.output + "/alignments/" + taxon + "_" + values[0] + ".fasta"
            with open(outfile, "w") as handle:
                handle.write(stdout)
            #Find the contigs we are interested in and write them into a temporary file, run a mafft alignment on that, and write mafft output to an outfile. 

            align = AlignIO.read(outfile, "fasta")
            summary_align = align_info.SummaryInfo(align)
            con_alignment = summary_align.iupac_consensus()
            record = SeqRecord(Seq(str(con_alignment)), id=values[0] + "_" + taxon + " |" + values[0], description="", name="")
            count = SeqIO.write(record, master_fasta, "fasta")
            #Take the output from the mafft alignment, and get the consensus alignment with IUPAC ambigious bases.
            #During this, rename the fasta header to reflect phyluce convention, and then append the new sequence to the 'master-fasta' file generated in phyluce 
            #during the 'get_fasta' script step. 

    os.remove("temp.file")
    dupes.close()

    added_loci = list(set(added_loci))
    return(added_loci, a)
###


def main(args):

    if not os.path.exists(args.output):
        os.makedirs(args.output)
        os.makedirs(args.output + '/alignments')
    else: 
        sys.exit("Output directory exists. Exiting to avoid overwriting files.")

    (added_loci, a) = generateconseq()

    config = open(args.conf,"r")

    header = config.readline()
    looping = False
    loci = []
    taxa = []
    for aline in config:
        aline = aline.rstrip('\n')
        if (re.search("[Loci]", aline) and looping == False):
            looping = True
        elif looping: 
            loci.append(aline)
        else: 
            taxa.append(aline)

    taxa = list(set(taxa))
    loci = list(set(loci))

    new_loci = set(loci) | set(added_loci)
    #unqiue list of loci (original and added): Union of the sets. 
    difference = len(new_loci)-len(loci)
    print("Found {} total loci, which is {} more than in the original {} file".format(len(new_loci), difference, args.conf))

    incomp = open(args.incomplete,"r")
    b = defaultdict(list)
    for aline in incomp:
        aline = aline.rstrip('\n')
        if (re.search("\[", aline)):
            taxon = aline.replace('[', '').replace(']', '')
        else: 
            b[taxon].append(aline)

    for t in b.items():
        A = (a[t[0]])
        B = (b[t[0]])

        missing_loci = set(B)-set(A)
        #difference the sets. 
        b[t[0]] = sorted(missing_loci)

    printoutincomplete(taxa, new_loci, b)


if __name__ == "__main__":
    args = get_args()
    main(args)
