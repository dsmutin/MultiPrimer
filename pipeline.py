#
# python pipeline.py -i /mnt/tank/scratch/dsmutin/misc/primers2primer/data/it.fna -o /mnt/tank/scratch/dsmutin/misc/primers2primer/test --primer3 /mnt/tank/scratch/dsmutin/tools/primer3/src/primer3_core

# Pipeline: ----

# 1. Initial set generation
#
# < evolutionary algorithm >
# 2. blastn
# 3. multimapping detection and filtering
# 4. primers matching
# 5*. mutations (if not final)
# </ evolutionary algorithm >
#
# 6. output

# 0. Imports ----

import os
import argparse
from Bio import SeqIO
import subprocess
import numpy as np
import pandas as pd
import re
import random

# 0. Functions
# Template generation


def primer_template(fasta_file,
                    PRIMER_PICK_PRIMER,
                    PRIMER_OPT_SIZE,
                    PRIMER_MIN_SIZE,
                    PRIMER_MAX_SIZE,
                    PRIMER_PRODUCT_SIZE_RANGE,
                    PRIMER_NUM_RETURN):
    # Читаем последовательности из fasta файла
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    output = []

    # Проходим по каждой последовательности
    for record in sequences:
        # Получаем имя файла без расширения и модификатора
        seq_id = os.path.basename(fasta_file).replace(".exon.mod.fna", "")

        # Формируем шаблон
        template = f"""SEQUENCE_ID={seq_id}_{record.id}
SEQUENCE_TEMPLATE={record.seq}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER={PRIMER_PICK_PRIMER}
PRIMER_PICK_RIGHT_PRIMER={PRIMER_PICK_PRIMER}
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE={PRIMER_OPT_SIZE}
PRIMER_MIN_SIZE={PRIMER_MIN_SIZE}
PRIMER_MAX_SIZE={PRIMER_MAX_SIZE}
PRIMER_PRODUCT_SIZE_RANGE={PRIMER_PRODUCT_SIZE_RANGE}
PRIMER_NUM_RETURN={PRIMER_NUM_RETURN}
PRIMER_EXPLAIN_FLAG=1
="""
        output.append(template)
    return "\n".join(output)

# Primer3 out parse to fasta


def parse_primer3_output(primer3_file):
    primers = []
    sequence_id = None

    with open(primer3_file, "r") as file:
        for line in file:
            line = line.strip()

            # Извлекаем SEQUENCE_ID
            if line.startswith("SEQUENCE_ID="):
                sequence_id = line.split("=")[1]

            # Извлекаем левый праймер
            if line.startswith("PRIMER_LEFT_") and "SEQUENCE" in line:
                primer_num = line.split("_")[2]
                primer_seq = line.split("=")[1]
                primers.append((sequence_id, primer_num, "LEFT", primer_seq))

            # Извлекаем правый праймер
            if line.startswith("PRIMER_RIGHT_") and "SEQUENCE" in line:
                primer_num = line.split("_")[2]
                primer_seq = line.split("=")[1]
                primers.append((sequence_id, primer_num, "RIGHT", primer_seq))

    return primers


def write_fasta(primers, output_file):
    fasta = open(output_file, "w")
    for primer in primers:
        sequence_id, primer_num, side, sequence = primer
        header = f">{sequence_id}_{primer_num}_{side}"
        fasta.write(f"{header}\n{sequence}\n")

    write_fasta(primers, output_file)

# Misc


def out_dir(iter):
    if args.output_tmp == "":
        return args.output + "/.tmp/" + str(iter) + "/"
    else:
        return args.output_tmp + str(iter) + "/"


def pairing(x):
    clear_primer = re.sub(r"(RIGHT)|(LEFT)", "", string=x)
    return [clear_primer+"RIGHT", clear_primer+"LEFT"]


def mutate_seq(x):
    indelrate=0.1
    nucleotide_code = "ATGCUWSMKRYBDHVN"
    rstate = random.random()
    if rstate <= args.mutation_rate:
        if rstate <= args.mutation_rate*indelrate:
            if rstate <= args.mutation_rate*indelrate/2:
                return random.choice(nucleotide_code)+x
            return ""
        return random.choice(nucleotide_code)
    else:
        return x

# 0. Argparsing ----
description = "Generation of primers based on fasta-files and blastn databases.\n\nTo use it, select one reference file to generate the initial primer set; blastn base to check primer universality and cut off multimapping; blastn bases to remove non-specific primers\n\nRequires primer3 and blastn pre-installed"

parser = argparse.ArgumentParser(description=description)

# Main
parser.add_argument("-i", "--input",
                    required=True,
                    help="Input FASTA file for generation. Primers are generated for different contigs separatly. Only gene-coding regions recommended (.fna)")

parser.add_argument("-tb", "--true_base",
                    required=True,
                    help="Input blastn database path for primer adjusting")

parser.add_argument("-fb", "--false_base",
                    required=True,
                    nargs="*",
                    help="Input blastn database path for non-specific testing. Wildcards are not accepted")

parser.add_argument("-c", "--contig_table",
                    required=True,
                    help=".tsv table with blast db information")

parser.add_argument("-o", "--output",
                    required=True,
                    help="Output path")

parser.add_argument("-t", "--threads",
                    required=False,
                    default="1",
                    help="number of threads")

parser.add_argument("-ot", "--output_tmp",
                    default="",
                    help="Output .tmp dicrectory path for calculations and data processing. .tmp in output directory as default")

# Evolutionary algoritm
parser.add_argument("-N", "--iterations",
                    default=5, type=int,
                    help="Maximum iterations of evolutionary algorithm. 100 by default")

parser.add_argument("-T", "--top",
                    default=10, type=int,
                    help="Top primers to mutate and use in next generation")

parser.add_argument("-M", "--mutation_rate",
                    default=0.05, type=int,
                    help="Mutation probability per position of primer")

parser.add_argument("-S", "--set_size",
                    default=10, type=int,
                    help="Size of mutated primers per primer")

parser.add_argument("-A", "--append",
                    default=True, type=bool,
                    help="Append best primers to array in evolutionary algoritm")

# Exec
parser.add_argument("--primer3",
                    required=False,
                    default="primer3",
                    help="primer3_core path or command to exec. 'primer3' as default")

parser.add_argument("--blastn",
                    required=False,
                    default="blastn",
                    help="blastn path or command to exec. 'blastn' as default")

parser.add_argument("--add_set",
                    required=False,
                    default=None,
                    nargs="*",
                    help="file to set of primers to append to initial primer3 generation. empty by default")

# Primer3 template
parser.add_argument("--PRIMER_PICK_PRIMER",
                    default=10,
                    help="primer3 template option. Number of primers to pick")

parser.add_argument("--PRIMER_NUM_RETURN",
                    default=10,
                    help="primer3 template option. initial set size per gene")

parser.add_argument("--PRIMER_OPT_SIZE",
                    default=25,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_MIN_SIZE",
                    default=15,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_MAX_SIZE",
                    default=30,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_PRODUCT_SIZE_RANGE",
                    default="100-1000",
                    help="primer3 template option. 2 values sepatated by '-'")

# Blastn template
parser.add_argument("--word_size",
                    default="11",
                    help="blastn template option")

parser.add_argument("--reward",
                    default="3",
                    help="blastn template option")

parser.add_argument("--penalty",
                    default="-3",
                    help="blastn template option")

parser.add_argument("--gapopen",
                    default="6",
                    help="blastn template option")

parser.add_argument("--gapextend",
                    default="3",
                    help="blastn template option")

parser.add_argument("--evalue",
                    default="1",
                    help="blastn template option")

# primer_check template
parser.add_argument("--max_mismatch",
                    default="5",
                    help="primer_check template option. maximum avialable mismatch")

parser.add_argument("--multimap_max",
                    default="1",
                    help="primer_check template option. maximum multimapped hits")

parser.add_argument("--negative_max",
                    default="0",
                    help="primer_check template option. maximum negative hits")

parser.add_argument("--min_ident",
                    default="70",
                    help="primer_check template option. minimal identity, percent")

args = parser.parse_args()


# 1. Initial set generation ----
print("\n---- MultiPrimer v.0.2 ----\n")
print("Arguments passed")

script_path = os.path.dirname(os.path.realpath(__file__)) + "/scripts/"
os.makedirs(out_dir(0), exist_ok=True)

# Make uniline fasta
uniline = "bash " + script_path + "uniline_fa.sh"
uniline += " -i " + args.input
uniline += " -o " + out_dir(0) + "input.fa"

subprocess.run(uniline, shell=True)
print("Input fasta parsed")

# Template generation
primer_temp = primer_template(
    out_dir(0) + "input.fa",
    args.PRIMER_PICK_PRIMER,
    args.PRIMER_OPT_SIZE,
    args.PRIMER_MIN_SIZE,
    args.PRIMER_MAX_SIZE,
    args.PRIMER_PRODUCT_SIZE_RANGE,
    args.PRIMER_NUM_RETURN)

template = open(out_dir(0)+"template", "w")
template.writelines(primer_temp)
template.close()

# Primer3 exec
primer3 = args.primer3 + " " + \
    out_dir(0) + "template" + " --output " + out_dir(0) + "output.p3"

subprocess.run(primer3, shell=True, executable="/bin/bash")
print("Primer3 done")

# Parse 2 fasta
primers = parse_primer3_output(out_dir(0) + "output.p3")
fasta = open(out_dir(0) + "output.fa", "w")
for primer in primers:
    sequence_id, primer_num, side, sequence = primer
    header = f">{sequence_id}_{primer_num}_{side}"
    fasta.write(f"{header}\n{sequence}\n")
fasta.close()

# Add primers 2 fasta
if args.add_set is not None:
    add_fasta = "cat " + args.add_set + " >> " + out_dir(0) + "output.fa"
    subprocess.run(add_fasta, shell=True, executable="/bin/bash")

# < evolutionary algorithm >

# blastn command
blastn = args.blastn + " -num_threads " + \
    args.threads + " -outfmt '6 qseqid sseqid evalue sstart send ppos mismatch' " + \
    " -word_size " + args.word_size + \
    " -reward " + args.reward + \
    " -penalty " + args.penalty + \
    " -gapopen " + args.gapopen + \
    " -gapextend " + args.gapextend + \
    " -evalue " + args.evalue

# primer_check command
primer_check = "bash " + script_path + "/primer_check.sh" + \
    " -p " + script_path + "/primer_filt.py" + \
    " -d " + args.contig_table + \
    " -m " + str(args.top) + \
    " -e " + args.max_mismatch + \
    " -i " + args.min_ident + \
    " -a " + args.multimap_max + \
    " -b " + args.negative_max

for iter in range(1, args.iterations+1):
    print("\nIteration", iter, "----")
    os.makedirs(out_dir(iter), exist_ok=True)

    # 2. blastn ----
    blastn_iter = blastn + " -query " + out_dir(iter-1) + "output.fa"

    # true base
    blastn_db = blastn_iter + " -db " + args.true_base + \
        " > " + out_dir(iter) + "positive_hits.tsv"
    subprocess.run(blastn_db, shell=True)

    print("Positive hits counted")

    # false bases
    for db_neg in args.false_base:
        blastn_db = blastn_iter + " -db " + db_neg + \
            " >> " + out_dir(iter) + "negative_hits.tsv"
        subprocess.run(blastn_db, shell=True)

    print("Negative hits counted")

    # 3. multimapping detection and filtering ----
    primer_check_iter = primer_check + \
        " -o " + out_dir(iter) + "clear_hits.tsv" + \
        " -r " + out_dir(iter) + "primer_check/" + \
        " -t " + out_dir(iter) + "positive_hits.tsv " + \
        out_dir(iter) + "negative_hits.tsv"

    subprocess.run(primer_check_iter, shell=True)

    # 4. primers matching ----
    try:
        primer_out = pd.read_table(out_dir(iter) + "clear_hits.tsv",
                                sep=' ', header=None)
    except:
        InterruptedError("Empty file after filtration, try to use other primer_check properties and review false databases")
    
    primer_vals = primer_out.iloc[:, 0].value_counts()

    primer_list = list(set(primer_out.iloc[:, 0]))

    paired_primers = []
    for _ in primer_list.copy():
        cp = pairing(_)
        paired_primers += [cp[0]] + [cp[1]]

    primer_list = list(set(paired_primers))
    primer_list_hash = [hash(_) for _ in primer_list]

    print("Maximum hits:", primer_vals.iloc[0])
    print("Mean hits:", round(sum(primer_vals)/len(primer_vals), 1))

    # grep in primers.fa from previous iter
    fasta = open(out_dir(iter-1) + "output.fa", "r")
    seqs = {}
    for iter_line, line in enumerate(fasta):
        if iter_line % 2 == 0:
            if np.isin(hash(line[1:-1]), primer_list_hash):
                keep = True
                line_name = line[1:-1]
            else:
                keep = False
        else:
            if keep:
                seqs[line_name] = line[:-1]
    fasta.close()

    # 5*. mutations ----
    if iter != args.iterations:  # (if not final)
        if args.append:
            seqs_mutated = seqs.copy()
        else:
            seqs_mutated = dict()
        for seqs_unique in seqs.keys():
            for seqs_iter in range(args.set_size):
                init_seq = seqs[seqs_unique]
                mutated_seq = init_seq
                while init_seq == mutated_seq:
                    mutated_seq = "".join([mutate_seq(_) for _ in init_seq])
                mseq="I" + str(iter)+"N"+str(seqs_iter)+"_"
                seqs_mutated[mseq+seqs_unique] = mutated_seq

        fasta = open(out_dir(iter)+"output.fa", "w")
        for fname in seqs_mutated.keys():
            fasta.write(">" + fname +"\n"+ seqs_mutated[fname]+"\n")
        fasta.close()
        print("Done")

# </ evolutionary algorithm >
# 6. output ----
fasta = open(args.output + "/output.fa", "w")
for fname in seqs.keys():
    fasta.write(">"+fname+"\n" + seqs[fname]+"\n")
fasta.close()
print("Done")
