import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Concatenate left and right primer with N")

parser.add_argument("-i", "--input",
                    required=True,
                    help="Primers fasta input. Pairs specified as 'LEFT' and 'RIGHT', all primers paired")

parser.add_argument("-t", "--tmp",
                    required=False, default="fasta_table.tmp",
                    help="tmp file for fasta-table")

parser.add_argument("-o", "--output",
                    required=True,
                    help="output")

parser.add_argument("--N",
                    required=False,
                    default=100, type=int,
                    help="Number of N to merge")

args = parser.parse_args()
script_path = os.path.dirname(os.path.realpath(__file__)) + "/"

# Main
fasta2table = "bash " + script_path +\
    "fasta2table.sh -i " + args.input +\
    " -o " + args.tmp

subprocess.run(fasta2table, shell=True)

fasta_table = open(args.tmp, "r")
output_fasta = open(args.output, "w")
sepN = "N"*args.N

# Iterate lines
for i, line in enumerate(fasta_table):
    dev = i % 4
    if dev == 0:
        inp_name = line[:-6]
    elif dev == 1:
        left = line[:-1] # remove \n
    elif dev == 2:
        # break if not exact
        if hash(inp_name) != hash(line[0:(len(line)-7)]):
            LookupError("Error on line "+str(i)+": check names in fasta file")
            break
    else:
        right = line
        output_fasta.write(inp_name+"\n"+left+sepN+right)

output_fasta.close()
fasta_table.close()

print("Primers merged")
