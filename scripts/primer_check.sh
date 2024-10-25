# Usage
# bash primers.sh -t TRUE_MAPS.tsv -o output.tsv -p path_to_primer_filter.py -d contigs [-r tmp_dir] -m max_out *FALSE_MAPS_1* *FALSE_MAPS_2* ...
# carefully: f can not accept wildcard. 

# contig table: contig_name \t genome_name
# prepare: 
# cat group_* | sed 's/.*>//g' | sed 's/ /\t/' | sed 's/ /_/' | sed 's/ .*//g' > contigs
# cat group_* | sed 's/:>/\t/g' | sed 's/ .*//g' | sed 's/.*\///g' > contigs

# bash primers.sh -t /mnt/tank/scratch/dsmutin/misc/primers2primer/data/primers/GCF_900099615.1.group_ac.blast.tsv -o /mnt/tank/scratch/dsmutin/misc/primers2primer/data/ac.tsv -p /mnt/tank/scratch/dsmutin/misc/primers2primer/primer_filt.py -d /mnt/tank/scratch/dsmutin/misc/primers2primer/contigs /mnt/tank/scratch/dsmutin/misc/primers2primer/data/primers/GCF_900088665.2.group_ac.blast.tsv /mnt/tank/scratch/dsmutin/misc/primers2primer/data/primers/GCF_900099615.1.group_ac.blast.tsv /mnt/tank/scratch/dsmutin/misc/primers2primer/data/primers/GCF_900243065.1.group_ac.blast.tsv

# read args
#shopt -s extglob 

while getopts t:o:p:d:r:m:e:i:a:b:*: flag
do
    case "${flag}" in
        t) true_file=${OPTARG};;
        o) output=${OPTARG};;
        p) path_to_filter_py=${OPTARG};;
        d) contig_table=${OPTARG};;
        r) tmp_dir=${OPTARG};;
        m) max_out=${OPTARG};;
        e) max_mismatch=${OPTARG};;
        i) min_ident=${OPTARG};;
        a) multimap_max=${OPTARG};;
        b) negative_max=${OPTARG};;
    esac
done
shift $(( OPTIND - 1 ))

if [ -z $max_mismatch ]; then
    max_mismatch=5
fi

if [ -z $multimap_max ]; then
    multimap_max=0
fi

if [ -z $negative_max ]; then
    negative_max=0
fi

if [ -z $min_ident ]; then
    min_ident=90 
fi

# functions
filter(){
    awk '$7 <= "'$max_mismatch'" && $6 >= "'$min_ident'"' $1
}

sort_unique() {
    negmax=$2
    filter $1 |\
        awk '{print $1}' |\
        sort |\
        uniq -c |\
        awk '{if ($1 >= "'$negmax'") print $2}'
        #-c | awk '{if ($1>2) print $2}'
}

count() {
    awk '{count[$1]++} END {for (word in count) print count[word], word}' $1 |\
    awk '{print $2 "\t" $1}'
}

# make tmp dir
if [ -z $tmp_dir ]; then
    tmp=./.tmp
else
    tmp=$tmp_dir
fi
mkdir -p $tmp

echo "Primer check initiated"

# count false mappers
for filter_files in "$@"; do
    bn=$(basename $filter_files)
    echo " - prepare false hits: ${bn%.*}"
    sort_unique $filter_files $negative_max >> $tmp/filter_files.tmp
done

echo " - count occurences"
#filter true mappers 
filter $true_file > $tmp/filtered.tmp

#unique hitters 
sort $tmp/filtered.tmp > $tmp/sorted_1.tmp
count $tmp/sorted_1.tmp | sort > $tmp/count.tmp
echo " --" $(wc -l < $tmp/count.tmp) primers sucsessfully mapped

#multimapping check 
echo " - multimapping check"

join $tmp/sorted_1.tmp $tmp/count.tmp | sort -k8 -nr > $tmp/amount_sorted.tmp

sort $tmp/filtered.tmp -k2 > $tmp/sorted_2.tmp
sort $contig_table -k2 > $tmp/sorted_contig.tmp
join -1 2 -2 2 $tmp/sorted_contig.tmp $tmp/sorted_2.tmp |\
    awk '{print $2 "---" $3}' > $tmp/sorted_with_genomes.tmp
count $tmp/sorted_with_genomes.tmp > $tmp/sorted_counts.tmp
sed 's/---/\t/g' $tmp/sorted_counts.tmp |\
    awk '{if ($3 > "'$multimap_max'") print $2}' >> $tmp/filter_files.tmp

# concat
sort $tmp/filter_files.tmp |\
    uniq > $tmp/remove.tmp

echo " - filtration"

# python exec
python $path_to_filter_py $tmp/amount_sorted.tmp $tmp/remove.tmp $output $max_out

# cleanup
if [ -z $tmp_dir ]; then
    echo "Clear ----"
    rm -r $tmp
fi