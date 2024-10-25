#create dbs and contigs table
prepare () {
    rm test/contigs
    rm -r test/blastn_base
    rm -r test/output
    
    bash prep_db.sh \
    -n test/blastn_base/true_base \
    -c test/contigs \
    -t test/.tmp \
    test/fasta_base/true_base/*
    
    bash prep_db.sh \
    -n test/blastn_base/false_base_1 \
    -c test/contigs \
    -t test/tmp \
    test/fasta_base/false_base_1/*
    
    bash prep_db.sh \
    -n test/blastn_base/false_base_2 \
    -c test/contigs \
    -t test/tmp \
    test/fasta_base/false_base_2/*
}

prepare

#exec
python pipeline.py \
    -i test/test.fna \
    -o test/output \
    -tb test/blastn_base/true_base \
    -fb test/blastn_base/false_base_1 \
    test/blastn_base/false_base_2 \
    -c test/contigs \
    --primer3 /mnt/tank/scratch/dsmutin/tools/primer3/src/primer3_core