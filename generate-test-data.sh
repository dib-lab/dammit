#!/usr/bin/env sh

DATA_DIR=dammit/tests/test-data
TEST_FILE=pom.20.fa
TEST_NAME=pom.20
TEST_PEP=pep.fa

#
# Execute
#

# Default pipeline and params
dammit run  --busco-group saccharomycetes_odb10  --n-threads 4 annotate $DATA_DIR/$TEST_FILE 

# Full pipeline, -o
#dammit run --pipeline full annotate $DATA_DIR/$TEST_FILE -o $TEST_NAME.dammit.full

# Global evalue cutoff
dammit run  --busco-group saccharomycetes_odb10  --n-threads 4 annotate --global-evalue 10.0 $DATA_DIR/$TEST_FILE -o $TEST_NAME.dammit.evalue10

# Passing a user database
dammit run  --busco-group saccharomycetes_odb10  --n-threads 4 --pipeline quick annotate $DATA_DIR/$TEST_FILE --user-database $DATA_DIR/$TEST_PEP -o $TEST_NAME.dammit.udb 

# multiple user databases
dammit run --n-threads 4 --busco-group saccharomycetes_odb10 --pipeline quick annotate --user-database $DATA_DIR/pep.fa --user-database $DATA_DIR/odb_subset.fa $DATA_DIR/$TEST_FILE -o $TEST_NAME.dammit.multi-udb

# Passing multiple musco groups
dammit run --n-threads 4 --pipeline quick --busco-group bacteria_odb10 --busco-group saccharomycetes_odb10 annotate -o pom.256.dammit.busco-multi $DATA_DIR/pom.256.fa 

# Using --regex-rename
dammit run --busco-group saccharomycetes_odb10 --pipeline quick annotate --regex-rename "(?P<name>^[a-zA-Z0-9\.]+)" -o pom.20.dammit.regex dammit/tests/test-data/pom.20.fa

# Rename with backmapping
dammit run --busco-group saccharomycetes_odb10 --pipeline quick annotate $DATA_DIR/$TEST_FILE --no-rename -o $TEST_NAME.dammit.norename


#
# Copy the data
#

cp $TEST_NAME.dammit/$TEST_NAME.dammit.fasta $DATA_DIR/ 
cp $TEST_NAME.dammit/$TEST_NAME.dammit.gff3 $DATA_DIR/

cp $TEST_NAME.dammit.evalue10/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.dammit.evalue10.fasta
cp $TEST_NAME.dammit.evalue10/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.dammit.evalue10.gff3

cp $TEST_NAME.dammit.udb/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.udb.dammit.fasta
cp $TEST_NAME.dammit.udb/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.udb.dammit.gff3

cp pom.256.dammit.busco-multi/pom.256.dammit.gff3 $DATA_DIR/pom.256.dammit.busco-multi.gff3

cp $TEST_NAME.dammit.norename/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.dammit.norename.gff3
cp $TEST_NAME.dammit.norename/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.dammit.norename.fasta

#cp $TEST_NAME.dammit.full/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.dammit.fasta.full
#cp $TEST_NAME.dammit.full/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.dammit.gff3.full

#cp $TEST_NAME.dammit.nr/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.dammit.fasta.nr
#cp $TEST_NAME.dammit.nr/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.dammit.gff3.nr

#cp $TEST_NAME.dammit.norename/$TEST_NAME.dammit.fasta $DATA_DIR/$TEST_NAME.dammit.fasta.norename
#cp $TEST_NAME.dammit.norename/$TEST_NAME.dammit.gff3 $DATA_DIR/$TEST_NAME.dammit.gff3.norename


