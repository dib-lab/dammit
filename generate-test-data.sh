#!/usr/bin/env sh

DATA_DIR=dammit/tests/test-data
TEST_FILE=pom.20.fa
TEST_PEP=pep.fa


dammit run annotate $DATA_DIR/$TEST_FILE 
dammit run annotate $DATA_DIR/$TEST_FILE --pipeline full -o $TEST_FILE.dammit.full 
#dammit annotate $DATA_DIR/$TEST_FILE --nr -o $TEST_FILE.dammit.nr 
dammit run annotate --evalue 10.0 $DATA_DIR/$TEST_FILE -o $TEST_FILE.dammit.evalue10 
dammit run annotate $DATA_DIR/$TEST_FILE --pipeline quick --user-databases $DATA_DIR/$TEST_PEP -o $TEST_FILE.dammit.udb 
dammit run annotate $DATA_DIR/$TEST_FILE --no-rename -o $TEST_FILE.dammit.norename

cp $TEST_FILE.dammit/$TEST_FILE.dammit.fasta $DATA_DIR/ 
cp $TEST_FILE.dammit/$TEST_FILE.dammit.gff3 $DATA_DIR/

cp $TEST_FILE.dammit.evalue10/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.evalue10
cp $TEST_FILE.dammit.evalue10/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.evalue10

cp $TEST_FILE.dammit.udb/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.udb
cp $TEST_FILE.dammit.udb/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.udb

cp $TEST_FILE.dammit.full/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.full
cp $TEST_FILE.dammit.full/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.full

#cp $TEST_FILE.dammit.nr/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.nr
#cp $TEST_FILE.dammit.nr/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.nr

cp $TEST_FILE.dammit.norename/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.norename
cp $TEST_FILE.dammit.norename/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.norename


