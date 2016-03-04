#!/usr/bin/env sh

DATA_DIR=dammit/tests/test-data
TEST_FILE=pom.single.fa

dammit annotate $DATA_DIR/$TEST_FILE
dammit annotate --evalue 10.0 $DATA_DIR/$TEST_FILE -o $TEST_FILE.dammit.evalue10

cp $TEST_FILE.dammit/$TEST_FILE.dammit.fasta $DATA_DIR/
cp $TEST_FILE.dammit/$TEST_FILE.dammit.gff3 $DATA_DIR/

cp $TEST_FILE.dammit.evalue10/$TEST_FILE.dammit.fasta $DATA_DIR/$TEST_FILE.dammit.fasta.evalue10
cp $TEST_FILE.dammit.evalue10/$TEST_FILE.dammit.gff3 $DATA_DIR/$TEST_FILE.dammit.gff3.evalue10
