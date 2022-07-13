#!/bin/bash

grep '>' $1 | sed 's/.*name=//g; s/;.*parent=/_/g; s/;.*//g; s/^/>/g' > formated_IDs.lst

fasta_file=$1
output=$(sed 's/.fasta//g; s/.fa//g; ' <<< "$1" | sed 's/$/_tidy.fa/')

#function renaming() {

python3 - << END
fasta= open('$fasta_file')
newnames= open('formated_IDs.lst')
newfasta= open('$output', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
END

rm formated_IDs.lst
