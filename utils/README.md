Steps:
1. make the executable using "make" command
2. run the executable: ./executable [input file] [output file]
3. take the output file from previous step and use the following command on it
   to produce a new output file:
    sed -e 's/\(.*\)/\L\1/' file > new file
4. then use the following command on that file to produce the final file:
    sed 's/  \+/ /g' file > final file
