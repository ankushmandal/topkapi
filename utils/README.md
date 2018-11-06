Steps:
1. make the executable using `make` command
2. run the executable: `./remv_bad_char [input file] [output file]`
3. take the output file from previous step and use the following command on it to produce the final file:
     ```sed '/^$/N;/^\n$/D' [output file] | sed 's/ \+/ /g' > [final file]```
