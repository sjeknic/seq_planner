# seq_planner

You may need to install tabulate for the script to work.

```
pip install tabulate
```

Copy and paste plasmid sequence into plasmid.txt.
Easiest is to remove existing plasmid.txt and make a new file.

```
rm plasmid.txt
vim plasmid.txt
```

Run script:

```
python find_plasmids.py
```

You will be asked to input the start and end of the range (in bp) and the script will find primers that can sequence that range.
The script assumes that the first base read from a primer will be 40 bp from the end of the primer (Buffer range). The script also assumes about 800 bases read from each primer (Read range).
In the future I want to change the functionality such that the script attempts multiple range values and allows the user to choose the one they want.

Currently, the script starts at the 5' end and finds only the read that uses the fewest primers. 
 
Primers are found using csv files from Benchling exports, therefore only primers made before the last export will be included. You can export data as an oligo spreadsheet using Benchling functions, however, if two folders are named the same (e.g. "My oligos"), only one will be included, so they must be exported one at a time and then moved to the same folder.

After a read, the user is given the option of removing oligos they do not want to use. The script will remove this oligo entirely and attempt to find a new read de novo, therefore mutliple oligos may change.

If the script fails to find a complete read, oligos that read part of the range will be presented, as well as the regions that the script expects will not be read. 
