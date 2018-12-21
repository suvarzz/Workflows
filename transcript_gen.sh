#!/bin/bash

# salmon location: ~/Softs/Salmon-latest_linux_x86_64/bin/

salmon index -t ~/Projects/Sources/Ensembl_r40/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa -i ~/Projects/Sources/transcripts_index --type quasi -k 31
