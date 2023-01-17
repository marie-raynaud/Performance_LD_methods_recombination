########################################################################
##################### LDHELMET - POST_TO_TEXT STEP #####################
########################################################################

#!/bin/bash

## Input file : output.post
## Output file : output.txt

time ./ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o output.txt output.post
