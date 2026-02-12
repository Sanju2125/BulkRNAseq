set -euo pipefail

ls 6.stringtie/*.gtf > 6.stringtie/mergelist.txt

stringtie --merge \
    -p 8 \
    -G 4.reference/gencode.v44.annotation.gtf \
    -o 6.stringtie/stringtie_merged.gtf \
    6.stringtie/mergelist.txt
