#!/bin/bash

# USAGE: ./extract_interpro.sh

echo "Extracting annotations from interproscan output file: $1"
a="$1"
b=$(echo $1 | perl -pe "s/(.+)\.([^\.]+)$/\1/")
echo "Extracting annotations from Pfam, TIGRFAM, ProSiteProfiles..."
grep -e 'Pfam' -e 'TIGRFAM' -e 'ProSiteProfiles' $1 | cut -f1,4-6,13-14 | sort | uniq | perl -pe "s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+)\t(.+)$/\1 \|\| \2/" | perl -pe "s/^([^\t]+\t[^\t]+\t[^\t]+\t)([^\t]+)\n$/\1\2\t\n/" | sort | uniq > "$b.interpro_filtered.parsed"

exit;
