# Orthocode
## Dependencies
OrthoCode requires:
- [`bc`](https://www.gnu.org/software/bc/manual/html_mono/bc.html)
- [`gnu-parallel`](https://www.gnu.org/software/parallel/)
- [ete3](https://etetoolkit.org/)
- [pandas](https://pandas.pydata.org/)
- [re](https://docs.python.org/3/library/re.html)

## What is it?
Orthocode is a set of tools to process orthogroup tables (initially from OrthoFinder), for downstream analysis, in particular to help with calculations of gene gains and losses taking phylogenetic context into account.

## What does it do?
1. `orthocounts.py` takes the files `Orthogroups.GeneCount.tsv` and `Orthogroups_UnassignedGenes.tsv` to produce `Orthogroups.FullGeneCount.tsv` and `Orthogroups.Presence.tsv` which contain, as their names imply, gene counts for the total number of orthogroups (including unassigned genes) and a matrix of the distribution of all orthogroups across the species included in the orthology analysis
2. `encode_orthogroups.py` takes the file `Orthogroups.Presence.tsv` and encodes each orthogroup using a unique binary number that can be back translated into a presence/absence matrix for easy interrogation (grep-ing) of specific groupings.

## Example

Let's suppose you have the following Orthogroups table

|Orthogroup|Species_1    |Species_2    |Species_3           |Species_4|
|----------|-------------|-------------|--------------------|---------|
|OG0001    |sp1_a1,sp1_a2|sp2_a1       |sp3_a1,sp3_a2,sp3_a3|         |
|OG0002    |             |sp2_b1,sp2_b1|sp3_b1              |sp4_b1   |
|OG0003    |sp1_c1       |sp2_c1       |sp3_c1              |         |
|OG0004    |             |             |                    |sp4_d1   |


The counts table would look like this:

|Orthogroup|Species_1|Species_2|Species_3|Species_4|Total|
|----------|---------|---------|---------|---------|-----|
|OG0001    |2        |1        |3        |0        |6    |
|OG0002    |0        |2        |1        |1        |4    |
|OG0003    |1        |1        |1        |0        |3    |
|OG0004    |0        |0        |0        |1        |1    |


The presence/absence matrix would look like this:

|Orthogroup|Species_1|Species_2|Species_3|Species_4|
|----------|---------|---------|---------|---------|
|OG0001    |1        |1        |1        |0        |
|OG0002    |0        |1        |1        |1        |
|OG0003    |1        |1        |1        |0        |
|OG0004    |0        |0        |0        |1        |


And the encoded matrix would look like hthis

|Orthogroup|Species_1|Species_2|Species_3|Species_4|Code|
|----------|---------|---------|---------|---------|----|
|OG0001    |1        |2        |4        |0        |7   |
|OG0002    |0        |2        |4        |8        |14  |
|OG0003    |1        |2        |4        |0        |7   |
|OG0004    |0        |0        |0        |8        |8   |

The encoded matrix can then be massively reduced, suppose you have a 100 column table, because the presence/absence of each orthogroup is encoded in the last column, you can drop 98 columns and still retain all the distribution of each orthogroup, such that:

```
$ cat species_list
Species_1
Species_2
Species_3
Species_4

$ decode.sh 7 species_list
7   1110

$ decode.sh 14 species_list
14	0111
```

So any orthogroup with code 7 will be present in species 1,2 and 3, but absent in species 4, in the same way that any orthogroup with code 14 will be present in species 2-4, but absent in species 1. Codes are unique and this is where we can reduce the amount of data to process and also do fast querying on the orthogroup tables.

In the example above, we can use simple command line tools to get the most abundant types of orthogroups (in terms of distribution across the species list)

```
$ cut -f6 Orthogroups.Encoded.tsv | tail -n+2 | sort -n | uniq -c | sort -n | tail
```

The encoded orthogroups table can then be used together with regular expressions to get orthogroups with specific distributions

```
$ full_code_list=$(tail -n+2 Orthogroups.Encoded.tsv | sort -n | uniq )
$ for code in ${full_code_list}
do
    decode.sh ${code} species_list
done > og_codes.tsv
```

Suppose you want orthogroups that are present in species 1, may or may not be present in species 2 and 3, and absent in species 4. For such a combination you can use the regular expression `1..0`:
```
1 -> present in species 1
. -> present or absent in species 2
. -> present or absent in species 3
0 -> absent in species 4
```

Such that:

```
$ code_list=$(grep -w "1..0" og_codes.tsv | cut -f1)
$ og_list=$(grep -w "${code_list}" Orthogroups.Encoded.tsv | cut -f1)
```