function usage(){
    echo "OrthoCode takes a Orthogroups tsv files produced by OrthoFinder 2+"
    echo
    echo "The files 'Orthogroups.GeneCount.tsv' and 'Orthogroups_UnassignedGenes.tsv'"
    echo "Are combined to produce 'Orthogroups.FullGeneCount.tsv' and 'Orthogroups.Presence.tsv'"
    echo "The resulting presence/absence matrix is encoded for fast querying of orthogroups"
    echo "The orthogroups are encoded from a list of species provided by the user or inferred from the tables"
    echo "The encoded orthogroups can then be used for downstream analyses"
    echo
    echo "Options:"
    echo "  --tsv_dir         -> Directory containing protein sequences in fasta format (mandatory)"
    echo "  --counts_file     -> Sequence similarity method to use in orthofinder (mandatory)"
    echo "  --unassigned_file -> Filters BLAST results based on total- or effective-length [TRUE]/FALSE"
    echo "  --species_list    -> Calculates LCR based effective length [TRUE]/FALSE"
    echo "  --threads         -> Number of CPU threads to use"
    echo
    echo "Notes:"
    echo "  --tsv_dir and --counts_file/--unassigned_file options are mutually exclusive"
    echo
    echo "Examples:"
    echo "  OrthoCode.sh --tsv_dir /path/to/my/orthogroups/folder --species_list species_list.txt --threads 16"
    echo "  OrthoCode.sh --tsv_dir /path/to/my/orthogroups/folder --threads 16"
    echo "  OrthoCode.sh --counts_file /path/to/counts.tsv --unassigned_file /path/to/counts.tsv --species_list species_list.txt --threads 16"
	}
export -f usage