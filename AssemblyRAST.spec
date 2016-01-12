/*
A KBase module: AssemblyRAST
This sample module contains one small method - filter_contigs.
*/

module AssemblyRAST {
    /*
        A string representing a ContigSet id.
    */
    typedef string contigset_id;

    /*
        A string representing a workspace name.
    */
    typedef string workspace_name;

    typedef structure {
        workspace_name workspace;
        contigset_id contigset_id;
        int min_length;
    } FilterContigsParams;

    /* 
        The workspace ID for a ContigSet data object.
        @id ws KBaseGenomes.ContigSet
    */
    typedef string ws_contigset_id;

    typedef structure {
        ws_contigset_id new_contigset_ref;
        int n_initial_contigs;
        int n_contigs_removed;
        int n_contigs_remaining;
    } FilterContigsResults;
	
    /*
        Filter contigs in a ContigSet by DNA length
    */
    funcdef filter_contigs(FilterContigsParams params) returns (FilterContigsResults) authentication required;
};