from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def run_blast(sequence_string):
    """
    Connects to the NCBI BLAST server to perform a nucleotide BLAST (blastn) search 
    against the nt database for a given sequence string.

    Args:
        sequence_string (str): The nucleotide sequence string to query.

    Returns:
        list: A list of up to 3 dictionaries containing match titles, e-values, and scores.
              Returns an empty list if the search fails or no matches are found.
    """

    print("\nConnecting to NCBI BLAST server. Please wait...")
    
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence_string)
        
        blast_record = NCBIXML.read(result_handle)
        
        matches = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                matches.append({
                    "title": alignment.title,
                    "e_value": hsp.expect,
                    "score": hsp.score
                })
        
        return matches[:3]
        
    except Exception as e:
        print(f"BLAST search could not be completed: {e}")
        return []