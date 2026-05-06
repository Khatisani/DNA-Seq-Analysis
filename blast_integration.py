from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def run_blast(sequence_string):

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