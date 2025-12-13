
def motif_search(sequence,  motif):

    sequence = sequence.upper()
    motif = motif.upper()
    
    motif_positions = []
    motif_length = len(motif)
    
    for i in range(len(sequence) - motif_length + 1):
        if sequence[i:i+motif_length] == motif:
            motif_positions.append(i + 1)
    
    return motif_positions