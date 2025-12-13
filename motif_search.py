
def motif_search(sequence,  motif):

    sequence = sequence.upper()
    motif = motif.upper()
    
    positions = []
    motif_length = len(motif)
        
    i = 0
    while True:
        i = sequence.find(motif, i)
        if i == -1:
            break
        positions.append(i + 1)
        i += 1 
        
    return positions