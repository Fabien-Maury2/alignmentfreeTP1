from loading import load_directory
from kmers import stream_kmers, kmer2str



def similarity(A, inter, B):
    simB = len(inter) / (len(B) + len(inter))
    simA = len(inter) / (len(A) + len(inter))
    return((simA, simB))


def jaccard(A, inter, B):
    r = len(inter) / (len(A) + len(B) + len(inter))
    return(r)

def shared_kmers(a,b,k): # returns a list containing 2 lists of all kmers of length k contained in file and b, respectively; and a list containing the shared elements of the 2 former
    ak = stream_kmers(a, k)
    bk = stream_kmers(b, k)
    inter = list(set(ak) & set(bk))
    return([ak, bk, inter])



if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---

            A, inter, B = shared_kmers(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
