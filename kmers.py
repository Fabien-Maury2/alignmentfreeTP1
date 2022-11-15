def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def nucl_compl(x):
  r = None
  if x == "A" or x == "a":
    r = "T"
  elif x == "T" or x == "t":
    r = "A"
  elif x == "C" or x == "c":
    r = "G"
  elif x == "G" or x == "g":
    r = "C"
  return(r)    

def binarize(seq):
  hash = ""
  for i in range (len(seq)):
    if seq[i] == "G" or seq[i] == "g":
      hash += "11"
    elif seq[i] == "T" or seq[i] == "t":
      hash += "10"
    elif seq[i] == "C" or seq[i] == "c":
      hash += "01"
    else:
      hash += "00"
  return(hash)

def encode(seq):
  if len(seq) > 1:
    hash1 = binarize(seq)
    hash1 = hash1[::-1]
    total1 = 0
    for i in range(len(hash1)):
      total1 += 2**i * int(hash1[i])
    cano = ""
    for i in range(len(seq)):
      cano += nucl_compl(seq[i])
    cano = cano[::-1]
    hash2 = binarize(cano)
    hash2 = hash2[::-1]
    total2 = 0
    for i in range(len(hash2)):
      total2 += 2**i * int(hash2[i])
    r = min(total1,total2)
  else:
    if seq == "C" or seq == "c":
      r = 1
    elif seq == "T" or seq == "t":
      r = 2
    elif seq == "G" or seq == "g":
      r = 3
    else:
      r = 0
  return(r)

def stream_kmers(seq, k):
  list_kmer = []
  kmer = 0
  for i in range(k-1):
    mask = (1 << (k-1)*2) - 1
    kmer = kmer << 2
    kmer += encode(seq[i])

  for nucl in seq[(k-1):]:
    kmer = kmer & mask
    kmer = kmer << 2
    kmer = kmer + encode(nucl)
    list_kmer.append(kmer)
  return(list_kmer)
