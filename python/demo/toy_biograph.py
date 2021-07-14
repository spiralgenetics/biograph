"""
Toy code to replicate and understand the
BioGraph Format creation and query process

How to use this document:
Start at the workflow() to see what operations are run at a high level.
Look at methods' documentation for a detailed description of each step.
The code is simplified and should be followable by anyone acquainted with python.
Comment lines are plentiful.
"""
MAXLINES = -1 # -1 for off, otherwise an even integer
LETTERS = ["A", "C", "G", "T"]

def make_pat_reads():
    """
    Reads used as an example in the patent
    """
    return ["ACTG",
            "CTGC",
            "CCGC",
            "ATGC"]

def make_reads():
    """
    Reads used for the WhitePaper diagrams
    Diploid Sample's genome:
    GTTACTCTGA--------TAGCCAT
    GTTACTCTGACGCTAACGTAGCCAT
    """
    ret = []
    read_len = 6
    s1 = "GTTACTCTGATAGCCAT"
    for i in range(len(s1) - read_len):
        ret.append(s1[i:i + read_len])
    s2 = "GTTACTCTGACGCTAACGTAGCCAT"
    for i in range(len(s2) - read_len):
        ret.append(s2[i:i + read_len])
    return ret

def make_suffixes(sequences):
    """
    Compute the suffixes for each sequence
    """
    ret = []
    for seq in sequences:
        for pos in range(len(seq)):
            ret.append(seq[pos:])
    return ret

def lex_sort(sequences):
    """
    Sort the sequences lexicographically
    """
    sequences.sort()
    return sequences

def rm_prefixes(sequences):
    """
    Remove sequences which are prefixes of another sequence.
    This process can be done by linearly scanning the list and removing any element which
    is a prefix of or identical to the next element
    """
    i = 0
    while i < len(sequences) - 1:
        if sequences[i + 1].startswith(sequences[i]):
            del(sequences[i])
        else:
            i += 1
    return sequences

def rotation(sequences):
    """
    Perform first letter rotation.
     - For each sequence, duplicate the sequence.
     - Move the first element to the end, separated by a delimiter which is
       lexicographically smaller than all other elements of the alphabet.
     - Postfix the other non-rotated copy with an element that is greater than all other
       elements of the alphabet.

    !NOTE! this is what the figures say. But I cannot recreate the tables unless i
    use a smaller element for the non-rotated ends than the delimited rotation...
    @Nils - Which is correct?
    """
    ret = []
    for i in sequences:
        ret.append(i + "$")
        ret.append(i[1:] + "#" + i[0])
    return ret

def make_table(sequences):
    """
    Make a table of all non-rotated elements.
     - Column 1: seq -> The element itself.
     - Column 2: rotate -> List all of the rotated elements which are between the
       current element and the next non-rotated element.
    Then delete all rotated elements in the seq column
     Remove final suffix from all non-rotated elements.
    """
    def find_rotated(pos, found=None):
        """
        Find the last base of all rotated elements in seq column from this
        position backwards to the previous non-rotated element
        """
        if found is None:
            found = []
        if pos < 0:
            return found
        if not sequences[pos].endswith("$"):
            found.append(sequences[pos][-1])
            return find_rotated(pos - 1, found[::-1])
        return found

    ret = []
    for i in range(len(sequences) - 1, -1, -1):
        cur_seq = sequences[i]
        if cur_seq.endswith("$"):
            my_rot = find_rotated(i-1)
            ret.append([cur_seq[:-1], my_rot])
    return ret[::-1]

def len_column(table):
    """
    Add length column containing the length of the original entry in the seq column.
    Insert a number column with numbered entries for each row.
    """
    for pos, i in enumerate(table):
        i.insert(0, pos)
        i.append(len(i[1]))
    return table

def offset_table(table):
    """
    Record the first entry number that begins with each letter.

    !NOTE!
    Need this edge case... If a alphabet letter does not appear in any entry, it's entry
    number should be set to the entry number of it's nearest successor.
    Is this just the sorted letter before? If there were no As, C becomes 0,
    if there were no Cs, Gs would become As (0) or Ts (some position?)
    """
    ret = {}
    for l in LETTERS:
        ret[l] = -1
    for i in table:
        if ret[i[1][0]] == -1:
            ret[i[1][0]] = i[0]
    return ret

def prefix_counts(table):
    """
    Compute prefix counts.
    In order to USE the table, an additional column per alphabet letter is needed. The
    columns logically represents the appearance count of that letter in the rotate column
    from all rows up to, but not including, the current.
    """
    running_row = [0] * len(LETTERS)
    for row in table:
        row.append(list(running_row))
        for pre in row[2]:
            running_row[LETTERS.index(pre)] += 1
    table.append([len(table), "", [], None, running_row])
    return table

def push_front(table, offset, X, Y, cur_context, base):
    """
    Push base to the end of the current context to find the new X/Y/context
    """
    print(f"Context [{X}, {Y}) '{cur_context}' -> pushing {base}")

    Xa = table[X][4][LETTERS.index(base)]
    Ya = table[Y][4][LETTERS.index(base)]

    Xb = Xa + offset[base]
    Yb = Ya + offset[base]

    cur_context = base + cur_context

    if table[Xb][3] < len(cur_context):
        Xb += 1

    # Halt the search if we create an invalid context
    if Xb >= Yb:
        return None, None, cur_context

    return Xb, Yb, cur_context

def sequence_search(table, offset, seq):
    """
    Output the step by step of searching for a sequence
    """
    print(f"Searching for '{seq}'")
    # Context start/end/sequence
    X = 0
    Y = len(table) - 1
    cur_context = ""

    # Must traverse read backwards
    next_pos = len(seq) - 1

    # if None, doesn't exist. if all the bases, found
    while X is not None and next_pos >= 0:
        cur_letter = seq[next_pos]
        X, Y, cur_context = push_front(table, offset, X, Y, cur_context, cur_letter)
        next_pos -= 1

    if X is None:
        print(f"Sequence does not exist")
    else:
        print(f"Found context [{X}, {Y}) '{cur_context}'")
    return X, Y, cur_context

def pop_back(table, offset, X, Y, cur_context): # pylint: disable=unused-argument
    """
    Drop the last character of a context and create a new context
    [0, 'ACTG', [], 4, [0, 0, 0, 0]]
    [1, 'ATGC', [], 4, [0, 0, 0, 0]]
    [2, 'CCGC', ['G'], 4, [0, 0, 0, 0]]
    [3, 'CGC', ['C'], 3, [0, 0, 1, 0]]
    [4, 'CTGC', ['A'], 4, [0, 1, 1, 0]]
    [5, 'GC', ['C', 'T'], 2, [1, 1, 1, 0]]
    [6, 'TGC', ['A', 'C'], 3, [1, 2, 1, 1]]
    [7, '', [], None, [2, 3, 1, 1]]
    """

    print(f"popping [{X}, {Y}) '{cur_context}'")
    base = cur_context[-1]
    curstring = cur_context[:-1]  # pylint: disable=unused-variable
    while table[X][4][LETTERS.index(base)] >= len(cur_context):
        X -= 1

    while table[Y][4][LETTERS.index(base)] >= len(cur_context):
        Y += 1

    return X, Y, cur_context

def output(data, maxlines=MAXLINES):
    """
    Helper to output a list of stuff.
    Useful for long lists in need of trucation.
    """
    if len(data) < maxlines or maxlines == -1:
        for i in data:
            print(i)
        return
    for i in data[:int(maxlines / 2)]:
        print(i)
    print("...")
    for i in data[-int(maxlines/2):]:
        print(i)

def format_seqset(table, offset):
    """
    Pretty format the seqset
    """
    tab = "\t"
    print(f"offsets\t\t\t\t{tab.join([str(x) for x in offset.values()])}")
    print(f"pos\tseq\trotate\tlength\t{tab.join(LETTERS)}")
    for row in table:
        print(f"{row[0]}\t{row[1]}\t{','.join(row[2])}\t{row[3]}\t{tab.join(str(x) for x in row[4])}")

def workflow(reads): # pylint: disable=redefined-outer-name
    """
    Goes through each step of creating a BioGraph file out of the reads
    Returns tuple of
        seqset (*fix and offset tables)
        readmap
    """
    print("reads\n----")
    output(reads)
    print("")

    print("suffixes\n----")
    sfx = make_suffixes(reads)
    output(sfx)
    print("suffixes sorted\n----")
    sfx = lex_sort(sfx)
    output(sfx)
    print("")

    print("prefixes removed\n----")
    pre = rm_prefixes(sfx)
    output(pre)
    print("")

    print("rotation\n----")
    rot = rotation(pre)
    output(rot)
    print("rotation sorted\n----")
    rot = lex_sort(rot)
    output(rot)
    print("")

    print("table\n----")
    fix = make_table(rot)
    output(fix)
    print("")

    print("length/number column\n----")
    fix = len_column(fix)
    output(fix)
    print("")

    print("offset table\n----")
    offset = offset_table(fix)
    print(offset)
    print("")

    # At this point, the second and third columns represent all the data that is needed
    # to find subsequences in the original list of sequences. This is all the need to
    # be stored. To help explain the methods of searching, we will leave the first
    # column in our figures, but it does not exist in implementation.

    print("prefix counts\n----")
    fix = prefix_counts(fix)
    output(fix)

    print("\nfinished making BioGraph seqset\n")
    format_seqset(fix, offset)
    print("")
    # In practice, we add all reads and their reverse complement to the BioGraph. Raw
    # reads have no absolute orientation (i.e. they can be from the direct or complement
    # strand relative to reference). By adding both, we're able to perform assembly
    # in a single direction and still detect all possible overlaps.

    # NOTE! Need to illustrate how a readmap is created also.
    # Though I don't think we need to do all the fancy hashing if it becomes too crazy
    # for simple list/dict data structures

    return fix, offset

def example_queries(prefix, offset):
    """
    Illustrate core operations on seqset (currently only works with pat reads)
    """
    # Lets look for a sequence that's not contained in the original sets of reads
    X, Y, curstring = sequence_search(prefix, offset, "GTGC")
    if X is not None: # we don't expect to find this
        print("Problem with GTGC")
    print("")

    # Lets look for a sequence that is contained in the original sets of reads
    X, Y, curstring = sequence_search(prefix, offset, "TGC")
    if X is None:
        print("Found no sequence for TGC.. problem")
    print("")

    # Now we'll pop from the context returned above in order to 'widen' it
    X, Y, curstring = pop_back(prefix, offset, X, Y, curstring)
    print(f"Popped context [{X}, {Y}) 'curstring'")

    # Find reads beginning with context
    # find_reads(readmap, X, Y, curstring)


if __name__ == '__main__':
    #reads = make_reads()
    reads = make_pat_reads()
    seqset = workflow(reads)
    example_queries(*seqset)
