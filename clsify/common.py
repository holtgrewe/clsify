def load_tsv(input_path):
    header = None
    records = []
    with open(input_path, "rt") as inputf:
        for line in inputf:
            if line.startswith("#"):
                continue
            arr = line.strip().split("\t")
            if not header:
                header = arr
            else:
                records.append(dict(zip(header, arr)))

    return header, records


def rev(seq):
    return list(reversed(seq))


def revcomp(seq):
    m = {"a": "t", "A": "T", "t": "a", "T": "A", "c": "g", "C": "G", "g": "c", "G": "C"}
    return "".join(reversed(list(map(lambda x: m.get(x, x), seq))))


def load_fasta(path):
    result = {}
    with open(path, "rt") as inputf:
        name = None
        seq_lines = []
        for line in inputf:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    result[name] = "".join(seq_lines)
                name = line[1:].split()[0]
                seq_lines = []
            elif name:  # ignore if first line is not ">"
                seq_lines.append(line)
        if name:
            result[name] = "".join(seq_lines)
    return result
