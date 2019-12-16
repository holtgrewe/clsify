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


def describe(*, pos, ref, alt, **kwargs):
    """Return description of variant."""
    if "-" in ref and "-" in alt:
        raise Exception("Cannot handle this yet!")
    elif "-" in ref:
        if ref[1] != "-":
            if len(ref) == 1:
                return "n.%ddel%sins%s" % (pos + 1, only_bases(ref[1:]), alt[1:])
            else:
                return "n.%d_%ddel%sins%s" % (
                    pos + 1,
                    pos + len(ref) - 1,
                    only_bases(ref[1:]),
                    alt[1:],
                )
        else:
            return "n.%d_%dins%s" % (pos, pos + 1, alt[1:])
    elif "-" in alt:
        if alt[1] != "-":
            if len(ref) == 1:
                return "n.%ddel%sins%s" % (pos + 1, ref[1:], only_bases(alt[1:]))
            else:
                return "n.%d_%ddel%sins%s" % (
                    pos + 1,
                    pos + len(ref) - 1,
                    ref[1:],
                    only_bases(alt[1:]),
                )
        else:
            if len(ref) == 2:
                return "n.%ddel%s" % (pos + 1, ref[1:])
            else:
                return "n.%d_%ddel%s" % (pos + 1, pos + len(ref) - 1, ref[1:])
    elif len(ref) == 1 and len(alt) == 1:
        return "n.%d%s>%s" % (pos, ref, alt)
    else:
        if len(ref) == 2:
            return "n.%ddel%s" % (pos + 1, ref[1:])
        else:
            return "n.%d_%ddel%s" % (pos + 1, pos + len(ref) - 1, ref[1:])


def only_bases(s):
    return "".join([x for x in s if x.upper() in "CGATN"])


def normalize_var(ref_seq, curr_var, offset):
    """Normalize the variant."""
    ref_seq = only_bases(ref_seq)
    ref = only_bases(curr_var["ref_bases"])
    alt = only_bases(curr_var["alt_bases"])
    orig_ref, orig_alt = ref, alt
    pos = curr_var["pos"] - 1 - offset
    delta = 0
    while ref[-1] == alt[-1] and pos > 0:
        delta += 1
        pos -= 1
        ref = ref_seq[pos] + ref
        alt = ref_seq[pos] + alt
        while ref[-1] == alt[-1] and len(ref) > 1 and len(alt) > 1:
            ref = ref[:-1]
            alt = alt[:-1]

    result = dict(curr_var)
    result["pos"] -= delta
    result["ref"] = ref if len(ref) >= len(alt) else ref + "-" * (len(alt) - len(ref))
    result["alt"] = alt if len(ref) <= len(alt) else alt + "-" * (len(ref) - len(alt))
    result["ref_bases"] = ref
    result["alt_bases"] = alt
    result["description"] = describe(**result)
    return result


def call_variants(ref, alt, offset=0):
    """Call variants from alignment (``ref`` and ``alt`` with gaps).

    The resulting positions are 1-based.
    """
    if len(ref) != len(alt):
        raise Exception("Invalid alignment row lengths: %d vs. %d", len(ref), len(alt))
    result = {}
    ref = ref.upper()
    alt = alt.upper()

    pos_ref = offset  # position in ref
    i = 0  # position in rows
    curr_var = {}
    while i < len(ref):
        if ref[i] == alt[i]:
            if curr_var:
                if curr_var["ref"].startswith("-") or curr_var["alt"].startswith("-"):
                    if curr_var["ali_pos"] > 1:  # left-pad indels if possible
                        curr_var = {
                            "pos": curr_var["pos"] - 1,
                            "ali_pos": curr_var["ali_pos"] - 1,
                            "ref": ref[curr_var["ali_pos"] - 2] + curr_var["ref"],
                            "alt": alt[curr_var["ali_pos"] - 2] + curr_var["alt"],
                        }
                    else:  # use base right of it, VCF-style
                        curr_var = {
                            "pos": curr_var["pos"],
                            "ali_pos": curr_var["ali_pos"],
                            "ref": curr_var["ref"]
                            + ref[curr_var["ali_pos"] + len(curr_var["ref"])],
                            "alt": curr_var["alt"]
                            + alt[curr_var["ali_pos"] + len(curr_var["alt"])],
                        }
                curr_var["ref_bases"] = only_bases(curr_var["ref"])
                curr_var["alt_bases"] = only_bases(curr_var["alt"])
                curr_var["description"] = describe(**curr_var)
                curr_var = normalize_var(ref, curr_var, offset)
                result[curr_var["pos"]] = curr_var
                curr_var = {}
        else:
            if curr_var:
                curr_var["ref"] = curr_var["ref"] + ref[i]
                curr_var["alt"] = curr_var["alt"] + alt[i]
            else:
                curr_var = {"pos": pos_ref + 1, "ali_pos": i + 1, "ref": ref[i], "alt": alt[i]}

        if ref[i] != "-":
            pos_ref += 1
        i += 1

    if curr_var:
        if curr_var["ref"].startswith("-") or curr_var["alt"].startswith("-"):
            if curr_var["ali_pos"] > 1:  # left-pad indels if possible
                curr_var = {
                    "pos": curr_var["pos"] - 1,
                    "ali_pos": curr_var["ali_pos"] - 1,
                    "ref": ref[curr_var["ali_pos"] - 2] + curr_var["ref"],
                    "alt": alt[curr_var["ali_pos"] - 2] + curr_var["alt"],
                }
            else:  # use base right of it, VCF-style
                curr_var = {
                    "pos": curr_var["ali_pos"],
                    "ali_pos": curr_var["ali_pos"],
                    "ref": curr_var["ref"] + ref[curr_var["ali_pos"] + len(curr_var["ref"])],
                    "alt": curr_var["alt"] + alt[curr_var["ali_pos"] + len(curr_var["alt"])],
                }
        curr_var["ref_bases"] = only_bases(curr_var["ref"])
        curr_var["alt_bases"] = only_bases(curr_var["alt"])
        curr_var = normalize_var(ref, curr_var, offset)
        curr_var["description"] = describe(**curr_var)
        result[curr_var["pos"]] = curr_var

    return result
