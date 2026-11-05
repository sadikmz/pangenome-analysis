"""
Burrows-Wheeler Transform demo.

encode() is fully implemented (naive rotation-matrix method, O(n^2 log n)
because of sorting n strings of length n -- fine for teaching, not for
production on megabase genome sequences).

decode() is left for you to implement -- see the TODO below.
"""

EOF = "$"  # sentinel: must sort before every other character in the alphabet


def encode(s: str) -> tuple[str, int]:
    """
    Compute the BWT of s.

    Returns:
        (bwt_string, original_index) -- original_index is the row in the
        sorted rotation matrix that corresponds to the un-rotated input.
        You need it (or the EOF marker) to decode unambiguously.
    """
    s = s + EOF
    n = len(s)

    # All rotations, e.g. "banana$" -> "banana$", "anana$b", "nana$ba", ...
    rotations = [s[i:] + s[:i] for i in range(n)]
    sorted_rotations = sorted(rotations)

    last_column = "".join(row[-1] for row in sorted_rotations)
    original_index = sorted_rotations.index(s)

    return last_column, original_index


def decode(bwt: str) -> str:
    """
    Invert the BWT back to the original string (including the trailing EOF).

    Two standard approaches, both valid -- pick one:

    1. Iterative table rebuild: repeatedly prepend `bwt` as a new column to
       a growing table and re-sort, n times. Simple to reason about,
       O(n^2 log n), mirrors the encode() logic above almost exactly.

    2. LF-mapping: build `first_column = sorted(bwt)`, then for each
       character in `bwt` compute its rank among identical characters seen
       so far (this is the "L" occurrence count). The same character at the
       same rank in `first_column` tells you which row you jump to next.
       Walk backward from the row containing EOF, prepending characters,
       for O(n) (or O(n log n) with a naive rank lookup).

    TODO(user): implement one of these. Approach 2 is the one real BWT-based
    tools (bzip2, FM-index) actually use, if you want the more realistic version.
    """
    raise NotImplementedError


if __name__ == "__main__":
    text = "banana"
    bwt, idx = encode(text)
    print(f"input:        {text}")
    print(f"bwt:          {bwt}")
    print(f"original_idx: {idx}")

    # Uncomment once decode() is implemented:
    # recovered = decode(bwt)
    # print(f"decoded:      {recovered}")
    # assert recovered == text + EOF
