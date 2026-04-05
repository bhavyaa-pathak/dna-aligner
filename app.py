from flask import Flask, render_template, request, jsonify
import re

app = Flask(__name__)

# ─── MOTIF PATTERNS ────────────────────────────────────────────────────────────
DNA_MOTIFS = {
    'TATA Box':          r'TATAAA',
    'Kozak Sequence':    r'[AG]CCATGG',
    'CpG Island':        r'CG',
    'EcoRI Site':        r'GAATTC',
    'BamHI Site':        r'GGATCC',
    'Start Codon':       r'ATG',
    'Stop Codons':       r'(?:TAA|TAG|TGA)',
    'Poly-A Signal':     r'AATAAA',
    'Splice Donor':      r'GT[AG]AGT',
}

PROTEIN_MOTIFS = {
    'N-Glycosylation':      r'N[^P][ST]',
    'RGD Motif':            r'RGD',
    'KDEL Retention':       r'KDEL',
    'Nuclear Localization': r'[KR]{3,}',
    'Casein Kinase II':     r'[ST].{2}[DE]',
    'Cysteine Pairs':       r'C.{1,8}C',
}


def parse_fasta(raw):
    """
    Parse raw input — skip FASTA header lines (starting with >)
    and return only the sequence.
    """
    lines = raw.strip().splitlines()
    seq_lines = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            continue  # skip accession/header line
        seq_lines.append(re.sub(r'[\s\d]', '', line))
    return ''.join(seq_lines).upper()


def detect_type(seq):
    dna_chars = set('ATGCN')
    if len(seq) == 0:
        return 'UNKNOWN'
    non_n = [c for c in seq if c != 'N']
    if not non_n:
        return 'DNA'
    if all(c in dna_chars for c in non_n):
        return 'DNA'
    return 'PROTEIN'


# ─── NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT ────────────────────────────────────────
def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    # limit for performance
    if n > 20000 or m > 20000:
        return None, None, None, "Sequences too long (max 20000 chars each for alignment)"
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i * gap
    for j in range(m + 1):
        dp[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

    # Traceback
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + score:
                a1.append(seq1[i-1])
                a2.append(seq2[j-1])
                i -= 1; j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            a1.append(seq1[i-1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(seq2[j-1])
            j -= 1

    a1 = ''.join(reversed(a1))
    a2 = ''.join(reversed(a2))

    # Identity
    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
    aligned_len = len(a1)
    identity = round(matches / aligned_len * 100, 2) if aligned_len else 0

    return a1, a2, identity, None


def find_motifs(seq, seq_type):
    motif_dict = DNA_MOTIFS if seq_type == 'DNA' else PROTEIN_MOTIFS
    found = {}
    for name, pattern in motif_dict.items():
        matches = [(m.start(), m.end(), m.group()) for m in re.finditer(pattern, seq)]
        if matches:
            found[name] = matches
    return found


def build_match_line(a1, a2):
    line = []
    for x, y in zip(a1, a2):
        if x == '-' or y == '-':
            line.append(' ')
        elif x == y:
            line.append('|')
        else:
            line.append('.')
    return ''.join(line)


# ─── ROUTES ────────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/align', methods=['POST'])
def align():
    data = request.json
    raw1 = data.get('seq1', '').strip()
    raw2 = data.get('seq2', '').strip()
    mode = data.get('mode', 'dna')       # 'dna' or 'protein'
    motif_on = data.get('motif', False)

    if not raw1 or not raw2:
        return jsonify({'error': 'Please enter both sequences.'}), 400

    seq1 = parse_fasta(raw1)
    seq2 = parse_fasta(raw2)

    if not seq1 or not seq2:
        return jsonify({'error': 'Could not parse sequences. Make sure they contain valid characters.'}), 400

    type1 = detect_type(seq1)
    type2 = detect_type(seq2)

    # Warn on type mismatch
    if mode == 'dna' and (type1 == 'PROTEIN' or type2 == 'PROTEIN'):
        return jsonify({'error': 'DNA mode selected but one or both sequences look like protein. Please switch to Protein mode.'}), 400
    if mode == 'protein' and (type1 == 'DNA' and type2 == 'DNA'):
        return jsonify({'error': 'Protein mode selected but sequences look like DNA. Please switch to DNA mode.'}), 400

    a1, a2, identity, err = needleman_wunsch(seq1, seq2)
    if err:
        return jsonify({'error': err}), 400

    match_line = build_match_line(a1, a2)
    matches = sum(1 for c in match_line if c == '|')
    gaps = sum(1 for c in match_line if c == ' ')
    mismatches = sum(1 for c in match_line if c == '.')

    response = {
        'seq1_clean': seq1,
        'seq2_clean': seq2,
        'aligned1': a1,
        'aligned2': a2,
        'match_line': match_line,
        'identity': identity,
        'length': len(a1),
        'matches': matches,
        'mismatches': mismatches,
        'gaps': gaps,
        'seq1_type': type1,
        'seq2_type': type2,
        'mode': mode,
    }

    if motif_on:
        seq_type = 'DNA' if mode == 'dna' else 'PROTEIN'
        m1 = find_motifs(seq1, seq_type)
        m2 = find_motifs(seq2, seq_type)
        response['motifs1'] = {k: [list(x) for x in v] for k, v in m1.items()}
        response['motifs2'] = {k: [list(x) for x in v] for k, v in m2.items()}

    return jsonify(response)


if __name__ == '__main__':
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
