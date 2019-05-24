import h5py
import logging
import os
import os.path as op
import sys

from flask import Flask
from flask import jsonify
from flask import request

import collections as col
from io import StringIO
from Bio import SeqIO
import mmap
import sqlite3

from flask_cors import CORS

__version__ = "0.1.0"

DEFAULT_MINIMIZER_K = 5
DEFAULT_MINIMIZER_W = 10

from skbio.alignment import StripedSmithWaterman

def calc_subst_mat():
    from Bio.SubsMat.MatrixInfo import blosum62
    import collections as col

    subst_mat = col.defaultdict(dict)

    for key, value in blosum62.items():
        subst_mat[key[0]][key[1]] = value
        subst_mat[key[1]][key[0]] = value
        subst_mat['*'][key[0]] = 1
        subst_mat[key[0]]['*'] = 1
        subst_mat[key[1]]['*'] = 1
        subst_mat['*'][key[1]] = 1

    subst_mat['*']['*'] = 1
    return subst_mat

subst_mat = calc_subst_mat()

def minimizers(sequence, k, w):
    queue = []
    minimizers = []
    
    for i in range(len(sequence) - k):
        subseq = sequence[i:i+k]

        queue += [subseq]
        if len(queue) > w + k - 1:
            queue = queue[-(w+k-1):]

        
        sorted_queue = sorted(queue)
        minimizer = sorted_queue[0]
        minimizers += [minimizer]
    return set(minimizers)

def kmer_ix_start_length(kmer):
    conn = sqlite3.connect(
        op.join(os.environ['BLAST_INDEX_LOCATION'],
            'fa.sqlite'))
    cursor = conn.cursor()
    cursor.execute(f'SELECT offset, length FROM kmer_index_offsets WHERE kmer="{kmer}"')
    row = cursor.fetchone()

    return row

def search(query, kmer_ix_list, results_frac=0.5):
    mins = [m for m in minimizers(query, 5, 10) if kmer_ix_start_length(m)]
    mins = sorted(mins,
        key=lambda x: kmer_ix_start_length(x)[1])
    seqs_found = col.defaultdict(int)
    MIN_TO_CHECK = 10
    num_to_check = max(
        int(len(mins)*results_frac),
        MIN_TO_CHECK
    )

    for _, minimizer in enumerate(mins[:num_to_check]):
        # print("minimizer:", minimizer)
        (ix_start, ixs_length) = kmer_ix_start_length(minimizer)
        ix_end = ix_start + ixs_length
        # print("minimizer", minimizer, ix_start, ix_end)

        seqs = kmer_ix_list[ix_start:ix_end]
        # print("seqs:", seqs)
        for seq in seqs:
            seqs_found[seq] += 1
            
    # print(seqs_found)
    return sorted(seqs_found.items(), key=lambda x: -x[1])

def get_offset(index):
    conn = sqlite3.connect(
        op.join(os.environ['BLAST_INDEX_LOCATION'],
            'fa.sqlite'))
    c = conn.cursor()
    c.execute('SELECT offset from offsets where ix == {};'.format(index))
    row = c.fetchone()
    return row[0]

def get_sequence_function(filename):
    '''
    Obtain the nth sequence from a fasta file
    
    Paramseters
    -----------
    filename: str
        The filename of the fasta file
    '''
    with open(filename, 'r+') as f:
        mm = mmap.mmap(f.fileno(), 0)
        # subtract one because offsets are 0-based and the input parameter
        # is one-based
        def get_sequence(number):
            start = get_offset(number)
            end = get_offset(number+1)

            # print('start:', start, 'end:', end, 'mm:', len(mm))
            # print('mm', mm[start:end])
            record = list(SeqIO.parse(StringIO(mm[start:end].decode('ascii')), "fasta"))[0]
            return record
        return get_sequence

get_sequence = get_sequence_function(
    op.join(os.environ['BLAST_INDEX_LOCATION'],
        'seq.fa'))

def create_app(test_config=None):
    """Create and configure an instance of the Flask application."""
    app = Flask(__name__, instance_relative_config=True)
    CORS(app)

    app.config.from_mapping(
        # a default secret that should be overridden by instance config
        SECRET_KEY="dev",
        # store the database in the instance folder
        DATABASE=os.path.join(app.instance_path, "flaskr.sqlite"),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.update(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    if 'BLAST_INDEX_LOCATION' not in os.environ:
        logging.critical("BLAST_INDEX_LOCATION needs to be defined as an environment variable")
        return None

    BIL = os.environ['BLAST_INDEX_LOCATION']

    kmer_ix_list_filename = op.join(BIL, 'sorted.merged.h5')
    kmer_ix_list = h5py.File(kmer_ix_list_filename, 'r')['ixs']


    @app.route("/api/v1/search", methods=['POST'])
    def hello():
        # query = 'MASTQNIVEEVQKMLDTYDTNKDGEITKAEAV'
        query = request.get_json()['searchString'];

        # print("request_json", request.get_json())
        results = search(query, kmer_ix_list)
        ssw = StripedSmithWaterman(query,
            protein=True, substitution_matrix=subst_mat)

        results_out = []

        for res in results[:50]:
            seq = get_sequence(res[0])

            desc = seq.description
            sequence = str(seq.seq)
            alignment = ssw(sequence)

            results_out += [{
                'description': seq.description,
                'aligned_query': alignment.aligned_query_sequence,
                'aligned_target': alignment.aligned_target_sequence,
                'target_seq': str(seq.seq),
                'target_description': str(seq.description),
                'score': int( alignment.optimal_alignment_score),
                'minimizer_matches': res[1]
            }]

        return jsonify({ "results": results_out})
    return app
    
app = create_app()
