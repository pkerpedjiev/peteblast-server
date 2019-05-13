import logging
import os
import os.path as op
import sys
from flask import Flask


import collections as col
from io import StringIO
from Bio import SeqIO
import mmap
import sqlite3

DEFAULT_MINIMIZER_K = 5
DEFAULT_MINIMIZER_W = 10

def minimizers(sequence, k, w):
    queue = []
    minimizers = []
    
    for i in range(len(sequence) - k):
        subseq = sequence[i:i+k]

        queue += [subseq]
        if len(queue) > w + k - 1:
            queue = queue[-(w+k-1):]

        
        sorted_queue = sorted(queue)
        #print("queue", queue)
        #print("sorte", sorted_queue)
        minimizer = sorted_queue[0]
        minimizers += [minimizer]
    return set(minimizers)

def search(query, pos_dict, kmer_ix_list, results_frac=0.5):
    mins = sorted(minimizers(query, 5, 10), key=lambda x: pos_dict[x][0])
    seqs_found = col.defaultdict(int)
    
    for i, minimizer in enumerate(mins[:int(len(mins)*results_frac)]):
        ix_start = pos_dict[minimizer][1] - pos_dict[minimizer][0]
        ix_end = pos_dict[minimizer][1]
        #print("minimizer", minimizer, ix_start, ix_end)

        seqs = kmer_ix_list[ix_start:ix_end]
        #print("seqs:", seqs)
        for seq in seqs:
            seqs_found[seq] += 1
            
    #print(seqs_found)
    return sorted(seqs_found.items(), key=lambda x: -x[1])


def create_app(test_config=None):
    """Create and configure an instance of the Flask application."""
    app = Flask(__name__, instance_relative_config=True)
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



    @app.route("/hello")
    def hello():
        conn = sqlite3.connect(
            op.join(os.environ['BLAST_INDEX_LOCATION'],
                'fa.sqlite'))
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM kmer_index_offsets WHERE kmer="AAAAM"')
        rows = cursor.fetchone()
        print("rows:", rows)
        return "Hello, World! " + os.environ['BLAST_INDEX_LOCATION']

    return app
