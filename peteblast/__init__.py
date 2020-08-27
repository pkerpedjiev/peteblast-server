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

__version__ = "0.1.1"

DEFAULT_MINIMIZER_K = 5
DEFAULT_MINIMIZER_W = 10

from skbio.alignment import StripedSmithWaterman
from elasticsearch import Elasticsearch
es = Elasticsearch(retry_on_timeout=True)

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

def results(query, max_results=30):    
#     print("query:", query)
#     t1 = time.time()
#     id_counts = search(query)
#     t2 = time.time()
    
#     print("time:", t2 - t1)
    mins = minimizers(query, k=5, w=10)
    es_query = {
        "query": {
            "bool": {
                "should": 
                    [{'match':{ "mins": m}} for m in mins]
                
            }
        }
    }
    
    res = es.search(index='seqs', body=es_query)
    ssw = StripedSmithWaterman(query, protein=True, substitution_matrix=subst_mat)
    results_out = []
    for hit in res['hits']['hits'][:max_results]:
        seq = hit['_source']
#         print("seq", seq)
        desc = seq['description']
        sequence = str(seq['seq'])
        alignment = ssw(sequence)
        match = 0
        total = 0
        for q, t in zip(alignment.aligned_query_sequence, alignment.aligned_target_sequence):
            if q == "-" or t == "-":
                continue
            if q == t:
                match += 1
            total += 1
        identity = match / total
        results_out += [
            {
                "description": seq['description'],
                "aligned_query": alignment.aligned_query_sequence,
                "aligned_target": alignment.aligned_target_sequence,
                "target_seq": str(seq['seq']),
                "score": int(alignment.optimal_alignment_score),
                "identity": "{:.2f}".format(identity),
            }
        ]
    return results_out

# get_sequence = get_sequence_function(
#     op.join(os.environ['BLAST_INDEX_LOCATION'],
#         'seq.fa'))

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


    @app.route("/api/v1/search/", methods=['POST'])
    def hello():
        # query = 'MASTQNIVEEVQKMLDTYDTNKDGEITKAEAV'
        query = request.get_json()['searchString'];

        # print("request_json", request.get_json())
        results_out = results(query)

        return jsonify({ "results": results_out})
    return app
    
app = create_app()
