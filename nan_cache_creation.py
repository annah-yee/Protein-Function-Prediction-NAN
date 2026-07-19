'''This file contains tools to create a cache holding information about
    each GO term'''

# neeed to have go-basic.obo downloaded first 

from goatools.obo_parser import GODag
import pickle


def create_DAG():
    godag = GODag("go-basic.obo")
    cache = {}

    for go_id, term in godag.items():
        cache[go_id] = {
            "name": term.name,
            "namespace": term.namespace,
            "parents": term.get_all_parents(),
            "children": term.get_all_children(),
            "depth": term.depth,
            "level": term.level
        }
    with open("go_cache.pkl", "wb") as f:
        pickle.dump(cache, f)

if __name__ == "__main__":
    create_DAG()