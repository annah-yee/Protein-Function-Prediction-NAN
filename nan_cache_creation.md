# Create Cache of GO Terms

Gotta do this so that we can hold information about each go term

## Import requirements

```python
from goatools.obo_parser import GODag
import pickle
```

## Create the cache

We're gonna save this stuff as a pickle object btw. 
Note running this will automatically reset the .pkl object (and effectively delete any additions from running `nan_freq` things). 

```python
def create_DAG():

    cache = {}

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
```
