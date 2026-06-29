import pickle

with open("go_cache.pkl", "rb") as f:
    cache = pickle.load(f)

# see the keys of the cache itself
print(type(cache))
print(list(cache.keys())[:5])

# see what one entry looks like
first_key = list(cache.keys())[0]
print(first_key)
print(cache[first_key])