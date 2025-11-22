import numpy as np
# TODO change memoization to 2D or !D flatten array of pascal triangle values
cache = {}
def Newton(n,k):
    if (n, k) in cache:
        return cache.get((n,k))
    else:
        val = 1 if n==k or k==0 else Newton(n-1, k-1) + Newton(n-1,k)
        cache[(n,k)] = val
        return val

def I(n):
    # Fixed range range(0, n) does go only to n-1 => range(0, n+1) <=> [0, n] n in N
    return np.sum(np.array([ 2.0**k*Newton(n, k) for k in range(0,n+1) ]))

# %%time
I_list = []
# Fixed iteration index, We need to start from 1 
for n in range(1, 31):
    I_list.append(I(n))

# Złożonośc obliczeniowa dla nie poprawionego alogrytnu to O(2^n), ponieważ w najgorszym przypadku, każde wywołanie funkcji rozgałezia sie na 2 kolejne wywołania.

# Alorytm po poprawieniu poprzez zastosowanie memoizacji oblicza każdą pare tylko raz, więc ilośc obliczęn funkcji sprowadza się do O(n*k) w najgorszym przypadku