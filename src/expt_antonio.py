import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    n = 10
    p1 = 0.75
    p2 = 1 - p1

    for a in [9]:
        max_f = -float('Infinity')
        max_i = -1
        max_p1 = -1
        for j in range(1000+1):
            p1 = 1.0*j/1000
            p2 = 1 - p1
            for i in range(a+1):
                try:
                    v1 = math.log(math.comb(n, i))
                    v2 = math.log(math.comb(n, a-i))
                    v3 = (n-a+2*i) * math.log(p1)
                    v4 = (n+a-2*i) * math.log(p2)
                except:
                    continue
                v = v1+v2+v3+v4

                if v > max_f:
                    max_f = v
                    max_i = i
                    max_p1 = p1
        print(a, max_i, max_p1, max_f)
