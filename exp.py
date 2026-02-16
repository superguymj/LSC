import os
import math

lines = []
with open('temp.txt', 'r') as f:
    lines = f.readlines()

n = int(math.sqrt(len(lines)))
for i in range(n):
    a = lines[i * n : i * n + n]
    a = list(map(int, a))
    for j in range(n):
        a[j] -= 1
    print(*a)
    
    
