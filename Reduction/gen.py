import random

n = random.randint(10, 12)
print(n)

lsc = [[(i + j) % n for j in range(n)] for i in range(n)]
random.shuffle(lsc)

_lsc = [[lsc[j][i] for j in range(n)] for i in range(n)]
lsc = _lsc
random.shuffle(lsc)


k = random.randint(2, 3)
for i in range(n):
    for j in range(n):
        if random.randint(0, k) == 0:
            lsc[i][j] = -1
        print(lsc[i][j], end=' ')
    print('')
    