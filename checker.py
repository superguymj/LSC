for i in range(10):
    lines = []
    with open("./data/{}.txt".format(i), 'r') as f:
        lines = f.readlines()
    
    n = int(lines[0])
    fixed = []
    for line in lines[1:]:
        x, y, w = map(int, line.strip().split())
        fixed.append((x, y, w))
    
    with open("./results/{}.txt".format(i), 'r') as f:
        lines = f.readlines()
        
    lsc = []
    for line in lines:
        lsc.append(list(map(int, line.strip().split())))
    flag = True
    for (x, y, w) in fixed:
        if lsc[x][y] != w:
            flag = False
    for p in lsc:
        for i in p:
            if i < 0 or i >= n:
                flag = False
    if flag == False:
        print("{}.txt".format(i))