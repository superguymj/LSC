import os
import math

def clean(fileName):
    lines = []
    with open(fileName, 'r') as f:
        lines = f.readlines()
        
    n = 0
    fixed = []
    for line in lines:
        if line[0] == 'p':
            n = int(line.strip().split(' ')[-2])
            n = int(math.sqrt(n))
        if line[0] == 'f':
            s = line.strip().split(' ')
            if len(s) == 3:
                x, y = map(int, s[-2:])
                x -= 1
                y -= 1
                i = x // n
                j = x % n
                fixed.append((i, j ,y))
    
    return (n, fixed)


def deep_traverse(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            # 拼接完整路径
            file_path = os.path.join(root, file)
            n, fixed = clean(file_path)
            
            with open('./data/' + file, 'w') as f:
                f.write(str(n) + '\n')
                for (i, j, c) in fixed:
                    f.write("{} {} {}\n".format(i, j, c))
            print(file + " success...")

# 使用示例
deep_traverse("./dataset/")  # 主文件夹路径