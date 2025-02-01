import os


lines = []
with open('temp.txt', 'r') as f:
    lines = f.readlines()

a = 0.0
# for i in range(0, 60, 2):
#     line = lines[i]
#     b = float(line.split()[-1][:-1])
#     a += b
    
for line in lines:
    a += float(line)
    
print(a / 30)
