import os
path = r'data'
f = file('merge.sdf', 'wb')
i = 0
for filename in os.listdir(path):
    name = path + '\\' + filename
    f0 = open(name, 'rb')
    while True:
        line = f0.readline()
        if line[:4] == '$$$$':
            i += 1
        if line != '':
            f.write(line)
        else:
            break
print i
