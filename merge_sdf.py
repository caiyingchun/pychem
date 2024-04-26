import os
import glob
f = open('merged.sdf', 'w')
num_confs = 0
for filename in glob.glob('sdf/*.sdf'):
    print(filename)
    f0 = open(filename, 'r')
    while True:
        line = f0.readline()
        if '$$$$' in line:
            num_confs += 1
        if line:
            f.write(line)
        else:
            break
f.close()
#os.system('rm -rf sdf/')
print(f'There are {num_confs} conformations.')
