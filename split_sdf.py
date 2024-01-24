fname = 'CNS_compound'
split_number= 2 #number of molecules in each file
number_of_sdfs = split_number
i = 0; j = 0; k = 0
with open(f'{fname}.sdf') as f:
    for line in f:
        if "$$$$" in line:
            k += 1
print(f'There are {k} molecules in {fname}.sdf')

f2 = open(f'{fname}_{j}.sdf', 'w')
with open(f'{fname}.sdf') as f:
    for line in f:
        f2.write(line)
        if "$$$$" in line:
            i += 1
        if i == k:
            f2.close()
            break
        if i == number_of_sdfs:
            number_of_sdfs += split_number 
            f2.close()
            j += 1
            f2 = open(f'{fname}_{j}.sdf', 'w')
