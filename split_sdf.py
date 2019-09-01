f = 'your_file_name'
split_number= 1000 #number of molecules in each file
number_of_sdfs = split_number
i = 0; j = 0
f2 = open(f + '_' + str(j) + '.sdf', 'w')
for line in open(f + '.sdf'):
    f2.write(line)
    if line[:4] == "$$$$":
        i += 1
    if i > number_of_sdfs:
        number_of_sdfs += split_number 
        f2.close()
        j += 1
        f2 = open(f + '_' + str(j) + '.sdf','w')
print(i)
