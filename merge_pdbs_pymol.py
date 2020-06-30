import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI


from pymol import cmd, finish_launching
import glob

def merge_pdb(pdb1, pdb2):
    #finish_launching()
    cmd.load('pdbs/' + pdb1)
    cmd.load('pdbs/' + pdb2)
    path = './docking_decoy/'
    merged_name = pdb1.split('.')[0] + '_' + pdb2.split('.')[0] + '.pdb'
    cmd.set('valence_mode', 1)
    cmd.save(path + merged_name, state = -1)
    #cmd.quit()
    cmd.delete('all')

paths = glob.glob('pdbs/*')
#paths.sort()
names = []
for path in paths:
    name = path.split('/')[-1]
    names.append(name)

print(names)
for name in names[1:]:
    print(name)
    merge_pdb(names[0], name)