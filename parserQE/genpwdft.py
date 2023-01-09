from munch import Munch

'''
    Generator of PWDFT inputs using info of
    quantum espresso outputs
'''
def load_default(title='My Title', name='pwdft_name', charge=0, mult=1):
    data = Munch

    data.title = title
    data.start = name
    data.charge = charge
    data.mult = mult

    return data


def genInputPWDFT(qe, data=None, fname='input.nw'):
    if data is None:
        data = load_default()
    try:
        with open(fname, "w+") as fout:
            print('# Input generator', file=fout)
            print('title \"{0}\"'.format(data.title), file=fout)
            print('\necho\n', file=fout)
            print('start {0}\n'.format(data.start), file=fout)
            print('memory  1900 mb\n', file=fout)
            print('permanent_dir ./perm', file=fout)
            print('scratch_dir   ./perm', file=fout)
            print('\ncharge {0}\n'.format(data.charge), file=fout)
            print('geometry noautosym nocenter noautoz', file=fout)
            print('system crystal', file=fout)
            print('   lattice_vectors', file=fout)
            print('     {0: 10.6} {1: 10.6} {2: 10.6}'.format(qe.chem.cell[0][0], 
                        qe.chem.cell[0][1], qe.chem.cell[0][2]), file=fout)
            print('     {0: 10.6} {1: 10.6} {2: 10.6}'.format(qe.chem.cell[1][0],
                        qe.chem.cell[1][1], qe.chem.cell[1][2]), file=fout)
            print('     {0: 10.6} {1: 10.6} {2: 10.6}'.format(qe.chem.cell[2][0],
                        qe.chem.cell[2][1], qe.chem.cell[2][2]), file=fout)
            print('end', file=fout)
            for i in range(qe.chem.natoms):
                print(' {0:3} {1: 10.6f} {2: 10.6f} {3: 10.6f}'.format(qe.chem.symb[i],
                qe.chem.coors[i][0],
                qe.chem.coors[i][1],
                qe.chem.coors[i][2]), file=fout)
            print('end\n', file=fout)
            #
            # Block of NWPW
            #
            print('nwpw', file=fout)
            print('   monkhorst-pack 1 1 1', file=fout)
            print('   cutoff {0}'.format(qe.chem.ecut/2.), file=fout)
            print('   mult {0}'.format(data.mult), file=fout)
            print('   xc {0}'.format(qe.chem.funct), file=fout)
            print('   lmbfgs grassman\n', file=fout)
            print('   2d-hcurve', file=fout)
            print(' end\n',file=fout)
            #
            # Extra parameters block
            #
            print('set nwpw:kbpp_ray     .true.', file=fout)
            print('set nwpw:kbpp_filter  .true.', file=fout)
            print('\ntask pspw energy ignore\n', file=fout)

    except IOError as error:
        print('Error to open file: {0}'.format(fname))
