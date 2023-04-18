import numpy as np

def read_theta_grid_file_into_dictionary(filename):
    def read_four_floats(file, length):
        var1 = np.zeros(length)
        var2 = np.zeros(length)
        var3 = np.zeros(length)
        var4 = np.zeros(length)
        for i in range(length):
            var1[i], var2[i], var3[i], var4[i] = [
                float(x) for x in file.readline().split()
            ]

        return var1, var2, var3, var4

    def read_three_floats(file, length):
        var1 = np.zeros(length)
        var2 = np.zeros(length)
        var3 = np.zeros(length)
        for i in range(length):
            var1[i], var2[i], var3[i] = [
                float(x) for x in file.readline().split()
            ]

        return var1, var2, var3

    try:
        with open(filename, 'r') as file:
            # Trapped lambdas
            junk_line = file.readline() #Expect this to be text nlambda
            nlambda = int(file.readline())
            junk_line = file.readline() #Expect this to be text lambda
            lambda_grid = np.zeros(nlambda)
            for i in range(nlambda):
                lambda_grid[i] = float(file.readline())
            bset = 1.0/lambda_grid

            # Theta grid size and scalar geometrical parameters
            junk_line = file.readline() #Expect this to be text ntgrid, nperiod, ntheta, drhodpsi, rmaj, shat, kxfac, q
            ntgrid, nperiod, ntheta, drhodpsi, rmaj, shat, kxfac, q = [
                type_func(x) for type_func, x in zip([int, int, int, float, float, float, float, float], file.readline().split())
            ]

            # Gbdrift, gradpar, grho, theta_grid
            junk_line = file.readline() #Expect this to be text gbdrift gradpar grho tgrid
            gbdrift, gradpar, grho, theta_grid = read_four_floats(file, 2*ntgrid+1)

            # Cvdrift, gds2, bmag, theta_grid
            junk_line = file.readline() #Expect this to be text cvdrift gds2 bmag tgrid
            cvdrift, gds2, bmag, theta_grid_tmp = read_four_floats(file, 2*ntgrid+1)

            # Sanity check
            if not np.allclose(theta_grid, theta_grid_tmp):
                print('Consistency error: Theta grid values during cvdrift read.')

            # gds21, gds22, theta_grid
            junk_line = file.readline() #Expect this to be text gds21 gds22 tgrid
            gds21, gds22, theta_grid_tmp = read_three_floats(file, 2*ntgrid+1)

            # Sanity check
            if not np.allclose(theta_grid, theta_grid_tmp):
                print('Consistency error: Theta grid values during gds21 read.')

            # cvdrift0, gbdrift0, theta_grid
            junk_line = file.readline() #Expect this to be text cvdrift0 gbdrift0 tgrid
            cvdrift0, gbdrift0, theta_grid_tmp = read_three_floats(file, 2*ntgrid+1)

            # Sanity check
            if not np.allclose(theta_grid, theta_grid_tmp):
                print('Consistency error: Theta grid values during cvdrift0 read.')

            # Rplot, Rprime, theta_grid
            junk_line = file.readline() #Expect this to be text Rplot Rprim tgrid

            if junk_line is not None and junk_line.strip() != '':
                rplot, rprime, theta_grid_tmp = read_three_floats(file, 2*ntgrid+1)

                # Sanity check
                if not np.allclose(theta_grid, theta_grid_tmp):
                    print('Consistency error: Theta grid values during Rplot read.')
            else:
                print('Missing data: Rplot/rprime not found setting to default values.')
                rplot = np.ones(2*ntgrid+1)
                rprime = np.zeros(2*ntgrid+1)

            # Zplot, Zprime, theta_grid
            junk_line = file.readline() #Expect this to be text Zplot Zprim tgrid

            if junk_line is not None  and junk_line.strip() != '':
                zplot, zprime, theta_grid_tmp = read_three_floats(file, 2*ntgrid+1)

                # Sanity check
                if not np.allclose(theta_grid, theta_grid_tmp):
                    print('Consistency error: Theta grid values during Zplot read.')
            else:
                print('Missing data: Zplot/zprime not found setting to default values.')
                zplot = np.ones(2*ntgrid+1)
                zprime = np.zeros(2*ntgrid+1)

            # Aplot, Aprime, theta_grid
            junk_line = file.readline() #Expect this to be text Aplot Aprim tgrid

            if junk_line is not None and junk_line.strip() != '':
                aplot, aprime, theta_grid_tmp = read_three_floats(file, 2*ntgrid+1)

                # Sanity check
                if not np.allclose(theta_grid, theta_grid_tmp):
                    print('Consistency error: Theta grid values during Aplot read.')
            else:
                print('Missing data: Aplot/aprime not found setting to default values.')
                aplot = np.ones(2*ntgrid+1)
                aprime = np.zeros(2*ntgrid+1)

            # Form result
            result = {
                'nbset' : nlambda,
                'nlambda' : nlambda,
                'lambda_grid' : lambda_grid,
                'bset' : bset,
                'ntgrid' : ntgrid,
                'nperiod' : nperiod,
                'ntheta' : ntheta,
                'drhodpsi' : drhodpsi,
                'rmaj' : rmaj,
                'shat' : shat,
                'kxfac' : kxfac,
                'q' : q,
                'gbdrift' : gbdrift,
                'gradpar' : gradpar,
                'grho' : grho,
                'theta_grid' : theta_grid,
                'cvdrift' : cvdrift,
                'gds2' : gds2,
                'bmag' : bmag,
                'gds21' : gds21,
                'gds22' : gds22,
                'cvdrift0' : cvdrift0,
                'gbdrift0' : gbdrift0,
                'rplot' : rplot,
                'rprime' : rprime,
                'zplot' : zplot,
                'zprime' : zprime,
                'aplot' : aplot,
                'aprime' : aprime
                }
            return result

    except FileNotFoundError:
        print(f"File '{filename}' not found.")

def write_theta_grid_data(filename, data):
    def write_four_floats(file, var1, var2, var3, var4):
        for i in range(len(var1)):
            file.write(' '.join([str(x[i]) for x in [var1, var2, var3, var4]]) + '\n')

    def write_three_floats(file, var1, var2, var3):
        for i in range(len(var1)):
            file.write(' '.join([str(x[i]) for x in [var1, var2, var3]]) + '\n')

    with open(filename, 'w') as file:
        file.write('nlambda' + '\n')
        file.write(str(data['nlambda']) + '\n')
        file.write('lambda' + '\n')
        for i in range(data['nlambda']):
            file.write(str(data['lambda_grid'][i]) + '\n')

        file.write('ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q' + '\n')
        file.write(' '.join([str(x) for x in [
            data['ntgrid'], data['nperiod'], data['ntheta'],
            data['drhodpsi'], data['rmaj'], data['shat'],
            data['kxfac'], data['q']]]) + '\n')

        file.write('gbdrift gradpar grho tgrid' + '\n')
        write_four_floats(file, data['gbdrift'], data['gradpar'],
                          data['grho'], data['theta_grid'])

        file.write('cvdrift gds2 bmag tgrid' + '\n')
        write_four_floats(file, data['cvdrift'], data['gds2'],
                          data['bmag'], data['theta_grid'])

        file.write('gds21 gds22 tgrid' + '\n')
        write_three_floats(file, data['gds21'], data['gds22'], data['theta_grid'])

        file.write('cvdrift0 gbdrift0 tgrid' + '\n')
        write_three_floats(file, data['cvdrift0'], data['gbdrift0'], data['theta_grid'])

        file.write('rplot rprime tgrid' + '\n')
        write_three_floats(file, data['rplot'], data['rprime'], data['theta_grid'])

        file.write('zplot zprime tgrid' + '\n')
        write_three_floats(file, data['zplot'], data['zprime'], data['theta_grid'])

        file.write('aplot aprime tgrid' + '\n')
        write_three_floats(file, data['aplot'], data['aprime'], data['theta_grid'])


def consistency_check_theta_grid_data(data, fix = False,
                                      min_bounce_tolerance = 1.0e-12, min_bounce_tolerance_large = 1.0e-8):

    passing = True

    # Check bset values are monotonically decreasing
    #
    # This is fixable -- we can just sort the data (and lambda_grid)
    if not all(np.diff(data['bset']) < 0):
        print('Bset values not monotonically decreasing :')
        print(data['bset'])
        if fix:
            data['bset'] = sort(data['bset'])[::-1]
            data['lambda_grid'] = sort(data['lambda_grid'])
        else:
            passing = False

    # Check bset/lambda values bounce on the grid.
    #
    # This is perhaps fixable if the difference isn't too large -- we
    # replace bset with the nearest bmag and then recalculate
    # lambda_grid
    for i in range(data['nbset']):
        min_bounce_value = np.abs(1 - data['bmag']*data['lambda_grid'][i]).min()
        if min_bounce_value > min_bounce_tolerance:
            min_loc = np.abs(data['bmag']-data['bset'][i]).argmin()
            print('Lambda value not bouncing on grid points:')
            print('Lambda :', data['lambda_grid'][i])
            print('Bset   :', data['bset'][i])
            print('Nearest bmag :', data['bmag'][min_loc])
            if fix:
                if min_bounce_value > min_bounce_tolerance_large:
                    print('Lambda grid correction would be too large -- not fixing')
                    passing = False
                else:
                    data['bset'][i] = data['bmag'][min_loc]
                    data['lambda_grid'][i] = 1/data['bset'][i]
            else:
                passing = False

    # Check bset is unique -- might want to move this after any potential auto-fixes to data
    if len(set(data['bset'])) < len(data['bset']):
        print('Duplicate bset values detected :')
        print(data['bset'])
        passing = False

    # Check the maximum magnetic field occurs at the end of the grid
    if data['bmag'].max() != data['bmag'][0]:
        print('First value of bmag does not match the maximum field:')
        print(data['bmag'])
        passing = False
    if data['bmag'].max() != data['bmag'][-1]:
        print('Last value of bmag does not match the maximum field:')
        print(data['bmag'])
        passing = False

    # Check grid sizes
    if (2*data['nperiod']-1) * data['ntheta'] + 1 != (2*data['ntgrid']+1):
        print('Inconsistent theta grid sizes:')
        print('nperiod :', data['nperiod'])
        print('ntheta  :', data['ntheta'])
        print('ntgrid  :', data['ntgrid'])
        print('First grid size :', (2*data['nperiod']-1) * data['ntheta'] + 1)
        print('Second grid size :', (2*data['ntgrid']+1))
        passing = False

    return passing

