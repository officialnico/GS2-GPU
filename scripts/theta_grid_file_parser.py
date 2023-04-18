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
            cvdrift,  gds2, bmag, theta_grid_tmp = read_four_floats(file, 2*ntgrid+1)

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
