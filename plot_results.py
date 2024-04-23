import numpy as np
import matplotlib.pyplot as plt

dpi=300
plt.style.use('plot_style.txt')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.markersize'] = 11 
plt.rcParams["figure.figsize"] = (6.4, 3.6)
 
if __name__ == '__main__':
    num_q = 3
    
    # data folder
    fdir = 'data/'
    # classical benchmark
    sa_tt_ksl = np.load(fdir + f'TT_{num_q}_spin_chain_sa_ksl.npy')
    sa_tt_taylor = np.load(fdir+f'TT_{num_q}_spin_chain_sa_Taylor'
                           + '.npy')
    sa_tt_chebyshev = np.load(fdir + f'TT_{num_q}_spin_chain_sa_'
                              + 'Chebyshev.npy')
    tt_time = np.load(fdir + f'TT_{num_q}_spin_chain_time.npy')

    # statevector results
    sa_statevector = np.load(fdir+f'{num_q}_spin_chain_SA_obs.npy')
    time = np.load(fdir + f'{num_q}_spin_chain_time.npy')

    # hadamard test results
    sa_corr = np.loadtxt(fdir + 'np_abs_correlation_'
                         'with_hadamard_test.csv', delimiter=';')
    sq_corr = np.loadtxt(fdir + 'np_sqrt_sum_squares_correlation_'
                         'with_hadamard_test.csv', delimiter=';')

    # plotting the data
    plt.plot(time, sa_tt_taylor, '-', label='TT-Taylor',
             linewidth=13)
    plt.plot(time, sa_tt_chebyshev, '-', label='TT-Chebyshev',
             linewidth=9)
    plt.plot(time, sa_tt_ksl, '-', label='TT-KSL',
             linewidth=5, color='gainsboro')
    plt.plot(time, sa_statevector, '--', label='Statevector')
    plt.plot(time, sa_corr, '*', label=r'Hadamard Test, $\mathfrak{Re}+i\mathfrak{Im}$')
    plt.xlabel(r'Time$')
    plt.ylabel('Absolute Value of Survival Amplitude, '
               r'$\left|\langle \psi | \psi \rangle \right|$')
    plt.xlim((min(time), max(time)))
    plt.yscale('log')
    plt.legend()
    plt.savefig('method_comparison.pdf', format='pdf', dpi=dpi)
    plt.savefig('method_comparison.png', format='png', dpi=dpi)
    plt.show()
