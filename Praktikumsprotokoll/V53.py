import pylab
import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
import numpy as np
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'lualatex',
    'pgf.preamble':r'\usepackage{unicode-math}\usepackage{siunitx}\usepackage{xfrac}',
})
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)
from scipy import stats
from scipy.stats import sem
from scipy.optimize import curve_fit
import scipy.constants as const

def rf(W_ex, W_lit):
    return noms(abs(((W_ex - W_lit)/W_lit)*100))

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

def lambda_G(a, n):
    return 2*a/n


# Formeln ---------------------------------------


# Literaturwerte ---------------------------------------  


# Importierung Messdaten ------------------------------

# Modenmessung
U_01 = 220
U_11 = 200
U_21 = 240
A_01 = 6.6
f_01 = 9.001e9

U_02 = 140
U_12 = 120
U_22 = 150
A_02 = 6.8
f_02 = 9.004e9

U_03 = 80
U_13 = 70
U_23 = 90
A_03 = 6.2
f_03 = 9.01e9



# Elektronische Abstimmung
U_e0 = 140
U_e1 = 120
U_e2 = 150
f_e0 = 9.004e9
f_e1 = 8.969e9
f_e2 = 9.044e9

f_breite = f_e0-f_e1
abstimm = f_breite/(U_e0-U_e1)
# Frequenz
f = 9e9
a_0 = ufloat(22.860, 0.046)*1e-03
l_1 = 60.3e-3
l_2 = 84.6e-3
lambda_g = 2*(l_2 - l_1)
f_micro = 3e08 * unp.sqrt((1/lambda_g)**2+(1/(2*a_0))**2)

v_phase = f_micro* lambda_g/const.c
print(v_phase)

# Dämpfung
P_SWR, l, P_d = np.genfromtxt('Daempfung.txt', unpack=True)

# SWR 1. Methode
d = np.array([3, 5, 7, 9])*1e-3
swr_1 = np.array([1.2, 1.65, 3.5, np.inf])

# SWR 2. Methode
d_1 = ufloat(64.3e-3,1e-3)
d_2 = ufloat(62.5e-3,1e-3)


min1 = 59.9e-3
min2 = 84.3e-3
lambda_g2 = 2*(min2-min1)
swr_2 = unp.sqrt(1+1/(unp.sin(np.pi *(d_1-d_2)/lambda_g2)**2))



# SWR 3. Methode
A_1 = ufloat(20,2)
A_2 = ufloat(46,2)
swr_3 = 10**((A_2-A_1)/20)

# Plot für Modenmessung ---------------------------------


def square(x, a, b , c):
    return a*x**2 + b*x + c


U_1 = np.array([U_11, U_01, U_21])

I_1 = np.array([0, A_01, 0])

par1 = np.polyfit(U_1, I_1, deg=2, cov=False)

a_U1 = par1[0]

b_U1 = par1[1]

c_U1 = par1[2]


U_2 = np.array([U_02, U_12, U_22])

I_2 = np.array([A_02, 0, 0])

par2 = np.polyfit(U_2, I_2, deg=2, cov=False)

a_U2 = par2[0]

b_U2 = par2[1]

c_U2 = par2[2]


U_3 = np.array([U_03, U_13, U_23])

I_3 = np.array([A_03, 0, 0])

par3 = np.polyfit(U_3, I_3, deg=2, cov=False)

a_U3 = par3[0]

b_U3 = par3[1]

c_U3 = par3[2]

fig, ax = plt.subplots(constrained_layout=True)

xplt = np.linspace(50, 260, 1000)

ax.plot(U_1, I_1, 'kx', label=r'Messdaten des 1. Modus')
ax.plot(xplt, square(xplt,*par1), 'r-', label=r'Fit an Messdaten des 1. Modus')
ax.plot(U_2, I_2, 'mx', label=r'Messdaten des 2. Modus')
ax.plot(xplt, square(xplt,*par2), 'b-', label=r'Fit an Messdaten des 2. Modus')
ax.plot(U_3, I_3, 'cx', label=r'Messdaten des 3. Modus')
ax.plot(xplt, square(xplt,*par3), 'g-', label=r'Fit an Messdaten des 3. Modus')
ax.set_xlabel(r'Reflektorspannung $U \:/\: \si{\volt}$')
ax.set_ylabel(r'Spannung am Detektor $U \:/\: \si{\volt}$')
plt.ylim(0, 10)
plt.legend(loc="best")
plt.grid()
plt.savefig("build/moden.pdf")
plt.close()


# Plot für Dämpfung ------------------------------

def dB(x,a,b):
    return 10**(a*x) + b

popt, pcov = curve_fit(
    dB,
    l,
    P_SWR,
    sigma=None,
    p0=[1,0]
    )

perr = np.sqrt(np.diag(pcov))
a1 = ufloat(popt[0], perr[0])
b1 = ufloat(popt[1], perr[1])


qopt, qcov = curve_fit(
    dB,
    l,
    P_d,
    sigma=None,
    p0=[1,0]
    )

qerr = np.sqrt(np.diag(qcov))
a2 = ufloat(qopt[0], qerr[0])
b2 = ufloat(qopt[1], qerr[1])

print(a1,a2)
print(b1,b2)

xplt = np.linspace(np.amin(l), np.amax(l), 1000)


fig, ax = plt.subplots(constrained_layout=True)

ax.plot(l, P_SWR, 'kx', label=r'An SWR-Meter gemessene Dämpfung')
ax.plot(l, P_d, 'mx', label=r'Dämpfung nach Eichkurve')
ax.plot(xplt, dB(xplt,*popt), 'r-', label=r'Fit an Messdaten')
ax.plot(xplt, dB(xplt,*qopt), 'b-', label=r'Fit an Eichkurve')
ax.set_xlabel(r'Schraubenstellung Dämpfungsglied $\:/\: \si{\milli\meter}$')
ax.set_ylabel(r'Dämpfung $\:/\: \si{\decibel}$')
plt.xticks([3.2e-3,3.4e-3,3.6e-3,3.8e-3,4.0e-3], [r"$3.2$",r"$3.4$",r"$3.6$",r"$3.8$",r"$4.0$"])
plt.legend(loc="best")
plt.grid()
plt.savefig("build/Dämpfung.pdf")
plt.close()

# Plot für Diskussion Schraubenoffset
b_eich = b2 - dB(3.2e-3, *qopt)
rf_eich = rf(b_eich, b1)
print('\n',b1, b_eich, rf_eich)

fig, ax = plt.subplots(constrained_layout=True)

ax.plot(l, P_SWR, 'kx', label=r'An SWR-Meter gemessene Dämpfung')
ax.plot(l, P_d, 'mx', label=r'Dämpfung nach Eichkurve')
ax.plot(xplt, dB(xplt,*popt), 'r-', label=r'Fit an Messdaten')
ax.plot(xplt, dB(xplt,*qopt), 'b-', label=r'Fit an Eichkurve')
ax.plot(xplt, dB(xplt,*qopt)-dB(3.2e-3,*qopt), 'g-', label=r'Fit mit Abzug von Schraubenoffset')
ax.set_xlabel(r'Schraubenstellung Dämpfungsglied $\:/\: \si{\milli\meter}$')
ax.set_ylabel(r'Dämpfung $\:/\: \si{\decibel}$')
plt.xticks([3.2e-3,3.4e-3,3.6e-3,3.8e-3,4.0e-3], [r"$3.2$",r"$3.4$",r"$3.6$",r"$3.8$",r"$4.0$"])
plt.legend(loc="best")
plt.grid()
plt.savefig("build/Dämpfung_offset.pdf")
plt.close()



 
# Erstellung von Tabellen

# Tabelle für  Modenmessung #######################################################################

first_row = ['Reflektorspannung' ,' ' ,' ' ,'Amplitude' , 'Frequenz']
second_row = ['$U_0\;/\;\si{\\volt} $' ,'$U_1\;/\;\si{\\volt}$', '$U_2\;/\;\si{\\volt}$', '$A_0\;/\;\si{\\volt}$', '$f_0\;/\;\si{\giga\hertz}$']
third_row = [U_01, U_11, U_21, A_01, f_01*1e-09]
fourth_row = [U_02, U_12, U_22, A_02, f_02*1e-09]
fifth_row = [U_03, U_13, U_23, A_03, f_03*1e-09]


tab_header = r'''
\begin{tabular}{c c S[table-format=3.3] S[table-format=3.3] S[table-format=2.2] }     %# X = Column Format, ex. 1.2
    \toprule
    { } & { } & {1. Modus} & {2. Modus} & {3. Modus} \\
    \cmidrule(lr{0.5em}){1-5}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0} & {1} & {2:3.3f} & {3:3.3f} & {4:2.2f} \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_moden.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row, fourth_row, fifth_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

# Tabelle für elektronische Abstimmung ###################################################################

first_row = ['Reflektorspannung $\;/\;\si{\\volt}$', 'Frequenz $\;/\;\si{\giga\hertz}$']
second_row = [U_e0, f_e0*1e-09]
third_row = [U_e1, f_e1*1e-09]
fourth_row = [U_e2, f_e2*1e-09]

tab_header = r'''
\begin{tabular}{c S[table-format=3.3] S[table-format=3.3] S[table-format=3.3] }     %# X = Column Format, ex. 1.2
    \toprule
    { } & {a)} & {b)} & {c)}  \\
    \cmidrule(lr{0.5em}){1-4}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0} & {1:3.3f} & {2:3.3f} & {3:3.3f}  \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_abstimmung.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row, fourth_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

# Tabelle für Frequenz- und Wellenlängenbestimmung ########################################################


first_row = [f*1e-09]
second_row = [l_1*1e03]
third_row = [l_2*1e03]
fourth_row = [lambda_g*1e03]
fifth_row = [a_0.n*1e03] 
sixth_row = [a_0.s*1e03]
seventh_row =[f_micro.n*1e-09]
eigth_row =[f_micro.s*1e-09]

tab_header = r'''
\begin{tabular}{S[table-format=1.1] S[table-format=2.1] S[table-format=2.1] S[table-format=2.2] S[table-format=2.3]@{${}\pm{}$}S[table-format=1.3] S[table-format=2.3]@{${}\pm{}$}S[table-format=1.3]}     %# X = Column Format, ex. 1.2
    \toprule
    {$f_{\text{Res}}\;/\;\si{\giga\hertz}$} & {1. Min$\;/\;\si{\milli\meter}$} & {2. Min$\;/\;\si{\milli\meter}$} & {$\lambda_g\;/\;\si{\milli\meter}$} & \multicolumn{2}{c}{a$\;/\;\si{\milli\meter}$} & \multicolumn{2}{c}{$f_{\text{mik}}\;/\;\si{\giga\hertz}$}  \\
    \cmidrule(lr{0.5em}){1-8}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0:1.1f} & {1:2.1f} & {2:2.1f} & {3:2.2f} & {4:2.3f} & {5:1.3f} & {6:2.3f} & {7:1.3f} \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_frequenz.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row, fourth_row, fifth_row, sixth_row, seventh_row, eigth_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

# Tabelle für Dämpfung #####################################################################

tab_header = r'''
\begin{tabular}{S[table-format=2.0] S[table-format=1.2] S[table-format=2.0] }     %# X = Column Format, ex. 1.2
    \toprule
    {SWR-Meter Ausschlag $\;/\;\si{\decibel}$} & {Mikrometereinstellung $\;/\;\si{\milli\meter} $} & {Dämpfung $\;/\;\si{\decibel}$} \\
    \cmidrule(lr{0.5em}){1-3}

'''

tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0:2.0f} & {1:1.2f} & {2:2.0f}  \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_daempfung.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(P_SWR, l*1e03, P_d):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

 #Tabelle für 1. Methode ##########################################################################

first_row = ['SWR']
second_row = [swr_1[0]]
third_row = [swr_1[1]]
fourth_row = [swr_1[2]]
fifth_row = [swr_1[3]]


tab_header = r'''
\begin{tabular}{c S[table-format=1.1] S[table-format=2.1] S[table-format=1.1] c}     %# X = Column Format, ex. 1.2
    \toprule
    {} & \multicolumn{4}{c}{Sondentiefe am Gleitschrauben-} \\
    {} &  \multicolumn{4}{c}{transformator in $\si{\milli\meter}$} \\
    \cmidrule(lr{0.5em}){1-5}
    {} & {3} & {5} & {7} & {9}   \\
    \cmidrule(lr{0.5em}){1-5}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0} & {1:1.1f} & {2:2.1f} & {3:1.1f} & {4:2.1f}   \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_methode1.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row, fourth_row, fifth_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',').replace('inf', '\infty') )
        a.write('\n')
    a.write(tab_footer)

 #Tabelle für 2. Methode ##########################################################################

first_row = [d_1.n*1e03]
second_row = [d_2.n*1e03]
third_row = [min1*1e03]
fourth_row = [min2*1e03]
fifth_row = [lambda_g2*1e03]
sixth_row = [noms(swr_2)]


tab_header = r'''
\begin{tabular}{S[table-format=2.1] S[table-format=2.1] S[table-format=2.1] S[table-format=2.1] S[table-format=2.2] S[table-format=1.3]}     %# X = Column Format, ex. 1.2
    \toprule
    {$d_1\;/\;\si{\milli\meter}$} & {$d_2\;/\;\si{\milli\meter}$} & {1.Min $\;/\;\si{\milli\meter}$} & {2.Min $\;/\;\si{\milli\meter}$} & {$\lambda_g\;/\;\si{\milli\meter}$} & {SWR}   \\
    \cmidrule(lr{0.5em}){1-6}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0:2.1f} & {1:2.1f} & {2:2.1f} & {3:2.1f} & {4:2.2f} & {5:1.3f}   \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_methode2.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row, fourth_row, fifth_row, sixth_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

# Tabelle für 3. Methode ##########################################################################

first_row = [A_1.n]
second_row = [A_2.n]
third_row = [swr_3.n]


tab_header = r'''
\begin{tabular}{S[table-format=2.0] S[table-format=2.0] S[table-format=2.3]}     %# X = Column Format, ex. 1.2
    \toprule
    {$A_1\;/\;\si{\deci\bel}$} & {$A_2\;/\;\si{\deci\bel}$} & {SWR}   \\
    \cmidrule(lr{0.5em}){1-3}

'''
tab_footer = r'''    
    \bottomrule
\end{tabular}
'''

_row_temp = r'         {0:2.0f} & {1:2.0f} & {2:2.3f}   \\'     # X = Column No. from 0, Y = Column Format

with open('build/table_methode3.tex', 'w') as a:
    a.write(tab_header)
    for row in zip(first_row, second_row, third_row):
        a.write(_row_temp.format(*row).replace('nan',' ').replace('.',',') )
        a.write('\n')
    a.write(tab_footer)

# Tabelle für 3. Methode ##########################################################################

#first_row = 

# tex file ##################################################################################


# tex file for f_breite

with open('build/f_breite.tex', 'w') as f:
    f.write(make_SI(abs(f_breite)*1e-06, r'\mega\hertz', figures=0))


# tex file for abstimm

with open('build/abstimm.tex', 'w') as f:
    f.write(make_SI(abstimm*1e-06, r'\mega\hertz\per\volt', figures=0))


# tex files for fit parameters

with open('build/a1.tex', 'w') as f:
    f.write(make_SI(a1, r''))


with open('build/a2.tex', 'w') as f:
    f.write(make_SI(a2, r''))


with open('build/b1.tex', 'w') as f:
    f.write(make_SI(b1, r'\decibel'))


with open('build/b2.tex', 'w') as f:
    f.write(make_SI(b2, r'\decibel'))


with open('build/a_U1.tex', 'w') as f:
    f.write(make_SI(a_U1, r'\volt\tothe{-1}', figures=2))

with open('build/b_U1.tex', 'w') as f:
    f.write(make_SI(b_U1, r'', figures=2))

with open('build/c_U1.tex', 'w') as f:
    f.write(make_SI(c_U1, r'\volt', figures=2))

with open('build/a_U2.tex', 'w') as f:
    f.write(make_SI(a_U2, r'\volt\tothe{-1}', figures=2))

with open('build/b_U2.tex', 'w') as f:
    f.write(make_SI(b_U2, r'', figures=2))

with open('build/c_U2.tex', 'w') as f:
    f.write(make_SI(c_U2, r'\volt', figures=2))

with open('build/a_U3.tex', 'w') as f:
    f.write(make_SI(a_U3, r'\volt\tothe{-1}', figures=2))

with open('build/b_U3.tex', 'w') as f:
    f.write(make_SI(b_U3, r'', figures=2))

with open('build/c_U3.tex', 'w') as f:
    f.write(make_SI(c_U3, r'\volt', figures=2))

with open('build/SWR_2.tex', 'w') as f:
    f.write(make_SI(ufloat(noms(swr_2),stds(swr_2)), r''))

with open('build/SWR_3.tex', 'w') as f:
    f.write(make_SI(swr_3, r''))

# tex file for v_phase

with open('build/v_phase.tex', 'w') as f:
    f.write(make_SI(v_phase, r'c', figures=2))

# tex file for b_eich,korrigiert and its rf

with open('build/b_eich.tex', 'w') as f:
    f.write(make_SI(b_eich, r'\decibel'))

with open('build/rf_eich.tex', 'w') as f:
    f.write(make_SI(rf_eich, r'\percent', figures=1))


# Berechnung von Grenzwellenlängen verschiedener Ordnung


with open('build/lambda_g_1.tex', 'w') as f:
    f.write(make_SI(lambda_G(a_0, 1)*1e03, r'\milli\meter', figures=1))

with open('build/lambda_g_2.tex', 'w') as f:
    f.write(make_SI(lambda_G(a_0, 2)*1e03, r'\milli\meter', figures=1))

with open('build/lambda_g_3.tex', 'w') as f:
    f.write(make_SI(lambda_G(a_0, 3)*1e03, r'\milli\meter', figures=1))

