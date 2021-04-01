# a func. to generate plots for Si.py example

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits import mplot3d
from matplotlib.colors import LightSource
import seaborn as sns

def generate_figs():
    sns.set ()
    sns.set_context ( "paper", font_scale = 2, rc = {"lines.linewidth": 4} )
    sns.set_style ( "ticks", {"xtick.major.size": 2, "ytick.major.size": 2} )
    fig_1 = plt.figure ( figsize = (6.5, 4.5) )
    ax_1 = fig_1.add_subplot ( 111 )
    ax_1.plot ( Tmp[0], h[0], 'o', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_1.yaxis.set_major_formatter ( FormatStrFormatter ( '%.2f' ) )
    ax_1.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_1.tick_params ( axis = "x", labelsize = 16 )
    ax_1.set_ylabel ( 'Band gap (eV)', fontsize = 16, labelpad = 10 )
    ax_1.tick_params ( axis = "y", labelsize = 16 )
    fig_1.tight_layout ()
    fig_1 = plt.figure ( figsize = (6.5, 4.5) )
    ax_1 = fig_1.add_subplot ( 111 )
    ax_1.plot ( Tmp[0], alpha[0], 'o', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_1.yaxis.set_major_formatter ( FormatStrFormatter ( '%.2f' ) )
    ax_1.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_1.tick_params ( axis = "x", labelsize = 16 )
    ax_1.set_ylabel ( 'Nonparabolic term (eV$^{-1}$)', fontsize = 16, labelpad = 10 )
    ax_1.tick_params ( axis = "y", labelsize = 16 )
    fig_1.tight_layout ()
    fig_3 = plt.figure ( figsize = (6.5, 4.5) )
    ax_3 = fig_3.add_subplot ( 111 )
    ax_3.plot ( e[0], DoS[0], 'None', linestyle = '-', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_3.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_3.set_xlabel ( 'Energy (eV)', fontsize = 16, labelpad = 10 )
    ax_3.tick_params ( axis = "x", labelsize = 16 )
    ax_3.set_ylabel ( 'Density of state (#/eV/m$^3$)', fontsize = 16, labelpad = 10 )
    ax_3.tick_params ( axis = "y", labelsize = 16 )
    ax_3.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_3.tight_layout ()
    fig_4 = plt.figure ( figsize = (6.5, 4.5) )
    ax_4 = fig_4.add_subplot ( 111 )
    ax_4.plot ( e[0], -1 * gVel[0] * 1e-5, 'None', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_4.yaxis.set_major_formatter ( FormatStrFormatter ( '%.0f' ) )
    ax_4.set_xlabel ( 'Energy (eV)', fontsize = 16, labelpad = 10 )
    ax_4.tick_params ( axis = "x", labelsize = 16 )
    ax_4.set_ylabel ( 'Group velocity (x10$^5$ m/s)', fontsize = 16, labelpad = 10 )
    fig_4.tight_layout ()
    fig_5 = plt.figure ( figsize = (6.5, 4.5) )
    ax_5 = fig_5.add_subplot ( 111 )
    ax_5.plot ( band[1::, ], 'None', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_5.yaxis.set_major_formatter ( FormatStrFormatter ( '%.0f' ) )
    ax_5.set_xlabel ( 'BZ tour', fontsize = 16, labelpad = 10 )
    ax_5.tick_params ( axis = "x", labelsize = 16 )
    ax_5.set_ylabel ( 'Energy (eV)', fontsize = 16 )
    ax_5.tick_params ( axis = "y", labelsize = 16 )
    ax_5.set_xticks ( [0, 199, 399, 599, 799] )
    ax_5.set_xticklabels ( ["W", "L", "$\Gamma$", "X", "W"] )
    fig_5.tight_layout ()
    fig_6 = plt.figure ( figsize = (6.5, 4.5) )
    ax_6 = fig_6.add_subplot ( 111 )
    ax_6.plot ( e[0], dis_no_inc[0], 'None', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_6.plot ( e[0], dis[0], 'None', linestyle = '-', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_6.plot ( e[0], dis_no_inc[-1], 'None', linestyle = '-.', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_6.plot ( e[0], dis[-1], 'None', linestyle = '-.', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_6.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_6.set_xlabel ( 'Energy (eV)', fontsize = 16, labelpad = 10 )
    ax_6.tick_params ( axis = "x", labelsize = 16 )
    ax_6.set_ylabel ( 'Fermi distribution', fontsize = 16, labelpad = 10 )
    ax_6.tick_params ( axis = "y", labelsize = 16 )
    fig_6.tight_layout ()
    fig_7 = plt.figure ( figsize = (6.5, 4.5) )
    ax_7 = fig_7.add_subplot ( 111 )
    ax_7.plot ( e[0], dfdE_no_inc[0], 'None', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_7.plot ( e[0], dfdE[0], 'None', linestyle = '-', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_7.plot ( e[0], dfdE_no_inc[-1], 'None', linestyle = '-.', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_7.plot ( e[0], dfdE[-1], 'None', linestyle = '-.', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_7.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_7.set_xlabel ( 'Energy (eV)', fontsize = 16, labelpad = 10 )
    ax_7.tick_params ( axis = "x", labelsize = 16 )
    ax_7.set_ylabel ( 'Fermi window (eV$^{-1}$)', fontsize = 16 )
    ax_7.tick_params ( axis = "y", labelsize = 16 )
    fig_7.tight_layout ()
    fig_8 = plt.figure ( figsize = (6.5, 4.5) )
    ax_8 = fig_8.add_subplot ( 111 )
    ax_8.plot ( Tmp[0], cc_no_inc[0], 'o', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_8.plot ( Tmp[0], JD_n_no_inc[0], 'o', linestyle = '-.', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_8.plot ( Tmp[0], cc_sc_no_inc[0], 'o', linestyle = '-.', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_8.plot ( Tmp[0], cc[0], 'o', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_8.plot ( Tmp[0], JD_n[0], 'o', linestyle = '-', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_8.plot ( Tmp[0], cc_sc[0], 'o', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_8.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_8.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_8.tick_params ( axis = "x", labelsize = 16 )
    ax_8.set_ylabel ( 'Carrier concentration( #/m$^{3}$)', fontsize = 16 )
    ax_8.tick_params ( axis = "y", labelsize = 16 )
    ax_8.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_8.tight_layout ()
    fig_9 = plt.figure ( figsize = (6.5, 4.5) )
    ax_9 = fig_9.add_subplot ( 111 )
    ax_9.plot ( Tmp[0], cc_no_inc[0], 'D', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_9.plot ( Tmp[0], cc_sc_no_inc[0], 'D', linestyle = '-', color = 'maroon',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'maroon',
                markeredgewidth = 1 )
    ax_9.plot ( Tmp[0], cc[0], 'o', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_9.plot ( Tmp[0], cc_sc[0], 'o', linestyle = '-', color = 'steelblue',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'steelblue',
                markeredgewidth = 1 )
    ax_9.plot ( Tmp[0], cc_direction_down[0], 'p', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_9.plot ( Tmp[0], cc_sc_direction_down[0], 'p', linestyle = '-', color = 'indigo',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'indigo',
                markeredgewidth = 1 )
    ax_9.plot ( Tmp[0], cc_1pct[0], 's', linestyle = 'None', color = 'black',
                markersize = 12, linewidth = 1.5,
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 1, zorder = 0 )
    ax_9.plot ( Tmp[0], cc_sc_1pct[0], 's', linestyle = '-', color = 'olive',
                markersize = 6, linewidth = 1.5,
                markerfacecolor = 'white',
                markeredgecolor = 'olive',
                markeredgewidth = 1 )
    ax_9.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_9.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_9.tick_params ( axis = "x", labelsize = 16 )
    ax_9.set_ylabel ( 'Carrier concentration( #/m$^{3}$)', fontsize = 16 )
    ax_9.tick_params ( axis = "y", labelsize = 16 )
    ax_9.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_9.tight_layout ()
    np.savetxt ( 'fig3_cc_no_inc', cc_no_inc[0] )
    np.savetxt ( 'fig3_cc_sc_no_inc', cc_sc_no_inc[0] )
    np.savetxt ( 'fig3_cc', cc[0] )
    np.savetxt ( 'fig3_cc_sc', cc_sc[0] )
    np.savetxt ( 'fig3_cc_direction_down', cc_direction_down[0] )
    np.savetxt ( 'fig3_cc_sc_direction_down', cc_sc_direction_down[0] )
    np.savetxt ( 'fig3_cc_1pct', cc_1pct[0] )
    np.savetxt ( 'fig3_cc_sc_1pct', cc_sc_1pct[0] )
    fig_10 = plt.figure ( figsize = (6.5, 4.5) )
    ax_10 = fig_10.add_subplot ( 111 )
    ax_10.plot ( Tmp[0], fermi_no_inc[0], 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1, zorder = 10 )
    ax_10.plot ( Tmp[0], JD_f_no_inc[0], 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_10.plot ( Tmp[0], fermi[0], 'o', linestyle = '-.', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1, zorder = 10 )
    ax_10.plot ( Tmp[0], JD_f[0], 'o', linestyle = '-.', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_10.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_10.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_10.tick_params ( axis = "x", labelsize = 16 )
    ax_10.set_ylabel ( 'E$_f$ (eV)', fontsize = 16 )
    ax_10.tick_params ( axis = "y", labelsize = 16 )
    ax_10.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_10.tight_layout ()
    fig_11 = plt.figure ( figsize = (6.5, 4.5) )
    ax_11 = fig_11.add_subplot ( 111 )
    ax_11.plot ( Tmp[0], fermi_no_inc[0], 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1, zorder = 10 )
    ax_11.plot ( Tmp[0], fermi[0], 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1, zorder = 10 )
    ax_11.plot ( Tmp[0], fermi_direction_down[0], 'o', linestyle = '-', color = 'indigo',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'indigo',
                 markeredgewidth = 1, zorder = 10 )
    ax_11.plot ( Tmp[0], fermi_1pct[0], 'o', linestyle = '-', color = 'olive',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1, zorder = 10 )
    ax_11.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_11.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_11.tick_params ( axis = "x", labelsize = 16 )
    ax_11.set_ylabel ( 'E$_f$ (eV)', fontsize = 16 )
    ax_11.tick_params ( axis = "y", labelsize = 16 )
    ax_11.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_11.tight_layout ()
    np.savetxt ( 'fig3_fermi_no_inc', fermi_no_inc[0] )
    np.savetxt ( 'fig3_fermi', fermi[0] )
    np.savetxt ( 'fig3_fermi_direction_down', fermi_direction_down[0] )
    np.savetxt ( 'fig3_fermi_1pct', fermi_1pct[0] )
    fig_12 = plt.figure ( figsize = (6.5, 4.5) )
    ax_12 = fig_12.add_subplot ( 111 )
    ax_12.plot ( Tmp[0], m_CB_no_inc[0] / Si.me, 'o', linestyle = '-', color = 'maroon',
                 markersize = 12, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1, zorder = 10 )
    ax_12.plot ( Tmp[0], m_CB_inc[0] / Si.me, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1, zorder = 10 )
    ax_12.plot ( Tmp[0], m_CB_inc_direction_down[0] / Si.me, 'o', linestyle = '-', color = 'indigo',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'indigo',
                 markeredgewidth = 1, zorder = 10 )
    ax_12.plot ( Tmp[0], m_CB_inc_1pct[0] / Si.me, 'o', linestyle = '-', color = 'olive',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1, zorder = 10 )
    ax_12.yaxis.set_major_formatter ( ScalarFormatter () )
    ax_12.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_12.tick_params ( axis = "x", labelsize = 16 )
    ax_12.set_ylabel ( 'Effective mass (m$_e$)', fontsize = 16 )
    ax_12.tick_params ( axis = "y", labelsize = 16 )
    ax_12.ticklabel_format ( axis = "y", style = "sci", scilimits = None )
    fig_12.tight_layout ()
    fig_13 = plt.figure ( figsize = (6.5, 4.5) )
    ax_13 = fig_13.add_subplot ( 111 )
    ax_13.plot ( Tmp[0], Coeff_no_inc[0] * 1e-5, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_13.plot ( Tmp[0], Coeff_no_np[0] * 1e-5, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_13.plot ( Tmp[0], (1 - vfrac) * Coeff[0] * 1e-5, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_13.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 2] * 1e3 * 1e-5, '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_13.plot ( ExpData_SiCfrac_5pct_direction_up[:, 0], ExpData_SiCfrac_5pct_direction_up[:, 2] * 1e3 * 1e-5, '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_13.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_13.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_13.tick_params ( axis = "x", labelsize = 16 )
    ax_13.set_ylabel ( 'Conductivity(x10$^5$ S/m)', fontsize = 16 )
    ax_13.tick_params ( axis = "y", labelsize = 16 )
    fig_13.tight_layout ()
    fig_14 = plt.figure ( figsize = (6.5, 4.5) )
    ax_14 = fig_14.add_subplot ( 111 )
    ax_14.plot ( Tmp[0], Coeff_no_inc[0] * 1e-5, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_14.plot ( Tmp[0], Coeff_direction_down_no_np[0] * 1e-5, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_14.plot ( Tmp[0], (1 - vfrac) * Coeff_direction_down[0] * 1e-5, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_14.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 2] * 1e3 * 1e-5, '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_14.plot ( ExpData_SiCfrac_5pct_direction_down[:, 0], ExpData_SiCfrac_5pct_direction_down[:, 2] * 1e3 * 1e-5,
                 '-o', linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_14.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_14.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_14.tick_params ( axis = "x", labelsize = 16 )
    ax_14.set_ylabel ( 'Conductivity(x10$^5$ S/m)', fontsize = 16 )
    ax_14.tick_params ( axis = "y", labelsize = 16 )
    fig_14.tight_layout ()
    fig_15 = plt.figure ( figsize = (6.5, 4.5) )
    ax_15 = fig_15.add_subplot ( 111 )
    ax_15.plot ( Tmp[0], Coeff_no_inc[0] * 1e-5, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_15.plot ( Tmp[0], Coeff_no_np_1pct[0] * 1e-5, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_15.plot ( Tmp[0], (1 - vfrac / 5) * Coeff_1pct[0] * 1e-5, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_15.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 2] * 1e3 * 1e-5, '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_15.plot ( ExpData_SiCfrac_1pct_direction_up[:, 0], ExpData_SiCfrac_1pct_direction_up[:, 2] * 1e3 * 1e-5, '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_15.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_15.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_15.tick_params ( axis = "x", labelsize = 16 )
    ax_15.set_ylabel ( 'Conductivity(x10$^5$ S/m)', fontsize = 16 )
    ax_15.tick_params ( axis = "y", labelsize = 16 )
    fig_15.tight_layout ()
    fig_16 = plt.figure ( figsize = (6.5, 4.5) )
    ax_16 = fig_16.add_subplot ( 111 )
    ax_16.plot ( Tmp[0], Coeff_no_inc[1] * 1e6, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_16.plot ( Tmp[0], Coeff_no_np[1] * 1e6, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_16.plot ( Tmp[0], Coeff[1] * 1e6, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_16.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 3], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_16.plot ( ExpData_SiCfrac_5pct_direction_up[:, 0], ExpData_SiCfrac_5pct_direction_up[:, 3], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_16.yaxis.set_major_formatter ( FormatStrFormatter ( '%.0f' ) )
    ax_16.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_16.tick_params ( axis = "x", labelsize = 16 )
    ax_16.set_ylabel ( 'Seebeck($\mu$V/K)', fontsize = 16 )
    ax_16.tick_params ( axis = "y", labelsize = 16 )
    fig_16.tight_layout ()
    fig_17 = plt.figure ( figsize = (6.5, 4.5) )
    ax_17 = fig_17.add_subplot ( 111 )
    ax_17.plot ( Tmp[0], Coeff_no_inc[1] * 1e6, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_17.plot ( Tmp[0], Coeff_direction_down_no_np[1] * 1e6, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_17.plot ( Tmp[0], Coeff_direction_down[1] * 1e6, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_17.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 3], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_17.plot ( ExpData_SiCfrac_5pct_direction_down[:, 0], ExpData_SiCfrac_5pct_direction_down[:, 3], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_17.yaxis.set_major_formatter ( FormatStrFormatter ( '%.0f' ) )
    ax_17.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_17.tick_params ( axis = "x", labelsize = 16 )
    ax_17.set_ylabel ( 'Seebeck($\mu$V/K)', fontsize = 16 )
    ax_17.tick_params ( axis = "y", labelsize = 16 )
    fig_17.tight_layout ()
    fig_18 = plt.figure ( figsize = (6.5, 4.5) )
    ax_18 = fig_18.add_subplot ( 111 )
    ax_18.plot ( Tmp[0], Coeff_no_inc[1] * 1e6, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_18.plot ( Tmp[0], Coeff_no_np_1pct[1] * 1e6, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_18.plot ( Tmp[0], Coeff_1pct[1] * 1e6, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_18.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 3], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_18.plot ( ExpData_SiCfrac_1pct_direction_up[:, 0], ExpData_SiCfrac_1pct_direction_up[:, 3], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_18.yaxis.set_major_formatter ( FormatStrFormatter ( '%.0f' ) )
    ax_18.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_18.tick_params ( axis = "x", labelsize = 16 )
    ax_18.set_ylabel ( 'Seebeck($\mu$V/K)', fontsize = 16 )
    ax_18.tick_params ( axis = "y", labelsize = 16 )
    fig_18.tight_layout ()
    fig_19 = plt.figure ( figsize = (6.5, 4.5) )
    ax_19 = fig_19.add_subplot ( 111 )
    ax_19.plot ( Tmp[0], Coeff_no_inc[2] * 1e3, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_19.plot ( Tmp[0], Coeff_no_np[2] * 1e3, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_19.plot ( Tmp[0], (1 - vfrac) * Coeff[2] * 1e3, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_19.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 4], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_19.plot ( ExpData_SiCfrac_5pct_direction_up[:, 0], ExpData_SiCfrac_5pct_direction_up[:, 4], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_19.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_19.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_19.tick_params ( axis = "x", labelsize = 16 )
    ax_19.set_ylabel ( 'Power factor(mW/mK$^2$)', fontsize = 16, labelpad = 10 )
    ax_19.tick_params ( axis = "y", labelsize = 16 )
    fig_19.tight_layout ()
    np.savetxt ( 'fig5_PF_no_np', Coeff_no_np[2] * 1e3 )
    np.savetxt ( 'fig5_PF', (1 - vfrac) * Coeff[2] * 1e3 )
    np.savetxt ( 'fig5_Temp_ExpData_SiCfrac_5pct_direction_up', ExpData_SiCfrac_5pct_direction_up[:, 0] )
    np.savetxt ( 'fig5_PF_ExpData_SiCfrac_5pct_direction_up', ExpData_SiCfrac_5pct_direction_up[:, 4] )
    fig_20 = plt.figure ( figsize = (6.5, 4.5) )
    ax_20 = fig_20.add_subplot ( 111 )
    ax_20.plot ( Tmp[0], Coeff_no_inc[2] * 1e3, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_20.plot ( Tmp[0], Coeff_direction_down_no_np[2] * 1e3, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_20.plot ( Tmp[0], (1 - vfrac) * Coeff_direction_down[2] * 1e3, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_20.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 4], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_20.plot ( ExpData_SiCfrac_5pct_direction_down[:, 0], ExpData_SiCfrac_5pct_direction_down[:, 4], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    np.savetxt ( 'fig5_PF_direction_down_no_np', Coeff_direction_down_no_np[2] * 1e3 )
    np.savetxt ( 'fig5_PF_direction_down', (1 - vfrac) * Coeff_direction_down[2] * 1e3 )
    np.savetxt ( 'fig5_Temp_ExpData_SiCfrac_5pct_direction_down', ExpData_SiCfrac_5pct_direction_down[:, 0] )
    np.savetxt ( 'fig5_PF_ExpData_SiCfrac_5pct_direction_down', ExpData_SiCfrac_5pct_direction_down[:, 4] )
    ax_20.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_20.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_20.tick_params ( axis = "x", labelsize = 16 )
    ax_20.set_ylabel ( 'Power factor(mW/mK$^2$)', fontsize = 16, labelpad = 10 )
    ax_20.tick_params ( axis = "y", labelsize = 16 )
    fig_20.tight_layout ()
    fig_21 = plt.figure ( figsize = (6.5, 4.5) )
    ax_21 = fig_21.add_subplot ( 111 )
    ax_21.plot ( Tmp[0], Coeff_no_inc[2] * 1e3, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_21.plot ( Tmp[0], Coeff_no_np_1pct[2] * 1e3, 'o', linestyle = '-', color = 'tan',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'tan',
                 markeredgewidth = 1 )
    ax_21.plot ( Tmp[0], (1 - vfrac / 5) * Coeff_1pct[2] * 1e3, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_21.plot ( ExpData_SiCfra_0pct_direction_up[:, 0], ExpData_SiCfra_0pct_direction_up[:, 4], '-o',
                 linestyle = 'None', color = 'black',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'black',
                 markeredgewidth = 1 )
    ax_21.plot ( ExpData_SiCfrac_1pct_direction_up[:, 0], ExpData_SiCfrac_1pct_direction_up[:, 4], '-o',
                 linestyle = 'None', color = 'olive',
                 markersize = 5, linewidth = 4,
                 markerfacecolor = 'white',
                 markeredgecolor = 'olive',
                 markeredgewidth = 1 )
    ax_21.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_21.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_21.tick_params ( axis = "x", labelsize = 16 )
    ax_21.set_ylabel ( 'Power factor(mW/mK$^2$)', fontsize = 16, labelpad = 10 )
    ax_21.tick_params ( axis = "y", labelsize = 16 )
    fig_21.tight_layout ()
    fig_22 = plt.figure ( figsize = (6.5, 4.5) )
    ax_22 = fig_22.add_subplot ( 111 )
    ax_22.plot ( Tmp[0], Coeff_no_inc[3], 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_22.plot ( Tmp[0], Coeff[3], 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_22.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_22.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_22.tick_params ( axis = "x", labelsize = 16 )
    ax_22.set_ylabel ( '$\kappa_e$(W/mK)', fontsize = 16, labelpad = 10 )
    ax_22.tick_params ( axis = "y", labelsize = 16 )
    fig_22.tight_layout ()
    fig_23 = plt.figure ( figsize = (6.5, 4.5) )
    ax_23 = fig_23.add_subplot ( 111 )
    ax_23.plot ( Tmp[0], Coeff_no_inc[4], 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_23.plot ( Tmp[0], Coeff[4], 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_23.yaxis.set_major_formatter ( FormatStrFormatter ( '%.2f' ) )
    ax_23.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_23.tick_params ( axis = "x", labelsize = 16 )
    ax_23.set_ylabel ( '$\Delta_1$(eV)', fontsize = 16, labelpad = 10 )
    ax_23.tick_params ( axis = "y", labelsize = 16 )
    fig_23.tight_layout ()
    fig_24 = plt.figure ( figsize = (6.5, 4.5) )
    ax_24 = fig_24.add_subplot ( 111 )
    ax_24.plot ( Tmp[0], Coeff_no_inc[5], 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_24.plot ( Tmp[0], Coeff[5], 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_24.yaxis.set_major_formatter ( FormatStrFormatter ( '%.2f' ) )
    ax_24.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_24.tick_params ( axis = "x", labelsize = 16 )
    ax_24.set_ylabel ( '$\Delta_2$([eV]$^2$)', fontsize = 16, labelpad = 10 )
    ax_24.tick_params ( axis = "y", labelsize = 16 )
    fig_24.tight_layout ()
    fig_25 = plt.figure ( figsize = (6.5, 4.5) )
    ax_25 = fig_25.add_subplot ( 111 )
    ax_25.plot ( Tmp[0], Coeff_no_inc[6] * 1e8, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_25.plot ( Tmp[0], Coeff[6] * 1e8, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_25.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_25.set_xlabel ( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_25.tick_params ( axis = "x", labelsize = 16 )
    ax_25.set_ylabel ( 'Lorenz number (x$10^{-8}$[V/K]$^2$)', fontsize = 16, labelpad = 10 )
    ax_25.tick_params ( axis = "y", labelsize = 16 )
    fig_25.tight_layout ()
    fig_26 = plt.figure ( figsize = (6.5, 4.5) )
    ax_26 = fig_26.add_subplot ( 111 )
    ax_26.plot ( Tmp[0], LD_nondegenerate_no_inc[0] * 1e9, 'o', linestyle = '-.', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_26.plot ( Tmp[0], LD_no_inc[0] * 1e9, 'o', linestyle = '-', color = 'maroon',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'maroon',
                 markeredgewidth = 1 )
    ax_26.plot ( Tmp[0], LD_nondegenerate[0] * 1e9, 'o', linestyle = '-.', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_26.plot ( Tmp[0], LD[0] * 1e9, 'o', linestyle = '-', color = 'steelblue',
                 markersize = 6, linewidth = 1.5,
                 markerfacecolor = 'white',
                 markeredgecolor = 'steelblue',
                 markeredgewidth = 1 )
    ax_26.yaxis.set_major_formatter ( FormatStrFormatter ( '%.1f' ) )
    ax_26.set_xlabel( 'Temperature (K)', fontsize = 16, labelpad = 10 )
    ax_26.tick_params( axis = "x", labelsize = 16 )
    ax_26.set_ylabel( 'Debyle length (nm)', fontsize = 16, labelpad = 10 )
    ax_26.tick_params( axis = "y", labelsize = 16 )
    fig_26.tight_layout()
    plt.show()
