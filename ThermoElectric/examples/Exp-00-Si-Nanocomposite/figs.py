

# import ThermoElectric as TE


TE.visualize(array_to_plot = [engr, e_density], linestyle = '-', marker = None, pic_name = 'density_state.png',
              x_label = 'Energy (eV)', y_label = 'Density of State (#/m$^3$)')

TE.visualize(array_to_plot = [engr, -1*grp_velocity*1e-5], linestyle = '-', marker = None, pic_name = 'group_velocity.png', x_label = 'Energy (eV)', y_label = 'Velocity (x10$^5$m/s)')

TE.visualize(array_to_plot = [engr, fermi_dist_0pct_heating[0]], linestyle = '-', marker = None, pic_name = 'fermi_dist.png', x_label = 'Energy (eV)', y_label = 'Fermi Dist.')

TE.visualize(array_to_plot = [engr, fermi_dist_0pct_heating[1]], linestyle = '-', marker = None, pic_name = 'fermi_window.png', x_label = 'Energy (eV)', y_label = 'Fermi Window')

TE.visualize(array_to_plot = [tmpr, JD_0pct_heating[1, :][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'fermi_approx_0pct_cooling.png', x_label = 'Temperature (K)', y_label = 'Fermi Level (eV)')

TE.visualize(array_to_plot = [tmpr, fermi_0pct_heating[1, :][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'fermi_level_0pct.png', x_label = 'Temperature (K)', y_label = 'Fermi Level (eV)')

TE.visualize(array_to_plot = [tmpr, fermi_1pct_heating[1, :][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'fermi_level_1pct.png', x_label = 'Temperature (K)', y_label = 'Fermi Level (eV)')

TE.visualize(array_to_plot = [tmpr, fermi_5pct_heating[1, :][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'fermi_level_5pct_heating.png', x_label = 'Temperature (K)', y_label = 'Fermi Level (eV)')

TE.visualize(array_to_plot = [tmpr, fermi_5pct_cooling[1, :][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'fermi_level_5pct_cooling.png', x_label = 'Temperature (K)', y_label = 'Fermi Level (eV)')


TE.visualize(array_to_plot = [tmpr, prop_0pct['Electrical_conductivity'][np.newaxis, :]*1e-5], linestyle = '-', marker = 'o', pic_name = 'conductivity_0pct.png', x_label = 'Temperature (K)', y_label = 'Conductivity (x10$^5$ S/m)')


TE.visualize(array_to_plot = [tmpr, prop_0pct['Seebeck_coefficient'][np.newaxis, :] * 1e6], linestyle = '-', marker = 'o', pic_name = 'seebeck_0pct.png', x_label = 'Temperature (K)', y_label = 'Seebeck($\mu$V/K)')

TE.visualize(array_to_plot = [tmpr, prop_0pct['Power_factor'][np.newaxis, :] * 1e3], linestyle = '-', marker = 'o', pic_name = 'power_factor_0pct.png', x_label = 'Temperature (K)', y_label = 'Power Factor (mW/mK$^2$)')

TE.visualize(array_to_plot = [tmpr, prop_0pct['Thermal_conductivity'][np.newaxis, :]], linestyle = '-', marker = 'o', pic_name = 'thermal_conductivity_0pct.png', x_label = 'Temperature (K)', y_label = 'Thermal Conduct.(W/mK)')

TE.visualize(array_to_plot = [tmpr, prop_0pct['Lorenz_number'][np.newaxis, :] * 1e8], linestyle = '-', marker = 'o', pic_name = 'lorenz_number_0pct.png', x_label = 'Temperature (K)', y_label = 'Lorenz Number (x$10^{-8}$[V/K]$^2$)')
