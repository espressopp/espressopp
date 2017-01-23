#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


'''
*****************************
units - convert to real units
*****************************

Espresso++ returns temperature, energy, pressure, box length etc. in dimensionless units. Usually user should take care about real length, energy, mass and charge units.  This python class is a helper in order to simplify the conversion which is based on basic units.  However, user always should use it carefully for complicated systems.

Currently it is implemented for SI units. Make sure that you are using
length in [nm]
energy in [kJ/mol]
mass in   [amu]
q in      [e]

and it will return you
pressure in     [bar]
temperature in  [K]
time in         [ps]
density in      [kg/m^3]
  
Example:



'''

import espressopp
import math

kB  = 1.3806488 * pow(10,-23) # m^2 * kg * s^-2 * K^-1
Na  = 6.0221413 * pow(10, 23) # mol^-1
amu = 1.6605389 #* pow(10,-27)

class Real_Units:
  def __init__(self, _length, _energy, _mass, _charge):
    self.length_factor = _length
    self.energy_factor = _energy
    self.mass_factor   = _mass
    self.charge_factor = _charge
    
    self.pressure_factor     = self.energy_factor / pow(self.length_factor, 3)
    self.temperature_factor  = self.energy_factor / (kB * Na) * 1000
    self.time_factor         = self.length_factor * math.sqrt( self.mass_factor / self.energy_factor)
    self.density_factor      = self.mass_factor * amu / pow(self.length_factor, 3)
  
  def length(self, dl_length):
    return dl_length * self.length_factor
  
  def energy(self, dl_energy):
    return dl_energy * self.energy_factor
  
  def mass(self, dl_mass):
    return dl_mass * self.mass_factor
  
  def charge(self, dl_charge):
    return dl_charge * self.charge_factor
  
  def pressure(self, dl_pressure):
    return dl_pressure * self.pressure_factor

  def temperature(self, dl_temperature):
    return dl_temperature * self.temperature_factor
  
  def time(self, dl_time):
    return dl_time * self.time_factor

  def density(self, dl_density):
    return dl_density * self.density_factor

  # the other way arround
  def dl_length(self, dl_length):
    return dl_length / self.length_factor
  
  def dl_energy(self, energy):
    return energy / self.energy_factor
  
  def dl_mass(self, mass):
    return mass / self.mass_factor
  
  def dl_charge(self, charge):
    return charge / self.charge_factor
  
  def dl_pressure(self, pressure):
    return pressure / self.pressure_factor

  def dl_temperature(self, temperature):
    return temperature / self.temperature_factor
  
  def dl_time(self, time):
    return time / self.time_factor

  def dl_density(self, density):
    return density / self.density_factor
