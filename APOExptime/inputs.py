p# -*- coding: utf-8 -*-

from synphot.models import BlackBody1D
from synphot import units
from os.path import dirname, abspath
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum
from astropy.io import ascii
from scipy import interpolate


class Instrument:
    """Instantiates an object of the Instrument class. This object will contain all the attributes of the selected
    instrument

    Args:
        Instr_name(str): Name of Instrument to be used.

    Attributes:
        efficiency(object): UnivariateInterpolatedSpline of instrument efficiency.
        readout_noise(float): Value of instrument readout noise
        filter_num(int): Number of filters for instrument
        gain(float): Gain of instrument
        """

    def __init__(self, Instr_name):
        """

        :param Instr_name:
        """
        para = ascii.read('../data/apo3_5m/' + Instr_name + "/" + Instr_name + '_param.data')

        efficiency = ascii.read('../data/apo3_5m/' + Instr_name + "/" + Instr_name + '_qe.data')

        for row in para['Filters']:
            filt = ascii.read('../data/apo3_5m/' + Instr_name + "/" + row)
            data_name = row.split('.dat')[0]
            filt_wavelength = filt["col1"]
            filt_throughput = filt["col2"]
            filt_interpolated = interpolate.InterpolatedUnivariateSpline(
                filt_wavelength, filt_throughput)
            setattr(Instrument, data_name, filt_interpolated)

            filt_range = [filt["col1"][0], filt["col1"][len(filt["col1"]) - 1]]
            range_name = row.split('filter.dat')[0] + 'range'
            setattr(Instrument, range_name, filt_range)

        qefficiency_wavelength = efficiency["col1"] * 10  # multiplied by 10 to turn to angstroms
        qefficiency_percent = efficiency["col2"] / 100  # divided by 100 to turn into decimal

        efficiency_interpolated = interpolate.InterpolatedUnivariateSpline(
            qefficiency_wavelength, qefficiency_percent)

        if para['isImager'][0] == 0:
            self.dispersion = para['Dispersion']

        self.efficiency = efficiency_interpolated
        self.readout_noise = para['readoutnoise[electrons]'][0]
        self.filter_num = para['FilterNum'][0]
        self.gain = para['gain'][0]
        self.scale = para['plate_scale[arcsec/pix]'][0]

class Observation:
    """Creates object for an observation given a certain telescope, instrument, sky conditions, and target"""

    def __init__(self, target, sky, instrument, telescope=None):
        """

        :param target:
        :param sky:
        :param instrument:
        :param telescope:
        """
        # telescope_transm = telescope.transmission
        self.telescope_area = (175 ** 2) * np.pi
        self.source = target.F_lambda
        self.skySED = sky.sky_emission
        self.skyTransmission = sky.sky_transmission
        self.seeing = sky.seeing
        self.rdnoise = instrument.readout_noise

        self.counts(self.source, instrument)
        self.skycounts(self.skySED, instrument)
        self.Npix = self.Npix(instrument)

    def Npix(self, instrument):
        """

        :param instrument:
        :return:
        """
        Npix = np.pi * ((self.seeing / 2) ** 2) / (instrument.scale ** 2)
        return Npix

    def skycounts(self, sky, instrument):
        att = dir(instrument)
        self.detector_qe = instrument.efficiency
        for row in att:
            if row.find('filter') > 0:
                filter_profile = getattr(instrument, row)
                integrate_range = getattr(instrument, row.replace('filter', 'range'))
                interpolationrange = range(integrate_range[0], integrate_range[1])
                s_integrade = sky_integradeInterpolate([sky, self.detector_qe, self.skyTransmission, filter_profile],
                                                       interpolationrange)
                s_prime = self.telescope_area * ((self.seeing / 2) ** 2) * s_integrade.integral(integrate_range[0],
                                                                                                integrate_range[1])
                count_name = row.replace('_filter', '') + '_skycountrate'
                setattr(Observation, count_name, s_prime)

    def counts(self, source, instrument):
        att = dir(instrument)
        self.detector_qe = instrument.efficiency
        for row in att:
            if row.find('filter') > 0:
                filter_profile = getattr(instrument, row)
                integrate_range = getattr(instrument, row.replace('filter', 'range'))
                interpolationrange = range(integrate_range[0], integrate_range[1])
                h = 6.626 * 10 ** (-27)  # ergs*s
                c = 2.9979 * 10 ** (18)  # A/s
                s_integrade = s_integradeInterpolate([source, self.detector_qe, self.skyTransmission, filter_profile],
                                                     interpolationrange)

                s_prime = self.telescope_area * (1 / (h * c)) * s_integrade.integral(integrate_range[0],
                                                                                     integrate_range[1])
                count_name = row.replace('_filter', '') + '_sourcecountrate'
                setattr(Observation, count_name, s_prime)

    def SNfromTime(self, exptime):
        """

        :param exptime:
        :return:
        """
        att = dir(self)
        returnList = []
        for row in att:
            if row.find('sourcecountrate') > 0:
                sourceCount = getattr(self, row)
                skyCount = getattr(self, row.replace('source', 'sky'))

                SN = (sourceCount * exptime) / np.sqrt(sourceCount * exptime + skyCount * exptime
                                                       + self.Npix * self.rdnoise ** 2)
                SNname = row.replace('sourcecountrate', 'SN')
                returnList.append([SN, SNname])

                setattr(Observation, SNname, SN)
        return returnList
        # PLOT SHIT HERE

    def TimefromSN(self, SN):
        att = dir(self)
        returnList = []
        for row in att:
            if row.find('sourcecountrate') > 0:
                sourceCount = getattr(self, row)
                skyCount = getattr(self, row.replace('source', 'sky'))

                t = (1. / (2. * sourceCount ** 2)) * (SN ** 2 * (sourceCount + skyCount) + np.sqrt(SN ** 4 * (
                        sourceCount + skyCount) ** 2 + 4. * self.Npix * (sourceCount * SN * self.rdnoise) ** 2))

                Tname = row.replace('sourcecountrate', 'time')
                returnList.append([t, Tname])

                setattr(Observation, Tname, t)
        return returnList
    #         # PLOT SHIT HERE


class Sky:
    """Object representing the sky flux and transmission.

    :param lunar_phase: Defines the lunar phase (0 is new and 1 is full). Defaults to 0.
    :type lunar_phase: float, optional
    :param seeing: Defines the seeing (in arcseconds). Defaults to 1.
    :type seeing: float, optional
    :param airmass: Defines the airmass. Defaults to 1.
    :type airmass: float, optional

    """

    def __init__(self,
                 lunar_phase=0,
                 seeing=1,
                 airmass=1
                 ):
        """Constructor method for the Sky class.

        """
        self.lunar_phase = lunar_phase
        self.airmass = airmass
        self.seeing = seeing
        self.sky_transmission = self.transmission()
        self.sky_emission = self.emission()

    def transmission(self):
        """Method to obtain the transmission of the sky.

        :return: The transmission of the sky interpolated.
        :rtype: Interpolated object
        """

        #Get the appropriate data file based on the airmass.
        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif 1.75 > self.airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif 1.75 <= self.airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        path_to_dir = dirname(abspath(__file__)) + '/data/Sky/'

        transmission = np.loadtxt(path_to_dir + trans_file)
        sky_transmission = interpolate.InterpolatedUnivariateSpline(
            transmission[:, 0] * 10, transmission[:, 1])

        return sky_transmission

    def emission(self):
        """Method to obtain the emission from the sky.

        :return: The sky emission of the sky interpolated.
        :rtype: Interpolated object
        """

        #Get the appropriate transmission file. Based on the lunar phase.
        if self.lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif 0.25 <= self.lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif self.lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        path_to_dir = dirname(abspath(__file__)) + '/data/Sky/'

        emission = np.loadtxt(path_to_dir + emission_file)

        sky_emission = interpolate.InterpolatedUnivariateSpline(
            emission[:, 0] * 10, (emission[:, 1] * 0.0001 * 10E-4))

        return sky_emission


class Target:
    """Object that represents the target which you wish to observe.

    :param magnitude: The magnitude of the object which you are observing.
    :type magnitude: float
    :param mag_system: The magnitude system used to measure the above magnitude.
    :type mag_system: str
    :param filter_range: The range (in Angstroms) which you wish to observe in.
    :type filter_range: tuple
    :param sed: The spectral energy distribution of the target star. Defaults to None.
    :type sed: Interpolated object, optional
    :param temp: The temperature of the target star (in Kelvin). Defaults to 5778.
    :type temp: float, optional

    """

    def __init__(self,
                 magnitude,
                 mag_system,
                 filter_range,
                 sed=None,
                 temp=5778
                 ):
        """Constructor for the target class.

        """

        #Define the magnitude system (as inputted by the user).
        if mag_system == 'VEGAMAG':
            sys = units.VEGAMAG
        elif mag_system == 'stmag':
            sys = u.STmag
        elif mag_system == 'abnu':
            sys = u.ABmag

        vega = SourceSpectrum.from_vega()

        self.magnitude = magnitude
        self.mag_system = mag_system
        self.SED = sed
        self.temp = temp
        self.inputFlux = units.convert_flux(filter_range, magnitude * sys, units.FLAM, vegaspec=vega)
        self.range = filter_range
        self.F_lambda = self.starF_lambda()

    def __str__(self):
        """Print magnitude and magnitude system of the star.

        :return: The magnitude and magnitude system of the star.
        :rtype: str

        """

        return "You are observing a {} magnitude star using the {} system."\
            .format(self.magnitude,self.mag_system)

    def f_lambda(self):
        """Method that calculates f_lambda beside on a black body spectrum, scaled to Vega.

        :return: The flux per wavelength (cgs) of the target star
        :rtype: Interpolated object

        """
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)
        sp_new = sp / np.mean(sp(self.range * u.AA, flux_unit=units.FLAM) / self.inputFlux)
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)
        f_lambda = interpolate.InterpolatedUnivariateSpline(range(1000, 30000), x)

        return f_lambda

def s_integradeInterpolate(functions, interpolation_range):
    """Takes multiple interpolation objects and multiplies them together then reinterpolates to output a single
    interpolation object :param functions: :param interpolation_range: :return:
    """
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range) * x

    return interpolate.InterpolatedUnivariateSpline(interpolation_range, (x * interpolation_range))


def sky_integradeInterpolate(functions, interpolation_range):
    """Takes multiple interpolation objects and multiplies them together then reinterpolates to output a single
    interpolation object :param functions: :param interpolation_range: :return:
    """
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range) * x

    return interpolate.InterpolatedUnivariateSpline(interpolation_range, x)
