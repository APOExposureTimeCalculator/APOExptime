# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:13:34 2019
@author: Alexander, Bryson, Roger, Hassan, Manny
"""
from synphot.models import BlackBody1D
from synphot import units
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum
from astropy.io import ascii
from scipy import interpolate


class Sky:
    """Object to contain attributes for sky conditions, given an airmass, seeing, and moon phase"""

    def __init__(self, lunar_phase=0, seeing=1, airmass=1):
        """

        :param lunar_phase:
        :param seeing:
        :param airmass:
        """
        self.lunar_phase = lunar_phase
        self.airmass = airmass
        self.seeing = seeing
        self.sky_transmission = self.transmission()
        self.sky_emission = self.emission()

    def transmission(self):
        """

        :return:
        """
        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif 1.75 > self.airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif 1.75 <= self.airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        transmission = np.loadtxt('../data/sky/' + trans_file)
        sky_transmission = interpolate.InterpolatedUnivariateSpline(
            transmission[:, 0] * 10, transmission[:, 1])
        return sky_transmission

    def emission(self):
        """

        :return:
        """
        if self.lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif 0.25 <= self.lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif self.lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        emission = np.loadtxt('../data/sky/' + emission_file)
        sky_emission = interpolate.InterpolatedUnivariateSpline(
            emission[:, 0] * 10, (emission[:, 1] * 0.0001 * 10E-4))
        return sky_emission


class Target:
    """Object to contain attributes of a target given an apparent magnitude, filter and mag system for given mag,
    and effective temprature, """

    def __init__(self, mag, magsystem, filtRange, sed=None, temp=5778, location=None):
        """

        :param mag:
        :param magsystem:
        :param filtRange:
        :param sed:
        :param temp:
        :param location:
        """
        if magsystem == 'VEGAMAG':
            sys = units.VEGAMAG
        elif magsystem == 'stmag':
            sys = u.STmag
        elif magsystem == 'abnu':
            sys = u.ABmag

        vega = SourceSpectrum.from_vega()

        self.mag = mag
        self.SED = sed
        self.temp = temp
        self.inputFlux = units.convert_flux(filtRange, mag * sys, units.FLAM, vegaspec=vega)
        self.range = filtRange
        self.F_lambda = self.starF_lambda()

    def starF_lambda(self):
        """

        """
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)
        sp_new = sp / np.mean(sp(self.range * u.AA, flux_unit=units.FLAM) / self.inputFlux)
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)
        F_lambda = interpolate.InterpolatedUnivariateSpline(range(1000, 30000), x)
        return F_lambda


class Observation:
    """Creates object for an observation given an  certain telescope, instrument, sky conditions, and target"""

    def __init__(self, target, sky, instrument, telescope=None):
        """

        :param target:
        :param sky:
        :param instrument:
        :param telescope:
        """
        # telescope_transm = telescope.transmission
        self.detector_qe = instrument.efficiency
        self.telescope_area = (175 ** 2) * np.pi
        self.source = target.F_lambda
        self.skySED = sky.sky_emission
        self.skyTransmission = sky.sky_transmission
        self.seeing = sky.seeing
        self.rdnoise = instrument.readout_noise
        self.isImager = instrument.isImager

        self.counts(self.source, instrument)
        self.skycounts(self.skySED, instrument)
        self.Npix = self.Npix(instrument)

    def Npix(self, instrument):
        """

        :param instrument:
        :return:
        """
        if self.isImager == 1:
            Npix = np.pi * ((self.seeing / 2) ** 2) / (instrument.scale ** 2)
        else:
            Npix = instrument.Npix_lam(range(instrument.range[0], instrument.range[1]))
        return Npix

    def skycounts(self, sky, instrument):
        att = dir(instrument)
        if self.isImager == 1:
            for row in att:
                if row.find('filter') > 0:
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument, row.replace('filter', 'range'))
                    interpolationrange = range(integrate_range[0], integrate_range[1])
                    s_integrade = s_integradeInterpolate([sky, self.detector_qe, self.skyTransmission, filter_profile],
                                                         interpolationrange)
                    s_prime_dlam = [self.telescope_area * ((self.seeing / 2) ** 2) * s_integrade[1], s_integrade[0]]
                    s_prime = np.trapz(s_prime_dlam[0], s_prime_dlam[1])
                    count_name = row.replace('_filter', '') + '_skycountrate'
                    setattr(Observation, count_name, s_prime)
        else:
            interpolationrange = range(instrument.range[0], instrument.range[1])
            h = 6.626 * 10 ** (-27)  # ergs*s
            c = 2.9979 * 10 ** (18)  # A/s
            s_integrade = s_integradeInterpolate([sky, self.detector_qe, self.skyTransmission],
                                                 interpolationrange)

            self.sky_prime_dlam = [(self.telescope_area * ((self.seeing / 2) ** 2) * .8 * s_integrade[1]),
                                   s_integrade[0]]

    def counts(self, source, instrument):
        att = dir(instrument)
        if self.isImager == 1:
            for row in att:
                if row.find('filter') > 0:
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument, row.replace('filter', 'range'))
                    interpolationrange = range(integrate_range[0], integrate_range[1])
                    h = 6.626 * 10 ** (-27)  # ergs*s
                    c = 2.9979 * 10 ** (18)  # A/s
                    s_integrade = s_integradeInterpolate(
                        [source, self.detector_qe, self.skyTransmission, filter_profile],
                        interpolationrange)

                    s_prime_dlam = [(self.telescope_area * (1 / (h * c)) * s_integrade[1] * interpolationrange),
                                    s_integrade[0]]
                    s_prime = np.trapz(s_prime_dlam[0], s_prime_dlam[1])
                    count_name = row.replace('_filter', '') + '_sourcecountrate'
                    setattr(Observation, count_name, s_prime)
        else:
            interpolationrange = range(instrument.range[0], instrument.range[1])
            h = 6.626 * 10 ** (-27)  # ergs*s
            c = 2.9979 * 10 ** (18)  # A/s
            s_integrade = s_integradeInterpolate([source, self.detector_qe, self.skyTransmission],
                                                 interpolationrange)

            self.s_prime_dlam = [(self.telescope_area * (1 / (h * c)) * .8 * s_integrade[1] * interpolationrange),
                                 s_integrade[0]]

    def SNfromTime(self, exptime):
        """

        :param exptime:
        :return:
        """
        att = dir(self)
        returnList = []
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    SN = (sourceCount * exptime) / np.sqrt(sourceCount * exptime + skyCount * exptime
                                                           + self.Npix * self.rdnoise ** 2)
                    SNname = row.replace('sourcecountrate', 'SN')
                    returnList.append([SN, SNname])

                    setattr(Observation, SNname, SN)

        else:
            SN_d_lam = (self.s_prime_dlam[0] * exptime) / np.sqrt(
                self.s_prime_dlam[0] * exptime + self.sky_prime_dlam[0] * exptime
                + (self.Npix * self.rdnoise ** 2))
            returnList = [np.array(self.s_prime_dlam[1]), SN_d_lam]

        return returnList
        # PLOT SHIT HERE

    def TimefromSN(self, SN):
        att = dir(self)
        returnList = []
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    t = (1. / (2. * sourceCount ** 2)) * (SN ** 2 * (sourceCount + skyCount) + np.sqrt(SN ** 4 * (
                            sourceCount + skyCount) ** 2 + 4. * self.Npix * (sourceCount * SN * self.rdnoise) ** 2))

                    Tname = row.replace('sourcecountrate', 'time')
                    returnList.append([t, Tname])

                    setattr(Observation, Tname, t)
        else:
            t_d_lam = (1. / (2. * self.s_prime_dlam[0] ** 2)) * (
                    SN ** 2 * (self.s_prime_dlam[0] + self.sky_prime_dlam[0]) + np.sqrt(SN ** 4 * (
                    self.s_prime_dlam[0] + self.sky_prime_dlam[0]) ** 2 + 4. * self.Npix * (self.s_prime_dlam[
                                                                                                0] * SN * self.rdnoise) ** 2))
            returnList = [np.array(self.s_prime_dlam[1]), t_d_lam]

        return returnList
    #         # PLOT SHIT HERE


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
        if para['isImager'][0] > 0:
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
            disp = ascii.read('../data/apo3_5m/' + Instr_name + "/" + Instr_name + '_disp.data')
            self.Npix_lam = interpolate.InterpolatedUnivariateSpline(disp['col2'], (disp['col1'] ** (-1)))
            self.range = [para['rangeMin'][0], para['rangeMax'][0]]

        self.efficiency = efficiency_interpolated
        self.readout_noise = para['readoutnoise[electrons]'][0]
        self.filter_num = para['FilterNum'][0]
        self.gain = para['gain'][0]
        self.scale = para['plate_scale[arcsec/pix]'][0]
        self.isImager = para['isImager'][0]


def s_integradeInterpolate(functions, interpolation_range):
    """Takes multiple interpolation objects and multiplies them together then reinterpolates to output a single
    interpolation object
    :param functions:
    :param interpolation_range:
    :return:


    """
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range) * x

    return [interpolation_range, x]
