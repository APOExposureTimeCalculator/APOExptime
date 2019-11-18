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
            emission[:, 0] * 10, (emission[:, 1] * 1E-8))
        return sky_emission


class Target:
    """Object to contain attributes of a target given an apparent magnitude, filter and mag system for given mag,
    and effective temprature, """

    def __init__(self, mag, magsystem, filtRange, sed=None, temp=5778):
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
        self.gain = instrument.gain

        self.counts(self.source, instrument)
        self.skycounts(self.skySED, instrument)
        self.Npix = self.loadinst(instrument)

    def loadinst(self, instrument):
        """

        :param instrument:
        :return:
        """
        if self.isImager == 1:
            Npix = np.pi * ((self.seeing / 2) ** 2) / (instrument.scale ** 2)
        else:

            self.disp_name = instrument.disp_name
            self.chan_name = instrument.chan_name
            spec_width = instrument.spec_width
            self.spec_range = instrument.spec_range
            if spec_width == 'slit':
                spec_width = self.seeing/instrument.scale

            Npix = spec_width*instrument.Npix_lam(range(instrument.range[0], instrument.range[1]))
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
                    s_prime_dlam = [self.telescope_area * (np.pi * (self.seeing / 2) ** 2) * s_integrade[1],
                                    s_integrade[0]]
                    s_prime = np.trapz(s_prime_dlam[0], s_prime_dlam[1])
                    count_name = row.replace('_filter', '') + '_skycountrate'
                    setattr(Observation, count_name, s_prime)
        else:
            sky_prime_dlam = []
            for i, row in enumerate(instrument.disp_name):
                interpolationrange = range(instrument.spec_range[i][0], instrument.spec_range[i][1])
                s_integrade = s_integradeInterpolate([sky, self.detector_qe[i], self.skyTransmission],
                                                     interpolationrange)

                sky_prime_dlam.append([(self.telescope_area * (np.pi * (self.seeing / 2) ** 2) * s_integrade[1]),
                                       s_integrade[0]])
                # note, A is diffrent for spectrograph and needs to be changed
            self.sky_prime_dlam = sky_prime_dlam

    def counts(self, source, instrument):
        att = dir(instrument)
        h = 6.626 * 10 ** (-27)  # ergs*s
        c = 2.9979 * 10 ** (18)  # A/s
        if self.isImager == 1:
            for row in att:
                if row.find('filter') > 0:
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument, row.replace('filter', 'range'))
                    interpolationrange = range(integrate_range[0], integrate_range[1])

                    s_integrade = s_integradeInterpolate(
                        [source, self.detector_qe, self.skyTransmission, filter_profile],
                        interpolationrange)

                    s_prime_dlam = [(self.telescope_area * (1 / (h * c)) * s_integrade[1] * interpolationrange),
                                    s_integrade[0]]
                    s_prime = np.trapz(s_prime_dlam[0], s_prime_dlam[1])
                    count_name = row.replace('_filter', '') + '_sourcecountrate'
                    setattr(Observation, count_name, s_prime)
        else:
            s_prime_dlam = []
            for i, row in enumerate(instrument.disp_name):
                interpolationrange = range(instrument.spec_range[i][0], instrument.spec_range[i][1])

                s_integrade = s_integradeInterpolate([source, self.detector_qe[i], self.skyTransmission],
                                                     interpolationrange)

                s_prime_dlam.append([(self.telescope_area * (1 / (h * c)) * s_integrade[1] * interpolationrange),
                                     s_integrade[0]])
            self.s_prime_dlam = s_prime_dlam
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
                                                           + self.Npix * (self.gain * self.rdnoise) ** 2)
                    SNname = row.replace('sourcecountrate', 'SN')
                    returnList.append([SN, SNname])

                    setattr(Observation, SNname, SN)

        else:
            for i, row in enumerate(self.disp_name):
                SN_d_lam = (self.s_prime_dlam[i][0] * exptime) / np.sqrt(
                    self.s_prime_dlam[i][0] * exptime + self.sky_prime_dlam[i][0] * exptime
                    + (self.Npix[i] * (self.gain[i] * self.rdnoise[i]) ** 2))
                returnList.append([np.array(self.s_prime_dlam[i][1]), SN_d_lam, self.chan_name[i], row])

        self.SN = returnList
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
                            sourceCount + skyCount) ** 2 + 4. * self.Npix * (sourceCount * SN * (
                        self.gain*self.rdnoise)) ** 2))

                    Tname = row.replace('sourcecountrate', 'time')
                    returnList.append([t, Tname])

                    setattr(Observation, Tname, t)
        else:
            for i, row in enumerate(self.disp_name):
                t_d_lam = (1. / (2. * self.s_prime_dlam[i][0] ** 2)) * (
                        SN ** 2 * (self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) + np.sqrt(SN ** 4 * (
                        self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) ** 2 + 4. * self.Npix[i] * (self.s_prime_dlam[i][
                                                                                                    0] * SN * (
                                                                                                        self.gain[i] * self.rdnoise[i])) ** 2))
                returnList.append([np.array(self.s_prime_dlam[1]), t_d_lam, self.chan_name[i], row])

        self.Time = returnList
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

    def __init__(self, Instr_name, Telescope_name='apo3_5m'):

        para = ascii.read('../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_param.data')
        self.isImager = para['isImager'][0]
        chan_name = []
        disp_name = []
        gain = []
        efficiency = []
        scale = []
        readout_noise = []
        spec_range = []
        spec_width = []
        Npix_lam = []

        if para['isImager'][0] > 0:
            for row in para['Filters']:
                filt = ascii.read('../data/' + Telescope_name + '/' + Instr_name + "/" + row)
                data_name = row.split('.dat')[0]
                filt_wavelength = filt["col1"]
                filt_throughput = filt["col2"]
                filt_interpolated = interpolate.InterpolatedUnivariateSpline(
                    filt_wavelength, filt_throughput)
                setattr(Instrument, data_name, filt_interpolated)

                filt_range = [filt["col1"][0], filt["col1"][len(filt["col1"]) - 1]]
                range_name = row.split('filter.dat')[0] + 'range'
                setattr(Instrument, range_name, filt_range)
            efficiency = ascii.read('../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_qe.data')

            qefficiency_wavelength = efficiency["col1"] * 10  # multiplied by 10 to turn to angstroms
            qefficiency_percent = efficiency["col2"] / 100  # divided by 100 to turn into decimal

            self.efficiency = interpolate.InterpolatedUnivariateSpline(
                qefficiency_wavelength, qefficiency_percent)
            self.readout_noise = para['readoutnoise[electrons]'][0]
            self.filter_num = para['FilterNum'][0]
            self.gain = para['gain'][0]
            self.scale = para['plate_scale[arcsec/pix]'][0]

        if para['isImager'][0] == 0:

            for i, chan in enumerate(para['channels']):
                effic = ascii.read(
                    '../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_' + chan + '_qe.data')
                qefficiency_wavelength = effic["col1"] * 10  # multiplied by 10 to turn to angstroms
                qefficiency_percent = effic["col2"] / 100  # divided by 100 to turn into decimal

                for i2, disp in enumerate(para[chan + '_dispersions']):
                    disp_name.append(disp)
                    chan_name.append(chan)
                    spec_width.append(para['Width'][i])
                    gain.append(para['gain'][i])
                    scale.append(para['plate_scale[arcsec/pix]'][i])
                    readout_noise.append(para['readoutnoise[electrons]'][i])
                    efficiency.append(interpolate.InterpolatedUnivariateSpline(
                        qefficiency_wavelength, qefficiency_percent))

                    dispersion = ascii.read(
                        '../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_' + disp)
                    dispersion_interplate = interpolate.InterpolatedUnivariateSpline(dispersion['col2'],
                                                                                     (dispersion['col1'] ** (-1)))
                    spec_range.append([dispersion['col2'].min(), para['col2'].max()])

                    Npix_lam.append(dispersion_interplate)

            self.disp_name = disp_name
            self.chan_name = chan_name
            self.spec_width = spec_width
            self.gain = gain
            self.scale = scale
            self.readout_noise = readout_noise
            self.efficiency = efficiency
            self.spec_range = spec_range
            self.Npix_lam = Npix_lam


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
