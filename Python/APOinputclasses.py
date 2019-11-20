# -*- coding: utf-8 -*-

from synphot.models import BlackBody1D
from synphot import units
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum
from astropy.io import ascii
from scipy import interpolate


class Sky:
    """Object that represents the sky.

    The `Sky` class is used to compute the transmission and emission of the sky. The transmission of the sky is inferred
    from the airmass, and the emission of the sky is based on the lunar phase.

    Parameters
    ----------
    lunar_phase : float, optional
        Floating point that represents the lunar phase where 0 means new moon and 1 is a full moon. Defaults to 0.

    seeing : float, optional
        The seeing. This parameter is used to calculate the background area in the S/N ratio equation.  Defaults to 1

    airmass : float, optional
        The airmass of the target. Large observatories typically have an airmass ~ 1. Defaults to 1.

    Attributes
    ----------

    lunar_phase : float
        The phase of the moon. 0 is a new moon and 1 is a full moon.

    airmass : float
        The airmass. This parameter is related to the altitude of the target.

    seeing : float
        The seeing parameter. For large aperature telescopes, this is typically 1 arcsecond.

    sky_transmission : Interpolated Object
        The transmission of the sky interpolated over a predetermined wavelength range.

    sky_emission : Interpolated Object
        The emission of the sky interpolated over a predetermined wavelength range

    """

    def __init__(self, lunar_phase=0, seeing=1, airmass=1):

        self.lunar_phase = lunar_phase
        self.airmass = airmass
        self.seeing = seeing

        self.sky_transmission = self.transmission()
        self.sky_emission = self.emission()

    def transmission(self):
        """Determine the transmission of the sky.

        The transmission of the sky is determined by the seeing. The package includes data files which read the
        appropriate transmission file based on the airmass.

        Returns
        -------
        sky_transmission : Interpolated Object
            The transmission of the sky interpolated over a given wavelength range specified in the data files.
        """

        # Find the appropriate airmass file.
        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif 1.75 > self.airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif 1.75 <= self.airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        # Load the data file
        transmission = np.loadtxt('../data/Sky/' + trans_file)

        # Interpolate the transmission
        sky_transmission = interpolate.InterpolatedUnivariateSpline(
            transmission[:, 0] * 10, transmission[:, 1])

        # Return the interpolated transmission.
        return sky_transmission

    def emission(self):
        """Determines the emission of the sky.

        The emission of the sky is primarily based on the lunar phase. This method computes the emission (photon flux)
        of the sky per wavelength based on the ``lunar_phase`` parameter.

        Returns
        -------
        sky_emission : Interpolated Object
            The emission of the sky interpolated over a given wavelength range specified in the data files.
        """

        # Find the appropriate date files.
        if self.lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif 0.25 <= self.lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif self.lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        # Load the data files
        emission = np.loadtxt('../data/sky/' + emission_file)

        # Interpolate
        sky_emission = interpolate.InterpolatedUnivariateSpline(
            emission[:, 0] * 10, (emission[:, 1] * 1E-8))

        # Return the interpolated emission
        return sky_emission


class Target:
    """This object represents the target star which you wish to compute an exposure time for.

    This class is intended to be used for unresolved or point source objects (i.e., stars) and we do not recommend using
    it for extended objects. The class can compute the spectrum of your target by taking the temperature and scaling a
    black body spectrum to match the specified magnitude.

    Parameters
    ----------
    mag : float
        The magnitude of the target object.

    magsystem : str
        The magnitude system used in the `mag` parameter (i.e., VEGAMAG).

    filt_range : tuple
        The wavelength range of the filter you wish to observe in.

    sed : arr, optional
        The spectral energy distribution of the target object. Defaults to None.

    temp : float, optional
        The temperature (K) of the target object which is used to compute a black body spectrum. Defaults to 5778.

    Attributes
    ----------
    mag : float
        The magnitude of the target object.

    magsystem : str
        The magnitude system used in the `mag` parameter (i.e., VEGAMAG).

    filt_range : tuple
        The wavelength range of the filter you wish to observe in.

    sed : arr, optional
        The spectral energy distribution of the target object.

    temp : float, optional
        The temperature (K) of the target object which is used to compute a black body spectrum.

    """

    def __init__(self, mag, magsystem, filt_range, sed=None, temp=5778):

        # Define the magnitude system.
        if magsystem.lower() == 'vegamag':
            sys = units.VEGAMAG
        elif magsystem.lower() == 'stmag':
            sys = u.STmag
        elif magsystem.lower() == 'abnu':
            sys = u.ABmag

        # Get Vega's spectrum.
        vega = SourceSpectrum.from_vega()

        # Set attributes.
        self.mag = mag
        self.SED = sed
        self.temp = temp
        self.inputFlux = units.convert_flux(filt_range, mag * sys, units.FLAM, vegaspec=vega)
        self.range = filt_range
        self.F_lambda = self.starF_lambda()

    def starF_lambda(self):
        """Compute the wavelength flux of the target object.

        This method creates a black body spectrum of temperature ``temp`` and scaled that spectrum to match the flux of
        a ``mag`` magnitude object.

        Returns
        --------
        F_lambda : Interpolated Object
            The spectrum of the star interpolated from 1000 A to 30000 A.
        """

        # Get the black body spectrum of an object at temperature "temp".
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)

        # Scale that black body to match the flux of a "mag" magnitude star.
        sp_new = sp / np.mean(sp(self.range * u.AA, flux_unit=units.FLAM) / self.inputFlux)
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)

        # Interpolate the flux.
        F_lambda = interpolate.InterpolatedUnivariateSpline(range(1000, 30000), x)

        # Return the interpolated flux.
        return F_lambda


class Observation:
    """Creates object for an observation given a certain telescope, instrument, sky conditions, and target.

    This object takes in the three classes specified above to compute a variety of things such as the signal/noise, the
    count rate from the source, the count rate of the sky, etc.

    Parameters
    ----------
    target : Object
        The ``APOinputclasses.Target`` class.

    sky : Object
        The ``APOinputclasses.Sky`` class.

    instrument : Object
        The ``APOinputclasses.Instrument`` class.

    Attributes
    ----------
    detector_qe : Interpolated Object
        The quantum efficiency of the detector.

    telescope_area : float
        The light collecting area of the telescope (in cgs).

    source : Interpolated Object
        The flux of the target object interpolated.

    skySED : Interpolated Object
        The emission of the sky interpolated.

    skyTransmission : Interpolated Object
        The transmission of the sky interpolated.

    seeing : float
        The seeing.

    rdnoise : float
        The readout noise of the instrument.

    isImager : bool
        True if the object is an imager. False if it is a spectrograph.

    gain : float
        The gain of the instrument.
    """

    def __init__(self, target, sky, instrument, telescope=None):

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
        self.Npix = self.Npix(instrument)

    def Npix(self, instrument):
        """The number of pixels covered by the source and sky.

        The number of pixels is used to compute the area covered by the sky on the detector as well as the amount of
        pixels that contributed to the readout noise. This method takes the seeing and the plate scale of the instrument
        to compute Npix.

        Parameters
        ----------
        instrument : object
            The ``APOinputclasses.Instrument`` class.

        Returns
        -------
        Npix : float
            The number of pixels.

        """

        # Determine whether the instrument is an imager or a spectrograph.
        if self.isImager == 1:
            Npix = np.pi * ((self.seeing / 2) ** 2) / (instrument.scale ** 2)
        else:
            Npix = instrument.Npix_lam(range(instrument.range[0], instrument.range[1]))

        # Return Npix
        return Npix

    def skycounts(self, sky, instrument):
        """Computes the amount of counts received by the sky.

        Parameters
        ----------
        sky : object
            The ``APOinputclasses.Sky`` class.

        instrument : object
            The ``APOinputclasses.Instrument`` class.

        """
        att = dir(instrument)
        if self.isImager == 1:
            for row in att:
                if row.find('filter') > 0:
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument, row.replace('filter', 'range'))
                    interpolationrange = range(integrate_range[0], integrate_range[1])
                    s_integrade = s_integradeInterpolate([sky, self.detector_qe, self.skyTransmission, filter_profile],
                                                         interpolationrange)
                    s_prime_dlam = [self.telescope_area * (np.pi*(self.seeing / 2) ** 2) * s_integrade[1], s_integrade[0]]
                    s_prime = np.trapz(s_prime_dlam[0], s_prime_dlam[1])
                    count_name = row.replace('_filter', '') + '_skycountrate'
                    setattr(Observation, count_name, s_prime)
        else:
            interpolationrange = range(instrument.range[0], instrument.range[1])
            h = 6.626 * 10 ** (-27)  # ergs*s
            c = 2.9979 * 10 ** (18)  # A/s
            s_integrade = s_integradeInterpolate([sky, self.detector_qe, self.skyTransmission],
                                                 interpolationrange)

            self.sky_prime_dlam = [(self.telescope_area * (np.pi*(self.seeing / 2) ** 2) * s_integrade[1]),
                                   s_integrade[0]]

    def counts(self, source, instrument):
        """The counts received from the source.

        Parameters
        -----------
        source : Interpolated Object
            The wavelength flux received from the source.

        instrument : object
            The ``APOinputclasses.Instrument`` class.
        """
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

            self.s_prime_dlam = [(self.telescope_area * (1 / (h * c)) * s_integrade[1] * interpolationrange),
                                         s_integrade[0]]

    def SNfromTime(self, exptime):
        """Computes the signal to noise ratio for a given exposure time.

        Parameters
        ----------
        exptime : float
            The exposure time for which you wish to compute the signal to noise ratio.

        Returns
        --------
        returnList : arr

        """

        att = dir(self)
        returnList = []
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    SN = (sourceCount * exptime) / np.sqrt(sourceCount * exptime + skyCount * exptime
                                                           + self.Npix * (self.gain*self.rdnoise) ** 2)
                    SNname = row.replace('sourcecountrate', 'SN')
                    returnList.append([SN, SNname])

                    setattr(Observation, SNname, SN)

        else:
            SN_d_lam = (self.s_prime_dlam[0] * exptime) / np.sqrt(
                self.s_prime_dlam[0] * exptime + self.sky_prime_dlam[0] * exptime
                + (self.Npix * (self.gain*self.rdnoise) ** 2))
            returnList = [np.array(self.s_prime_dlam[1]), SN_d_lam]

        return returnList
        # PLOT SHIT HERE

    def TimefromSN(self, SN):
        """Computes the exposure time need to achieve a desired signal to noise ratio.

        Parameters
        ----------
        SN : float
            The desired signal to noise ratio.

        Returns
        --------
        returnList : arr

        """

        att = dir(self)
        returnList = []
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    t = (1. / (2. * sourceCount ** 2)) * (SN ** 2 * (sourceCount + skyCount) + np.sqrt(SN ** 4 * (
                            sourceCount + skyCount) ** 2 + 4. * self.Npix * (sourceCount * SN * (self.gainself.rdnoise)) ** 2))

                    Tname = row.replace('sourcecountrate', 'time')
                    returnList.append([t, Tname])

                    setattr(Observation, Tname, t)
        else:
            t_d_lam = (1. / (2. * self.s_prime_dlam[0] ** 2)) * (
                    SN ** 2 * (self.s_prime_dlam[0] + self.sky_prime_dlam[0]) + np.sqrt(SN ** 4 * (
                    self.s_prime_dlam[0] + self.sky_prime_dlam[0]) ** 2 + 4. * self.Npix * (self.s_prime_dlam[
                                                                                                0] * SN * (self.gain*self.rdnoise)) ** 2))
            returnList = [np.array(self.s_prime_dlam[1]), t_d_lam]

        return returnList
    #         # PLOT SHIT HERE


class Instrument:
    """Object that represents the instrument used to observe.

    It is important to note that since this exposure time calculator is designed for the Astrophysical Research
    Consortium (ARC) 3.5m telescope, the list of instruments available is exclusive to this telescope. That list is::
        * ARCTIC        (Imager)
        * AGILE         (Imager)
        * ARCES         (Spectrograph)
        * DIS           (Spectrograph)
        * TRIPLESEC     (Spectrograph)
        * NICFPS        (Spectrograph)

    Parameters
    -----------
    Instr_name : (str)
        The name of the instrument used.

    Attributes
    ----------
    efficiency: Interpolated Object
        UnivariateInterpolatedSpline of the instrument efficiency.

    readout_noise : float
        Value of instrument readout noise.

    filter_num : int
        Number of filters for the instrument.

    gain : float
        Gain of the instrument.

    scale : float
        The plate scale of the instrument.
    """

    def __init__(self, Instr_name):

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
            spec_width = para['Width']
            disp = ascii.read('../data/apo3_5m/' + Instr_name + "/" + Instr_name + '_disp.data')
            self.Npix_lam = interpolate.InterpolatedUnivariateSpline(disp['col2'], (spec_width*disp['col1'] ** (-1)))
            self.range = [para['rangeMin'][0], para['rangeMax'][0]]

        self.efficiency = efficiency_interpolated
        self.readout_noise = para['readoutnoise[electrons]'][0]
        self.filter_num = para['FilterNum'][0]
        self.gain = para['gain'][0]
        self.scale = para['plate_scale[arcsec/pix]'][0]
        self.isImager = para['isImager'][0]


def s_integradeInterpolate(functions, interpolation_range):
    """The integrand of the count equation.

    This objects takes in all of the interpolated objects that goes into the count equation (sky transmission, sky
    emission, telescope throughput, instrument effiency, etc.) and multiplies them together. It then re-interpolates it
    based on the specified interpolation range.

    Parameters
    ----------
    functions : arr-like
        List of interpolated objects to go into the count equation.

    interpolation_range : tuple
        The range that wish you to interpolate over.

    Returns
    -------
    interpolation_range, x : tuple
        Tuple where the first element is the interpolation range and the second element is the re-interpolated object.


    """
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range) * x

    return [interpolation_range, x]
