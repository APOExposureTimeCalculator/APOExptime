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
        The airmass of the target. Max airmass handled is 3 Defaults to 1.

    Attributes
    ----------

    lunar_phase : float
        The phase of the moon. 0 is a new moon and 1 is a full moon.

    airmass : float
        The airmass. This parameter is related to the altitude of the target.

    seeing : float
        The seeing parameter. For large aperature telescopes, this is typically 1 arcsecond. Defaults to 1 arcsecond

    sky_transmission : Interpolated Object
        The transmission of the sky interpolated from 3000A - 30000A.

    sky_emission : Interpolated Object
        The emission of the sky interpolated from 3000A - 30000A

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

    magsystem : str The magnitude system used in the `mag` parameter. The 3 options available are 'VEGAMAG', 'stmag',
    and 'abnu'. String IS case sensitive

    filt_range : tuple The wavelength range of the filter you wish to observe in. Default is wavelength range
    corresponding to the Johnson V band

    sed : arr, optional
        Optional ability to enter your own spectral energy distribution of the target object. Defaults to None.

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
        1 if the object is an imager. 0 if it is a spectrograph.

    gain : float
        The gain of the instrument.
    """

    def __init__(self, target, sky, instrument):

        # TODO Redo value loading from objects to be more streamline and the same for spectrograph and imager

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

        self.Npix = self.Npix(instrument)
        self.counts(self.source, instrument)
        self.skycounts(self.skySED, instrument)

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
        # TODO Remove spectrograph insturment loading from here. Do after instrument values are outputed in
        #  consistant manner

        # Determine whether the instrument is an imager or a spectrograph.
        if self.isImager == 1:
            Npix = np.pi * ((self.seeing / 2) ** 2) / (instrument.scale ** 2)
        else:
            Npix = []
            self.disp_name = instrument.disp_name
            self.chan_name = instrument.chan_name
            spec_width = instrument.spec_width
            self.spec_range = instrument.spec_range
            self.disp_efficiency = instrument.disp_efficiency
            if spec_width[0] == 'slit':
                spec_width = self.seeing / instrument.scale[0]

            for i, row in enumerate(instrument.disp_name):
                Npix.append(spec_width * instrument.Npix_lam[i](
                    range(instrument.spec_range[i][0], instrument.spec_range[i][1])))
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
                    interpolationrange = range(int(integrate_range[0]), int(integrate_range[1]))
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
                s_integrade = s_integradeInterpolate(
                    [sky, self.detector_qe[i], self.skyTransmission, self.disp_efficiency[i]],
                    interpolationrange)

                sky_prime_dlam.append([(self.telescope_area * (self.seeing * 1.1) * s_integrade[1]),
                                       s_integrade[0]])
            self.sky_prime_dlam = sky_prime_dlam

    def counts(self, source, instrument):
        """The counts received from the source.

        Parameters
        -----------
        source : Interpolated Object
            The wavelength flux received from the source.

        instrument : object
            The ``APOinputclasses.Instrument`` class.
        """
        # TODO make spectrograph and instrument count calcuating the same function, only changing the output to be an
        #  extra integrate step for mager

        att = dir(instrument)
        h = 6.626 * 10 ** (-27)  # ergs*s
        c = 2.9979 * 10 ** (18)  # A/s
        if self.isImager == 1:
            for row in att:
                if row.find('filter') > 0:
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument, row.replace('filter', 'range'))
                    interpolationrange = range(int(integrate_range[0]), int(integrate_range[1]))

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

                s_integrade = s_integradeInterpolate(
                    [source, self.detector_qe[i], self.skyTransmission, self.disp_efficiency[i]],
                    interpolationrange)

                s_prime_dlam.append([(self.telescope_area * (1 / (h * c)) * s_integrade[1] * interpolationrange),
                                     s_integrade[0]])
            self.s_prime_dlam = s_prime_dlam

    def SNfromTime(self, exptime):
        """Computes the signal to noise ratio for a given exposure time.

        Parameters
        ----------
        exptime : float
            The exposure time for which you wish to compute the signal to noise ratio.

        Returns
        --------
        returnList : bytearray
            Array containing signal to noise and filter/dispersion names for each filter/dispersion of instrument

        """
        # TODO make process for imager and spectrograph the same
        att = dir(self)
        returnList = [1]
        self.exptime = exptime
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    SN = (sourceCount * exptime) / np.sqrt(sourceCount * exptime + skyCount * exptime
                                                           + self.Npix * (self.rdnoise) ** 2)
                    SNname = row.replace('sourcecountrate', 'SN')
                    returnList.append([SN, SNname])

                    setattr(Observation, SNname, SN)

        else:
            for i, row in enumerate(self.disp_name):
                SN_d_lam = (self.s_prime_dlam[i][0] * exptime) / np.sqrt(
                    self.s_prime_dlam[i][0] * exptime + self.sky_prime_dlam[i][0] * exptime
                    + (self.Npix[i] * (self.rdnoise[i]) ** 2))
                returnList.append([np.array(self.s_prime_dlam[i][1]), SN_d_lam, self.chan_name[i], row])

        self.SN = returnList
        return returnList

    def TimefromSN(self, SN):
        """Computes the exposure time need to achieve a desired signal to noise ratio.

        Parameters
        ----------
        SN : float
            The desired signal to noise ratio.

        Returns
        --------
        returnList : bytearray
            Array containing Exposure time and filter/dispersion names for each filter/dispersion of instrument

        """

        att = dir(self)
        returnList = [0]
        if self.isImager == 1:
            for row in att:
                if row.find('sourcecountrate') > 0:
                    sourceCount = getattr(self, row)
                    skyCount = getattr(self, row.replace('source', 'sky'))

                    t = (1. / (2. * sourceCount ** 2)) * (SN ** 2 * (sourceCount + skyCount) + np.sqrt(SN ** 4 * (
                            sourceCount + skyCount) ** 2 + 4. * self.Npix * (sourceCount * SN * (
                        self.rdnoise)) ** 2))

                    Tname = row.replace('sourcecountrate', 'time')
                    returnList.append([t, Tname])

                    setattr(Observation, Tname, t)
        else:
            for i, row in enumerate(self.disp_name):
                t_d_lam = (1. / (2. * self.s_prime_dlam[i][0] ** 2)) * (
                        SN ** 2 * (self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) + np.sqrt(SN ** 4 * (
                        self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) ** 2 + 4. * self.Npix[i] * (
                            self.s_prime_dlam[i][0] * SN * (self.rdnoise[i])) ** 2))
                returnList.append([np.array(self.s_prime_dlam[i][1]), t_d_lam, self.chan_name[i], row])

        self.Time = returnList
        return returnList


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
        disp_efficiency = []

        if para['isImager'][0] > 0:
            for row in para['Filters']:
                filt = ascii.read('../data/' + Telescope_name + '/' + Instr_name + "/" + row)
                data_name = row.split('.dat')[0]
                filt_wavelength = filt["col1"]
                filt_throughput = filt["col2"]
                filt_interpolated = interpolate.InterpolatedUnivariateSpline(
                    filt_wavelength, filt_throughput)
                setattr(Instrument, data_name, filt_interpolated)
                # TODO change output so it's an array and not an attribute for every filter
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
                if str(chan) != '--':
                    effic = ascii.read(
                        '../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_' + chan + '_qe.data')
                    qefficiency_wavelength = effic["col1"] * 10  # multiplied by 10 to turn to angstroms
                    qefficiency_percent = effic["col2"] / 100  # divided by 100 to turn into decimal

                    for disp in para[chan + '_dispersions']:
                        if str(disp) != '--':
                            disp_name.append(disp)
                            chan_name.append(chan)
                            spec_width.append(para['Width'][i])
                            gain.append(para['gain'][i])
                            scale.append(para['plate_scale[arcsec/pix]'][i])
                            readout_noise.append(para['readoutnoise[electrons]'][i])

                            efficiency.append(interpolate.InterpolatedUnivariateSpline(
                                qefficiency_wavelength, qefficiency_percent))

                            disp_effic = ascii.read(
                                '../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_' + disp + '_effic.data')
                            disp_efficiency.append(interpolate.InterpolatedUnivariateSpline(disp_effic['col1'],
                                                                                            (disp_effic['col2'] / 100)))

                            dispersion = ascii.read(
                                '../data/' + Telescope_name + '/' + Instr_name + "/" + Instr_name + '_' + disp + '_disp.data')
                            dispersion_interplate = interpolate.InterpolatedUnivariateSpline(dispersion['col2'],
                                                                                             (dispersion['col1'] ** (
                                                                                                 -1)))
                            spec_range.append([dispersion['col2'].min(), dispersion['col2'].max()])

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
            self.disp_efficiency = disp_efficiency
    # TODO clean instrument attributes so that imager and spectrograph use same outputs. Should make it more concise
    #  and easier to follow. Possiably make more outputs arrays


def s_integradeInterpolate(functions, interpolation_range):
    """The integrand of the count equation.

    This objects takes in all of the interpolated objects that goes into the count equation (sky transmission,
    sky emission, telescope throughput, instrument effiency, etc.) and multiplies them together. It then Outputs an
    array with the values of the product vs wavelength

    Parameters
    ----------
    functions : arr-like
        List of interpolated objects to go into the count equation.

    interpolation_range : tuple
        The range that wish you to interpolate over.

    Returns ------- interpolation_range, x : tuple Tuple where the first element is the interpolation range and the
        second element is the product array from the multiplication.


    """
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range) * x

    return [interpolation_range, x]
