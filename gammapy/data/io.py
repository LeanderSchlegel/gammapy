# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings
from astropy.io import fits
from astropy.table import Table
from gammapy.data import EventListMetaData, EventList, ObservationTable
from gammapy.utils.scripts import make_path
from gammapy.utils.metadata import CreatorMetaData

from .metadata import ObservationMetaData


class ObservationTableReader:
    """Reader class for ObservationTable
    Based on class EventListReader!!!

    Format specification: gadf v.0.2 and v.0.3

    Parameters
    ----------
    hdu : str
        Name of observation-index HDU. Deafult is OBS-INDEX.
    checksum : bool
        If True checks both DATASUM and CHECKSUM cards in the file headers. Default is False.
    """

    def __init__(self, hdu="OBS-INDEX", checksum=False):
        self.hdu = hdu
        self.checksum = checksum

    @staticmethod
    def from_hdu(obs_hdu, format):
        """Create ObservationTable from HDU with specific format."""
        table = Table.read(obs_hdu)
        meta = ObservationMetaData.from_header(table.meta, format)
        return ObservationTable(table=table, meta=meta)

    @staticmethod
    def identify_format_from_hduclass(obs_hdu):
        """Identify format for HDU header keywords."""
        hduclass = obs_hdu.header.get("HDUCLASS", "unknown")
        return hduclass.lower()

    def read(self, filename, format="gadf03"):
        """Read ObservationTable from file.

        Parameters
        ----------
        filename : `pathlib.Path`, str
            Filename
        format : {"gadf02"}, {"gadf03"} optional
            format of the ObservationTable. Default is 'gadf03'.
            If None, will try to guess from header.
        """
        filename = make_path(filename)

        with fits.open(filename) as hdulist:
            obs_hdu = hdulist[self.hdu]

            if self.checksum:
                if obs_hdu.verify_checksum() != 1:
                    warnings.warn(
                        f"Checksum verification failed for HDU {self.hdu} of {filename}.",
                        UserWarning,
                    )

            if format is None:
                format = self.identify_format_from_hduclass(obs_hdu)

            if format == "gadf02" or format == "gadf03":
                return self.from_gadf_hdu(obs_hdu)
            else:
                raise ValueError(f"Unknown format :{format}")


class EventListReader:
    """Reader class for EventList.

    Format specification: :ref:`gadf:iact-events`

    Parameters
    ----------
    hdu : str
        Name of events HDU. Default is "EVENTS".
    checksum : bool
        If True checks both DATASUM and CHECKSUM cards in the file headers. Default is False.
    """

    def __init__(self, hdu="EVENTS", checksum=False):
        self.hdu = hdu
        self.checksum = checksum

    @staticmethod
    def from_gadf_hdu(events_hdu):
        """Create EventList from gadf HDU."""
        table = Table.read(events_hdu)
        meta = EventListMetaData.from_header(table.meta)
        return EventList(table=table, meta=meta)

    @staticmethod
    def identify_format_from_hduclass(events_hdu):
        """Identify format for HDU header keywords."""
        hduclass = events_hdu.header.get("HDUCLASS", "unknown")
        return hduclass.lower()

    def read(self, filename, format="gadf"):
        """Read EventList from file.

        Parameters
        ----------
        filename : `pathlib.Path`, str
            Filename
        format : {"gadf"}, optional
            format of the EventList. Default is 'gadf'.
            If None, will try to guess from header.
        """
        filename = make_path(filename)

        with fits.open(filename) as hdulist:
            events_hdu = hdulist[self.hdu]

            if self.checksum:
                if events_hdu.verify_checksum() != 1:
                    warnings.warn(
                        f"Checksum verification failed for HDU {self.hdu} of {filename}.",
                        UserWarning,
                    )

            if format is None:
                format = self.identify_format_from_hduclass(events_hdu)

            if format == "gadf" or format == "ogip":
                return self.from_gadf_hdu(events_hdu)
            else:
                raise ValueError(f"Unknown format :{format}")


class EventListWriter:
    """Writer class for EventList."""

    def __init__(self):
        pass

    @staticmethod
    def _to_gadf_table_hdu(event_list):
        """Convert input event list to a `~astropy.io.fits.BinTableHDU` according gadf."""
        bin_table = fits.BinTableHDU(event_list.table, name="EVENTS")

        # A priori don't change creator information
        if event_list.meta.creation is None:
            event_list.meta.creation = CreatorMetaData()
        else:
            event_list.meta.creation.update_time()

        bin_table.header.update(event_list.meta.to_header())
        return bin_table

    def to_hdu(self, event_list, format="gadf"):
        """
        Convert input event list to a `~astropy.io.fits.BinTableHDU` according to format.

        Parameters
        ----------
        format : str, optional
            Output format, currently only "gadf" is supported. Default is "gadf".

        Returns
        -------
        hdu : `astropy.io.fits.BinTableHDU`
            EventList converted to FITS representation.
        """
        if format != "gadf":
            raise ValueError(f"Only the 'gadf' format supported, got {format}")

        return self._to_gadf_table_hdu(event_list)
