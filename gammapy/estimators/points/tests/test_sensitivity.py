# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
from numpy.testing import assert_allclose
from gammapy.datasets import SpectrumDataset, SpectrumDatasetOnOff, Datasets
from gammapy.estimators import (
    FluxPoints,
    SensitivityEstimator,
    ParameterSensitivityEstimator,
)
from gammapy.irf import EDispKernelMap
from gammapy.maps import MapAxis, RegionNDMap
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel


@pytest.fixture()
def spectrum_dataset():
    e_true = MapAxis.from_energy_bounds("1 TeV", "10 TeV", nbin=20, name="energy_true")
    e_reco = MapAxis.from_energy_bounds("1 TeV", "10 TeV", nbin=4)

    background = RegionNDMap.create(region="icrs;circle(0, 0, 0.1)", axes=[e_reco])
    background.data += 3600
    background.data[0] *= 1e3
    background.data[-1] *= 1e-3
    edisp = EDispKernelMap.from_diagonal_response(
        energy_axis_true=e_true, energy_axis=e_reco, geom=background.geom
    )

    exposure = RegionNDMap.create(
        region="icrs;circle(0, 0, 0.1)", axes=[e_true], unit="m2 h", data=1e6
    )

    spectrum_dataset = SpectrumDataset(
        name="test", exposure=exposure, edisp=edisp, background=background
    )
    return spectrum_dataset


def test_cta_sensitivity_estimator(spectrum_dataset):
    geom = spectrum_dataset.background.geom

    dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=spectrum_dataset,
        acceptance=RegionNDMap.from_geom(geom=geom, data=1),
        acceptance_off=RegionNDMap.from_geom(geom=geom, data=5),
    )

    sens = SensitivityEstimator(gamma_min=25, bkg_syst_fraction=0.075)
    table = sens.run(dataset_on_off)

    assert len(table) == 4
    assert table.colnames == [
        "e_ref",
        "e_min",
        "e_max",
        "e2dnde",
        "excess",
        "background",
        "criterion",
    ]
    assert table["e_ref"].unit == "TeV"
    assert table["e2dnde"].unit == "erg / (cm2 s)"

    row = table[0]
    assert_allclose(row["e_ref"], 1.33352, rtol=1e-3)
    assert_allclose(row["e2dnde"], 2.74559e-08, rtol=1e-3)
    assert_allclose(row["excess"], 270000, rtol=1e-3)
    assert_allclose(row["background"], 3.6e06, rtol=1e-3)
    assert row["criterion"] == "bkg"

    row = table[1]
    assert_allclose(row["e_ref"], 2.37137, rtol=1e-3)
    assert_allclose(row["e2dnde"], 6.04795e-11, rtol=1e-3)
    assert_allclose(row["excess"], 334.454, rtol=1e-3)
    assert_allclose(row["background"], 3600, rtol=1e-3)
    assert row["criterion"] == "significance"

    row = table[3]
    assert_allclose(row["e_ref"], 7.49894, rtol=1e-3)
    assert_allclose(row["e2dnde"], 1.42959e-11, rtol=1e-3)
    assert_allclose(row["excess"], 25, rtol=1e-3)
    assert_allclose(row["background"], 3.6, rtol=1e-3)
    assert row["criterion"] == "gamma"


def test_integral_estimation(spectrum_dataset):
    dataset = spectrum_dataset.to_image()
    geom = dataset.background.geom

    dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=dataset,
        acceptance=RegionNDMap.from_geom(geom=geom, data=1),
        acceptance_off=RegionNDMap.from_geom(geom=geom, data=5),
    )

    sens = SensitivityEstimator(gamma_min=25, bkg_syst_fraction=0.075)
    table = sens.run(dataset_on_off)
    flux_points = FluxPoints.from_table(
        table, sed_type="e2dnde", reference_model=sens.spectral_model
    )

    assert_allclose(table["excess"].data.squeeze(), 270540, rtol=1e-3)
    assert_allclose(flux_points.flux.data.squeeze(), 7.52e-9, rtol=1e-3)


def test_parameter_sensitivity_estimator(spectrum_dataset):
    spectral_model = PowerLawSpectralModel()
    default_value = spectral_model.amplitude.value

    spectrum_dataset.models = SkyModel(spectral_model=spectral_model)
    datasets = Datasets(spectrum_dataset)

    estimator = ParameterSensitivityEstimator(spectral_model.amplitude, 0, rtol=1e-2)

    value = estimator.run(datasets)
    assert_allclose(value, 4.67553e-12, rtol=1e-2)

    assert_allclose(spectral_model.amplitude.value, default_value, rtol=1e-2)


def test_warning_for_bad_model(caplog, spectrum_dataset):
    geom = spectrum_dataset.background.geom
    spectral_model = PowerLawSpectralModel()
    spectral_model.amplitude.value = -3e-14
    spectrum_dataset.models = SkyModel(spectral_model=spectral_model)

    dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=spectrum_dataset,
        acceptance=RegionNDMap.from_geom(geom=geom, data=1),
        acceptance_off=RegionNDMap.from_geom(geom=geom, data=5),
    )

    estimator = SensitivityEstimator(
        gamma_min=25, bkg_syst_fraction=0.075, spectral_model=spectral_model
    )
    estimator.run(dataset_on_off)
    assert "WARNING" in [_.levelname for _ in caplog.records]
    assert (
        "Spectral model predicts negative flux. Results of estimator should be interpreted with caution"
        in [_.message for _ in caplog.records]
    )
