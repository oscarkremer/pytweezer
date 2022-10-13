import pytest
import numpy as np
from pytweezer.utils.polarizability import clausius_mossoti, filter_coupled_dipole, lattice_dispersion_relation

def test_clausius_mossoti_polarizability():
    decimal = 3
    spacing = 1.0
    index = np.array([[1, 2]])
    alpha = clausius_mossoti(spacing, index)
    expected = np.array([[0], [0.1194]])
    np.testing.assert_array_almost_equal(alpha, expected, decimal=decimal, err_msg='Clausius Mossoti Polarizability Error')


def test_filter_coupled_dipole_polarizability():
    decimal = 3
    spacing = 1.0
    index = np.array([[1, 2]])
    alpha = filter_coupled_dipole(spacing, index)
    expected = np.array([[0], [0.000029131505380 - 0.003023300651752j]])
    np.testing.assert_array_almost_equal(alpha, expected, decimal=decimal, err_msg='Filter Coupled Dipole Polarizability Error')


def test_lattice_dispersion_relation_multiple_args():
    decimal = 3
    spacing = 1.0;
    index = np.array([1, 2])
    kvec = np.array([0, 0, 1])
    E = np.array([1, 0, 0])
    alpha = lattice_dispersion_relation(spacing, index, kvec=kvec, E=E, k=2*np.pi)
    expected = np.array([[0], [-0.0014 + 0.0057j]])
    np.testing.assert_array_almost_equal(alpha, expected, decimal=decimal, err_msg='Lattice Dispersion Relation Polarizability \
        Error for Multiple Parameters')


def test_lattice_dispersion_relation_polarizability():
    decimal = 3
    spacing = 1.0
    index = np.array([1, 2])
    alpha = lattice_dispersion_relation(spacing, index)
    expected = np.array([[0], [-0.002627924391542 + 0.004518926661470j]])
    np.testing.assert_array_almost_equal(alpha, expected, decimal=decimal, err_msg='Lattice Dispersion Relation Polarizability Error')


def test_lattice_dispersion_relation_args():
    decimal = 3
    spacing = 1.0;
    index = np.array([1, 2])
    kvec = np.array([0, 0, 1])
    E = np.array([1, 0])
    with pytest.raises(ValueError):
        alpha = lattice_dispersion_relation(spacing, index, kvec=kvec, E=E, k=2*np.pi)
    