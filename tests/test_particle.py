import unittest
import numpy as np

from pysingfel.radiationDamageMPI import *
from pysingfel.particle import symmpdb, merge_chains, parse_symmetry_operations, apply_symmetry_operations, symmetrize_chain

from Bio import PDB


class particleTests(unittest.TestCase):

    def setUp(self):
        parser = PDB.PDBParser()
        self.__structure = parser.get_structure('structure', 'test_files/1uf2.pdb')

    def test_readPDB(self):
        p = Particle()
        p.readPDB('test_files/5g3x.pdb', ff='WK')

        self.assertTrue(isinstance(p.ffTable, np.ndarray))


    def test_merge_chains(self):
        """ Test that merging a number of pdb chains into one chain works. """

        structure = self.__structure

        chains = [c for c in structure.get_chains()]
        print len(chains)

        print len([a for a in structure.get_atoms()])
        merged_chain_atoms_list = merge_chains(chains)

        self.assertEqual( len([a for a in structure.get_atoms()]), len(merged_chain_atoms_list) )


    def test_calFromPDB(self):
        p = Particle()
        p.readPDB('test_files/5g3x.pdb', ff='WK')

        # generate and apply rotation to particle
        quaternion = getRandomRotation('z')
        rotateParticle(quaternion, p)

        geomFile = 'test_files/s2e.geom'
        beamFile = 'test_files/s2e_readPDB.beam'

        det = Detector(geomFile)  # read geom file
        beam = Beam(beamFile)  # read beam file (note that beam info must exist for reading pdb)
        det.init_dp(beam)  # initialize diffraction pattern

        detector_intensity = calculate_molecularFormFactorSq(p, det)
        # Using LCLS AMO beamline fluence: 8.6480e+05 J/cm^2 (assuming pulse energy 4.1 mJ)
        detector_intensity *= det.solidAngle * det.PolarCorr * 8.6480e9 / 1.6e-19
        detector_counts = convert_to_poisson(detector_intensity)

        self.assertTrue(isinstance(detector_counts, np.ndarray))

if __name__ == '__main__':
    unittest.main()
