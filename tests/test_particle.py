import unittest
import numpy
import copy

from pysingfel.radiationDamageMPI import *
from pysingfel.particle import symmpdb, merge_chains, parse_symmetry_operations, apply_symmetry_operations, transrotate_chain

from Bio import PDB


class particleTests(unittest.TestCase):

    def setUp(self):
        self.__1uf2_fname = 'test_files/1uf2.pdb'
        self.__2bfu_fname = 'test_files/2bfu.pdb'

        parser = PDB.PDBParser()
        self.__2bfu = parser.get_structure('structure', self.__2bfu_fname)
        self.__1uf2 = parser.get_structure('structure', self.__1uf2_fname)

    def test_readPDB(self):
        p = Particle()
        p.readPDB('test_files/5g3x.pdb', ff='WK')

        self.assertTrue(isinstance(p.ffTable, numpy.ndarray))

    def test_symm_pdb(self):
        """ Test the symmetrization of a pdb entry. """
        fname = self.__2bfu_fname

        symmetrized_sorted_atoms = symmpdb(fname)

        print symmetrized_sorted_atoms[0]


    def test_merge_chains(self):
        """ Test that merging a number of pdb chains into one chain works. """

        structure = self.__1uf2

        chains = [c for c in structure.get_chains()]
        print len(chains)

        print len([a for a in structure.get_atoms()])
        merged_chain_atoms_list = merge_chains(chains)

        self.assertEqual( len([a for a in structure.get_atoms()]), len(merged_chain_atoms_list) )

    def test_parse_symmetry_operations(self):
        """ Test that we are able to read the symmetry operations from a given pdb file. """
        fname = self.__1uf2_fname
        symmetries, translations = parse_symmetry_operations(fname)

        self.assertIsInstance( symmetries, dict )
        self.assertIsInstance( translations, dict )

        self.assertEqual( symmetries.keys(), [('A', 'B', 'P', 'C', 'D', 'Q', 'E', 'F', 'R', 'G', 'H', 'S', 'I', 'J', 'T', 'K')] )

        # Check first two symmetry operations.
        self.assertEqual( symmetries.values()[0][0][0], 1.0 )
        self.assertEqual( symmetries.values()[0][0][1], 0.0 )
        self.assertEqual( symmetries.values()[0][0][2], 0.0 )

        self.assertEqual( symmetries.values()[0][1][0], 0.0 )
        self.assertEqual( symmetries.values()[0][1][1], 1.0 )
        self.assertEqual( symmetries.values()[0][1][2], 0.0 )

        self.assertEqual( symmetries.values()[0][2][0], 0.0 )
        self.assertEqual( symmetries.values()[0][2][1], 0.0 )
        self.assertEqual( symmetries.values()[0][2][2], 1.0 )

        self.assertEqual( symmetries.values()[0][3][0],  0.500000 )
        self.assertEqual( symmetries.values()[0][3][1], -0.309017 )
        self.assertEqual( symmetries.values()[0][3][2], -0.809017 )

        self.assertEqual( symmetries.values()[0][4][0], -0.309017 )
        self.assertEqual( symmetries.values()[0][4][1],  0.809017 )
        self.assertEqual( symmetries.values()[0][4][2], -0.500000 )

        self.assertEqual( symmetries.values()[0][5][0],  0.809017 )
        self.assertEqual( symmetries.values()[0][5][1],  0.500000 )
        self.assertEqual( symmetries.values()[0][5][2],  0.309017 )

        # Check all translation operations.
        self.assertEqual( max(translations.values()[0]), 0.0 )
        self.assertEqual( min(translations.values()[0]), 0.0 )

    def test_symmetrize_chain(self):
        """ Test the symmetrization of a given chain. """
        chain = self.__1uf2.get_chains().next()
        cloned_chain = copy.deepcopy( chain )

        symmetry = numpy.array([[1.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0]])

        translation = numpy.array([0.0, 0.0, 0.0])

        transrotate_chain( cloned_chain, symmetry, translation)
        symmetrized_atoms = cloned_chain.get_atoms()

        for atom in chain.get_atoms():
            transformed_atom = symmetrized_atoms.next()
            self.assertEqual( (atom.get_vector() - transformed_atom.get_vector()).norm(), 0.0 )

        ### Point reflection
        chain = self.__1uf2.get_chains().next()
        cloned_chain = copy.deepcopy( chain )
        symmetry = numpy.array([[-1.0, 0.0, 0.0],
                                [0.0, -1.0, 0.0],
                                [0.0, 0.0, -1.0]])

        transrotate_chain( cloned_chain, symmetry, translation)
        symmetrized_atoms = cloned_chain.get_atoms()

        for atom in chain.get_atoms():
            transformed_atom = symmetrized_atoms.next()
            self.assertEqual( (atom.get_vector() + transformed_atom.get_vector()).norm(), 0.0 )

        #### Rotation about z axis by 45 degrees
        chain = self.__1uf2.get_chains().next()
        cloned_chain = copy.deepcopy( chain )

        reference_angle = math.pi/4.
        c, s = math.cos( reference_angle ), math.sin( reference_angle )

        rotation = numpy.array([[c, -s, 0.0],
                                [s,  c, 0.0],
                                [0.0, 0.0, 1.0]])

        transrotate_chain( cloned_chain, rotation, translation)
        symmetrized_atoms = cloned_chain.get_atoms()

        for atom in chain.get_atoms():
            vector = atom.get_vector().get_array()
            transformed_atom = symmetrized_atoms.next()
            Tvector = transformed_atom.get_vector().get_array()
            vector[2] = 0.0
            Tvector[2] = 0.0
            self.assertAlmostEqual( math.acos( vector.dot(Tvector) / numpy.linalg.norm(vector)/numpy.linalg.norm(Tvector)) , reference_angle)

    def test_apply_symmetry_operations(self):
        """ Test that applying symmetry operations to a given chain works. """

        fname = self.__2bfu_fname
        structure = self.__2bfu
        number_of_atoms_in_asymmetric_unit = len([a for a in structure.get_atoms()])
        chains = structure.get_chains()

        symmetries, translations = parse_symmetry_operations(fname)
        number_of_operations = len(symmetries.values()[0])/3

        atoms = []
        for chain in chains:
            atoms += apply_symmetry_operations(chain, symmetries, translations)

        self.assertEqual( len(atoms), number_of_atoms_in_asymmetric_unit * number_of_operations )


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

        self.assertTrue(isinstance(detector_counts, numpy.ndarray))

if __name__ == '__main__':
    unittest.main()
