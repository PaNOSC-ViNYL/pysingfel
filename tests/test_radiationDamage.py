import unittest
import os, shutil

from pysingfel.radiationDamage import *


class radiationDamageTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Setting up the test class."""
        os.mkdir("diffr_out")

    @classmethod
    def tearDownClass(cls):
        """ Tear down the test class."""
        shutil.rmtree("diffr_out")

    def test_diffraction_calculation(self):
        """
        Test diffraction calculation using just one process.
        """
        diffraction_parameters = {'uniform_rotation': True,
                                  'calculate_Compton': False,
                                  'slice_interval': 100,
                                  'number_of_slices': 100,
                                  'pmi_start_ID': 1,
                                  'pmi_stop_ID': 1,
                                  'number_of_diffraction_patterns': 1,
                                  'beam_parameter_file': 'test_files/s2e.beam',
                                  'beam_geometry_file': 'test_files/s2e.geom',
                                  'number_of_MPI_processes': 1}

        uniform_rotation = diffraction_parameters['uniform_rotation']
        calculate_Compton = diffraction_parameters['calculate_Compton']
        slice_interval = diffraction_parameters['slice_interval']
        number_of_slices = diffraction_parameters['number_of_slices']
        pmi_start_ID = diffraction_parameters['pmi_start_ID']
        pmi_stop_ID = diffraction_parameters['pmi_stop_ID']
        number_of_diffraction_patterns = diffraction_parameters['number_of_diffraction_patterns']
        beam_parameter_file = diffraction_parameters['beam_parameter_file']
        beam_geometry_file = diffraction_parameters['beam_geometry_file']

        parameters = {'inputDir': 'pmi',
                      'outputDir': 'diffr_out',
                      'beamFile': beam_parameter_file,
                      'geomFile': beam_geometry_file,
                      'uniformRotation': uniform_rotation,
                      'calculateCompton': calculate_Compton,
                      'sliceInterval': slice_interval,
                      'numSlices': number_of_slices,
                      'pmiStartID': pmi_start_ID,
                      'pmiEndID': pmi_stop_ID,
                      'numDP': number_of_diffraction_patterns,
                      'rotationAxis': 'xyz'}

        diffract(parameters)

        # Check expected files exist.
        self.assertTrue(os.path.isdir(os.path.abspath('diffr_out')))
        self.assertIn('diffr_out_0000001.h5', os.listdir(os.path.abspath('diffr_out')))

if __name__ == '__main__':
    unittest.main()

