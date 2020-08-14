import unittest
import fluordynamics as fd
import os

class ImportTest(unittest.TestCase):

    def testTrajectory_rkappa(self):
        traj = fd.fluorburst.Trajectory.from_file(rkappa_filename=os.path.join(fd._TEST_DIR, 'testdata', 'test_rkappa1.dat'))
        self.assertEqual(traj.time[0], 0)

    def testTrajectory_dyecoords(self):
        traj = fd.fluorburst.Trajectory.from_file(rkappa_filename=os.path.join(fd._TEST_DIR, 'testdata', 'test_rkappa1.dat'), 
                                              don_coords_filename=os.path.join(fd._TEST_DIR, 'testdata', 'test_donorcoords1.xvg'),
                                              acc_coords_filename=os.path.join(fd._TEST_DIR, 'testdata', 'test_acceptorcoords1.xvg'))
        self.assertEqual(traj.time[0], 0)


    def testSpecies(self):
        species = fd.fluorburst.Species(name='all', probability=1, filelist_rkappa=[os.path.join(fd._TEST_DIR, 'testdata', 'test_rkappa1.dat'), os.path.join(fd._TEST_DIR, 'testdata', 'test_rkappa2.dat')])
        self.assertEqual(species.name, 'all')
        self.assertEqual(species.trajectories[0].time[0], 0)
        self.assertAlmostEqual(species.trajectories[0].weight, 0.4, places=3)

    def testEnsemble(self):
        ensemble = fd.fluorburst.Ensemble(os.path.join(fd._TEST_DIR, 'testdata'), {'species': {'name': ["all"], "unix_pattern_rkappa": ["*.dat"], "probability": [1]}})
        self.assertEqual(ensemble.species[0].trajectories[0].time[0], 0)


class IntegrationTest(unittest.TestCase):

    def setUp(self):
        self.parameters = fd.fluorburst.readParameters(os.path.join(fd._TEST_DIR, 'test_parameters.json'))

    def test_singleCore_withoutAnisotropy(self):
        self.parameters['sampling']['multiprocessing'] = False
        self.experiment = fd.fluorburst.Experiment(os.path.join(fd._TEST_DIR, 'testdata'), self.parameters, binwidth=0.025, verbose=False, compute_anisotropy=False)

    def test_multiCore_withoutAnisotropy(self):
        self.parameters['sampling']['multiprocessing'] = True
        self.experiment = fd.fluorburst.Experiment(os.path.join(fd._TEST_DIR, 'testdata'), self.parameters, binwidth=0.025, verbose=False, compute_anisotropy=False)

    # def test_multiCore_withAnisotropy(self):
    #     self.parameters['sampling']['multiprocessing'] = True
    #     self.experiment = fd.fluorburst.Experiment(os.path.join(fd._TEST_DIR, 'testdata'), self.parameters, binwidth=0.025, verbose=False, compute_anisotropy=True)


if __name__ == "__main__":
    unittest.main()
