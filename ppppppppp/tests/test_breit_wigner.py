from unittest import TestCase

from ppppppppp.massShape import BreitWigner


class BreitWignerTest(TestCase):

    def test_condition_daughter_masses(self):
        """Test that the condition for daughter masses raises an exception, in init.."""
        kwargs = {'name': 'test', 'mass': 1, 'width': 1, 'spin': 1,
                  'mother_mass': 1, 'bachelor_mass': 1, 'daughter_mass1': 1, 'daughter_mass2': 1,
                  'rr': 1, 'rd': 1}

        self.assertRaises(Exception, BreitWigner, **kwargs)

    def test_condition_spin(self):
        """Test that the condition for spin raises an exception."""
        kwargs = {'name': 'test', 'mass': 4, 'width': 1, 'spin': 3,
                  'mother_mass': 1, 'bachelor_mass': 1, 'daughter_mass1': 1, 'daughter_mass2': 1,
                  'rr': 1, 'rd': 1}

        self.assertRaises(Exception, BreitWigner, **kwargs)

    def test_eval_condition_daughter_masses(self):
        """Test that the condition for daughter masses raises an exception, in eval."""

        kwargs = {'name': 'test', 'mass': 4, 'width': 1, 'spin': 1,
                  'mother_mass': 1, 'bachelor_mass': 1, 'daughter_mass1': 1, 'daughter_mass2': 1,
                  'rr': 1, 'rd': 1}

        bw = BreitWigner(**kwargs)

        s12 = 1
        self.assertRaises(Exception, bw.eval, s12)

    def test_condition_mother_bachelor(self):
        """Test that the condition for mother and bachelor mass is fulfilled."""

        kwargs = {'name': 'test', 'mass': 4, 'width': 1, 'spin': 1,
                  'mother_mass': 1, 'bachelor_mass': 1, 'daughter_mass1': 1, 'daughter_mass2': 1,
                  'rr': 1, 'rd': 1}

        bw = BreitWigner(**kwargs)

        s12 = 4
        output = bw.eval(s12)

        self.assertEqual(output, (1+1j))

    def test_eval(self):
        """Test that the eval output/pipeline for arbitrary input values."""
        kwargs = {'name': 'test', 'mass': 4, 'width': 1, 'spin': 1,
                  'mother_mass': 17, 'bachelor_mass': 0.1,
                  'daughter_mass1': 1, 'daughter_mass2': 1,
                  'rr': 1, 'rd': 1}

        bw = BreitWigner(**kwargs)

        s12 = 4
        output = bw.eval(s12)

        self.assertEqual(output, (0.01974374728564802+0j))
