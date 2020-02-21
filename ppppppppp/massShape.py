# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Classes implementing the resonance amplitudes.

The variable naming always refers to a decay of a mother M into three final-state
particles 1,2,3. So mMother is the mass of the decaying particle and
fs_masses is a vector containing the three masses of the final state
particles {m1,m2,m3}.

Now we model the process not as direct decay [M -> 1 + 2 + 3], but
through two two-particle decays: [M -> Bachelor + Isobar] and [Isobar ->
daughter1 + daughter2], where the Isobar (an intermediate 2-body
resonance) can either be formed by the final-state particles (12), (13),
or (23). Therefore the daughter masses and the bachelor mass are all
elements of fs_masses, but the assignment depends on which two
particles form the Isobar.
"""


import numpy as np


