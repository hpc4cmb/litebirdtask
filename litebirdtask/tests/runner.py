# Copyright (c) 2015-2023 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import os
import shutil
import sys
import unittest

from toast.mpi import MPI, use_mpi
from toast.tests.mpi import MPITestRunner

from . import ops_sim_observe as test_ops_sim_observe


def test(name=None, verbosity=2):
    # We run tests with COMM_WORLD if available
    comm = None
    rank = 0
    if use_mpi:
        comm = MPI.COMM_WORLD
        rank = comm.rank

    outdir = "litebirdtask_test_output"

    if rank == 0:
        outdir = os.path.abspath(outdir)
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.makedirs(outdir)

    if comm is not None:
        outdir = comm.bcast(outdir, root=0)

    # Run python tests.

    loader = unittest.TestLoader()
    mpirunner = MPITestRunner(verbosity=verbosity, warnings="ignore")
    suite = unittest.TestSuite()

    if name is None:
        suite.addTest(loader.loadTestsFromModule(test_ops_sim_observe))
    else:
        modname = "litebirdtask.tests.{}".format(name)
        if modname not in sys.modules:
            result = f"'{name}' is not a valid test.  Try"
            for name in sys.modules:
                if name.startswith("litebirdtask.tests."):
                    result += '\n  - "{}"'.format(name.replace("litebirdtask.tests.", ""))
            result += "\n"
            raise RuntimeError(result)
        suite.addTest(loader.loadTestsFromModule(sys.modules[modname]))

    ret = 0
    _ret = mpirunner.run(suite)
    if not _ret.wasSuccessful():
        ret += 1

    if comm is not None:
        ret = comm.allreduce(ret, op=MPI.SUM)

    if ret > 0:
        print(f"{ret} Processes had failures")
        sys.exit(6)

    return ret
