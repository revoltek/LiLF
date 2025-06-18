"""
This code is based on MockPyStep.py from https://dp3.readthedocs.io/en/latest/steps/PythonStep.html
Polarization conversion implemented by Jurjen de Jong.

To use this script in DPPP/DP3, you can save this file in the same folder as from where you run your DPPP command.
Or you save this python file somewhere you like and run:
export PYTHONPATH=/somewhere/you/like:$PYTHONPATH

"""
from subprocess import check_output
import re
import numpy as np
import shutil
import sys

#hacky way to figure out the DPPP/DP3 version (important to run this script properly)
DP3name = shutil.which('DP3')
if not DP3name:
    DP3name = shutil.which('DPPP')
try:
    rgx = r'[0-9]+(\.[0-9]+)+'
    grep_version_string = str(check_output(DP3name+' --version', shell=True), 'utf-8')
    DP3_VERSION = float(re.search(rgx, grep_version_string).group()[0:3])
except AttributeError:
    print('WARNING: grep of DP3 version failed.')
    DP3_VERSION=0.0 # set to default

if DP3_VERSION > 5.3:
    from dp3 import Fields

try:
    from dppp import DPStep as Step
except:
    if DP3_VERSION >= 6:
        from dp3.pydp3 import Step
    else:
        from dp3 import Step

class PolConv(Step):
    """
    Convert UV data polarization.
    lin2circ --> convert from linear to circular UV data
    circ2lin --> convert from circular to linear UV data

    ----------------------------------------------------------------
    Text from: George Heald (2010) and Reinout van Weeren (2014) from lin2circ.py code

    Assuming that the definitions follow the IAU ones, which are
    described by Hamaker & Bregman (1996, A&AS, 117, 161).

    In particular, we use the coordinate transformation given in Section 3,

              1     /  1   +i  \
    C_A =  -------  |          |
           sqrt(2)  \  1   -i  /


    The Hermitian conjugate of this is,

               1     /   1   1  \
    C+_A =  -------  |          |
            sqrt(2)  \  -i  +i  /

    So V_RL = C_A * V_XY * C+_A

    where V_XY is the visibilities in linear coordinates,
    V_RL is the visibilities in circular coordinates.

    Since the transformation matrices are hermitian, C_A ^-1 = C+_A and C+_A^-1 = C_A.
    So, we have:

    V_XY = C+_A*V_RL*C_A
    ----------------------------------------------------------------
    """

    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.
        Set fetch_weights to True if the weights need to be read.
        Similarly for fetch_uvw.

        Args:
          parset: Parameter set for the entire pipeline
          prefix: Prefix for this step."
        """

        super().__init__()

        try:
            self.lin2circ = bool(parset.getInt(prefix + "lin2circ"))
        except RuntimeError:
            self.lin2circ = False
        except AttributeError:
            # DP3 Python bindings have been renamed.
            try:
                self.lin2circ = bool(parset.get_int(prefix + "lin2circ"))
            except RuntimeError:
                self.lin2circ = False

        try:
            self.circ2lin = bool(parset.getInt(prefix + "circ2lin"))
        except RuntimeError:
            self.circ2lin=False
        except AttributeError:
            # DP3 Python bindings have been renamed.
            try:
                self.circ2lin = bool(parset.get_int(prefix + "circ2lin"))
            except RuntimeError:
                self.circ2lin=False

        if self.lin2circ and self.circ2lin:
            sys.exit("Cannot do both lin2circ and circ2lin."
                     "\nChoose:"
                     "\npystep.circ2lin=1 pystep.lin2circ=0"
                     "\nOR"
                     "\npystep.circ2lin=0 pystep.lin2circ=0")

        elif not self.lin2circ and not self.circ2lin:
            sys.exit("Cannot do both lin2circ and circ2lin."
                     "\nChoose:"
                     "\npystep.circ2lin=1 pystep.lin2circ=0"
                     "\nOR"
                     "\npystep.circ2lin=0 pystep.lin2circ=0")

        self.fetch_uvw = True

    def get_required_fields(self):
        if DP3_VERSION>5.3:
            return (Fields.DATA | Fields.FLAGS | Fields.WEIGHTS | Fields.UVW)
        else:
            pass

    def get_provided_fields(self):
        if DP3_VERSION>5.3:
            return Fields()
        else:
            pass

    def update_info(self, dpinfo):
        """
        Process metadata. This will be called before any call to process.

        Args:
          dpinfo: DPInfo object with all metadata, see docs in pydp3.cc
        """

        super().update_info(dpinfo)

        # Make sure data is read
        self.info().set_need_vis_data()

        # Make sure data will be written
        self.info().set_write_data()

    def show(self):
        """Print a summary of the step and its settings"""

        print("\nPolConv")
        if self.lin2circ:
            print("\nConverting UV data polarization from linear to circular\n")
        elif self.circ2lin:
            print("\nConverting UV data polarization from circular to linear\n")


    def process(self, dpbuffer):
        """
        Process one time slot of data. This function MUST call process_next_step.

        Args:
          dpbuffer: DPBuffer object which can contain data, flags and weights
                    for one time slot.
        """

        data = np.array(dpbuffer.get_data(), copy=False, dtype=np.complex64)


        if self.circ2lin:
            """
            circ2lin
            XX =   RR +  RL +  LR +  LL
            XY =  iRR - iRL + iLR - iLL
            YX = -iRR - iRL + iLR + iLL
            YY =   RR -  RL -  LR +  LL
            """

            newdata = np.transpose(np.array([
                0.5 * (data[:, :, 0] + data[:, :, 1] + data[:, :, 2] + data[:, :, 3]),
                0.5 * (1j * data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] - 1j * data[:, :, 3]),
                0.5 * (-1j * data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] + 1j * data[:, :, 3]),
                0.5 * (data[:, :, 0] - data[:, :, 1] - data[:, :, 2] + data[:, :, 3])]),
                (1, 2, 0))

        elif self.lin2circ:
            """
            lin2circ
            RR = XX - iXY + iYX + YY
            RL = XX + iXY + iYX - YY
            LR = XX - iXY - iYX - YY
            LL = XX + iXY - iYX + YY
            """

            newdata = np.transpose(np.array([
                0.5 * (data[:, :, 0] - 1j * data[:, :, 1] + 1j * data[:, :, 2] + data[:, :, 3]),
                0.5 * (data[:, :, 0] + 1j * data[:, :, 1] + 1j * data[:, :, 2] - data[:, :, 3]),
                0.5 * (data[:, :, 0] - 1j * data[:, :, 1] - 1j * data[:, :, 2] - data[:, :, 3]),
                0.5 * (data[:, :, 0] + 1j * data[:, :, 1] - 1j * data[:, :, 2] + data[:, :, 3])]),
                (1, 2, 0))

        data *= 0 # trick to change the UV data
        data += newdata

        # Send processed data to the next step
        if DP3_VERSION>5.3:
            next_step = self.get_next_step()
            if next_step is not None:
                next_step.process(dpbuffer)
        else:
            self.process_next_step(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """

        if self.circ2lin:
            print('Converted UV data from circular (RR,RL,LR,LL) to linear polarization (XX,XY,YX,YY)')
        elif self.lin2circ:
            print('Converted UV data from linear (XX,XY,YX,YY) to circular polarization (RR,RL,LR,LL)')
        pass