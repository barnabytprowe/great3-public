"""@file test_offsets.py
Module containing some basic routines to check that the offsets assumed by Melanie's code match the
offsets contained in truth catalogues.
"""
import os
import sys
import numpy as np
sys.path.append("..")
import pipelines
sys.path.append(os.path.join("..", ".."))
import great3sims.mapper
sys.path.append(os.path.join("..", "..", "public-scripts"))
import branch_id
sys.path.append(os.path.join("..", "..", "server", "great3"))
import evaluate

TRUTH_DIR = "/Users/browe/great3/truth"
NFIELDS = 10
NSUBFIELDS = 200

def melanie_identifier(subfield_index):
    """Get the 9-digit identifier integer according to Melanie's rules in presubmission.py.
    """
    return

def print_offsets_times_70(experiment, obs_type, truth_dir=TRUTH_DIR):
    """Prints the offsets multiplied by 70 for comparison with branch_id x_offset and y_offset
    values as per Melanie's suggestion on Issue #16.
    """
    subfield_index, x_offset_deg, y_offset_deg = evaluate.get_generate_variable_offsets(
        experiment, obs_type, truth_dir=truth_dir)
    print "70. * x_offset ="
    print x_offset_deg * 70.
    print "70. * y_offset ="
    print y_offset_deg * 70.
    return

if __name__ == "__main__":

    experiment = "control"

    # First test the offsets in branch_id match those in the offset_truth files
    for obs_type in ("ground", "space"):

        subfield_index_truth, x_offset_truth, y_offset_truth = \
            evaluate.get_generate_variable_offsets(experiment, obs_type, truth_dir=TRUTH_DIR)
        x_offset_presub = branch_id.x_offset[experiment+"-"+obs_type+"-variable"]
        y_offset_presub = branch_id.y_offset[experiment+"-"+obs_type+"-variable"]
        # Assert that the truth * 70 = presub version
        np.testing.assert_array_equal(
            (70. * x_offset_truth).astype(int), np.asarray(x_offset_presub),
            err_msg="Truth x_offsets do not match those in public-scripts.branch_id")
        np.testing.assert_array_equal(
            (70. * y_offset_truth).astype(int), np.asarray(y_offset_presub),
            err_msg="Truth x_offsets do not match those in public-scripts.branch_id")

        # Then try testing the whole x-y using the (no longer) hacked presubmission
        pipelines.build_submission(
            "im3shape-1", experiment, obs_type, "variable",
            submission_dir=os.path.join("..", "submissions"),
            presubmission_exec=os.path.join("..", "..", "public-scripts", "presubmission.py"))

        # Then build the truth catalogues saving the x and y
        field_index, theta, map_E, map_B, maperr = evaluate.get_generate_variable_truth(
            experiment, obs_type, truth_dir=TRUTH_DIR,
            corr2_params=os.path.join("..", "..", "server", "great3", "corr2.params"),
            make_plots=False, output_xy_prefix="./cats/truth_xy_"+experiment+"-"+obs_type)

        try:
            xysubfile = "./cats/test_xy_"+experiment+"-"+obs_type+"-variable"+(
                "-sub%03d" % 0)+".asc"
            np.loadtxt(xysubfile)
        except IOError as err:
            # Now the hacked presubmission has moved all of these cataloges to ".", and we want to
            # move them to ./cats, so do this
            print "Moving cats to ./cats"
            import subprocess
            import glob
            call_list = (["mv",]+glob.glob("test*.asc")+["./cats"])
            print call_list
            retcode = subprocess.call(["mv",]+glob.glob("truth*.asc")+["./cats"])

        # Then check each subfield
        for isub in range(NSUBFIELDS):
            xysubfile = "./cats/test_xy_"+experiment+"-"+obs_type+"-variable"+(
                "-sub%03d" % isub)+".asc"
            xytruthfile = "./cats/truth_xy_"+experiment+"-"+obs_type+("-sub%03d" % isub)+".asc"
            print "Comparing "+xysubfile+" to "+xytruthfile
            subdata = np.loadtxt(xysubfile)
            truthdata = np.loadtxt(xytruthfile)
            try:
                np.testing.assert_array_almost_equal(
                    subdata, truthdata, decimal=4,
                    err_msg="Barney's truth and presubmission-generated x,y arrays do not match "+
                    "at 3 d.p. for subfield = "+str(isub))
            except AssertionError as err:
                print str(err)
                pass

    #print_offsets_times_70("control", "ground")
    #print_offsets_times_70("control", "space")
    #print "All tests pass so far..."
