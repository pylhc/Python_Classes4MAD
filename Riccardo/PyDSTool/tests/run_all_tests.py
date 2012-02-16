#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      14 Sep 2006
#$Revision: 1.0.3 $
"""
    Run all PyDSTool tests and examples

    Robert Clewley, March 2006.
"""
import os, time

# set this option to False if you have no external compiler present
do_external=True

# Compatible with no external compiler present
flist = ['Symbolic_test', 'interp_vode_test', 'interp_pcwc', 'vode_event_test',
    'mapsystem_test', 'impfn_test', 'expfn_test', 'objectdelete_test',
    'PointInfo_test', 'saveload_test', 'HH_model', 'HH_model_testbounds',
    'vode_withJac_test', 'vode_withJac_Symbolic_test', 
    'HH_loaded', 'fingermodel_auxvartest', 'fingermodel_vode',
    'numeric_diff_test', 'Symbolic_Diff_test', 'IF_model_test',
    'IF_squarespike_model', 'IF_delaynet', 'imp_load_test',
    'joe_pest', 'pest_test1', 'pest_test2', 'pest_test3',
    'ModelSpec_test', 'SLIP_2D_maps', 'PyCont_Brusselator',
    'PyCont_LPNeuron', 'PyCont_Catalytic', 'PyCont_ABCReaction',
    'PyCont_DiscPredPrey']

# Compatible with external compiler, if supported
dopri_list = ['imprecise_event_test', 'interp_dopri_test', 'dopri_event_test',
    'HH_model_Cintegrator', 'HH_loaded_Cintegrator',
    'IF_delaynet_syn', 'pest_test3_Cintegrator', 'CIN',
    'HH_model_Cintegrator_testbounds', 'Dopri_backwards_test']

radau_list = ['test_hybrid_extinputs', 'SLIP_2D_pdc',
              'freefinger_noforce_radau', 'sloppycell_example']

auto_list = ['PyCont_MorrisLecar_TypeII', 'PyCont_HindmarshRose']



def test(flist, errstr=""):
    for f in flist:
        fname = f+'.py'
        print "\n***** Testing script %s ****************************\n"%fname
        try:
            e=os.system('python '+fname)
        except:
            print "\n***** Testing failed on test file %s"%fname
            raise
        else:
            if e==0 or e==3:
                time.sleep(2)
                print "\n"
            else:
                print "\n***** Testing failed on test file %s.py"%f
                print errstr
                raise RuntimeError

# ---------------------------------------------------------------------------

print "***** Running all tests in order...\n"
print "Note: Depending on your settings, you may have to close matplotlib windows by hand in order to continue to the next test script\n"

test(flist)

if do_external:
    print "\n***** Now running tests that use external compilers...\n"
    test(dopri_list, "Dopri ODEsystems do not seem to work with your setup.")
    test(auto_list, "The interface to AUTO does not seem to work with your setup.")
    test(radau_list, "The Radau ODEsystem does not seem to work with your setup.")

print "\n***** All testing appears to have been successful"
