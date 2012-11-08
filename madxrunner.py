"""
This class helps to run MADX.
It just calls MADX with the given file or string as input.
It can either be loaded as a module or called standalone.

version 2 20101019 (tbach)
 - simplified calling with default arguments
 - added possibility to change stdout

version 1 20101019 (tbach)
 - init
"""

import subprocess
import optparse
import sys

globalMadxPath = "/afs/cern.ch/group/si/slap/bin/madx"

def runForInputFile(filePath, madxPath=globalMadxPath, stdout = None):
    """
    Run MADX with the given file as input, returns the errorcode from MADX as int.
    Arguments:
     - filepath: Path to the input file
     - madxPath: Optional argument to the madx binary, it uses the default madx path given by the global variable globalMadxPath.
     - stdout: Optional argument to specify what should happen with messages to stdout.
               Some of the Possible values:
                   - None: default behaviour, messages are printed to stdout
                   - open(os.devnull, "w"): Suppress any std messages
                   - open(path): Write to file located at path
    """
    fileOptions = open(filePath, "r")
    stringInput = fileOptions.read()
    returnCode = runForInputString(stringInput, madxPath, stdout)
    return returnCode

def runForInputString(stringInput, madxPath=globalMadxPath, stdout = None):
    """
    Run MADX with the given string as input, returns the errorcode from MADX as int.
    Arguments:
     - stringInput: input string
     - madxPath: Optional argument to the madx binary, it uses the default madx path given by the global variable globalMadxPath.
     - stdout: Optional argument to specify what should happen with messages to stdout.
               Some of the Possible values:
                   - None: default behaviour, messages are printed to stdout
                   - open(os.devnull, "w"): Suppress any std messages
                   - open(path): Write to file located at path
    """
    process = subprocess.Popen(madxPath, shell=False, stdin=subprocess.PIPE, stdout=stdout)
    process.communicate(stringInput)
    returnCode = process.wait()
    return returnCode
    
    
if __name__ == "__main__":
    parser = optparse.OptionParser(usage="python %prog -f madx_input_file [other options]", version="%prog 2")
    parser.add_option("-f", "--file",
                    help="Filename or path to file with MADX input",
                    default="./", dest="file")
    parser.add_option("-m", "--madxpath",
                    help="Filename or path to MADX binary, default is: " + globalMadxPath,
                    default=globalMadxPath, dest="madxpath")
    (options, args) = parser.parse_args()
    globalMadxPath = options.madxpath
    
    myReturnCode = runForInputFile(options.file)
    sys.exit(myReturnCode)
