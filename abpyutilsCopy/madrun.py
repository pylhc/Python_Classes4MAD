import re
import sys
import os
import pexpect

class madrun:
  """
  r=m.madrun();
  r.s('a=4;value,a^2;')
  r.p('a=4;value,a^2;')
  """
  prompt=['\r\nX: ==>\r\n', pexpect.TIMEOUT]
  def __init__(self,prg='madx'):
    self.mad=pexpect.spawn(prg)
    self.mad.expect(self.prompt)
    print self.mad.before;
    self.mad.sendline('option,-echo;')
    self.mad.expect(self.prompt)


  def s(self,string):
    """
    Send command to madx and return string
    print r.s('call,file=job.sample.madx;')
    """
    self.mad.sendline(string)
    self.mad.expect(self.prompt)
    return self.mad.before

  def p(self,string):
    """
    Send command to madx and print string
    r.p('call,file=job.sample.madx;')
    """
    print self.s(string)

  def i(self):
    """
    Enter in madx command line
    r.i()
    """
    readline.parse_and_bind("tab: complete")
    fn1=home+'/.history_madx'
    fn2=home+'/.ipython/history'
    os.system("touch " + fn1)
    os.system("touch " + fn2)
    print "Type Ctrl-D to enter ipython"
    cmd=''
    readline.clear_history()
    readline.read_history_file(fn1)
    print 'X: ==> '
    try:
      while( 1):
        cmd=cmd+re.sub('\!.*','',raw_input())
        if cmd.endswith(';'):
          self.p(cmd)
          cmd=''
          print 'X: ==> '
    except EOFError:
      pass
    readline.write_history_file(fn1)
    readline.clear_history()
    readline.read_history_file(fn2)
    print "Type madx() to enter MadX"

