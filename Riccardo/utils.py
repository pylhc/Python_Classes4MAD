import os

def parsedate(s):
  cmd="date -d '%s' " % s
  cmd+="+'%Y-%m-%d %H:%M:%S'"
  date=os.popen(cmd).readline()
  return date.strip()

def today():
  return parsedate('today')
