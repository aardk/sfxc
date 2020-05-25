#!/usr/bin/env python
from multiprocessing.pool import ThreadPool
import subprocess
from functools import partial
from argparse import ArgumentParser
import getpass
import os
import re

def parse_ranksfile(ranksfile):
  lines = open(ranksfile, 'r').readlines()
  r = re.compile(r'(rank )(\d+)(\s*=\s*)(\S+)( slot)')
  ranks = {}
  for line in lines:
    m = r.match(line)
    if m != None:
      g = m.groups()
      rank = int(g[1])
      host = g[3]
      try:
        ranks[host].append(rank)
      except KeyError:
        ranks[host] = [rank]
  return ranks

def parse_options():
  description = "Let all processes in SFXC job dump their state, and fetch the results"
  parser = ArgumentParser(description=description)
  parser.add_argument("ranksfile", help = "The ranksfile for the job for which the state should be dumped")
  parser.add_argument("-d", "--output-dir", default = ".", help = "Directory to which information shouuld be copied")
  current_user = getpass.getuser()
  parser.add_argument("-u", "--user", default = current_user, help = "Owner of SFXC processes for which state should be dumped, default = current_user (" + current_user + ")")
  parser.add_argument("-i", "--identity_file", help = "SSH key to use")
  nproc = 20
  parser.add_argument("-n", "--nproc", type = int, default = nproc, help = "Number of simultaneous processes, default = %d"%nproc)
  args = parser.parse_args()
  if not os.path.isdir(args.output_dir):
    parser.error('Invalid output directory: "' + args.output_dir + '"')
  if args.nproc <= 0:
    parser.error('Invalid number of processes: ' + str(args.nproc))
  ranks = parse_ranksfile(args.ranksfile)
  return ranks, args

def create_logs(host, user, identity_file):
  cmd = "ssh"
  if identity_file:
    cmd += " -i %s"%identity_file
  cmd += " %s@%s killall -u %s -s USR1 sfxc"%(user, host, user)
  proc = subprocess.Popen(cmd.split())
  proc.communicate()

def get_logs(host, ranks, user, identity_file, target_dir):
  cmd = "scp"
  if identity_file:
    cmd += " -i %s"%identity_file
  for rank in ranks[host]:
    cmd += " %s@%s:/tmp/sfxc-%s-rank%d.json"%(user, host, user, rank)
  cmd += " %s"%target_dir
  proc = subprocess.Popen(cmd.split())
  proc.communicate()

if __name__ == "__main__":
  ranks, args = parse_options()
  pool = ThreadPool(processes=args.nproc)
  pool.map(partial(create_logs, 
                   user=args.user, 
                   identity_file=args.identity_file), 
           ranks.keys())
  pool.map(partial(get_logs, 
                   ranks=ranks, 
                   user=args.user, 
                   identity_file=args.identity_file, 
                   target_dir=args.output_dir), 
           ranks.keys())
