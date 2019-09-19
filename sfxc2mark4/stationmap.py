import json
# Default one letter codes, copied from DiFX2Mark4
default_codes = {}
default_codes["Ai"] = "A"
default_codes["Bd"] = "B"
default_codes["Sh"] = "C"
default_codes["13"] = "D"
default_codes["Wf"] = "E"
default_codes["Eb"] = "F"
default_codes["Ef"] = "F"
default_codes["Gb"] = "G"
default_codes["Ho"] = "H"
default_codes["Ma"] = "I"
default_codes["Cc"] = "J"
default_codes["Kk"] = "K"
default_codes["Mc"] = "M"
default_codes["Md"] = "M"
default_codes["Ny"] = "N"
default_codes["Kb"] = "O"
default_codes["Oh"] = "P"
default_codes["Tc"] = "Q"
default_codes["Zc"] = "R"
default_codes["Nt"] = "S"
default_codes["Ts"] = "T"
default_codes["Ur"] = "U"
default_codes["Wz"] = "V"
default_codes["On"] = "X"
default_codes["O8"] = "X"
default_codes["Yb"] = "Y"
default_codes["Ys"] = "Y"
default_codes["Mh"] = "Z"
default_codes["Ap"] = "a"
default_codes["Br"] = "b"
default_codes["Cm"] = "c"
default_codes["Cn"] = "d"
default_codes["Fd"] = "f"
default_codes["Hn"] = "h"
default_codes["Jc"] = "j"
default_codes["Kp"] = "k"
default_codes["La"] = "l"
default_codes["Mk"] = "m"
default_codes["Nl"] = "n"
default_codes["Ov"] = "o"
default_codes["Pt"] = "p"
default_codes["Qb"] = "q"
default_codes["Ro"] = "r"
default_codes["Sc"] = "s"
default_codes["Ti"] = "t"
default_codes["Ur"] = "u"
default_codes["Pv"] = "v"
default_codes["Wb"] = "w"
default_codes["Y"] = "y"

one_letter_codes = default_codes

def create_one_letter_mapping(vex, codefile=None):
  global one_letter_codes
  if codefile == None:
    one_letter_codes = default_codes
  else:
    fp = open(codefile, 'r')
    # json.load produces unicode, but we need plain strings
    codes = json.load(fp)
    one_letter_codes = {v[0].encode('utf-8'): v[1].encode('utf-8') for v in codes.items()}
    fp.close()

  # Make a list of spare codes, these will be used for stations
  # which we haven't explicitly assigned a one letter code to. 
  # NB: When we have used up all spare codes, we will use the codes for 
  # stations which are not in the current experiment. 
  allcodes = "abcdefghijklmnopqrstuvwxyz"
  allcodes = set(allcodes + allcodes.upper())
  sparecodes = list(allcodes - set(one_letter_codes.values()))

  stations = []
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    stations.append(st)
  unused = list(set(one_letter_codes.keys()) - set(stations))
  unused = [one_letter_codes[x] for x in unused]
  spares = unused + sparecodes
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    try:
      code = one_letter_codes[st]
    except KeyError:
      code = spares.pop()
      print 'Warning: station {} did not have a default one letter code, '\
            'assigned code: "{}"'.format(st, code)
    one_letter_codes[st] = code
