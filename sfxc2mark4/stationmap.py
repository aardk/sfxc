station_codes = {}
sparecodes = ["L", "W", "e", "g", "i", "x", "z"]
station_codes["Ai"] = "A"
station_codes["Bd"] = "B"
station_codes["C"] = "Sh"
station_codes["13"] = "D"
station_codes["Wf"] = "E"
station_codes["Eb"] = "F"
station_codes["Ef"] = "F"
station_codes["Gb"] = "G"
station_codes["Ho"] = "H"
station_codes["Ma"] = "I"
station_codes["Cc"] = "J"
station_codes["Kk"] = "K"
station_codes["Mc"] = "M"
station_codes["Md"] = "M"
station_codes["Ny"] = "N"
station_codes["Kb"] = "O"
station_codes["Oh"] = "P"
station_codes["Tc"] = "Q"
station_codes["Zc"] = "R"
station_codes["Nt"] = "S"
station_codes["Ts"] = "T"
station_codes["Ur"] = "U"
station_codes["Wz"] = "V"
station_codes["On"] = "X"
station_codes["O8"] = "X"
station_codes["Yb"] = "Y"
station_codes["Ys"] = "Y"
station_codes["Mh"] = "Z"
station_codes["Ap"] = "a"
station_codes["Br"] = "b"
station_codes["Cm"] = "c"
station_codes["Cn"] = "d"
station_codes["Fd"] = "f"
station_codes["Hn"] = "h"
station_codes["Jc"] = "j"
station_codes["Kp"] = "k"
station_codes["La"] = "l"
station_codes["Mk"] = "m"
station_codes["Nl"] = "n"
station_codes["Ov"] = "o"
station_codes["Pt"] = "p"
station_codes["Qb"] = "q"
station_codes["Ro"] = "r"
station_codes["Sc"] = "s"
station_codes["Ti"] = "t"
station_codes["Ur"] = "u"
station_codes["Pv"] = "v"
station_codes["Wb"] = "w"
station_codes["Y"] = "y"

def create_one_letter_mapping(vex):
  # There are 7 unused station codes, if those are used up we use codes from
  # stations which are not in the current experiment
  stations = []
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    stations.append(st)
  unused = list(set(station_codes.keys()) - set(stations))
  unused = [station_codes[x] for x in unused]
  spares = unused + sparecodes
  stations = {}
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    try:
      code = station_codes[st]
    except KeyError:
      code = spares.pop()
    stations[st] = code
  return stations
