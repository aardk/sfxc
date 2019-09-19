#!/usr/bin/env python
import create_type1
import create_type3
import stationmap

vex, ctrl, rootid, ovexfile, codefile = create_type1.parse_args()
stationmap.create_one_letter_mapping(vex, codefile)
create_type1.process_job(vex, ctrl, rootid, ovexfile)
create_type3.process_job(vex, ctrl, rootid)
