#!/usr/bin/env python

import create_type1
import create_type3

vex, ctrl, rootid = create_type1.parse_args()
create_type1.process_job(vex, ctrl, rootid)
create_type3.process_job(vex, ctrl, rootid)
