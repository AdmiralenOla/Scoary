#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

if (sys.version_info >= (3,0)):
	sys.stdout.write("You seem to be using Python 3. Scoary was written for Python 2.7.x. Try 'python2.7 Scoary.py [arguments]'\n")
	sys.exit(1)
elif (sys.version_info < (2,7)):
	sys.stdout.write("Warning: Scoary was written for Python 2.7.x. Your version is " + str(sys.version) + "\n")
	sys.stdout.write("Please make sure to use Python 2.7 to avoid errors.\n")
	
import Scoary_methods

if __name__ == '__main__':
	Scoary_methods.main()
