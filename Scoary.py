#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Scoary - Microbial Pan-GWAS. Associates genes in Roary output with phenotypes.
# Copyright (C) 2016  Ola Brynildsrud (ola.brynildsrud@fhi.no)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
