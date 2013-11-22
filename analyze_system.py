import sys_AND
import sys_OR
import sys_NOT
"""
print "AND gate:"
for i in xrange(10):
	print "0.9", 0.1 * i, sys_AND.analyze_system(100000, 10000, 0.9, 0.1 * i, [ "SoutQ" ])[0]
for i in xrange(10):
	pass
	#print "0.1", 0.1 * i, sys_AND.analyze_system(10000, 1000, 0.9, 0.1 * i, [ "SoutQ" ])[0]
"""
"""
print sys_AND.analyze_system(6000, 1000, 0.1, 0.9, [ "SoutQ" ])
print sys_AND.analyze_system(6000, 1000, 0.9, 0.1, [ "SoutQ" ])
print sys_AND.analyze_system(6000, 1000, 0.9, 0.9, [ "SoutQ" ])
print
print "OR gate:"
print sys_OR.analyze_system(6000, 1000, 0.1, 0.1, [ "SoutQ" ])
print sys_OR.analyze_system(6000, 1000, 0.1, 0.9, [ "SoutQ" ])
print sys_OR.analyze_system(6000, 1000, 0.9, 0.1, [ "SoutQ" ])
print sys_OR.analyze_system(6000, 1000, 0.9, 0.9, [ "SoutQ" ])
"""
"""
print "OR gate:"
for i in xrange(10):
	for j in xrange(10):
		print 0.1 * i, 0.1 * j, sys_OR.analyze_system(4*60*60, 1000, 0.1 * i, 0.1 * j, [ "SoutQ" ])[0], " ",
	print
"""
"""
print "AND gate:"
print "",
for j in xrange(10):
	print 0.1 * j, "",
print
for i in xrange(10):
	print 0.1 * i, "",
	for j in xrange(10):
		print sys_AND.analyze_system(4*60*60, 1000, 0.1 * i, 0.1 * j, [ "SoutQ" ])[0],
	print
"""
"""
AND gate:
[5.5694645018831359]
[33.912661214380613]
[33.912661203073547]
[111.96244132318184]

OR gate:
[9.9170786510815923]
[138.24692674600129]
[138.24692723677998]
[149.72267537264725]

real	0m11.724s
user	0m11.705s
sys	0m0.023s

"""
print "NOT gate:"
for i in xrange(10):
	print 0.1 * i , sys_NOT.analyze_system(6000, 1000, 0.1 * i, [ "SoutQ" ])
