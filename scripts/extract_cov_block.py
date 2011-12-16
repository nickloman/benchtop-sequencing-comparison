import sys

dump = False
for ln in sys.stdin:
	if ln.startswith('>>Per base sequence quality'):
		dump = True
		continue
	if dump and ln.startswith('>>END_MODULE'):
		break	
	if dump:
		print ln, 
	
